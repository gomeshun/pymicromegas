from collections.abc import Mapping
import ctypes
import os
from pathlib import Path


BRIDGE_C = r"""
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../include/micromegas.h"
#include "../include/micromegas_aux.h"
#include "lib/pmodel.h"

int pymicromegas_assign_values(int n, const char **names, const double *values)
{
    int err = 0;
    for (int i = 0; i < n; ++i) {
        err = assignVal(names[i], values[i]);
        if (err) return i + 1;
    }
    return 0;
}

int pymicromegas_sort_odd_particles(char *name, int name_size)
{
    char cdm_name[128] = "";
    int err = sortOddParticles(cdm_name);
    if (name && name_size > 0) {
        strncpy(name, cdm_name, name_size - 1);
        name[name_size - 1] = '\0';
    }
    return err;
}

double pymicromegas_dark_omega(double *xf, int fast, double beps, int *err)
{
    return darkOmega(xf, fast, beps, err);
}

double pymicromegas_dark_omega2(int fast, double beps)
{
    return darkOmega2(fast, beps);
}

double pymicromegas_vsigma(double temperature, double beps, int fast)
{
    return vSigma(temperature, beps, fast);
}

int pymicromegas_load_heff_geff(const char *fname)
{
    return loadHeffGeff((char *)fname);
}

double pymicromegas_find_val(const char *name, int *err)
{
    double value = 0.0;
    int status = findVal((char *)name, &value);
    if (err) *err = status;
    return value;
}
"""


BRIDGE_MAKEFILE = r"""
.PHONY: libs

AllFlags = ../CalcHEP_src/FlagsForMake
ifeq (,$(wildcard $(AllFlags) ))
$(error File $(AllFlags) is absent. Compile micrOMEGAs first)
endif
include ../CalcHEP_src/FlagsForMake

cLib = $(CALCHEP)/lib
SSS = $(wildcard lib/*.a) ../lib/micromegas.a $(cLib)/dynamic_me.a ../lib/micromegas.a \
 work/work_aux.a $(wildcard lib/*.a) $(cLib)/sqme_aux.$(SO) $(cLib)/libSLHAplus.a \
 $(cLib)/num_c.a $(cLib)/serv.a $(cLib)/ntools.a $(LX11)

ifneq ($(LHAPDFPATH),)
  SSS += -L$(LHAPDFPATH) -lLHAPDF $(cLib)/dummy.a
  DLSET = export LD_RUN_PATH=$(LHAPDFPATH);
else
  SSS += $(cLib)/dummy.a
  DLSET =
endif

pymicromegas_native.so: pymicromegas_native.c libs work/bin
	$(DLSET) $(CC) $(CFLAGS) -fPIC -shared -o $@ pymicromegas_native.c $(SSS) $(lDL) -lm -lpthread

libs:
	$(MAKE) -C work
	$(MAKE) -C lib
	$(MAKE) -C ../sources

work/bin:
	ln -s $(shell pwd)/../CalcHEP_src/bin $(shell pwd)/work/bin
"""


class NativeMicrOmegas:
    """ctypes bridge for a compiled micrOMEGAs project."""

    source_name = "pymicromegas_native.c"
    makefile_name = "pymicromegas_native.mk"
    library_name = "pymicromegas_native.so"

    def __init__(self, project_path, runner=None):
        self.project_path = Path(project_path)
        self.runner = runner
        self.library_path = self.project_path / self.library_name
        self._cdll = None

    @property
    def cdll(self):
        if self._cdll is None:
            self._cdll = ctypes.CDLL(str(self.library_path))
            self._configure_api()
        return self._cdll

    def build(self, force=False):
        self._write_if_changed(self.project_path / self.source_name, BRIDGE_C)
        self._write_if_changed(self.project_path / self.makefile_name, BRIDGE_MAKEFILE)
        if force or not self.library_path.exists():
            if self.runner is None:
                raise RuntimeError("No project runner is available to compile native bridge.")
            process = self.runner(
                ["make", "-f", self.makefile_name, self.library_name],
                shell=False,
                verbose=True,
            )
            if process.returncode:
                raise RuntimeError(process.stdout)
        self._cdll = None
        return self

    def function(self, name, restype=None, argtypes=None):
        """Return a raw ctypes function from the native library."""
        func = getattr(self.cdll, name)
        if restype is not None:
            func.restype = restype
        if argtypes is not None:
            func.argtypes = argtypes
        return func

    def assign_values(self, parameters):
        """Assign values from a dict-like object or pandas.Series."""
        names = list(parameters.keys())
        if isinstance(parameters, Mapping):
            raw_values = parameters.values()
        else:
            raw_values = parameters.values
        values = [float(value) for value in raw_values]
        name_array = (ctypes.c_char_p * len(names))(
            *[name.encode() for name in names]
        )
        value_array = (ctypes.c_double * len(values))(*values)
        err = self.cdll.pymicromegas_assign_values(
            len(names), name_array, value_array
        )
        if err:
            raise RuntimeError(f"micrOMEGAs rejected parameter '{names[err - 1]}'.")

    def sort_odd_particles(self):
        name = ctypes.create_string_buffer(128)
        err = self.cdll.pymicromegas_sort_odd_particles(name, len(name))
        if err:
            raise RuntimeError(f"micrOMEGAs failed to sort odd particles: {err}")
        return name.value.decode()

    def load_heff_geff(self, path):
        err = self.cdll.pymicromegas_load_heff_geff(os.fsencode(path))
        # micrOMEGAs loadHeffGeff returns a positive line count on success.
        if err <= 0:
            raise RuntimeError(f"micrOMEGAs could not load DOF file: {path}")
        return err

    def find_val(self, name):
        err = ctypes.c_int()
        value = self.cdll.pymicromegas_find_val(name.encode(), ctypes.byref(err))
        if err.value:
            raise RuntimeError(f"micrOMEGAs could not find variable '{name}'.")
        return value

    def dark_omega(self, parameters=None, dof_fname=None, fast=1, beps=1.0e-4):
        if parameters is not None:
            self.assign_values(parameters)
        if dof_fname is not None:
            self.load_heff_geff(dof_fname)
        self.sort_odd_particles()
        xf = ctypes.c_double()
        err = ctypes.c_int()
        omega = self.cdll.pymicromegas_dark_omega(
            ctypes.byref(xf), int(fast), float(beps), ctypes.byref(err)
        )
        return {"Xf": xf.value, "Omega": omega, "err": err.value}

    def dark_omega2(self, parameters=None, fast=1, beps=1.0e-4):
        if parameters is not None:
            self.assign_values(parameters)
        self.sort_odd_particles()
        return self.cdll.pymicromegas_dark_omega2(int(fast), float(beps))

    def v_sigma(self, temperature, beps=1.0e-5, fast=1):
        return self.cdll.pymicromegas_vsigma(
            float(temperature), float(beps), int(fast)
        )

    def _configure_api(self):
        lib = self._cdll
        lib.pymicromegas_assign_values.argtypes = [
            ctypes.c_int,
            ctypes.POINTER(ctypes.c_char_p),
            ctypes.POINTER(ctypes.c_double),
        ]
        lib.pymicromegas_assign_values.restype = ctypes.c_int
        lib.pymicromegas_sort_odd_particles.argtypes = [
            ctypes.c_char_p,
            ctypes.c_int,
        ]
        lib.pymicromegas_sort_odd_particles.restype = ctypes.c_int
        lib.pymicromegas_dark_omega.argtypes = [
            ctypes.POINTER(ctypes.c_double),
            ctypes.c_int,
            ctypes.c_double,
            ctypes.POINTER(ctypes.c_int),
        ]
        lib.pymicromegas_dark_omega.restype = ctypes.c_double
        lib.pymicromegas_dark_omega2.argtypes = [ctypes.c_int, ctypes.c_double]
        lib.pymicromegas_dark_omega2.restype = ctypes.c_double
        lib.pymicromegas_vsigma.argtypes = [
            ctypes.c_double,
            ctypes.c_double,
            ctypes.c_int,
        ]
        lib.pymicromegas_vsigma.restype = ctypes.c_double
        lib.pymicromegas_load_heff_geff.argtypes = [ctypes.c_char_p]
        lib.pymicromegas_load_heff_geff.restype = ctypes.c_int
        lib.pymicromegas_find_val.argtypes = [
            ctypes.c_char_p,
            ctypes.POINTER(ctypes.c_int),
        ]
        lib.pymicromegas_find_val.restype = ctypes.c_double

    @staticmethod
    def _write_if_changed(path, content):
        if path.exists() and path.read_text() == content:
            return
        path.write_text(content)
