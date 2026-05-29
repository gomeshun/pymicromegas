from __future__ import annotations

import ctypes
import tempfile
from pathlib import Path

from .base import MicrOmegasProcess


class RelicDensityProcesses(MicrOmegasProcess):
    def _prepare(self, parameters=None):
        if parameters is not None:
            self.model.assign(parameters)
        self.model.set_gauge()
        self.model.clear_feeble_list()
        self.model.sort_odd_particles()

    def get_infl_decay(self, h0, gamma):
        trh = ctypes.c_double()
        tmax = ctypes.c_double()
        aend = ctypes.c_double()
        err = self.lib.pymicromegas_get_infl_decay(
            float(h0), float(gamma), ctypes.byref(trh), ctypes.byref(tmax), ctypes.byref(aend)
        )
        return {"TRH": trh.value, "Tmax": tmax.value, "aEnd": aend.value, "err": err}

    def get_infl_decay_plus(self, hsm, hi, gamma, alpha, omega):
        trh = ctypes.c_double()
        tmax = ctypes.c_double()
        aend = ctypes.c_double()
        err = self.lib.pymicromegas_get_infl_decay_plus(
            float(hsm),
            float(hi),
            float(gamma),
            float(alpha),
            float(omega),
            ctypes.byref(trh),
            ctypes.byref(tmax),
            ctypes.byref(aend),
        )
        return {"TRH": trh.value, "Tmax": tmax.value, "aEnd": aend.value, "err": err}

    def print_thermal_sets(self):
        return self.model._read_bridge_text(self.lib.pymicromegas_print_thermal_sets)

    def freeze_out_channels(self, xf, cut=0.01, beps=1e-4, percent=False):
        tmp_path = None
        try:
            with tempfile.NamedTemporaryFile("w+", delete=False) as tmp_file:
                tmp_path = Path(tmp_file.name)
            total = self.lib.pymicromegas_print_channels(
                str(tmp_path).encode("UTF-8"), float(xf), float(cut), float(beps), int(bool(percent))
            )
            if total < 0:
                raise RuntimeError("Failed to write freeze-out channels.")
            return {"total": total, "text": tmp_path.read_text()}
        finally:
            if tmp_path is not None:
                tmp_path.unlink(missing_ok=True)

    def dark_omega(self, parameters=None, dof_fname=None, fast=1, beps=1e-4):
        if parameters is not None:
            self.model.assign(parameters)
        if dof_fname is not None:
            self.model.load_heff_geff(dof_fname)
        self.model.set_gauge()
        self.model.clear_feeble_list()
        self.model.sort_odd_particles()
        xf = ctypes.c_double()
        err = ctypes.c_int()
        omega = self.lib.pymicromegas_dark_omega(
            ctypes.byref(xf), int(fast), float(beps), ctypes.byref(err)
        )
        return {"Xf": xf.value, "Omega": omega, "err": err.value}

    def dark_omega2(self, parameters=None, fast=1, beps=1e-4):
        self._prepare(parameters=parameters)
        err = ctypes.c_int()
        omega = self.lib.pymicromegas_dark_omega2(int(fast), float(beps), ctypes.byref(err))
        if err.value:
            raise RuntimeError(f"darkOmega2 failed with error code {err.value}.")
        return omega

    def dark_omega_tr(self, parameters=None, tr=1e10, yr=0.0, fast=1, beps=1e-4):
        self._prepare(parameters=parameters)
        err = ctypes.c_int()
        omega = self.lib.pymicromegas_dark_omega_tr(
            float(tr), float(yr), int(fast), float(beps), ctypes.byref(err)
        )
        return {"Omega": omega, "TR": float(tr), "YR": float(yr), "err": err.value}

    def dark_omega_fo(self, parameters=None, fast=1, beps=1e-4):
        self._prepare(parameters=parameters)
        xf = ctypes.c_double()
        err = ctypes.c_int()
        omega = self.lib.pymicromegas_dark_omega_fo(
            ctypes.byref(xf), int(fast), float(beps), ctypes.byref(err)
        )
        return {"Xf": xf.value, "Omega": omega, "err": err.value}

    def dark_omega2_tr(self, parameters=None, tr=1e10, y1r=0.0, y2r=0.0, fast=1, beps=1e-4):
        self._prepare(parameters=parameters)
        err = ctypes.c_int()
        omega = self.lib.pymicromegas_dark_omega2_tr(
            float(tr), float(y1r), float(y2r), int(fast), float(beps), ctypes.byref(err)
        )
        return {"Omega": omega, "TR": float(tr), "Y1R": float(y1r), "Y2R": float(y2r), "err": err.value}

    def dark_omega_n(self, parameters=None, fast=1, beps=1e-4):
        self._prepare(parameters=parameters)
        err = ctypes.c_int()
        omega = self.lib.pymicromegas_dark_omega_n(int(fast), float(beps), ctypes.byref(err))
        return {"Omega": omega, "err": err.value}

    def dark_omega_n_tr(self, y_initial, parameters=None, tr=1e10, fast=1, beps=1e-4):
        self._prepare(parameters=parameters)
        values = [float(value) for value in y_initial]
        n_cdm = max(int(self.model.n_cdm), 0)
        if n_cdm and len(values) == n_cdm:
            values = [0.0] + values
        array_size = max(len(values), n_cdm + 1)
        y_array = (ctypes.c_double * array_size)()
        for index, value in enumerate(values):
            y_array[index] = value
        err = ctypes.c_int()
        omega = self.lib.pymicromegas_dark_omega_n_tr(
            float(tr), y_array, int(fast), float(beps), ctypes.byref(err)
        )
        result = {"Omega": omega, "TR": float(tr), "Y": list(y_array), "err": err.value}
        if n_cdm:
            result["Y_by_sector"] = list(y_array)[1 : n_cdm + 1]
        return result

    def dark_omega_infl(self, branching, parameters=None, beps=1e-4):
        self._prepare(parameters=parameters)
        tfo = ctypes.c_double()
        err = ctypes.c_int()
        omega = self.lib.pymicromegas_dark_omega_infl(
            float(branching), float(beps), ctypes.byref(tfo), ctypes.byref(err)
        )
        return {"Omega": omega, "Tfo": tfo.value, "branching": float(branching), "err": err.value}