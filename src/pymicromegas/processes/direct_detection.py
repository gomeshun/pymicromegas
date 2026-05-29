from __future__ import annotations

import ctypes

from .base import MicrOmegasProcess


class DirectDetectionProcesses(MicrOmegasProcess):
    def nucleon_amplitudes(self, wimp):
        proton_scalar = ctypes.c_double()
        proton_axial = ctypes.c_double()
        neutron_scalar = ctypes.c_double()
        neutron_axial = ctypes.c_double()
        function = self.model.function(
            "nucleonAmplitudes",
            restype=ctypes.c_int,
            argtypes=[ctypes.c_char_p, ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double)],
        )
        err = function(
            self.model._encode_string(wimp),
            ctypes.byref(proton_scalar),
            ctypes.byref(proton_axial),
            ctypes.byref(neutron_scalar),
            ctypes.byref(neutron_axial),
        )
        return {
            "err": err,
            "proton_scalar": proton_scalar.value,
            "proton_axial": proton_axial.value,
            "neutron_scalar": neutron_scalar.value,
            "neutron_axial": neutron_axial.value,
        }

    def dnde_recoil(self, energy, recoil_spectrum):
        function = self.model.function(
            "dNdERecoil",
            restype=ctypes.c_double,
            argtypes=[ctypes.c_double, ctypes.POINTER(ctypes.c_double)],
        )
        return function(float(energy), self.model._as_double_array(recoil_spectrum, size=self.model.RECOIL_SIZE))