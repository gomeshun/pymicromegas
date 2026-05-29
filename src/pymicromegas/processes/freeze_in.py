from __future__ import annotations

import ctypes
import tempfile
from pathlib import Path

from .base import MicrOmegasProcess


class FreezeInProcesses(MicrOmegasProcess):
    def v_sigma_plus23(self, process, temperature):
        err = ctypes.c_int()
        value = self.lib.pymicromegas_v_sigma_plus23(
            self.model._encode_string(process), float(temperature), ctypes.byref(err)
        )
        return {"value": value, "err": err.value}

    def v_sigma_plus24(self, process, temperature):
        err = ctypes.c_int()
        value = self.lib.pymicromegas_v_sigma_plus24(
            self.model._encode_string(process), float(temperature), ctypes.byref(err)
        )
        return {"value": value, "err": err.value}

    def y_freeze_in_22(self, process, t0, tr=1e10, plot_dydt=False):
        err = ctypes.c_int()
        value = self.lib.pymicromegas_y_freeze_in22(
            self.model._encode_string(process), float(t0), float(tr), int(bool(plot_dydt)), ctypes.byref(err)
        )
        return {"Y": value, "T0": float(t0), "TR": float(tr), "process": str(process), "err": err.value}

    def prepare_freeze_in(self, parameters=None, particle_name=None, sector=1, reset_feeble=True):
        if parameters is not None:
            self.model.assign(parameters)
        self.model.set_gauge()
        sorted_particle_name = self.model.sort_odd_particles()
        if particle_name is None and sorted_particle_name:
            particle_name = sorted_particle_name
        particle_name = self.model._resolve_particle_name(particle_name, sector=sector)
        if reset_feeble:
            self.model.clear_feeble_list()
        self.model.to_feeble_list(particle_name)
        return particle_name

    def freeze_in_channels(self, cut=0.0, percent=False):
        tmp_path = None
        try:
            with tempfile.NamedTemporaryFile("w+", delete=False) as tmp_file:
                tmp_path = Path(tmp_file.name)
            err = self.lib.pymicromegas_print_channels_fi(
                str(tmp_path).encode("UTF-8"), float(cut), int(bool(percent))
            )
            if err:
                raise RuntimeError("Failed to write freeze-in channels.")
            return self.parse_freeze_in_channels(tmp_path.read_text())
        finally:
            if tmp_path is not None:
                tmp_path.unlink(missing_ok=True)

    @staticmethod
    def parse_freeze_in_channels(text):
        channels = []
        for line in text.splitlines():
            terms = line.split()
            if len(terms) < 5 or terms[3] != "->":
                continue
            channels.append(
                {
                    "weight": float(terms[0]),
                    "in": (terms[1].rstrip(","), terms[2]),
                    "out": tuple(term.strip() for term in " ".join(terms[4:]).split(",")),
                }
            )
        return channels

    def dark_omega_freeze_in(
        self,
        parameters=None,
        particle_name=None,
        sector=1,
        tr=1e10,
        reset_feeble=True,
        channels=False,
        channel_cut=0.0,
        channel_percent=False,
    ):
        particle_name = self.prepare_freeze_in(
            parameters=parameters,
            particle_name=particle_name,
            sector=sector,
            reset_feeble=reset_feeble,
        )
        err = ctypes.c_int()
        omega = self.lib.pymicromegas_dark_omega_fi(
            float(tr), particle_name.encode("UTF-8"), ctypes.byref(err)
        )
        result = {"Omega": omega, "TR": float(tr), "particle": particle_name, "err": err.value}
        if channels:
            result["channels"] = self.freeze_in_channels(cut=channel_cut, percent=channel_percent)
        return result

    def dark_omega_freeze_in_22(
        self,
        process,
        parameters=None,
        particle_name=None,
        sector=1,
        tr=1e10,
        reset_feeble=True,
    ):
        particle_name = self.prepare_freeze_in(
            parameters=parameters,
            particle_name=particle_name,
            sector=sector,
            reset_feeble=reset_feeble,
        )
        err = ctypes.c_int()
        omega = self.lib.pymicromegas_dark_omega_fi22(
            float(tr), str(process).encode("UTF-8"), particle_name.encode("UTF-8"), ctypes.byref(err)
        )
        return {
            "Omega": omega,
            "TR": float(tr),
            "particle": particle_name,
            "process": str(process),
            "err": err.value,
        }

    def dark_omega_freeze_in_decay(
        self,
        bath_particle,
        parameters=None,
        particle_name=None,
        sector=1,
        tr=1e10,
        reset_feeble=True,
    ):
        particle_name = self.prepare_freeze_in(
            parameters=parameters,
            particle_name=particle_name,
            sector=sector,
            reset_feeble=reset_feeble,
        )
        omega = self.lib.pymicromegas_dark_omega_fi_decay(
            float(tr), str(bath_particle).encode("UTF-8"), particle_name.encode("UTF-8")
        )
        return {
            "Omega": omega,
            "TR": float(tr),
            "particle": particle_name,
            "bath_particle": str(bath_particle),
        }