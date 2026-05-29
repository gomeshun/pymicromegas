from __future__ import annotations

import ctypes

from .base import MicrOmegasProcess


class IndirectDetectionProcesses(MicrOmegasProcess):
    def calc_spectrum(self, key=0):
        arrays = [self.model._new_double_array() for _ in range(6)]
        err = ctypes.c_int()
        double_pointer = ctypes.POINTER(ctypes.c_double)
        function = self.model.function(
            "calcSpectrum",
            restype=ctypes.c_double,
            argtypes=[ctypes.c_int, double_pointer, double_pointer, double_pointer, double_pointer, double_pointer, double_pointer, ctypes.POINTER(ctypes.c_int)],
        )
        sigma_v = function(int(key), *arrays, ctypes.byref(err))
        names = ["gamma", "positron", "antiproton", "nu_e", "nu_mu", "nu_tau"]
        result = {"sigma_v": sigma_v, "err": err.value}
        result.update({name: self.model._numpy_from_array(array) for name, array in zip(names, arrays)})
        return result

    def calc_spectrum_plus(self, process, out_particle):
        spectrum = self.model._new_double_array()
        err = ctypes.c_int()
        function = self.model.function(
            "calcSpectrumPlus",
            restype=ctypes.c_double,
            argtypes=[ctypes.c_char_p, ctypes.c_int, ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_int)],
        )
        value = function(self.model._encode_string(process), int(out_particle), spectrum, ctypes.byref(err))
        return {"value": value, "err": err.value, "spectrum": self.model._numpy_from_array(spectrum)}

    def decay_spectrum(self, particle_name, out_particle):
        spectrum = self.model._new_double_array()
        function = self.model.function(
            "decaySpectrum",
            restype=ctypes.c_int,
            argtypes=[ctypes.c_char_p, ctypes.c_int, ctypes.POINTER(ctypes.c_double)],
        )
        err = function(self.model._encode_string(particle_name), int(out_particle), spectrum)
        return {"err": err, "spectrum": self.model._numpy_from_array(spectrum)}

    def basic_spectra(self, mass, pdg, out_particle, uncertainty=False):
        spectrum = self.model._new_double_array()
        function_name = "spectraUncertainty" if uncertainty else "basicSpectra"
        function = self.model.function(
            function_name,
            restype=ctypes.c_int,
            argtypes=[ctypes.c_double, ctypes.c_int, ctypes.c_int, ctypes.POINTER(ctypes.c_double)],
        )
        err = function(float(mass), int(pdg), int(out_particle), spectrum)
        return {"err": err, "spectrum": self.model._numpy_from_array(spectrum)}

    def planck_cmb(self, v_sigma, gamma_spectrum, electron_spectrum, old=False):
        function = self.model.function(
            "PlanckCMB_old" if old else "PlanckCMB",
            restype=ctypes.c_double,
            argtypes=[ctypes.c_double, ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double)],
        )
        return function(float(v_sigma), self.model._as_double_array(gamma_spectrum), self.model._as_double_array(electron_spectrum))

    def dwarf_signal(self, v_sigma, proton_spectrum):
        function = self.model.function(
            "DwarfSignal",
            restype=ctypes.c_double,
            argtypes=[ctypes.c_double, ctypes.POINTER(ctypes.c_double)],
        )
        return function(float(v_sigma), self.model._as_double_array(proton_spectrum))

    def x_interp(self, x, spectrum):
        function = self.model.function(
            "xInterp",
            restype=ctypes.c_double,
            argtypes=[ctypes.c_double, ctypes.POINTER(ctypes.c_double)],
        )
        return function(float(x), self.model._as_double_array(spectrum))

    def z_interp(self, z, spectrum):
        function = self.model.function(
            "zInterp",
            restype=ctypes.c_double,
            argtypes=[ctypes.c_double, ctypes.POINTER(ctypes.c_double)],
        )
        return function(float(z), self.model._as_double_array(spectrum))

    def spectr_info(self, e_min, spectrum):
        e_total = ctypes.c_double()
        function = self.model.function(
            "spectrInfo",
            restype=ctypes.c_double,
            argtypes=[ctypes.c_double, ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double)],
        )
        n_total = function(float(e_min), self.model._as_double_array(spectrum), ctypes.byref(e_total))
        return {"Ntotal": n_total, "Etotal": e_total.value}

    def spectr_int(self, e_min, e_max, spectrum):
        function = self.model.function(
            "spectrInt",
            restype=ctypes.c_double,
            argtypes=[ctypes.c_double, ctypes.c_double, ctypes.POINTER(ctypes.c_double)],
        )
        return function(float(e_min), float(e_max), self.model._as_double_array(spectrum))

    def gamma_flux_tab(self, fi, dfi, sigma_v, spectrum):
        observed = self.model._new_double_array()
        function = self.model.function(
            "gammaFluxTab",
            restype=None,
            argtypes=[ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double)],
        )
        function(float(fi), float(dfi), float(sigma_v), self.model._as_double_array(spectrum), observed)
        return self.model._numpy_from_array(observed)

    def gamma_flux_tab_gc(self, longitude, latitude, d_longitude, d_latitude, sigma_v, spectrum):
        observed = self.model._new_double_array()
        function = self.model.function(
            "gammaFluxTabGC",
            restype=None,
            argtypes=[ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double)],
        )
        function(
            float(longitude),
            float(latitude),
            float(d_longitude),
            float(d_latitude),
            float(sigma_v),
            self.model._as_double_array(spectrum),
            observed,
        )
        return self.model._numpy_from_array(observed)

    def solar_modulation(self, phi, mass, spectrum):
        output = self.model._new_double_array()
        function = self.model.function(
            "solarModulation",
            restype=None,
            argtypes=[ctypes.c_double, ctypes.c_double, ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double)],
        )
        function(float(phi), float(mass), self.model._as_double_array(spectrum), output)
        return self.model._numpy_from_array(output)

    def pbar_background_tab(self, emax):
        spectrum = self.model._new_double_array()
        function = self.model.function(
            "pBarBackgroundTab",
            restype=None,
            argtypes=[ctypes.c_double, ctypes.POINTER(ctypes.c_double)],
        )
        function(float(emax), spectrum)
        return self.model._numpy_from_array(spectrum)

    def positron_flux(self, energy, sigma_v, spectrum):
        function = self.model.function(
            "posiFlux",
            restype=ctypes.c_double,
            argtypes=[ctypes.c_double, ctypes.c_double, ctypes.POINTER(ctypes.c_double)],
        )
        return function(float(energy), float(sigma_v), self.model._as_double_array(spectrum))

    def positron_flux_tab(self, e_min, sigma_v, spectrum):
        output = self.model._new_double_array()
        function = self.model.function(
            "posiFluxTab",
            restype=None,
            argtypes=[ctypes.c_double, ctypes.c_double, ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double)],
        )
        function(float(e_min), float(sigma_v), self.model._as_double_array(spectrum), output)
        return self.model._numpy_from_array(output)

    def pbar_flux_tab(self, e_min, sigma_v, spectrum):
        output = self.model._new_double_array()
        function = self.model.function(
            "pbarFluxTab",
            restype=None,
            argtypes=[ctypes.c_double, ctypes.c_double, ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double)],
        )
        function(float(e_min), float(sigma_v), self.model._as_double_array(spectrum), output)
        return self.model._numpy_from_array(output)

    def basic_nu_spectra(self, for_sun, mass, pdg, polarization=0):
        neutrino = self.model._new_double_array()
        antineutrino = self.model._new_double_array()
        function = self.model.function(
            "basicNuSpectra",
            restype=ctypes.c_int,
            argtypes=[ctypes.c_int, ctypes.c_double, ctypes.c_int, ctypes.c_int, ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double)],
        )
        err = function(int(for_sun), float(mass), int(pdg), int(polarization), neutrino, antineutrino)
        return {"err": err, "nu": self.model._numpy_from_array(neutrino), "nu_bar": self.model._numpy_from_array(antineutrino)}

    def muon_contained(self, neutrino_spectrum, antineutrino_spectrum, rho):
        muon = self.model._new_double_array()
        function = self.model.function(
            "muonContained",
            restype=None,
            argtypes=[ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.c_double, ctypes.POINTER(ctypes.c_double)],
        )
        function(self.model._as_double_array(neutrino_spectrum), self.model._as_double_array(antineutrino_spectrum), float(rho), muon)
        return self.model._numpy_from_array(muon)

    def muon_upward(self, neutrino_spectrum, antineutrino_spectrum):
        muon = self.model._new_double_array()
        function = self.model.function(
            "muonUpward",
            restype=None,
            argtypes=[ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double)],
        )
        function(self.model._as_double_array(neutrino_spectrum), self.model._as_double_array(antineutrino_spectrum), muon)
        return self.model._numpy_from_array(muon)