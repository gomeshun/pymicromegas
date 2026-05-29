import math
import subprocess
import sys
import unittest
import uuid
from pathlib import Path

from pymicromegas import MICROMEGAS_DIR, MICROMEGAS_VERSION, MicrOmegas, Project, PyMicrOmegas


SINGLETDM_VALID_PARAMETERS = {
    "Q": 100.0,
    "Mh": 125.0,
    "laS": 0.2,
    "laSH": 0.1,
    "Mdm1": 50.0,
}

SINGLETDM_FREEZEIN_PARAMETERS = {
    "Q": 100.0,
    "Mh": 120.0,
    "laS": 0.0,
    "laSH": 1e-11,
    "Mdm1": 70.0,
}


def assert_valid_omega(testcase, result):
    testcase.assertIn("Omega", result)
    testcase.assertTrue(math.isfinite(result["Omega"]))
    testcase.assertGreaterEqual(result["Omega"], 0.0)


class ImportBehaviorTest(unittest.TestCase):
    def test_import_has_no_build_side_effect_output(self):
        process = subprocess.run(
            [sys.executable, "-c", "import pymicromegas; print('imported')"],
            check=True,
            capture_output=True,
            text=True,
        )
        self.assertEqual(process.stdout.strip(), "imported")


class VendoredSourceTest(unittest.TestCase):
    def test_micromegas_source_version_is_7_1(self):
        self.assertEqual(MICROMEGAS_VERSION, "7.1")
        self.assertEqual(Path(MICROMEGAS_DIR).name, "micromegas_7.1")
        self.assertTrue((Path(MICROMEGAS_DIR) / "man" / "manual_7.1.tex").is_file())


class MicrOmegasIntegrationTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.interface = PyMicrOmegas()
        cls.interface.ensure_micromegas_built()
        cls.micromegas_path = Path(cls.interface.path)

    def test_default_model_project_calc_omega(self):
        project = Project("SingletDM")
        project.compile()
        result = project.calc_omega(SINGLETDM_VALID_PARAMETERS)
        assert_valid_omega(self, result)

    def test_custom_project_from_mdl_files_calc_omega(self):
        project_name = f"pymg_unittest_{uuid.uuid4().hex[:8]}"
        model_dir = self.micromegas_path / "SingletDM" / "work" / "models"
        mdl_paths = sorted(model_dir.glob("*.mdl"))
        self.assertTrue(mdl_paths)

        try:
            self.interface.create_newproject(project_name)
            project = Project(project_name)
            project.load_mdl_files(mdl_paths)
            project.compile()
            result = project.calc_omega(SINGLETDM_VALID_PARAMETERS)
            assert_valid_omega(self, result)
        finally:
            if self.interface.project_exists(project_name):
                self.interface.remove_project(project_name)

    def test_default_model_ctypes_dark_omega(self):
        model = MicrOmegas("SingletDM")
        result = model.dark_omega(SINGLETDM_VALID_PARAMETERS)
        assert_valid_omega(self, result)
        self.assertEqual(result["err"], 0)
        self.assertIsNotNone(model.lib.pymicromegas_cdm1())
        self.assertGreater(model.lib.pymicromegas_mcdm1(), 0.0)
        self.assertAlmostEqual(model.lib.pymicromegas_cdm_fraction(1), 1.0)

    def test_default_model_ctypes_dark_omega_freeze_in(self):
        model = MicrOmegas("SingletDM")
        result = model.dark_omega_freeze_in(SINGLETDM_FREEZEIN_PARAMETERS, channels=True)
        assert_valid_omega(self, result)
        self.assertEqual(result["err"], 0)
        self.assertEqual(result["particle"], "~x1")
        self.assertTrue(model.is_feeble("~x1"))
        self.assertTrue(result["channels"])

    def test_default_model_ctypes_header_wrappers(self):
        model = MicrOmegas("SingletDM")
        model.assign(SINGLETDM_VALID_PARAMETERS)
        model.sort_odd_particles()

        self.assertEqual(str(model.pdg_name(25)), "h")
        self.assertGreater(int(model.particle_number("~x1")), 0)
        self.assertAlmostEqual(float(model.particle_mass("~x1")), SINGLETDM_VALID_PARAMETERS["Mdm1"])
        self.assertGreater(float(model.h_eff(1.0)), 0.0)
        self.assertGreater(float(model.g_eff(1.0)), 0.0)
        self.assertGreater(float(model.hubble(1.0)), 0.0)

        next_odd = model.next_odd(0)
        self.assertEqual(next_odd["name"], "~x1")
        self.assertAlmostEqual(next_odd["mass"], SINGLETDM_VALID_PARAMETERS["Mdm1"])

        lep = model.lsp_nlsp_lep()
        self.assertIn("excluded", lep)
        self.assertIn("cross_section_limit", lep)
        self.assertIn("~x1", model.print_masses())

        freeze_out = model.dark_omega_fo(SINGLETDM_VALID_PARAMETERS)
        assert_valid_omega(self, freeze_out)
        self.assertEqual(freeze_out["err"], 0)


if __name__ == "__main__":
    unittest.main()