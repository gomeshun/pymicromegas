import math
import subprocess
import sys
import unittest
import uuid
from pathlib import Path

from pymicromegas import MicrOmegas, Project, PyMicrOmegas


SINGLETDM_VALID_PARAMETERS = {
    "Q": 100.0,
    "Mh": 125.0,
    "laS": 0.2,
    "laSH": 0.1,
    "Mdm1": 50.0,
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
            project = self.interface.create_newproject(project_name, return_project=True)
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


if __name__ == "__main__":
    unittest.main()