import os
import subprocess
from pathlib import Path

from setuptools import setup
from setuptools.command.build_py import build_py as _build_py
from wheel.bdist_wheel import bdist_wheel as _bdist_wheel


class build_py(_build_py):
    def run(self):
        super().run()
        micromegas_dir = Path(self.build_lib) / "pymicromegas" / "micromegas_5.0.8"
        if not micromegas_dir.is_dir():
            raise RuntimeError(f"micrOMEGAs source tree was not copied: {micromegas_dir}")

        env = os.environ.copy()
        env["PYMICROMEGAS_MICROPATH"] = str(micromegas_dir)
        subprocess.run(["make"], cwd=micromegas_dir, env=env, check=True)

        build_marker = micromegas_dir / "include" / "microPath.h"
        flags_file = micromegas_dir / "CalcHEP_src" / "FlagsForMake"
        library_file = micromegas_dir / "lib" / "micromegas.a"
        missing = [path for path in (build_marker, flags_file, library_file) if not path.is_file()]
        if missing:
            missing_paths = ", ".join(str(path) for path in missing)
            raise RuntimeError(f"micrOMEGAs build did not produce expected files: {missing_paths}")


class bdist_wheel(_bdist_wheel):
    def finalize_options(self):
        super().finalize_options()
        self.root_is_pure = False


setup(cmdclass={"build_py": build_py, "bdist_wheel": bdist_wheel})