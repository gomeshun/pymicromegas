# pymicromegas

Python interface for micrOMEGAs (unofficial).

The primary API is the `MicrOmegas` class. It creates or loads a micrOMEGAs
project, builds a small shared-library bridge, and calls micrOMEGAs through
`ctypes` without going through a generated command-line executable for each
calculation.

The lower-level `Project` and `PyMicrOmegas` classes are still available for
compatibility with earlier pymicromegas workflows and for users who want to work
directly with micrOMEGAs project directories.

## Source Version

This package currently vendors micrOMEGAs 7.1 from Zenodo record `20267206`:

- source tree: `src/pymicromegas/micromegas_7.1/`
- archive checksum: `md5:aa066ac8d712a9c5eca4134a168f1f15`
- local integration patch: `sources/omega.c` uses C99 `isfinite()` instead of
  the obsolete `finite()` call for macOS/clang compatibility

The vendored source tree is included as package data and is built by the Python
build backend when the package is installed. Importing `pymicromegas` itself does
not run `make`.

## License

This repository is distributed under the GNU General Public License v3.0; see
the `LICENSE` file for the full text.

The package vendors micrOMEGAs and associated third-party components under
`src/pymicromegas/micromegas_7.1/`. Those components remain subject to their
respective upstream license notices.

## Install

From a source checkout:

```bash
uv add <path-or-url>
```

or:

```bash
pip install .
```

## Quick Start

Use `MicrOmegas` for normal Python calculations:

```python
from pymicromegas import MicrOmegas

parameters = {
    "Q": 100.0,
    "Mh": 125.0,
    "laS": 0.2,
    "laSH": 0.1,
    "Mdm1": 50.0,
}

model = MicrOmegas("SingletDM")
omega = model.dark_omega(parameters)

print(omega["Omega"])
print(omega["Xf"])
```

Freeze-in calculations are also available through `MicrOmegas`:

```python
from pymicromegas import MicrOmegas

parameters = {
  "Q": 100.0,
  "Mh": 120.0,
  "laS": 0.0,
  "laSH": 1e-11,
  "Mdm1": 70.0,
}

model = MicrOmegas("SingletDM")
omega_fi = model.dark_omega_freeze_in(parameters, channels=True)

print(omega_fi["Omega"])
print(omega_fi["channels"][:3])
```

For a custom model, pass the CalcHEP `.mdl` files when constructing the class:

```python
from pymicromegas import MicrOmegas

model = MicrOmegas(
    "my_model",
    mdl_paths=[
        "path/to/vars1.mdl",
        "path/to/func1.mdl",
        "path/to/prtcls1.mdl",
        "path/to/lgrng1.mdl",
    ],
)

result = model.dark_omega({"parname1": 1.0, "parname2": 10.0})
```

`MicrOmegas.lib` exposes the raw `ctypes.CDLL` object. Use
`MicrOmegas.function(name, restype, argtypes)` or
`MicrOmegas.call(name, *args, restype=..., argtypes=...)` when no convenience
method exists.

## Main Classes

### `MicrOmegas`

- Main Python-facing class for calculations.
- Builds a generated shared library in the micrOMEGAs project directory.
- Groups major physics calculations under process delegates:
  `model.relic_density`, `model.freeze_in`, `model.indirect_detection`, and
  `model.direct_detection`. The legacy direct method names on `MicrOmegas`
  remain available and delegate to these process objects.
- Provides convenience methods such as `assign`, `find_value`, `dark_omega`,
  `dark_omega2`, `dark_omega_freeze_in`, `dark_omega_freeze_in_22`,
  `dark_omega_freeze_in_decay`, `dark_omega_tr`, `dark_omega_fo`,
  `dark_omega2_tr`, `dark_omega_n`, and `dark_omega_infl`.
- Adds Pythonic wrappers for simple public `micromegas.h` APIs, including
  particle lookups, LEP/Z constraints, thermodynamic functions, equilibrium
  abundances, cross-section helpers, FreezeIn helper yields, halo/profile
  helpers, spectrum-table utilities, neutrino tables, and scalar/nucleon form
  factor helpers. Check `model.available_header_methods` for the functions
  exported by the current project library.
- Converts common C output patterns into Python values: text printers return
  strings, scalar output pointers return dictionaries, and `NZ=250` spectra are
  returned as `numpy.ndarray` objects.
- Exposes v7 N-component metadata through bridge helpers such as
  `pymicromegas_cdm_name`, `pymicromegas_cdm_mass`, and
  `pymicromegas_cdm_fraction` on `model.lib`.

### `Project`

- Compatibility wrapper around `make` and `./pymicromegas_main ...` in a
  micrOMEGAs project directory.
- Useful when you want the raw stdout produced by a micrOMEGAs-style executable.
- Default bundled model directories and user-created projects share the same API.

```python
from pymicromegas import Project

project = Project("SingletDM")
project.compile()
result = project.calc_omega(parameters)
```

The generated executable receives arguments as:

```text
./pymicromegas_main <flags> <n parameters> <DOF file or None> <parameter name 1> ... <parameter name n> <parameter value 1> ... <parameter value n>
```

### `PyMicrOmegas`

- Lower-level helper for `newProject`, `make`, and `make clean` in the vendored
  micrOMEGAs directory.
- Ensures the micrOMEGAs source tree is built when project compilation needs it.

## Runtime Flags

`Project.run` accepts flag names from `pymicromegas.FLAGS`, including:

```python
[
    "MASSES_INFO",
    "CONSTRAINTS",
    "MONOJET",
    "HIGGSBOUNDS",
    "HIGGSSIGNALS",
    "LILITH",
    "SMODELS",
    "OMEGA",
    "FREEZEIN",
    "INDIRECT_DETECTION",
    "RESET_FORMFACTORS",
    "CDM_NUCLEON",
    "CDM_NUCLEUS",
    "NEUTRINO",
    "DECAYS",
    "CROSS_SECTIONS",
    "SHOWPLOTS",
    "CLEAN",
]
```

The `MicrOmegas` class is preferred for new code. Use `Project` flags when you
specifically need the command-line runner behavior.

## Tests

The integration tests use the standard-library `unittest` module. They compile
and run both the bundled `SingletDM` model and a newly-created project loaded
from `.mdl` files:

```bash
uv run python -m unittest discover -s tests
```
