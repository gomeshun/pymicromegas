# pymicromegas
Python interface of micromegas (unofficial)

The Python package uses a `src/` layout.  The vendored micrOMEGAs source tree is
packaged under `src/pymicromegas/micromegas_5.0.8/` and is built by the Python
build backend when the package is installed with tools such as `pip` or `uv`.
Importing `pymicromegas` itself does not run `make`.

## Features
pymicromegas:

- directly receives your model parameters as Python `dict` (does not generate internal `.par` files)
- can switch on/off by `flags` parameter

If you don't like Python, instead, you can use `main.c/cpp` in
`src/pymicromegas/` just as normal `main.c/cpp` of micromegas.
These modified `main` files receive arguments like:
```
./pymicromegas_main <integer to define flags> <n: number of parameters> <DOF file or None> <parameter name 1> ... <parameter name n> <parameter value 1> ... <parameter value n>
```


# How to use it

Install pymicromegas with `pip install .` or `uv add <path-or-url>`. The
micrOMEGAs source tree is included as package data, so it can also be inspected
in `src/pymicromegas/micromegas_5.0.8/` from a source checkout. Then

```python
from pymicromegas import PyMicrOmegas

interf = PyMicrOmegas()
project = interf.create_newproject("test")
project.load_mdl_files(["the", "list of", "your", ".mdl file", "paths"])
project.compile()  # after compiling, you can reload it as Project(project_name).


args = {
  "parname1" : 1.0  # parameter values
  "parname2" : 10   
  "parname3" : 1e-8
}

##########
# flags: defined in pymicromegas.FLAGS. 
# 
# List of available flags:
# ['MASSES_INFO','CONSTRAINTS','MONOJET','HIGGSBOUNDS', 'HIGGSSIGNALS', 'LILITH', 'SMODELS', 'OMEGA', 'FREEZEIN', 'INDIRECT_DETECTION', 'RESET_FORMFACTORS', 'CDM_NUCLEON', 'CDM_NUCLEUS', 'NEUTRINO', 'DECAYS', 'CROSS_SECTIONS', 'SHOWPLOTS', 'CLEAN']
##########
flags = ["MASSES_INFO","OMEGA"]

#process = project.run(args,flags)  # return subprocess.CompletedProcess
#print(process.stdout)  # print the output text of micromegas

output_dict = project(args,flags)  # directly return parsed output (at present, relic density only)
output_dict = project.calc_omega(args)  # directly return parsed output about relic density with channels

print(output_dict)  
```

# Class

## `MicrOmegas`
- Integrated project and calculation wrapper that builds a generated shared
  library in the micrOMEGAs project directory and loads it with `ctypes`.
- The generated bridge is small and is linked with the same project `Makefile`
  command, so future
  micrOMEGAs makefile changes are reused as much as possible.
- `MicrOmegas.lib` exposes the raw `ctypes.CDLL` object.  Use
  `MicrOmegas.function(name, restype, argtypes)` or
  `MicrOmegas.call(name, *args, restype=..., argtypes=...)` to call linked
  micrOMEGAs functions directly when no convenience method exists.

```python
from pymicromegas import MicrOmegas

mo = MicrOmegas("test", mdl_paths=["your_model_1.mdl", "your_model_2.mdl"])

parameters = {
  "parname1": 1.0,
  "parname2": 10.0,
}

omega = mo.dark_omega(parameters)
print(omega["Omega"])

# Direct ctypes access to linked micrOMEGAs/project symbols is also available.
mass = mo.find_value("Mcdm")
raw_lib = mo.lib
```

## `PyMicrOmegas`
- wrapper class of doing `newProject`, `make`, `make clean` in the micromegas directory.
- micrOMEGAs is built during package installation.  In editable/source-tree
  development, project compilation also verifies the build artifacts and runs
  `make` if they are missing.

If you want to modify micromegas, 
1. clean
1. modify 
1. make again
    
## `Project`
  - wrapper class of `make`, `./pymicromegas_main ...`, in project directories.
  - default micrOMEGAs model directories and user-created projects are both
    callable through the same API.
  - `Project.__call__` to directly return parsed micromegas outputs (callable object, used as if it is like a function. See the previous example.)

# Tests

The integration tests are written with the standard-library `unittest` module.
They compile and run both a bundled default model and a newly-created project
loaded from `.mdl` files:

```bash
uv run python -m unittest discover -s tests
```


# TODO
- etc...
