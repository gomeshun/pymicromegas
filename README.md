# pymicromegas
Python interface of micromegas.

## Features
pymicromegas:

- directly reseives your model parameters as Python `dict` (does not generate internal `.par` files)
- can switch on/off by `flags` parameter

If you don't like Python, instead, you can use `main.c/cpp` in `/pymicromegas/` just as normal `main.c/cpp` of micromegas.
These modified `main` files receive arguments like:
```
./main <integer to define flags> <n: number of parameters> <parmeter name 1> ... <parmeter name n> <parmeter value 1> ... <parmeter value n>
```


# How to use it

`git clone` to download pymicromegas. The micrOMEGAs source tree is included
as `micromegas_5.0.8/`, so it can be inspected and built directly. Then

```python
from pymicromegas import PyMicrOmegas

interf = PyMicrOmegas()
project = interf.create_newproject("test")
project.load_mdl_files(["the", "list of", "your", ".mdl file", "paths"])
project.compile()  # once a project is compiled, you can directly call the compiled project as Project(project_name).


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

outout_dict = project(args,flags)  # directly return parsed output (at present, relic density only)
output_dict = project.calc_omega(args)  # directly return parsed output about relic density with channels

print(output_dict)  
```

# Class

## `MicrOmegas`
- Integrated project and calculation wrapper that builds a generated shared
  library in the micrOMEGAs project directory and loads it with `ctypes`.
- It does not patch the upstream micrOMEGAs C sources.  The generated bridge is
  small and is linked with the same project `Makefile` command, so future
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
- When pymicromegas is imported for the first time, it installs (make) the
  checked-in `micromegas_5.0.8/` source tree if it has not been built yet.

If you want to modify micromegas, 
1. clean
1. modify 
1. make again
    
## `Project`
  - wrapper class of `make`, `./main ...`, in project directories.
  - `Project.__call__` to directly return parsed micromegas outputs (callable object, used as if it is like a function. See the previous example.)


# TODO
- etc...
