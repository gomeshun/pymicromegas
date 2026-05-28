# pymicromegas
Python interface of micromegas.

## Features
pymicromegas:

- directly receives your model parameters as Python `dict` (does not generate internal `.par` files)
- can switch on/off by `flags` parameter

If you don't like Python, instead, you can use `main.c/cpp` in `/pymicromegas/` just as normal `main.c/cpp` of micromegas.
These modified `main` files receive arguments like:
```
./main <integer to define flags> <n: number of parameters> <parmeter name 1> ... <parmeter name n> <parmeter value 1> ... <parmeter value n>
```


# How to use it

`git clone` to download pymicromegas. Then 

```python
from pymicromegas import MicrOmegas

micromegas = MicrOmegas("test", create=True)
micromegas.load_mdl_files(["the", "list of", "your", ".mdl file", "paths"])
micromegas.compile()  # once compiled, load it again with MicrOmegas(project_name) using the default create=False.


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

#process = micromegas.run(args, flags)  # return subprocess.CompletedProcess
#print(process.stdout)  # print the output text of micromegas

parsed_dict = micromegas(args, flags)  # directly return parsed output (at present, relic density only)
omega_dict = micromegas.calc_omega(args)  # directly return parsed output about relic density with channels

print(omega_dict)  
```

The older `PyMicrOmegas` and `Project` classes are still available for
backwards compatibility.

## ctypes interface

`MicrOmegas` can also call selected micrOMEGAs functions through a small
project-local shared library loaded with `ctypes`.  The generated C file lives
inside the micrOMEGAs project directory and links against the normal micrOMEGAs
libraries, so the bundled upstream C sources do not need to be patched.

```python
micromegas.assign_values(args)
print(micromegas.sort_odd_particles())
print(micromegas.dark_omega(args))
print(micromegas.find_val("parname1"))

# Advanced users can call exported micrOMEGAs symbols directly.
assign_val = micromegas.function("assignVal")
```

# Class

## `MicrOmegas`
- standalone wrapper class for micrOMEGAs management and project operations.
- use `MicrOmegas(project_name)` to load an existing project, or
  `MicrOmegas(project_name, create=True)` to create and load a new project.

## `PyMicrOmegas`
- wrapper class of doing `newProject`, `make`, `make clean` in the micromegas directory.
- When pymicromegas imported for the first time, it unzip `miccromegas_5.0.8.tgz` and install (make) it

If you want to modify micromegas, 
1. clean
1. modify 
1. make again
    
## `Project`
  - wrapper class of `make`, `./main ...`, in project directories.
  - `Project.__call__` to directly return parsed micromegas outputs (callable object, used as if it is like a function. See the previous example.)


# TODO
- etc...
