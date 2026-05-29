import os
import ctypes
import shutil
import subprocess
import sys
from pathlib import Path
import numpy as np
from pandas import Series


PACKAGE_DIR = Path(__file__).resolve().parent
MICROMEGAS_DIR = PACKAGE_DIR / "micromegas_5.0.8"
PYTHON_MAIN_C = PACKAGE_DIR / "main.c"
PYTHON_MAIN_CPP = PACKAGE_DIR / "main.cpp"
PYTHON_MAIN_SOURCE = "pymicromegas_main.c"
PYTHON_MAIN_CPP_SOURCE = "pymicromegas_main.cpp"
PYTHON_MAIN_EXECUTABLE = "pymicromegas_main"
USER_PROJECT_MARKER = ".pymicromegas-project"

dir_micromegas = str(MICROMEGAS_DIR) + os.sep



FLAGS = {
    "MASSES_INFO"        : (1 << 0),
    "CONSTRAINTS"        : (1 << 1),
    "MONOJET"            : (1 << 2),
    "HIGGSBOUNDS"        : (1 << 3),
    "HIGGSSIGNALS"       : (1 << 4),
    "LILITH"             : (1 << 5),
    "SMODELS"            : (1 << 6),
    "OMEGA"              : (1 << 7),
    "FREEZEIN"           : (1 << 8),
    "INDIRECT_DETECTION" : (1 << 9),
    "RESET_FORMFACTORS"  : (1 << 10),
    "CDM_NUCLEON"        : (1 << 11),
    "CDM_NUCLEUS"        : (1 << 12),
    "NEUTRINO"           : (1 << 13),
    "DECAYS"             : (1 << 14),
    "CROSS_SECTIONS"     : (1 << 15),
    "SHOWPLOTS"          : (1 << 16),
    "CLEAN"              : (1 << 17)
}


#def flag_to_int(flag_name,bool_flag):
#    if type(bool_flag) is not bool: raise TypeError("invalid input for {}".format(key))
#    if bool_flag:
#        return FLAGS[flag_name]
#    else:
#        return 0

    
#def flags_to_int(dict_flags):
#    int_flags = [flag_to_int(flag_name,bool_flag) for flag_name,bool_flag in dict_flags.items() ] 
#    return sum(int_flags)


######## utils ########
def to_abspath(path):
    return str(Path(path).resolve())


def micromegas_env():
    env = os.environ.copy()
    env["PYMICROMEGAS_MICROPATH"] = str(MICROMEGAS_DIR)
    return env
    
    
def run_bash(command,shell=True,stdout=subprocess.PIPE,encoding="UTF-8",check=True,input=None,cwd=None,verbose=True,env=None,stderr=subprocess.PIPE):
    if verbose: print(command)
    process = subprocess.run(
        command,
        shell=shell,
        stdout=stdout,
        stderr=stderr,
        encoding=encoding,
        check=check,
        input=input,
        cwd=cwd,
        env=micromegas_env() if env is None else env,
    )
    return process

def is_valid_project_name(project_name):
    """
    Check that given project_name is valid.
    Valid project names consist of [a-z]/[A-Z]/[",", ".", "-", "_"].
    """
    return project_name.replace(",","").replace(".","").replace("-","").replace("_","").isalnum()


def get_keys(dict_like_object):
    return dict_like_object.keys()

def get_values(dict_like_object):
    if isinstance(dict_like_object,Series): 
        return dict_like_object.values
    else:
        return dict_like_object.values()


def ensure_project_workdirs(project_path):
    work_path = Path(project_path) / "work"
    for dirname in ("tmp", "results", "so_generated", "lanhep"):
        (work_path / dirname).mkdir(parents=True, exist_ok=True)

######## class definitions ########
class PyMicrOmegas:    
    
    def __init__(self,verbose=False):
        self.path = dir_micromegas
        self.verbose = verbose
        os.environ["PYMICROMEGAS_MICROPATH"] = str(MICROMEGAS_DIR)
        
        if not os.path.isdir(self.path):
            raise RuntimeError(f"micromegas directory is missing: {self.path}")
    
    
    def run_bash(self,command,shell=True,stdout=subprocess.PIPE,encoding="UTF-8",check=False,input=None,verbose=True):
        return run_bash(command,shell=shell,stdout=stdout,encoding=encoding,check=check,input=input,cwd=self.path,verbose=verbose)
    
    def is_micromegas_built(self):
        path = Path(self.path)
        return (
            (path / "include" / "microPath.h").is_file()
            and (path / "CalcHEP_src" / "FlagsForMake").is_file()
            and (path / "lib" / "micromegas.a").is_file()
        )
    
    def ensure_micromegas_built(self):
        if not self.is_micromegas_built():
            return self.compile_micromegas()
        return None
    
    
    def compile_micromegas(self):
        print("Compiling micromegas...")
        return self.run_bash("make",check=True)

    
    def clean_micromegas(self):
        print("Cleaning micromegas...")
        return self.run_bash("make clean",input="y\n")
    
    
    def project_exists(self,project_name):
        return os.path.isdir(self.path + project_name)
    
    def install_python_main(self,project_name):
        if not self.project_exists(project_name): 
            raise RuntimeError("project '{}' does not exist yet.".format(project_name))
        project_path = Path(self.path) / project_name
        ensure_project_workdirs(project_path)
        shutil.copy2(PYTHON_MAIN_C, project_path / PYTHON_MAIN_SOURCE)
        shutil.copy2(PYTHON_MAIN_CPP, project_path / PYTHON_MAIN_CPP_SOURCE)
        return project_path / PYTHON_MAIN_SOURCE
       
        
    def load_modified_main(self,project_name,return_project=False):
        print("Loading modified main files...")
        self.install_python_main(project_name)
        
        if return_project:
            return Project(project_name)
        else: 
            return None
    
    def create_newproject(self,project_name,return_project=False):
        print("Creating new project...")
        if not is_valid_project_name(project_name): raise RuntimeError(f"'{project_name}' is not valid project name.")
        if self.project_exists(project_name): raise RuntimeError("project '{}' exists already.".format(project_name))
            
        commands = [ "./newProject {0}" ]
        process = self.run_bash("\n".join(commands).format(project_name),check=True)
        project_path = Path(self.path) / project_name
        (project_path / USER_PROJECT_MARKER).write_text("created by pymicromegas\n")
        self.install_python_main(project_name)
         
        if return_project:
            return Project(project_name)
        else: 
            return process
        
    
    def load_project(self,project_name):
        if self.project_exists(project_name):
            return Project(project_name)
        else:
            print(f"Project {project_name} is not created yet. Start creating...")
            return self.create_newproject(project_name,return_project=True)
        
    
    def remove_project(self,project_name):
        print("Removing project...")
        if not is_valid_project_name(project_name): raise RuntimeError(f"{project_name} is not valid project name.")
        if not self.project_exists(project_name): raise RuntimeError(f"project '{project_name}' does not exists.")
        
        project = Project(project_name)
        if not project.is_user_defined_project: raise RuntimeError(f"{project_name} is not created by users. remove_project cannot remove default projects.")
        shutil.rmtree(Path(self.path) / project_name)
        return None
    
#    def load_mdl(self,project_name,mdl_paths):
#        if self.project_exists(project_name): raise RuntimeError("project '{}' exists already.".format(project_name))
#        for path in mdl_paths:
#            if path[-4:] != ".mdl": raise RuntimeError("{} is not .mdl file".format(path))
#            commands = [ "cd {0}","cp main=main.cpp"]
    
    
    
class Project:
    
    
    def __init__(self,project_name):
        self.interface = PyMicrOmegas()
        if not self.interface.project_exists(project_name): raise RuntimeError("Project {} does not exist yet. Create it by PyMicrOmegas.create_newproject.".format(project_name))
        self.project_name = project_name
        self.path = str(Path(dir_micromegas) / project_name) + os.sep
        self.models_path = str(Path(self.path) / "work" / "models") + os.sep
        self.main_source = PYTHON_MAIN_SOURCE
        self.main_executable = PYTHON_MAIN_EXECUTABLE
    
    
    def run_bash(self,command,shell=True,stdout=subprocess.PIPE,encoding="UTF-8",check=False,input=None,verbose=True):
        '''
        run './main' and return the output string.
        '''
        return run_bash(command,shell=shell,stdout=stdout,encoding=encoding,check=check,input=input,cwd=self.path,verbose=verbose)
    
    
    def load_mdl_files(self,mdl_paths):
        print("Loading .mdl files...")
        if type(mdl_paths) not in [list, tuple]: raise RuntimeError("input arguments must be list or tuple.")
        if len(mdl_paths)==0: raise RuntimeError("input argument is empty list or tuple.")
        for mdl_path in mdl_paths:
            if not os.path.isfile(mdl_path): raise RuntimeError("{} does not existing file".format(mdl_path))
            abs_mdl_path = to_abspath(mdl_path)
            print(f"Loading {abs_mdl_path}...")
            shutil.copy2(abs_mdl_path, self.models_path)
        return None
    
    def install_python_main(self):
        return self.interface.install_python_main(self.project_name)
    
        
    def compile(self,main=None):
        #process = subprocess.run("bash install_project.sh " + project_name,shell=True,stdout=subprocess.PIPE)
        self.interface.ensure_micromegas_built()
        ensure_project_workdirs(self.path)
        if main is None:
            self.install_python_main()
            main = self.main_source
            self.main_executable = PYTHON_MAIN_EXECUTABLE
        else:
            self.main_executable = Path(main).stem
        return self.run_bash("make main={}".format(main),check=True)
           
      
    def clean(self):
        return self.run_bash("make clean")
    
    
    @property
    def vars(self):
        vars = np.loadtxt(self.models_path+"vars1.mdl",skiprows=3,delimiter="|",dtype=str,comments="=")  # .mdl file sometimes contain "=======..." line
        return vars
    
    @property
    def is_user_defined_project(self):
        return os.path.exists(self.path+USER_PROJECT_MARKER) or os.path.exists(self.path+"main_original.c")
    
    def ensure_executable(self):
        executable = Path(self.path) / self.main_executable
        if not executable.exists():
            self.compile()
        return executable
    
    
    def run(self,dict_parameters,flags=None,dof_fname=None):
        """
        kwargs: see 'PyMicrOmegas.FlAGS.keys()'. 
        
        dict_parameters: dict or pandas.Series
        """
        
        #int_flags = flags_to_int(kwargs)
        int_flags = 0 if (flags is None) else sum( FLAGS[key] for key in flags )
        n_inputvals = len(dict_parameters)
        if dof_fname is None:
            dof_fname = "None"
        else:
            if not os.path.isfile(dof_fname): raise RuntimeError("{} does not existing file".format(dof_fname))
            dof_fname = to_abspath(dof_fname)
        
        par_names = " ".join(map(str,get_keys(dict_parameters)))
        par_vals  = " ".join(map(str,get_values(dict_parameters)))
        args = f"{int_flags} {n_inputvals} {dof_fname} {par_names} {par_vals}"
        
        executable = self.ensure_executable()
        return self.run_bash("./{} {}".format(executable.name, args),check=True,verbose=False)
    
    
    def parse_omega(self,micromegas_output,flags=None,with_channels=False):
        def parse_channel(channel_str):
            """
            parse channel string.
            If channel_str is not proper, return None
            """
            terms = channel_str.split()
            if (len(terms) != 5) or (len(terms[3]) <= 2) or (terms[3][:2] != "->"): return None
            
            branching = terms[0]
            ch_in = (terms[1],terms[2])
            ch_out = (terms[3][2:],terms[4])
            return {"Br":branching,"in":ch_in,"out":ch_out}
        
        if flags is None: raise RuntimeError("No flags specified")
        return_dict = {}
        if "OMEGA" in flags:
            floatnize = lambda key_val: (key_val[0],float(key_val[1]))
            lines = micromegas_output.split("\n")
            #ind = lines.index("==== Calculation of relic density =====")
            try:
                ind = lines.index("==== Calculation of relic density =====")
            except ValueError as e:
                raise ValueError(f"Cannot find the line \"==== Calculation of relic density =====\"\noutput:{micromegas_output}")
            key_vals = dict([floatnize(key_val.split("=")) for key_val in lines[ind+1].split()])
            return_dict.update(key_vals)
            
            #### parse channels ####
            if with_channels:
                channels = []
                for line in lines[ind+4:]:
                    ch = parse_channel(line)
                    if ch is not None:
                        channels.append(ch)
                    else: 
                        break
                return_dict["channels"] = channels
            
        return return_dict
        
    
    def __call__(self,dict_parameters,flags=None,dof_fname=None):
        output = self.run(dict_parameters,flags,dof_fname).stdout
        return self.parse_omega(output,flags)
    
    def calc_omega(self,dict_parameters,dof_fname=None):
        flags = ["OMEGA"]
        output = self.run(dict_parameters,flags,dof_fname=dof_fname).stdout
        return self.parse_omega(output,flags,with_channels=True)    


class MicrOmegas:
    """ctypes based micrOMEGAs project wrapper.

    The class creates or loads one micrOMEGAs project, builds a small generated
    shared-library bridge in that project, and exposes both convenience methods
    and the raw ``ctypes.CDLL`` object for direct calls to linked micrOMEGAs
    symbols.
    """

    BRIDGE_SOURCE = r"""
#include <string.h>

#include "../include/micromegas.h"
#include "../include/micromegas_aux.h"
#include "lib/pmodel.h"

int pymicromegas_assign_value(const char *name, double value)
{
    return assignVal((char *)name, value);
}

int pymicromegas_assign_values(int n_values, const char **names, const double *values)
{
    int err = 0;
    for(int i = 0; i < n_values; ++i)
    {
        err = assignVal((char *)names[i], values[i]);
        if(err) return i + 1;
    }
    return 0;
}

double pymicromegas_find_value(const char *name)
{
    return findValW((char *)name);
}

void pymicromegas_set_gauge(int force_ug, int vzdecay, int vwdecay)
{
    ForceUG = force_ug;
    VZdecay = vzdecay;
    VWdecay = vwdecay;
}

int pymicromegas_sort_odd_particles(char *cdm_name, int cdm_name_size)
{
    char local_name[64] = "";
    int err = sortOddParticles(local_name);
    if(!err && cdm_name && cdm_name_size > 0)
    {
        strncpy(cdm_name, local_name, (size_t)cdm_name_size - 1);
        cdm_name[cdm_name_size - 1] = '\0';
    }
    return err;
}

int pymicromegas_load_heff_geff(const char *path)
{
    return loadHeffGeff((char *)path);
}

double pymicromegas_dark_omega(double *xf, int fast, double beps, int *err)
{
    return darkOmega(xf, fast, beps, err);
}

double pymicromegas_dark_omega2(int fast, double beps)
{
    return darkOmega2(fast, beps);
}

const char *pymicromegas_cdm1(void)
{
    return CDM1;
}

const char *pymicromegas_cdm2(void)
{
    return CDM2;
}

double pymicromegas_mcdm(void)
{
    return Mcdm;
}

double pymicromegas_mcdm1(void)
{
    return Mcdm1;
}

double pymicromegas_mcdm2(void)
{
    return Mcdm2;
}
"""

    def __init__(self, project_name, mdl_paths=None, build=True, verbose=False):
        if not is_valid_project_name(project_name):
            raise RuntimeError(f"'{project_name}' is not a valid project name.")

        self.project_name = project_name
        self.verbose = verbose
        self.interface = PyMicrOmegas(verbose=verbose)
        self.interface.ensure_micromegas_built()
        self.micromegas_path = Path(self.interface.path)
        self.path = self.micromegas_path / project_name
        self.models_path = self.path / "work" / "models"
        self.bridge_source = self.path / "pymicromegas_bridge.c"
        self.library_path = self.path / ("libpymicromegas.dylib" if sys.platform == "darwin" else "libpymicromegas.so")
        self._lib = None

        self._ensure_project()
        if mdl_paths is not None:
            self.load_mdl_files(mdl_paths)
        if build:
            self.build_library()
            self.load_library()

    def _run(self, command, check=True, input=None):
        if self.verbose:
            print(command)
        return subprocess.run(
            command,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            encoding="UTF-8",
            check=check,
            input=input,
            cwd=self.path,
            env=micromegas_env(),
        )

    def _ensure_project(self):
        if self.path.is_dir():
            return
        process = subprocess.run(
            f"./newProject {self.project_name}",
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            encoding="UTF-8",
            check=True,
            cwd=self.micromegas_path,
            env=micromegas_env(),
        )
        (self.path / USER_PROJECT_MARKER).write_text("created by pymicromegas\n")
        if self.verbose:
            print(process.stdout)

    def load_mdl_files(self, mdl_paths):
        if not isinstance(mdl_paths, (list, tuple)):
            raise RuntimeError("mdl_paths must be a list or tuple.")
        if len(mdl_paths) == 0:
            raise RuntimeError("mdl_paths cannot be empty.")
        self.models_path.mkdir(parents=True, exist_ok=True)
        for mdl_path in mdl_paths:
            mdl_path = Path(mdl_path)
            if not mdl_path.is_file():
                raise RuntimeError(f"{mdl_path} is not an existing file")
            shutil.copy2(str(mdl_path.resolve()), str(self.models_path))

    @property
    def vars(self):
        return np.loadtxt(
            str(self.models_path / "vars1.mdl"),
            skiprows=3,
            delimiter="|",
            dtype=str,
            comments="=",
        )

    def write_bridge_source(self, force=False):
        if force or not self.bridge_source.exists() or self.bridge_source.read_text() != self.BRIDGE_SOURCE:
            self.bridge_source.write_text(self.BRIDGE_SOURCE)
        return self.bridge_source

    def _shared_link_command(self):
        make_result = self._run(f"make -n main={self.bridge_source.name}", check=True)
        executable_name = self.bridge_source.with_suffix("").name
        source_name = self.bridge_source.name
        shared_flags = "-dynamiclib -fPIC -Wl,-undefined,dynamic_lookup" if sys.platform == "darwin" else "-shared -fPIC -Wl,-export-dynamic"
        for line in make_result.stdout.splitlines():
            if source_name in line and f"-o {executable_name}" in line:
                return line.replace(
                    f"-o {executable_name}",
                    f"{shared_flags} -o {self.library_path.name}",
                    1,
                )
        raise RuntimeError(
            "Could not derive the micrOMEGAs link command. "
            f"make output was:\n{make_result.stdout}"
        )

    def build_library(self, force=False):
        ensure_project_workdirs(self.path)
        self.write_bridge_source(force=force)
        if force or not self.library_path.exists():
            self._run("make libs work/bin", check=True)
            self._run(self._shared_link_command(), check=True)
        return self.library_path

    def load_library(self, force=False):
        if force or self._lib is None:
            if not self.library_path.exists():
                self.build_library()
            self._lib = ctypes.CDLL(str(self.library_path))
            self._configure_bridge_functions()
        return self._lib

    @property
    def lib(self):
        return self.load_library()

    def _configure_bridge_functions(self):
        lib = self._lib
        if lib is None:
            raise RuntimeError("micrOMEGAs bridge library is not loaded.")
        lib.pymicromegas_assign_value.argtypes = [ctypes.c_char_p, ctypes.c_double]
        lib.pymicromegas_assign_value.restype = ctypes.c_int
        lib.pymicromegas_assign_values.argtypes = [
            ctypes.c_int,
            ctypes.POINTER(ctypes.c_char_p),
            ctypes.POINTER(ctypes.c_double),
        ]
        lib.pymicromegas_assign_values.restype = ctypes.c_int
        lib.pymicromegas_find_value.argtypes = [ctypes.c_char_p]
        lib.pymicromegas_find_value.restype = ctypes.c_double
        lib.pymicromegas_set_gauge.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.c_int]
        lib.pymicromegas_set_gauge.restype = None
        lib.pymicromegas_sort_odd_particles.argtypes = [ctypes.c_char_p, ctypes.c_int]
        lib.pymicromegas_sort_odd_particles.restype = ctypes.c_int
        lib.pymicromegas_load_heff_geff.argtypes = [ctypes.c_char_p]
        lib.pymicromegas_load_heff_geff.restype = ctypes.c_int
        lib.pymicromegas_dark_omega.argtypes = [
            ctypes.POINTER(ctypes.c_double),
            ctypes.c_int,
            ctypes.c_double,
            ctypes.POINTER(ctypes.c_int),
        ]
        lib.pymicromegas_dark_omega.restype = ctypes.c_double
        lib.pymicromegas_dark_omega2.argtypes = [ctypes.c_int, ctypes.c_double]
        lib.pymicromegas_dark_omega2.restype = ctypes.c_double
        lib.pymicromegas_cdm1.argtypes = []
        lib.pymicromegas_cdm1.restype = ctypes.c_char_p
        lib.pymicromegas_cdm2.argtypes = []
        lib.pymicromegas_cdm2.restype = ctypes.c_char_p
        lib.pymicromegas_mcdm.argtypes = []
        lib.pymicromegas_mcdm.restype = ctypes.c_double
        lib.pymicromegas_mcdm1.argtypes = []
        lib.pymicromegas_mcdm1.restype = ctypes.c_double
        lib.pymicromegas_mcdm2.argtypes = []
        lib.pymicromegas_mcdm2.restype = ctypes.c_double

    def assign(self, parameters):
        names = list(map(str, get_keys(parameters)))
        values = [float(value) for value in get_values(parameters)]
        name_array = (ctypes.c_char_p * len(names))(
            *[name.encode("UTF-8") for name in names]
        )
        value_array = (ctypes.c_double * len(values))(*values)
        err = self.lib.pymicromegas_assign_values(len(names), name_array, value_array)
        if err:
            raise RuntimeError(f"Could not assign parameter '{names[err - 1]}' (error code: {err}).")
        return None

    def find_value(self, name):
        return self.lib.pymicromegas_find_value(str(name).encode("UTF-8"))

    def set_gauge(self, force_ug=0, vzdecay=0, vwdecay=0):
        self.lib.pymicromegas_set_gauge(int(force_ug), int(vzdecay), int(vwdecay))

    def sort_odd_particles(self):
        cdm_name = ctypes.create_string_buffer(64)
        err = self.lib.pymicromegas_sort_odd_particles(cdm_name, len(cdm_name))
        if err:
            raise RuntimeError(f"Failed to sort odd particles (error code: {err}).")
        return cdm_name.value.decode("UTF-8")

    def load_heff_geff(self, dof_fname):
        dof_path = Path(dof_fname)
        if not dof_path.is_file():
            raise RuntimeError(f"{dof_fname} is not an existing file")
        err = self.lib.pymicromegas_load_heff_geff(str(dof_path.resolve()).encode("UTF-8"))
        if err < 0:
            raise RuntimeError("Failed to load Heff/Geff data: invalid file format")
        if err == 0:
            raise RuntimeError(f"Failed to load Heff/Geff data: cannot open file {dof_fname}")
        return err

    def dark_omega(self, parameters=None, dof_fname=None, fast=1, beps=1e-4):
        if parameters is not None:
            self.assign(parameters)
        if dof_fname is not None:
            self.load_heff_geff(dof_fname)
        self.set_gauge()
        self.sort_odd_particles()
        xf = ctypes.c_double()
        err = ctypes.c_int()
        omega = self.lib.pymicromegas_dark_omega(
            ctypes.byref(xf), int(fast), float(beps), ctypes.byref(err)
        )
        return {"Xf": xf.value, "Omega": omega, "err": err.value}

    def dark_omega2(self, parameters=None, fast=1, beps=1e-4):
        if parameters is not None:
            self.assign(parameters)
        self.set_gauge()
        self.sort_odd_particles()
        return self.lib.pymicromegas_dark_omega2(int(fast), float(beps))

    def function(self, name, restype=ctypes.c_double, argtypes=None):
        func = getattr(self.lib, name)
        func.restype = restype
        if argtypes is not None:
            func.argtypes = argtypes
        return func

    def call(self, name, *args, restype=ctypes.c_double, argtypes=None):
        return self.function(name, restype=restype, argtypes=argtypes)(*args)

    def __getattr__(self, name):
        if name.startswith("__"):
            raise AttributeError(name)
        try:
            return getattr(self.lib, name)
        except AttributeError as exc:
            raise AttributeError(name) from exc
