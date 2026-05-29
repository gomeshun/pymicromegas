import os
import ctypes
import hashlib
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any
import numpy as np
from pandas import Series

from .processes import (
    DirectDetectionProcesses,
    FreezeInProcesses,
    IndirectDetectionProcesses,
    RelicDensityProcesses,
)


PACKAGE_DIR = Path(__file__).resolve().parent
MICROMEGAS_VERSION = "7.1"
MICROMEGAS_PATH_ENV = "PYMICROMEGAS_MICROPATH"
DEFAULT_MICROMEGAS_DIR = PACKAGE_DIR / f"micromegas_{MICROMEGAS_VERSION}"
MICROMEGAS_DIR = Path(os.environ.get(MICROMEGAS_PATH_ENV, DEFAULT_MICROMEGAS_DIR)).expanduser().resolve()
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


def get_micromegas_dir():
    return Path(os.environ.get(MICROMEGAS_PATH_ENV, DEFAULT_MICROMEGAS_DIR)).expanduser().resolve()


def micromegas_env():
    env = os.environ.copy()
    env[MICROMEGAS_PATH_ENV] = str(get_micromegas_dir())
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
        self.micromegas_dir = get_micromegas_dir()
        self.path = str(self.micromegas_dir) + os.sep
        self.verbose = verbose
        os.environ[MICROMEGAS_PATH_ENV] = str(self.micromegas_dir)
        
        if not os.path.isdir(self.path):
            raise RuntimeError(
                "micrOMEGAs is not bundled with pymicromegas. "
                f"Download micrOMEGAs {MICROMEGAS_VERSION} from the official site or Zenodo, "
                f"unpack it, and set {MICROMEGAS_PATH_ENV} to that directory. "
                f"Missing: {self.path}"
            )
    
    
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
        self.path = str(Path(self.interface.path) / project_name) + os.sep
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

    SPECTRUM_SIZE = 250
    RECOIL_SIZE = 120
    ALL_DD_EXPERIMENTS = 0xFFFFFFF
    DIRECT_DETECTION_EXPERIMENTS = {
        "xenon1t_2018": 1,
        "darkside_2018": 2,
        "pico_2019": 4,
        "cresst_2019": 8,
        "lz5t_2024": 16,
    }

    PROCESS_DELEGATES = {
        "get_infl_decay": "relic_density",
        "get_infl_decay_plus": "relic_density",
        "print_thermal_sets": "relic_density",
        "freeze_out_channels": "relic_density",
        "dark_omega": "relic_density",
        "dark_omega2": "relic_density",
        "dark_omega_tr": "relic_density",
        "dark_omega_fo": "relic_density",
        "dark_omega2_tr": "relic_density",
        "dark_omega_n": "relic_density",
        "dark_omega_n_tr": "relic_density",
        "dark_omega_infl": "relic_density",
        "v_sigma_plus23": "freeze_in",
        "v_sigma_plus24": "freeze_in",
        "y_freeze_in_22": "freeze_in",
        "prepare_freeze_in": "freeze_in",
        "freeze_in_channels": "freeze_in",
        "parse_freeze_in_channels": "freeze_in",
        "dark_omega_freeze_in": "freeze_in",
        "dark_omega_freeze_in_22": "freeze_in",
        "dark_omega_freeze_in_decay": "freeze_in",
        "calc_spectrum": "indirect_detection",
        "calc_spectrum_plus": "indirect_detection",
        "decay_spectrum": "indirect_detection",
        "basic_spectra": "indirect_detection",
        "planck_cmb": "indirect_detection",
        "dwarf_signal": "indirect_detection",
        "x_interp": "indirect_detection",
        "z_interp": "indirect_detection",
        "spectr_info": "indirect_detection",
        "spectr_int": "indirect_detection",
        "gamma_flux_tab": "indirect_detection",
        "gamma_flux_tab_gc": "indirect_detection",
        "solar_modulation": "indirect_detection",
        "pbar_background_tab": "indirect_detection",
        "positron_flux": "indirect_detection",
        "positron_flux_tab": "indirect_detection",
        "pbar_flux_tab": "indirect_detection",
        "basic_nu_spectra": "indirect_detection",
        "muon_contained": "indirect_detection",
        "muon_upward": "indirect_detection",
        "nucleon_amplitudes": "direct_detection",
        "dnde_recoil": "direct_detection",
    }

    HEADER_FUNCTION_SIGNATURES = {
        "pNum": (ctypes.c_int, [ctypes.c_char_p]),
        "pMass": (ctypes.c_double, [ctypes.c_char_p]),
        "pdg2name": (ctypes.c_char_p, [ctypes.c_int]),
        "isSMP": (ctypes.c_int, [ctypes.c_int]),
        "readVar": (ctypes.c_int, [ctypes.c_char_p]),
        "Zinvisible": (ctypes.c_int, []),
        "ZpLimCMS": (ctypes.c_double, [ctypes.c_char_p]),
        "setPDT": (ctypes.c_int, [ctypes.c_char_p]),
        "setLHAPDF": (ctypes.c_int, [ctypes.c_char_p, ctypes.c_int]),
        "restorePDF": (ctypes.c_int, [ctypes.c_char_p]),
        "hCollider": (ctypes.c_double, [ctypes.c_double, ctypes.c_int, ctypes.c_int, ctypes.c_double, ctypes.c_double, ctypes.c_char_p, ctypes.c_char_p, ctypes.c_double, ctypes.c_int]),
        "monoJet": (ctypes.c_double, []),
        "convStrFun3": (ctypes.c_double, [ctypes.c_double, ctypes.c_double, ctypes.c_int, ctypes.c_int, ctypes.c_int]),
        "hEff": (ctypes.c_double, [ctypes.c_double]),
        "gEff": (ctypes.c_double, [ctypes.c_double]),
        "gEff2": (ctypes.c_double, [ctypes.c_double]),
        "T_s3": (ctypes.c_double, [ctypes.c_double]),
        "s3_T": (ctypes.c_double, [ctypes.c_double]),
        "h1eff": (ctypes.c_double, [ctypes.c_double, ctypes.c_double, ctypes.c_int]),
        "g1eff": (ctypes.c_double, [ctypes.c_double, ctypes.c_double, ctypes.c_int]),
        "p1eff": (ctypes.c_double, [ctypes.c_double, ctypes.c_double, ctypes.c_int]),
        "n1eff": (ctypes.c_double, [ctypes.c_double, ctypes.c_double, ctypes.c_int]),
        "Hubble": (ctypes.c_double, [ctypes.c_double]),
        "HubbleTime": (ctypes.c_double, [ctypes.c_double, ctypes.c_double]),
        "freeStreaming": (ctypes.c_double, [ctypes.c_double, ctypes.c_double, ctypes.c_double]),
        "hEffLnDiff": (ctypes.c_double, [ctypes.c_double]),
        "vSigmaA": (ctypes.c_double, [ctypes.c_double, ctypes.c_int, ctypes.c_double]),
        "vSigmaS": (ctypes.c_double, [ctypes.c_double, ctypes.c_int, ctypes.c_double]),
        "vSigmaMem": (ctypes.c_double, [ctypes.c_double]),
        "checkTE": (ctypes.c_double, [ctypes.c_int, ctypes.c_double, ctypes.c_int, ctypes.c_double]),
        "Ta": (ctypes.c_double, [ctypes.c_double]),
        "aT": (ctypes.c_double, [ctypes.c_double]),
        "Ha": (ctypes.c_double, [ctypes.c_double]),
        "rhoIa": (ctypes.c_double, [ctypes.c_double]),
        "rhoSMa": (ctypes.c_double, [ctypes.c_double]),
        "Za": (ctypes.c_double, [ctypes.c_double]),
        "ZaEq": (ctypes.c_double, [ctypes.c_double]),
        "Ya": (ctypes.c_double, [ctypes.c_double]),
        "YaEq": (ctypes.c_double, [ctypes.c_double]),
        "ZaN": (ctypes.c_double, [ctypes.c_double, ctypes.c_char_p]),
        "ZaNEq": (ctypes.c_double, [ctypes.c_double, ctypes.c_char_p]),
        "YaN": (ctypes.c_double, [ctypes.c_double, ctypes.c_char_p]),
        "YaNEq": (ctypes.c_double, [ctypes.c_double, ctypes.c_char_p]),
        "vs_aExt": (ctypes.c_double, [ctypes.c_double]),
        "vs_sExt": (ctypes.c_double, [ctypes.c_double]),
        "vs_3Ext": (ctypes.c_double, [ctypes.c_double]),
        "vs_4Ext": (ctypes.c_double, [ctypes.c_double]),
        "vs1120F": (ctypes.c_double, [ctypes.c_double]),
        "vs2200F": (ctypes.c_double, [ctypes.c_double]),
        "vs1100F": (ctypes.c_double, [ctypes.c_double]),
        "vs1210F": (ctypes.c_double, [ctypes.c_double]),
        "vs1122F": (ctypes.c_double, [ctypes.c_double]),
        "vs2211F": (ctypes.c_double, [ctypes.c_double]),
        "vs1110F": (ctypes.c_double, [ctypes.c_double]),
        "vs2220F": (ctypes.c_double, [ctypes.c_double]),
        "vs1112F": (ctypes.c_double, [ctypes.c_double]),
        "vs1222F": (ctypes.c_double, [ctypes.c_double]),
        "vs1220F": (ctypes.c_double, [ctypes.c_double]),
        "vs2210F": (ctypes.c_double, [ctypes.c_double]),
        "vs2221F": (ctypes.c_double, [ctypes.c_double]),
        "vs1211F": (ctypes.c_double, [ctypes.c_double]),
        "TCoeffF": (ctypes.c_double, [ctypes.c_double]),
        "Y1F": (ctypes.c_double, [ctypes.c_double]),
        "Y2F": (ctypes.c_double, [ctypes.c_double]),
        "YF": (ctypes.c_double, [ctypes.c_double]),
        "oneChannel": (ctypes.c_double, [ctypes.c_double, ctypes.c_double, ctypes.c_char_p, ctypes.c_char_p, ctypes.c_char_p, ctypes.c_char_p]),
        "Yeq": (ctypes.c_double, [ctypes.c_double]),
        "Yeq1": (ctypes.c_double, [ctypes.c_double]),
        "Yeq2": (ctypes.c_double, [ctypes.c_double]),
        "defThermalSet": (ctypes.c_int, [ctypes.c_int, ctypes.c_char_p]),
        "YdmNEq": (ctypes.c_double, [ctypes.c_double, ctypes.c_char_p]),
        "YdmN": (ctypes.c_double, [ctypes.c_double, ctypes.c_char_p]),
        "vSigmaN": (ctypes.c_double, [ctypes.c_double, ctypes.c_char_p]),
        "setFastBeps": (None, [ctypes.c_int, ctypes.c_double]),
        "decayAbundance": (ctypes.c_double, [ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_int, ctypes.c_double, ctypes.c_int]),
        "YFi": (ctypes.c_double, [ctypes.c_double]),
        "Zi": (ctypes.c_double, [ctypes.c_int]),
        "SYukawa": (ctypes.c_double, [ctypes.c_double, ctypes.c_double]),
        "SHulthen": (ctypes.c_double, [ctypes.c_double, ctypes.c_double]),
        "setClumpConst": (None, [ctypes.c_double, ctypes.c_double]),
        "rhoClumpsConst": (ctypes.c_double, [ctypes.c_double]),
        "HaloFactor": (ctypes.c_double, [ctypes.c_double, ctypes.c_double]),
        "gammaFlux": (ctypes.c_double, [ctypes.c_double, ctypes.c_double, ctypes.c_double]),
        "gammaFluxGC": (ctypes.c_double, [ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double]),
        "hProfileZhao": (ctypes.c_double, [ctypes.c_double]),
        "setProfileZhao": (None, [ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double]),
        "hProfileEinasto": (ctypes.c_double, [ctypes.c_double]),
        "setProfileEinasto": (None, [ctypes.c_double, ctypes.c_double]),
        "noClumps": (ctypes.c_double, [ctypes.c_double]),
        "pBarBackgroundFlux": (ctypes.c_double, [ctypes.c_double]),
        "pbarFlux": (ctypes.c_double, [ctypes.c_double, ctypes.c_double]),
        "FSRdNdE": (ctypes.c_double, [ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_int, ctypes.c_int]),
        "nuAttenuation": (ctypes.c_double, [ctypes.c_int, ctypes.c_double, ctypes.c_double]),
        "atmNuFlux": (ctypes.c_double, [ctypes.c_int, ctypes.c_double, ctypes.c_double]),
        "atmNuFluxI": (ctypes.c_double, [ctypes.c_int, ctypes.c_double, ctypes.c_double]),
        "calcScalarFF": (None, [ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double]),
        "calcScalarQuarkFF": (None, [ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double]),
        "SetFermi": (None, [ctypes.c_double, ctypes.c_double, ctypes.c_double]),
        "FermiFF": (ctypes.c_double, [ctypes.c_int, ctypes.c_double]),
        "Maxwell": (ctypes.c_double, [ctypes.c_double]),
        "SHMpp": (ctypes.c_double, [ctypes.c_double]),
        "setSpinDepFF": (ctypes.c_int, [ctypes.c_int, ctypes.c_int]),
        "MaxGapLim": (ctypes.c_double, [ctypes.c_double, ctypes.c_double]),
        "widthSMh": (ctypes.c_double, [ctypes.c_double]),
        "brSMhGG": (ctypes.c_double, [ctypes.c_double]),
        "brSMhAA": (ctypes.c_double, [ctypes.c_double]),
        "darkOmegaNu": (ctypes.c_double, [ctypes.c_double, ctypes.c_double, ctypes.c_int, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double]),
        "LasymmS": (ctypes.c_double, [ctypes.c_double]),
        "callSuperIsoSLHA": (ctypes.c_int, []),
        "Xe1TpEff0": (ctypes.c_double, [ctypes.c_double]),
        "Xe1TpEff1": (ctypes.c_double, [ctypes.c_double]),
        "Xe1TpEff2": (ctypes.c_double, [ctypes.c_double]),
        "XENON1T": (ctypes.c_double, [ctypes.c_double]),
        "DS50": (ctypes.c_double, [ctypes.c_double]),
        "DS50_noB": (ctypes.c_double, [ctypes.c_double]),
        "CRESST_III": (ctypes.c_double, [ctypes.c_double]),
        "CRESST_III_SDn": (ctypes.c_double, [ctypes.c_double]),
        "PICO60": (ctypes.c_double, [ctypes.c_double]),
        "PICO60_SDp": (ctypes.c_double, [ctypes.c_double]),
        "XENON1T_SDp": (ctypes.c_double, [ctypes.c_double]),
        "XENON1T_SDn": (ctypes.c_double, [ctypes.c_double]),
        "PandaX4T": (ctypes.c_double, [ctypes.c_double]),
        "PandaXS2": (ctypes.c_double, [ctypes.c_double]),
        "LZ5T": (ctypes.c_double, [ctypes.c_double]),
        "LZ5T_SDn": (ctypes.c_double, [ctypes.c_double]),
        "LZ5T_SDp": (ctypes.c_double, [ctypes.c_double]),
        "XENON_NT": (ctypes.c_double, [ctypes.c_double]),
        "Darwin_XLZD": (ctypes.c_double, [ctypes.c_double]),
        "Neutrino_FloorXeSI": (ctypes.c_double, [ctypes.c_double]),
    }

    HEADER_METHODS = {
        "particle_number": "pNum",
        "particle_mass": "pMass",
        "pdg_name": "pdg2name",
        "is_sm_particle": "isSMP",
        "read_var_file": "readVar",
        "z_invisible": "Zinvisible",
        "zp_lim_cms": "ZpLimCMS",
        "set_pdt": "setPDT",
        "set_lhapdf": "setLHAPDF",
        "restore_pdf": "restorePDF",
        "h_collider": "hCollider",
        "monojet": "monoJet",
        "conv_str_fun3": "convStrFun3",
        "h_eff": "hEff",
        "g_eff": "gEff",
        "g_eff2": "gEff2",
        "temperature_from_s3": "T_s3",
        "s3_from_temperature": "s3_T",
        "h1_eff": "h1eff",
        "g1_eff": "g1eff",
        "p1_eff": "p1eff",
        "n1_eff": "n1eff",
        "hubble": "Hubble",
        "hubble_time": "HubbleTime",
        "free_streaming": "freeStreaming",
        "h_eff_ln_diff": "hEffLnDiff",
        "v_sigma_a": "vSigmaA",
        "v_sigma_s": "vSigmaS",
        "v_sigma_mem": "vSigmaMem",
        "check_te": "checkTE",
        "t_of_a": "Ta",
        "a_of_t": "aT",
        "h_of_a": "Ha",
        "rho_inflaton_a": "rhoIa",
        "rho_sm_a": "rhoSMa",
        "z_of_a": "Za",
        "z_eq_of_a": "ZaEq",
        "y_of_a": "Ya",
        "y_eq_of_a": "YaEq",
        "z_n_of_a": "ZaN",
        "z_n_eq_of_a": "ZaNEq",
        "y_n_of_a": "YaN",
        "y_n_eq_of_a": "YaNEq",
        "vs_a_ext": "vs_aExt",
        "vs_s_ext": "vs_sExt",
        "vs_3_ext": "vs_3Ext",
        "vs_4_ext": "vs_4Ext",
        "vs1120_f": "vs1120F",
        "vs2200_f": "vs2200F",
        "vs1100_f": "vs1100F",
        "vs1210_f": "vs1210F",
        "vs1122_f": "vs1122F",
        "vs2211_f": "vs2211F",
        "vs1110_f": "vs1110F",
        "vs2220_f": "vs2220F",
        "vs1112_f": "vs1112F",
        "vs1222_f": "vs1222F",
        "vs1220_f": "vs1220F",
        "vs2210_f": "vs2210F",
        "vs2221_f": "vs2221F",
        "vs1211_f": "vs1211F",
        "t_coeff_f": "TCoeffF",
        "y1_f": "Y1F",
        "y2_f": "Y2F",
        "y_f": "YF",
        "one_channel": "oneChannel",
        "y_eq": "Yeq",
        "y_eq1": "Yeq1",
        "y_eq2": "Yeq2",
        "define_thermal_set": "defThermalSet",
        "y_dm_n_eq": "YdmNEq",
        "y_dm_n": "YdmN",
        "v_sigma_n": "vSigmaN",
        "set_fast_beps": "setFastBeps",
        "decay_abundance": "decayAbundance",
        "y_fi": "YFi",
        "zi": "Zi",
        "sommerfeld_yukawa": "SYukawa",
        "sommerfeld_hulthen": "SHulthen",
        "set_clump_const": "setClumpConst",
        "rho_clumps_const": "rhoClumpsConst",
        "halo_factor": "HaloFactor",
        "gamma_flux": "gammaFlux",
        "gamma_flux_gc": "gammaFluxGC",
        "halo_profile_zhao": "hProfileZhao",
        "set_profile_zhao": "setProfileZhao",
        "halo_profile_einasto": "hProfileEinasto",
        "set_profile_einasto": "setProfileEinasto",
        "no_clumps": "noClumps",
        "pbar_background_flux": "pBarBackgroundFlux",
        "pbar_flux": "pbarFlux",
        "fsr_dnde": "FSRdNdE",
        "nu_attenuation": "nuAttenuation",
        "atm_nu_flux": "atmNuFlux",
        "atm_nu_flux_i": "atmNuFluxI",
        "calc_scalar_ff": "calcScalarFF",
        "calc_scalar_quark_ff": "calcScalarQuarkFF",
        "set_fermi": "SetFermi",
        "fermi_ff": "FermiFF",
        "maxwell": "Maxwell",
        "shmpp": "SHMpp",
        "set_spin_dependent_ff": "setSpinDepFF",
        "max_gap_limit": "MaxGapLim",
        "sm_higgs_width": "widthSMh",
        "sm_higgs_br_gg": "brSMhGG",
        "sm_higgs_br_aa": "brSMhAA",
        "dark_omega_nu": "darkOmegaNu",
        "lepton_asymmetry_s": "LasymmS",
        "call_superiso_slha": "callSuperIsoSLHA",
        "xenon1t_efficiency0": "Xe1TpEff0",
        "xenon1t_efficiency1": "Xe1TpEff1",
        "xenon1t_efficiency2": "Xe1TpEff2",
        "xenon1t_limit": "XENON1T",
        "ds50_limit": "DS50",
        "ds50_no_background_limit": "DS50_noB",
        "cresst_iii_limit": "CRESST_III",
        "cresst_iii_sdn_limit": "CRESST_III_SDn",
        "pico60_limit": "PICO60",
        "pico60_sdp_limit": "PICO60_SDp",
        "xenon1t_sdp_limit": "XENON1T_SDp",
        "xenon1t_sdn_limit": "XENON1T_SDn",
        "pandax4t_limit": "PandaX4T",
        "pandaxs2_limit": "PandaXS2",
        "lz5t_limit": "LZ5T",
        "lz5t_sdn_limit": "LZ5T_SDn",
        "lz5t_sdp_limit": "LZ5T_SDp",
        "xenon_nt_limit": "XENON_NT",
        "darwin_xlzd_limit": "Darwin_XLZD",
        "neutrino_floor_xe_si": "Neutrino_FloorXeSI",
    }

    BRIDGE_SOURCE = r"""
#include <stdio.h>
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

int pymicromegas_print_vars(const char *path)
{
    FILE *file = fopen(path, "w");
    if(!file) return 1;
    printVar(file);
    fclose(file);
    return 0;
}

int pymicromegas_print_masses(const char *path, int sort)
{
    FILE *file = fopen(path, "w");
    if(!file) return 1;
    printMasses(file, sort);
    fclose(file);
    return 0;
}

int pymicromegas_print_higgs(const char *path)
{
    FILE *file = fopen(path, "w");
    if(!file) return 1;
    printHiggs(file);
    fclose(file);
    return 0;
}

double pymicromegas_particle_width(const char *name)
{
    txtList decays = NULL;
    double width = pWidth(name, &decays);
    if(decays) cleanTxtList(decays);
    return width;
}

int pymicromegas_decay_info(const char *name, const char *path, double *width)
{
    FILE *file = fopen(path, "w");
    if(!file) return 1;
    *width = decay2Info((char *)name, file);
    fclose(file);
    return 0;
}

const char *pymicromegas_next_odd(int number, double *mass)
{
    return nextOdd(number, mass);
}

int pymicromegas_lsp_nlsp_lep(double *cross_section_limit)
{
    return LspNlsp_LEP(cross_section_limit);
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

int pymicromegas_to_feeble_list(const char *name)
{
    return toFeebleList((char *)name);
}

int pymicromegas_is_feeble(const char *name)
{
    return isFeeble((char *)name);
}

int pymicromegas_n_feeble(void)
{
    return nFeeble;
}

double pymicromegas_dark_omega(double *xf, int fast, double beps, int *err)
{
    return darkOmega(xf, fast, beps, err);
}

double pymicromegas_dark_omega2(int fast, double beps, int *err)
{
    return darkOmega2((double)fast, beps, err);
}

double pymicromegas_dark_omega_fo(double *xf, int fast, double beps, int *err)
{
    return darkOmegaFO(xf, fast, beps, err);
}

double pymicromegas_dark_omega_tr(double tr, double yr, int fast, double beps, int *err)
{
    return darkOmegaTR(tr, yr, fast, beps, err);
}

double pymicromegas_dark_omega2_tr(double tr, double y1r, double y2r, int fast, double beps, int *err)
{
    return darkOmega2TR(tr, y1r, y2r, (double)fast, beps, err);
}

double pymicromegas_dark_omega_n(int fast, double beps, int *err)
{
    return darkOmegaN(fast, beps, err);
}

double pymicromegas_dark_omega_n_tr(double tr, double *y, int fast, double beps, int *err)
{
    return darkOmegaNTR(tr, y, fast, beps, err);
}

double pymicromegas_dark_omega_infl(double branching, double beps, double *tfo, int *err)
{
    return darkOmegaInfl(branching, beps, tfo, err);
}

int pymicromegas_get_infl_decay(double h0, double gamma, double *trh, double *tmax, double *aend)
{
    return getInflDecay(h0, gamma, trh, tmax, aend);
}

int pymicromegas_get_infl_decay_plus(double hsm, double hi, double gamma, double alpha, double omega, double *trh, double *tmax, double *aend)
{
    return getInflDecayPlus(hsm, hi, gamma, alpha, omega, trh, tmax, aend);
}

double pymicromegas_v_sigma_plus23(const char *process, double temperature, int *err)
{
    return vSigmaPlus23((char *)process, temperature, err);
}

double pymicromegas_v_sigma_plus24(const char *process, double temperature, int *err)
{
    return vSigmaPlus24((char *)process, temperature, err);
}

double pymicromegas_y_freeze_in22(const char *process, double t0, double tr, int plot_dydt, int *err)
{
    return YfreezeIn22((char *)process, t0, tr, plot_dydt, err);
}

int pymicromegas_print_thermal_sets(const char *path)
{
    FILE *file = fopen(path, "w");
    if(!file) return 1;
    FILE *old_stdout = stdout;
    stdout = file;
    printThermalSets();
    fflush(file);
    stdout = old_stdout;
    fclose(file);
    return 0;
}

double pymicromegas_print_channels(const char *path, double xf, double cut, double beps, int percent)
{
    FILE *file = fopen(path, "w");
    double result;
    if(!file) return -1.0;
    result = printChannels(xf, cut, beps, percent, file);
    fclose(file);
    return result;
}

double pymicromegas_dark_omega_fi(double tr, const char *name, int *err)
{
    return darkOmegaFi(tr, (char *)name, err);
}

double pymicromegas_dark_omega_fi22(double tr, const char *process, const char *name, int *err)
{
    return darkOmegaFi22(tr, (char *)process, (char *)name, err);
}

double pymicromegas_dark_omega_fi_decay(double tr, const char *bath_particle, const char *name)
{
    return darkOmegaFiDecay(tr, (char *)bath_particle, (char *)name);
}

int pymicromegas_print_channels_fi(const char *path, double cut, int percent)
{
    FILE *file = fopen(path, "w");
    if(!file) return 1;
    printChannelsFi(cut, percent, file);
    fclose(file);
    return 0;
}

const char *pymicromegas_cdm_name(int sector)
{
    if(CDM && sector >= 1 && sector <= Ncdm && CDM[sector]) return CDM[sector];
    if(sector == 1) return CDM1;
    if(sector == 2) return CDM2;
    return NULL;
}

double pymicromegas_cdm_mass(int sector)
{
    if(McdmN && sector >= 1 && sector <= Ncdm) return McdmN[sector];
    if(sector == 1) return Mcdm;
    return 0.0;
}

double pymicromegas_cdm_fraction(int sector)
{
    if(fracCDM && sector >= 1 && sector <= Ncdm) return fracCDM[sector];
    return sector == 1 ? 1.0 : 0.0;
}

const char *pymicromegas_cdm1(void)
{
    return pymicromegas_cdm_name(1);
}

const char *pymicromegas_cdm2(void)
{
    return pymicromegas_cdm_name(2);
}

double pymicromegas_mcdm(void)
{
    return Mcdm;
}

double pymicromegas_mcdm1(void)
{
    return pymicromegas_cdm_mass(1);
}

double pymicromegas_mcdm2(void)
{
    return pymicromegas_cdm_mass(2);
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
        bridge_tag = hashlib.sha256(self.BRIDGE_SOURCE.encode("UTF-8")).hexdigest()[:12]
        library_name = f"libpymicromegas_{bridge_tag}.dylib" if sys.platform == "darwin" else f"libpymicromegas_{bridge_tag}.so"
        self.library_path = self.path / library_name
        self._lib = None
        self.relic_density = RelicDensityProcesses(self)
        self.freeze_in = FreezeInProcesses(self)
        self.indirect_detection = IndirectDetectionProcesses(self)
        self.direct_detection = DirectDetectionProcesses(self)

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

    def bridge_source_is_stale(self):
        return not self.bridge_source.exists() or self.bridge_source.read_text() != self.BRIDGE_SOURCE

    def bridge_library_is_stale(self):
        if not self.library_path.exists():
            return True
        if self.bridge_source_is_stale():
            return True
        return self.bridge_source.exists() and self.bridge_source.stat().st_mtime > self.library_path.stat().st_mtime

    def write_bridge_source(self, force=False):
        if force or self.bridge_source_is_stale():
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
        bridge_library_was_stale = self.bridge_library_is_stale()
        self.write_bridge_source(force=force)
        if force or bridge_library_was_stale:
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
        lib.pymicromegas_print_vars.argtypes = [ctypes.c_char_p]
        lib.pymicromegas_print_vars.restype = ctypes.c_int
        lib.pymicromegas_print_masses.argtypes = [ctypes.c_char_p, ctypes.c_int]
        lib.pymicromegas_print_masses.restype = ctypes.c_int
        lib.pymicromegas_print_higgs.argtypes = [ctypes.c_char_p]
        lib.pymicromegas_print_higgs.restype = ctypes.c_int
        lib.pymicromegas_particle_width.argtypes = [ctypes.c_char_p]
        lib.pymicromegas_particle_width.restype = ctypes.c_double
        lib.pymicromegas_decay_info.argtypes = [
            ctypes.c_char_p,
            ctypes.c_char_p,
            ctypes.POINTER(ctypes.c_double),
        ]
        lib.pymicromegas_decay_info.restype = ctypes.c_int
        lib.pymicromegas_next_odd.argtypes = [ctypes.c_int, ctypes.POINTER(ctypes.c_double)]
        lib.pymicromegas_next_odd.restype = ctypes.c_char_p
        lib.pymicromegas_lsp_nlsp_lep.argtypes = [ctypes.POINTER(ctypes.c_double)]
        lib.pymicromegas_lsp_nlsp_lep.restype = ctypes.c_int
        lib.pymicromegas_set_gauge.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.c_int]
        lib.pymicromegas_set_gauge.restype = None
        lib.pymicromegas_sort_odd_particles.argtypes = [ctypes.c_char_p, ctypes.c_int]
        lib.pymicromegas_sort_odd_particles.restype = ctypes.c_int
        lib.pymicromegas_load_heff_geff.argtypes = [ctypes.c_char_p]
        lib.pymicromegas_load_heff_geff.restype = ctypes.c_int
        lib.pymicromegas_to_feeble_list.argtypes = [ctypes.c_char_p]
        lib.pymicromegas_to_feeble_list.restype = ctypes.c_int
        lib.pymicromegas_is_feeble.argtypes = [ctypes.c_char_p]
        lib.pymicromegas_is_feeble.restype = ctypes.c_int
        lib.pymicromegas_n_feeble.argtypes = []
        lib.pymicromegas_n_feeble.restype = ctypes.c_int
        lib.pymicromegas_dark_omega.argtypes = [
            ctypes.POINTER(ctypes.c_double),
            ctypes.c_int,
            ctypes.c_double,
            ctypes.POINTER(ctypes.c_int),
        ]
        lib.pymicromegas_dark_omega.restype = ctypes.c_double
        lib.pymicromegas_dark_omega2.argtypes = [ctypes.c_int, ctypes.c_double, ctypes.POINTER(ctypes.c_int)]
        lib.pymicromegas_dark_omega2.restype = ctypes.c_double
        lib.pymicromegas_dark_omega_fo.argtypes = [
            ctypes.POINTER(ctypes.c_double),
            ctypes.c_int,
            ctypes.c_double,
            ctypes.POINTER(ctypes.c_int),
        ]
        lib.pymicromegas_dark_omega_fo.restype = ctypes.c_double
        lib.pymicromegas_dark_omega_tr.argtypes = [
            ctypes.c_double,
            ctypes.c_double,
            ctypes.c_int,
            ctypes.c_double,
            ctypes.POINTER(ctypes.c_int),
        ]
        lib.pymicromegas_dark_omega_tr.restype = ctypes.c_double
        lib.pymicromegas_dark_omega2_tr.argtypes = [
            ctypes.c_double,
            ctypes.c_double,
            ctypes.c_double,
            ctypes.c_int,
            ctypes.c_double,
            ctypes.POINTER(ctypes.c_int),
        ]
        lib.pymicromegas_dark_omega2_tr.restype = ctypes.c_double
        lib.pymicromegas_dark_omega_n.argtypes = [ctypes.c_int, ctypes.c_double, ctypes.POINTER(ctypes.c_int)]
        lib.pymicromegas_dark_omega_n.restype = ctypes.c_double
        lib.pymicromegas_dark_omega_n_tr.argtypes = [
            ctypes.c_double,
            ctypes.POINTER(ctypes.c_double),
            ctypes.c_int,
            ctypes.c_double,
            ctypes.POINTER(ctypes.c_int),
        ]
        lib.pymicromegas_dark_omega_n_tr.restype = ctypes.c_double
        lib.pymicromegas_dark_omega_infl.argtypes = [
            ctypes.c_double,
            ctypes.c_double,
            ctypes.POINTER(ctypes.c_double),
            ctypes.POINTER(ctypes.c_int),
        ]
        lib.pymicromegas_dark_omega_infl.restype = ctypes.c_double
        lib.pymicromegas_get_infl_decay.argtypes = [
            ctypes.c_double,
            ctypes.c_double,
            ctypes.POINTER(ctypes.c_double),
            ctypes.POINTER(ctypes.c_double),
            ctypes.POINTER(ctypes.c_double),
        ]
        lib.pymicromegas_get_infl_decay.restype = ctypes.c_int
        lib.pymicromegas_get_infl_decay_plus.argtypes = [
            ctypes.c_double,
            ctypes.c_double,
            ctypes.c_double,
            ctypes.c_double,
            ctypes.c_double,
            ctypes.POINTER(ctypes.c_double),
            ctypes.POINTER(ctypes.c_double),
            ctypes.POINTER(ctypes.c_double),
        ]
        lib.pymicromegas_get_infl_decay_plus.restype = ctypes.c_int
        lib.pymicromegas_v_sigma_plus23.argtypes = [ctypes.c_char_p, ctypes.c_double, ctypes.POINTER(ctypes.c_int)]
        lib.pymicromegas_v_sigma_plus23.restype = ctypes.c_double
        lib.pymicromegas_v_sigma_plus24.argtypes = [ctypes.c_char_p, ctypes.c_double, ctypes.POINTER(ctypes.c_int)]
        lib.pymicromegas_v_sigma_plus24.restype = ctypes.c_double
        lib.pymicromegas_y_freeze_in22.argtypes = [
            ctypes.c_char_p,
            ctypes.c_double,
            ctypes.c_double,
            ctypes.c_int,
            ctypes.POINTER(ctypes.c_int),
        ]
        lib.pymicromegas_y_freeze_in22.restype = ctypes.c_double
        lib.pymicromegas_print_thermal_sets.argtypes = [ctypes.c_char_p]
        lib.pymicromegas_print_thermal_sets.restype = ctypes.c_int
        lib.pymicromegas_print_channels.argtypes = [
            ctypes.c_char_p,
            ctypes.c_double,
            ctypes.c_double,
            ctypes.c_double,
            ctypes.c_int,
        ]
        lib.pymicromegas_print_channels.restype = ctypes.c_double
        lib.pymicromegas_dark_omega_fi.argtypes = [
            ctypes.c_double,
            ctypes.c_char_p,
            ctypes.POINTER(ctypes.c_int),
        ]
        lib.pymicromegas_dark_omega_fi.restype = ctypes.c_double
        lib.pymicromegas_dark_omega_fi22.argtypes = [
            ctypes.c_double,
            ctypes.c_char_p,
            ctypes.c_char_p,
            ctypes.POINTER(ctypes.c_int),
        ]
        lib.pymicromegas_dark_omega_fi22.restype = ctypes.c_double
        lib.pymicromegas_dark_omega_fi_decay.argtypes = [
            ctypes.c_double,
            ctypes.c_char_p,
            ctypes.c_char_p,
        ]
        lib.pymicromegas_dark_omega_fi_decay.restype = ctypes.c_double
        lib.pymicromegas_print_channels_fi.argtypes = [ctypes.c_char_p, ctypes.c_double, ctypes.c_int]
        lib.pymicromegas_print_channels_fi.restype = ctypes.c_int
        lib.pymicromegas_cdm_name.argtypes = [ctypes.c_int]
        lib.pymicromegas_cdm_name.restype = ctypes.c_char_p
        lib.pymicromegas_cdm_mass.argtypes = [ctypes.c_int]
        lib.pymicromegas_cdm_mass.restype = ctypes.c_double
        lib.pymicromegas_cdm_fraction.argtypes = [ctypes.c_int]
        lib.pymicromegas_cdm_fraction.restype = ctypes.c_double
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
        self._available_header_functions = set()
        for name, (restype, argtypes) in self.HEADER_FUNCTION_SIGNATURES.items():
            try:
                function = getattr(lib, name)
            except AttributeError:
                continue
            function.restype = restype
            function.argtypes = argtypes
            self._available_header_functions.add(name)

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

    def _resolve_particle_name(self, particle_name=None, sector=1):
        if particle_name is not None:
            return str(particle_name)
        encoded_name = self.lib.pymicromegas_cdm_name(int(sector))
        if encoded_name is None:
            raise RuntimeError(f"Could not resolve CDM particle name for sector {sector}.")
        return encoded_name.decode("UTF-8")

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

    @staticmethod
    def _encode_string(value):
        if value is None:
            return None
        if isinstance(value, bytes):
            return value
        return str(value).encode("UTF-8")

    @classmethod
    def _coerce_argument(cls, argtype, value):
        if argtype is ctypes.c_char_p:
            return cls._encode_string(value)
        if argtype in (ctypes.c_int, ctypes.c_uint, ctypes.c_long, ctypes.c_ulong):
            return int(value)
        if argtype is ctypes.c_double:
            return float(value)
        return value

    @staticmethod
    def _decode_result(restype, value):
        if restype is ctypes.c_char_p:
            return None if value is None else value.decode("UTF-8")
        if restype is ctypes.c_int:
            return int(value)
        if restype is ctypes.c_double:
            return float(value)
        return value

    @property
    def available_header_functions(self):
        return sorted(getattr(self, "_available_header_functions", set()))

    @property
    def available_header_methods(self):
        available = set(getattr(self, "_available_header_functions", set()))
        return sorted(method for method, function in self.HEADER_METHODS.items() if function in available)

    def _call_header_function(self, name, *args):
        if name not in self.HEADER_FUNCTION_SIGNATURES:
            raise AttributeError(name)
        if name not in getattr(self, "_available_header_functions", set()):
            raise RuntimeError(f"micrOMEGAs symbol '{name}' is not available in this project library.")
        restype, argtypes = self.HEADER_FUNCTION_SIGNATURES[name]
        if len(args) != len(argtypes):
            raise TypeError(f"{name} expects {len(argtypes)} arguments, got {len(args)}.")
        function = getattr(self.lib, name)
        coerced_args = [self._coerce_argument(argtype, value) for argtype, value in zip(argtypes, args)]
        return self._decode_result(restype, function(*coerced_args))

    def particle_number(self, particle_name) -> int:
        return int(self.lib.pNum(self._encode_string(particle_name)))

    def particle_mass(self, particle_name) -> float:
        return float(self.lib.pMass(self._encode_string(particle_name)))

    def pdg_name(self, pdg) -> str | None:
        encoded_name = self.lib.pdg2name(int(pdg))
        return None if encoded_name is None else encoded_name.decode("UTF-8")

    def h_eff(self, temperature) -> float:
        return float(self.lib.hEff(float(temperature)))

    def g_eff(self, temperature) -> float:
        return float(self.lib.gEff(float(temperature)))

    def hubble(self, temperature) -> float:
        return float(self.lib.Hubble(float(temperature)))

    def _read_bridge_text(self, bridge_function, *args):
        tmp_path = None
        try:
            with tempfile.NamedTemporaryFile("w+", delete=False) as tmp_file:
                tmp_path = Path(tmp_file.name)
            err = bridge_function(str(tmp_path).encode("UTF-8"), *args)
            if isinstance(err, int) and err:
                raise RuntimeError("micrOMEGAs failed to write text output.")
            return tmp_path.read_text()
        finally:
            if tmp_path is not None:
                tmp_path.unlink(missing_ok=True)

    def print_vars(self):
        return self._read_bridge_text(self.lib.pymicromegas_print_vars)

    def print_masses(self, sort=True):
        return self._read_bridge_text(self.lib.pymicromegas_print_masses, int(bool(sort)))

    def print_higgs(self):
        return self._read_bridge_text(self.lib.pymicromegas_print_higgs)

    def particle_width(self, particle_name):
        return self.lib.pymicromegas_particle_width(self._encode_string(particle_name))

    def decay_info(self, particle_name):
        tmp_path = None
        width = ctypes.c_double()
        try:
            with tempfile.NamedTemporaryFile("w+", delete=False) as tmp_file:
                tmp_path = Path(tmp_file.name)
            err = self.lib.pymicromegas_decay_info(
                self._encode_string(particle_name), str(tmp_path).encode("UTF-8"), ctypes.byref(width)
            )
            if err:
                raise RuntimeError(f"Failed to calculate decay information for {particle_name}.")
            return {"width": width.value, "text": tmp_path.read_text()}
        finally:
            if tmp_path is not None:
                tmp_path.unlink(missing_ok=True)

    def next_odd(self, number):
        mass = ctypes.c_double()
        encoded_name = self.lib.pymicromegas_next_odd(int(number), ctypes.byref(mass))
        name = None if encoded_name is None else encoded_name.decode("UTF-8")
        return {"name": name, "mass": mass.value}

    def lsp_nlsp_lep(self):
        cross_section_limit = ctypes.c_double()
        excluded = self.lib.pymicromegas_lsp_nlsp_lep(ctypes.byref(cross_section_limit))
        return {"excluded": bool(excluded), "cross_section_limit": cross_section_limit.value}

    def to_feeble_list(self, particle_name):
        encoded_name = None if particle_name is None else str(particle_name).encode("UTF-8")
        err = self.lib.pymicromegas_to_feeble_list(encoded_name)
        if err:
            raise RuntimeError(f"Could not add '{particle_name}' to the feeble-particle list.")
        return None

    def clear_feeble_list(self):
        self.to_feeble_list(None)

    def is_feeble(self, particle_name):
        return bool(self.lib.pymicromegas_is_feeble(str(particle_name).encode("UTF-8")))

    @property
    def n_feeble(self):
        return self.lib.pymicromegas_n_feeble()

    @property
    def n_cdm(self):
        return ctypes.c_int.in_dll(self.lib, "Ncdm").value

    @classmethod
    def _new_double_array(cls, size=None):
        return (ctypes.c_double * int(size or cls.SPECTRUM_SIZE))()

    @staticmethod
    def _numpy_from_array(array, size=None):
        view = np.ctypeslib.as_array(array)
        if size is not None:
            view = view[:size]
        return view.copy()

    @classmethod
    def _as_double_array(cls, values, size=None):
        size = int(size or cls.SPECTRUM_SIZE)
        array = cls._new_double_array(size)
        for index, value in enumerate(values):
            if index >= size:
                break
            array[index] = float(value)
        return array

    def function(self, name, restype: object = ctypes.c_double, argtypes=None):
        func = getattr(self.lib, name)
        func.restype = restype
        if argtypes is not None:
            func.argtypes = argtypes
        return func

    def call(self, name, *args, restype: object = ctypes.c_double, argtypes=None):
        return self.function(name, restype=restype, argtypes=argtypes)(*args)

    def __getattr__(self, name) -> Any:
        if name.startswith("__"):
            raise AttributeError(name)
        delegate_name = self.PROCESS_DELEGATES.get(name)
        if delegate_name is not None:
            return getattr(getattr(self, delegate_name), name)
        if name in self.HEADER_METHODS:
            function_name = self.HEADER_METHODS[name]

            def wrapper(*args):
                return self._call_header_function(function_name, *args)

            wrapper.__name__ = name
            return wrapper
        try:
            return getattr(self.lib, name)
        except AttributeError as exc:
            raise AttributeError(name) from exc
