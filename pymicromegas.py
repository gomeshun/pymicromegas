import os
import subprocess
import sys
import warnings
from pathlib import Path 
import numpy as np
from pandas import Series

try:
    from ._native import NativeMicrOmegas
except ImportError:
    from _native import NativeMicrOmegas


dir_micromegas = os.path.dirname(__file__) + "/micromegas_5.0.8/"



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

NATIVE_OMEGA_FLAGS = ("OMEGA",)


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
    
    
def run_bash(command,shell=True,stdout=subprocess.PIPE,encoding="UTF-8",check=True,input=None,cwd=None,verbose=True):
    if verbose: print(command)
    process = subprocess.run(command,shell=shell,stdout=stdout,encoding=encoding,check=check,input=input,cwd=cwd)
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

######## class definitions ########
class PyMicrOmegas:    
    
    def __init__(self,verbose=False):
        self.path = dir_micromegas
        
        if not os.path.isdir(self.path):
            print("PyMicrOmegas: unzip micromegas...")
            run_bash("tar -xzf micromegas_5.0.8.tgz",cwd=os.path.dirname(__file__))

        if os.path.isfile(self.path + "include/microPath.h"): 
            pass
        
        else:
            print("PyMicrOmegas: micromegas are not installed yet. Start install automatically...")
            output = self.compile_micromegas()
            print("PyMicrOmegas: micromegas are installed.")
            if verbose: print(output)
            #raise RuntimeError("micromegas is not properly installed. Run 'make all' in micromegas_5.0.8 directory")
    
    
    def run_bash(self,command,shell=True,stdout=subprocess.PIPE,encoding="UTF-8",check=False,input=None,verbose=True):
        return run_bash(command,shell=shell,stdout=stdout,encoding=encoding,check=check,input=input,cwd=self.path,verbose=verbose)
    
    
    def compile_micromegas(self):
        print("Compiling micromegas...")
        return self.run_bash("make")

    
    def clean_micromegas(self):
        print("Cleaning micromegas...")
        return self.run_bash("make clean",input="y\n")
    
    
    def project_exists(self,project_name):
        return os.path.isdir(self.path + project_name)
       
        
    def load_modified_main(self,project_name,return_project=False):
        print("Loading modified main files...")
        if not self.project_exists(project_name): 
            raise RuntimeError("project '{}' does not exist yet.".format(project_name))
            
        commands = [ "mv {0}/main.c {0}/main_original.c", "mv {0}/main.cpp {0}/main_original.cpp", "cp ../main.c {0}/", "cp ../main.cpp {0}/"]
        process = self.run_bash("\n".join(commands).format(project_name))
        
        if return_project:
            return Project(project_name)
        else: 
            return process
    
    def create_newproject(self,project_name,return_project=False):
        print("Creating new project...")
        if not is_valid_project_name(project_name): raise RuntimeError(f"'{project_name}' is not valid project name.")
        if self.project_exists(project_name): raise RuntimeError("project '{}' exists already.".format(project_name))
            
        commands = [ "./newProject {0}" ]
        process = self.run_bash("\n".join(commands).format(project_name))
        process = self.load_modified_main(project_name)
         
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
        if ("/" in project_name) or () : raise RuntimeError(f"{project_name} is not valid project name.")  # "rm" is dangerous!!! We shold check project name is properly passed. 
        if not self.project_exists(project_name): raise RuntimeError(f"project '{project_name}' does not exists.")
        
        project = Project(project_name)
        if not project.is_user_defined_project: raise RuntimeError(f"{project_name} is not created by users. remove_project cannot remove default projects.")
        command = f"rm -r {project_name}"
        process = self.run_bash(command)
        return process
    
#    def load_mdl(self,project_name,mdl_paths):
#        if self.project_exists(project_name): raise RuntimeError("project '{}' exists already.".format(project_name))
#        for path in mdl_paths:
#            if path[-4:] != ".mdl": raise RuntimeError("{} is not .mdl file".format(path))
#            commands = [ "cd {0}","cp main=main.cpp"]
    
    
    
class Project:
    
    
    def __init__(self,project_name):
        if not PyMicrOmegas().project_exists(project_name): raise RuntimeError("Project {} does not exist yet. Create it by PyMicrOmegas.create_newproject.".format(project_name))
        self.project_name = project_name
        self.path = dir_micromegas + project_name + "/"
        self.models_path = self.path + "work/models/"
        self._native = None
    
    
    def run_bash(self,command,shell=True,stdout=subprocess.PIPE,encoding="UTF-8",check=False,input=None,verbose=True):
        '''
        run './main' and return the output string.
        '''
        return run_bash(command,shell=shell,stdout=stdout,encoding=encoding,check=check,input=input,cwd=self.path,verbose=verbose)
    
    
    def load_mdl_files(self,mdl_paths):
        print("Loading .mdl files...")
        if type(mdl_paths) not in [list, tuple]: raise RuntimeError("input arguments must be list or tuple.")
        if len(mdl_paths)==0: raise RuntimeError("input argument is empty list or tuple.")
        processes = []
        for mdl_path in mdl_paths:
            if not os.path.isfile(mdl_path): raise RuntimeError("{} does not existing file".format(mdl_path))
            abs_mdl_path = to_abspath(mdl_path)
            print(f"Loading {abs_mdl_path}...")
            process = self.run_bash("cp {} {}".format(abs_mdl_path,self.models_path))
            processes.append(process)
        return processes
    
        
    def compile(self,main="main.c"):
        #process = subprocess.run("bash install_project.sh " + project_name,shell=True,stdout=subprocess.PIPE)
        return self.run_bash("make main={}".format(main))


    def native(self,build=True,force=False):
        """
        Return the ctypes bridge for this project.

        The bridge is generated inside the project directory and links against
        the already-installed micrOMEGAs libraries, so upstream C sources are
        not patched.
        """
        if self._native is None:
            self._native = NativeMicrOmegas(self.path,runner=self.run_bash)
        if build:
            self._native.build(force=force)
        return self._native


    def compile_native(self,force=False):
        return self.native(build=True,force=force)
           
      
    def clean(self):
        return self.run_bash("make clean")
    
    
    @property
    def vars(self):
        vars = np.loadtxt(self.models_path+"vars1.mdl",skiprows=3,delimiter="|",dtype=str,comments="=")  # .mdl file sometimes contain "=======..." line
        return vars
    
    @property
    def is_user_defined_project(self):
        return os.path.exists(self.path+"main_original.c")
    
    
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
        
        args = (
            ["./main", str(int_flags), str(n_inputvals), dof_fname]
            + list(map(str,get_keys(dict_parameters)))
            + list(map(str,get_values(dict_parameters)))
        )

        return self.run_bash(args,shell=False,verbose=False)
    
    
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
            key_vals = dict([floatnize(key_val.split("=")) for key_val in lines[ind+1].split(" ")])
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
        flag_names = tuple(flags or ())
        # Relic-density-only calls can use the ctypes bridge; other flag
        # combinations still use the subprocess path to preserve stdout parsing.
        if flag_names == NATIVE_OMEGA_FLAGS:
            try:
                return_dict = self.native().dark_omega(dict_parameters,dof_fname=dof_fname)
                return {"Xf":return_dict["Xf"],"Omega":return_dict["Omega"]}
            except (OSError, RuntimeError, subprocess.SubprocessError) as exc:
                warnings.warn(f"Native micrOMEGAs call failed; falling back to subprocess: {exc}",RuntimeWarning)
        output = self.run(dict_parameters,flags,dof_fname).stdout
        return self.parse_omega(output,flags)
    
    def calc_omega(self,dict_parameters,dof_fname=None):
        flags = ["OMEGA"]
        output = self.run(dict_parameters,flags,dof_fname=dof_fname).stdout
        return self.parse_omega(output,flags,with_channels=True)    
