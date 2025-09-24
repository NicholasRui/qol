import qol.tools.formatter as formatter

import os
import shutil

class MesaPythonScript:
    """
    Keeps track of information related to Python scripts inserted
    within MESA workflow
    """
    def __init__(self, name, template, const_args=[], prereqs=[], products=[], data_path='data/'):
        """
        name: an informative name
        template: .py file to copy

        relative path (rel_path) will be "script_{name}.py"

        arrange arguments of script to take first all of the constant args,
        then all of the prereqs, then the products after
        fed into script using sys package
        """
        self.name = name
        self.data_path = data_path
        
        self.rel_path = f'script_{name}.py'
        self.template = template
    
        self.const_args = const_args # constant arguments which do not depend on output of other runs
        self.data_prereqs = prereqs # model files and other things which are required for this to work
        self.data_products = products # model files and other things which are saved by this python script

        self.args_in = prereqs + products + const_args # arguments in ./rn

        self.LOGS_dir = None
        self.photos_dir = None
        
    def rn_string(self):
        return f'ipython tasks/{self.rel_path} ' + ' '.join([formatter.to_fortran(arg_in) for arg_in in self.args_in])

    re_string = rn_string # restarting this is the same as running

    def save(self, run_path):
        abs_path = os.path.join(run_path, 'tasks/', self.rel_path)
        shutil.copy(self.template, abs_path)
