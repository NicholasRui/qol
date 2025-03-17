import qol.tools.formatter as formatter

import shutil

class MesaPythonScript:
    """
    Keeps track of information related to Python scripts inserted
    within MESA workflow
    """
    def __init__(self, rel_path, template, const_args=[], prereqs=[], products=[]):
        """
        template: .py file to copy
        prereqs: list of input files
        products: list of output files

        arrange arguments of script to take first all of the constant args,
        then all of the prereqs, then the products after
        fed into script using sys package
        """
        self.rel_path = rel_path
        self.template = template
    
        self.const_args = const_args # constant arguments which do not depend on output of other runs
        self.prereqs = prereqs # model files and other things which are required for this to work
        self.products = products # model files and other things which are saved by this inlist

        self.args_in = const_args + prereqs + products

    def rn_string(self):
        return f'ipython {self.rel_path} ' + ' '.join([formatter.to_fortran(arg_in) for arg_in in self.args_in])

    def save(self, run_path):
        abs_path = f'{run_path}/{self.rel_path}'
        shutil.copy(self.template, abs_path)