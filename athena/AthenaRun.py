import numpy as np

#from qol.athena.table.AthenaTable import AthenaTable
from qol.athena.read import read_tab

import glob
import os

# this is a continuing work in progress -- will continue to generalize this as it fits my use cases


class AthenaRun:
    """
    Holds information about a given Athena++ run.
    """
    def __init__(self, athinput_fname, run_path='.'):
        """
        athinput_fname: path to the athinput file used
        """
        self.athinput_fname = athinput_fname
        self.run_path = run_path
        self.athinput_args = {}

        self.athoutput_cache = {} # keys are (output_number, block_number, output_index)

        self.parse_athinput()

    def parse_athinput(self):
        """
        Read the athinput file and store all the arguments. An input like:
            <output1>
            file_type  = tab
        will be stored under the key 'output1.file_type'.

        Intended to be used as a helper for __init__.
        """
        # read the file
        with open(self.athinput_fname, 'r') as f:
            lines = f.readlines()

        # loop through all lines
        current_sect = ''

        for line in lines:
            # remove comments
            line = line.split('#')[0]

            # skip empty lines
            if len(line.split()) == 0:
                continue

            # recognize if the line is indicating a section
            if line.startswith('<') and line.endswith('>\n'):
                current_sect = line.strip('<>\n')
                continue

            # recognize if line has an argument in it
            if '=' in line:
                argname, argval = line.split('=', 1)
                argname, argval = argname.strip(), argval.strip() # get rid of leading/trailing whitespace

                # try to cast argval into appropriate type
                try:
                    argval = int(argval)
                except:
                    try:
                        argval = float(argval)
                    except:
                        pass
                
                self.athinput_args[f'{current_sect}.{argname}'] = argval

                continue

            # if you got here, throw an error because we have failed to interpret the line
            raise RuntimeError(f'Failed to interpret following line:\n  > {line}')
    
    def get_output_fnames(self, output_number, block_number=0):
        """
        return list of fnames within a given numbered output

        output_number: The number of the output, i.e., <output3> in athinput has output_number=3
        block_number: block number
        """
        if f'output{output_number}.file_type' not in self.athinput_args.keys():
            raise ValueError(f'Output {output_number} does not exist. (not specified in athinput file)')
        if self.athinput_args[f'output{output_number}.file_type'] != 'tab':
            raise NotImplementedError(f"Currently file_type must be 'tab': {self.athinput_args[f'output{output_number}.file_type']} not supported.")

        output_pattern = f"{self.athinput_args['job.problem_id']}.block{block_number}.out{output_number}.*.tab"
        output_pattern = os.path.join(self.run_path, output_pattern)

        output_fnames = glob.glob(output_pattern)
        output_fnames.sort()
        return output_fnames

    def get_output_indices(self, output_number, block_number=0):
        """
        return list of indices within a given numbered output

        output_number: The number of the output, i.e., <output3> in athinput has output_number=3
        block_number: block number
        """
        output_fnames = self.get_output_fnames(output_number, block_number)

        output_indices = [int(output_fname.split('.tab')[0].split('.')[-1]) for output_fname in output_fnames]

        return output_indices

    def get_output_tab(self, output_number, output_index, block_number=0, cache=True):
        """
        Retrieve AthenaTable version of output tab

        output_number: The number of the output, i.e., <output3> in athinput has output_number=3
        block_number: block number

        # TODO: make it so either output_number or output_index can be array-like, and return a list of tables
        """
        if f'output{output_number}.file_type' not in self.athinput_args.keys():
            raise ValueError(f'Output {output_number} does not exist. (not specified in athinput file)')
        if self.athinput_args[f'output{output_number}.file_type'] != 'tab':
            raise NotImplementedError(f"Currently file_type must be 'tab': {self.athinput_args[f'output{output_number}.file_type']} not supported.")
        
        # check if it's in cache and return if so
        if (output_number, block_number, output_index) in self.athoutput_cache.keys():
            return self.athoutput_cache[(output_number, block_number, output_index)]

        # retrieve table
        tab = read_tab(f"{self.athinput_args['job.problem_id']}.block{block_number}.out{output_number}.{str(output_index).zfill(5)}.tab")

        # cache if this is desired
        if cache:
            self.athoutput_cache[(output_number, block_number, output_index)] = tab
        
        return tab
