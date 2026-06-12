import numpy as np

#from qol.athena.table.AthenaTable import AthenaTable
from qol.athena.read import read_tab, read_vtk
import qol.athena.const as athconst

import glob
import os

# NOTE: this is a continuing work in progress -- will continue to generalize this as it fits my use cases


class AthenaRunOutput:
    """
    Holds information about a given Athena++ run.
    """
    def __init__(self, athinput_fname, run_path='.'): # TODO: connect this to AthInput and don't use this notation
        """
        athinput_fname: relative path to the athinput file used
        """
        self.athinput_fname = athinput_fname
        self.run_path = run_path
        self.athinput_args = athinput_args = {}

        self.athoutput_cache = {} # keys are (output_number, block_number, output_index)

        self.parse_athinput()

        ###############################################
        ## DEFINE CONVENIENT ATTRIBUTES IF AVAILABLE ##
        ###############################################

        # if T0, rho0, L0, mu0 are specified for implicit radiation units
        self.has_radiation_units = False
        if set(['radiation.T_unit', 'radiation.density_unit', 'radiation.length_unit', 'radiation.molecular_weight']).issubset(athinput_args.keys()):
            self.T0 = T0 = athinput_args['radiation.T_unit']
            self.rho0 = rho0 = athinput_args['radiation.density_unit']
            self.L0 = L0 = athinput_args['radiation.length_unit']
            self.mu0 = mu0 = athinput_args['radiation.molecular_weight']

            # derived units
            arad = athconst.radiation_aconst_cgs
            clight = 2.99792458e10

            self.Rgas = Rgas = 8.314462618e7 / mu0
            self.v0 = v0 = np.sqrt(Rgas * T0)
            self.press0 = press0 = rho0 * Rgas * T0
            self.t0 = t0 = L0 / v0

            self.E0 = E0 = arad * T0 ** 4
            self.F0 = F0 = E0 * clight

            self.has_radiation_units = True








    def parse_athinput(self):
        """
        Read the athinput file and store all the arguments. An input like:
            <output1>
            file_type  = tab
        will be stored under the key 'output1.file_type'.

        Intended to be used as a helper for __init__.
        """
        # read the file
        with open(os.path.join(self.run_path, self.athinput_fname), 'r') as f:
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
        
        ext = self.athinput_args[f'output{output_number}.file_type']
        assert ext in ['tab', 'vtk']
        
        output_pattern = f"{self.athinput_args['job.problem_id']}.block{block_number}.out{output_number}.*.{ext}"
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

        ext = self.athinput_args[f'output{output_number}.file_type']
        assert ext in ['tab', 'vtk']
        
        output_indices = [int(output_fname.split(f'.{ext}')[0].split('.')[-1]) for output_fname in output_fnames]

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
        
        ext = self.athinput_args[f'output{output_number}.file_type']
        assert ext in ['tab', 'vtk']
        
        # check if it's in cache and return if so
        if (output_number, block_number, output_index) in self.athoutput_cache.keys():
            return self.athoutput_cache[(output_number, block_number, output_index)]

        # retrieve table
        tab_fname = os.path.join(self.run_path, f"{self.athinput_args['job.problem_id']}.block{block_number}.out{output_number}.{str(output_index).zfill(5)}.{ext}")
        if ext == 'tab':
            tab = read_tab(tab_fname)
        else:
            tab = read_vtk(tab_fname)

        # cache if this is desired
        if cache:
            self.athoutput_cache[(output_number, block_number, output_index)] = tab
        
        return tab

    def get_output_time(self, output_number, output_index, block_number=0):
        """
        read 'time' from first line of output file
        """
        if f'output{output_number}.file_type' not in self.athinput_args.keys():
            raise ValueError(f'Output {output_number} does not exist. (not specified in athinput file)')
        if self.athinput_args[f'output{output_number}.file_type'] != 'tab':
            raise NotImplementedError(f"Currently file_type must be 'tab': {self.athinput_args[f'output{output_number}.file_type']} not supported.")
        
        tab_fname = os.path.join(self.run_path, f"{self.athinput_args['job.problem_id']}.block{block_number}.out{output_number}.{str(output_index).zfill(5)}.tab")
        with open(tab_fname, 'r') as f:
            line = f.readline()
        
        time = float(line.split('time=')[1].split('cycle=')[0].replace(' ', ''))
        return time

    def get_output_cycle(self, output_number, output_index, block_number=0):
        """
        read 'cycle' from first line of output file
        """
        if f'output{output_number}.file_type' not in self.athinput_args.keys():
            raise ValueError(f'Output {output_number} does not exist. (not specified in athinput file)')
        if self.athinput_args[f'output{output_number}.file_type'] != 'tab':
            raise NotImplementedError(f"Currently file_type must be 'tab': {self.athinput_args[f'output{output_number}.file_type']} not supported.")
        
        tab_fname = os.path.join(self.run_path, f"{self.athinput_args['job.problem_id']}.block{block_number}.out{output_number}.{str(output_index).zfill(5)}.tab")
        with open(tab_fname, 'r') as f:
            line = f.readline()
        
        cycle = int(line.split('cycle=')[1].split('variables=')[0].replace(' ', ''))
        return cycle
    
    def get_all_output_tabs(self, output_number, block_number=0, ext='tab'):
        iis = self.get_output_indices(output_number=output_number, block_number=block_number, ext=ext)
        tabs = [self.get_output_tab(output_number=output_number, output_index=ii, block_number=block_number) for ii in iis]
        return tabs

    def get_all_output_times(self, output_number, block_number=0):
        iis = self.get_output_indices(output_number=output_number, block_number=block_number, ext='tab')
        times = np.array([self.get_output_time(output_number=output_number, output_index=ii, block_number=block_number) for ii in iis])
        return times

    def get_all_output_cycles(self, output_number, block_number=0):
        iis = self.get_output_indices(output_number=output_number, block_number=block_number, ext='tab')
        cycles = np.array([self.get_output_cycle(output_number=output_number, output_index=ii, block_number=block_number) for ii in iis])
        return cycles






