import numpy as np

from qol.athena.table.AthenaTable import AthenaTable

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
        self.athinput_args = {}

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










