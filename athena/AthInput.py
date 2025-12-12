"""
Class which helps initialize an athinput file (input file for Athena++).
"""
import numpy as np

import warnings

import qol.athena.alias as alias


class AthArg:
    """
    Helper object which stores individual arguments, in the following way:

    <block>
    name = value   # comment
    """
    def __init__(self, block, name, value, comment=None):
        self.block = block
        self.name = name
        self.value = value
        self.comment = comment

        # also store string version of value, and write it in scientific notation if float
        self.value_str = f'{value:e}' if isinstance(value, float) else str(value)

class AthInput:
    """
    TODO: allow this to read in an athinput file
    TODO: bind [] indexing to argument
    """
    def __init__(self):
        self.athargs = []
        self.header = ''

    def add_arg(self, block, name, value, comment=None):
        """
        Add an argument
        """
        self.athargs.append(AthArg(block=block, name=name, value=value, comment=comment))
    
    # Aliases to call add_arg for common athinput blocks
    add_arg_comment = alias.add_arg_comment
    add_arg_job = alias.add_arg_job
    add_arg_mesh = alias.add_arg_mesh
    add_arg_time = alias.add_arg_time
    add_arg_problem = alias.add_arg_problem
    add_arg_hydro = alias.add_arg_hydro
    add_arg_radiation = alias.add_arg_radiation
    add_arg_output1 = alias.add_arg_output1
    add_arg_output2 = alias.add_arg_output2

    def add_header_line(self, line):
        """
        Add a comment line at the top of the document
        """
        assert '\n' not in line
        self.header += f'# {line}\n'

    def get_blocks(self):
        blocks = list(set([arg.block for arg in self.athargs]))

        # Sort it in the following way:
        #  1. elements in some named list
        #  2. output blocks
        #  3. everything else
        first_blocks = ['comment', 'job', 'mesh', 'time', 'problem']
        order_map = {block: ii for ii, block in enumerate(first_blocks)}

        sort_func = lambda block: (0, order_map[block]) if block in first_blocks else (1, block) if block.startswith('output') else (2, block)
        blocks = sorted(blocks, key=sort_func)

        return blocks


    def save(self, fname):
        """
        Save AthInput file at filename
        """
        blocks = self.get_blocks()

        outstr = ''
        for block in blocks:
            outstr += f'<{block}>\n'

            # get arguments within block
            # sort_func = lambda atharg: atharg.name
            # athargs_block = sorted([atharg for atharg in self.athargs if atharg.block == block], key=sort_func)
            athargs_block = [atharg for atharg in self.athargs if atharg.block == block]
            atharg_names = [atharg.name for atharg in athargs_block]
            atharg_comments = [atharg.comment for atharg in athargs_block]

            # add each argument to outstr, making sure equal and comment signs are aligned with each other
            max_len_name = max([len(name) for name in atharg_names])
            atharg_strs = [f'{atharg.name.ljust(max_len_name)} = {atharg.value_str}' for atharg in athargs_block]
            max_len_atharg_str = max([len(atharg_str) for atharg_str in atharg_strs])
            atharg_strs = [f'{atharg_str.ljust(max_len_atharg_str)}   # {atharg_comments[ii]}\n' if atharg_comments[ii] is not None else f'{atharg_str}\n' for ii, atharg_str in enumerate(atharg_strs)]

            outstr += ''.join(atharg_strs) + '\n'

        with open(fname, 'w') as f:
            f.write(self.header)
            f.write('\n')
            f.write(outstr)
