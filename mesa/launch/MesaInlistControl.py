import qol.helper.formatter as formatter

class MesaInlistControl:
    """
    Simple class which stores a single MESA control
    """
    def __init__(self, namelist, control, value, category=None, comment=None):
        """
        namelist: e.g., star_job, controls, etc.
        control: name of inlist option
        value: value of option
        category: optional -- allows grouping of options together under comment string given by this input
        comment: optional -- add trailing comment
        """
        self.namelist = namelist
        self.control = control
        self.value = value
        if category is None:
            self.category = ''
        else:    
            self.category = category
        self.comment = comment

    def inlist_string(self):
        """
        Return the string which should appear in inlist
        """
        fortran_value = formatter.to_fortran(self.value)
        
        if self.comment is not None:
            comment_str = f' ! {self.comment}'
        else:
            comment_str = ''

        return f'{self.control} = {fortran_value}{comment_str}'