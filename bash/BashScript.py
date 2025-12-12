"""
Object for writing bash scripts
"""
class BashScript:
    """
    Stores information for bash script
    """
    def __init__(self):
        pass

        self.tasks = []

    def add_task(self, task):
        """
        Add unix line indicating what to do
        """
        if task == '': # if empty line, just add a line break
            self.tasks.append('\n')
            return

        if task[-1] != '\n': # if task doesn't end in line break, add one
            task += '\n'
        
        self.tasks.append(task) # add task

    def save(self, abs_path):
        text = '#!/bin/bash\n\n'

        # add tasks
        text += ''.join(self.tasks)

        with open(abs_path, 'w') as f:
            f.write(text)
