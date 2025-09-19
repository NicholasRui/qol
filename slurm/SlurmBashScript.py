# Object which automatically writes bash scripts for submitting jobs


class SlurmBashScript:
    """
    Stores information for slurm bash script
    """
    def __init__(self,
                 job_name='job',
                 time='7-00:00:00', # specify walltime as a string for now
                 ntasks=1, nodes=1,
                 ntasks_per_node=1,
                 mem_per_cpu='10G', # specify as string for now
                 output='output.out', error='error.out', # absolute paths
                 mail_user=None, # email address
                 mail_type='BEGIN,FAIL,END' # conditions for emailing
                 ):
        """
        """
        self.job_name = job_name
        self.time = time
        self.ntasks = ntasks
        self.nodes = nodes
        self.ntasks_per_node = ntasks_per_node
        self.mem_per_cpu = mem_per_cpu
        self.output = output
        self.error = error
        self.mail_user = mail_user
        self.mail_type = mail_type

        self.tasks = []

    def save(self, abs_path):
        """
        save bash file at some absolute path abs_path
        """
        text = ''
        text += '#!/bin/bash\n'

        # TODO abstract this a bit
        text += f'#SBATCH --time={self.time}     # walltime\n' if self.time is not None else ''
        text += f'#SBATCH --ntasks={self.ntasks}     # number of processor cores (i.e. tasks)\n' if self.ntasks is not None else ''
        text += f'#SBATCH --nodes={self.nodes}     # number of nodes\n' if self.nodes is not None else ''
        text += f'#SBATCH --ntasks-per-node={self.ntasks_per_node}     # number of tasks per node\n' if self.ntasks_per_node is not None else ''
        text += f'#SBATCH --mem-per-cpu={self.mem_per_cpu}     # memory per CPU core\n' if self.mem_per_cpu is not None else ''
        text += f'#SBATCH -J "{self.job_name}"     # job name\n' if self.job_name is not None else ''
        text += f'#SBATCH --output={self.output}\n' if self.output is not None else ''
        text += f'#SBATCH --error={self.error}\n' if self.error is not None else ''
        text += f'#SBATCH --mail-user={self.mail_user}     # email address\n' if self.mail_user is not None else ''
        text += f'#SBATCH --mail-type={self.mail_type}     # when to email\n' if self.mail_type is not None and self.mail_user is not None else '' # only populate if mail_user also specified
        text += '\n'

        # add tasks
        text += ''.join(self.tasks)

        with open(abs_path, 'w') as f:
            f.write(text)

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

    def set_job_name(self, job_name):
        self.job_name = job_name
    
    def set_time(self, time):
        self.time = time
    
    def set_ntasks(self, ntasks):
        self.ntasks = ntasks
    
    def set_nodes(self, nodes):
        self.nodes = nodes
    
    def set_mem_per_cpu(self, mem_per_cpu):
        self.mem_per_cpu = mem_per_cpu
    
    def set_output(self, output):
        self.output = output
    
    def set_error(self, error):
        self.error = error
    
    def set_mail_user(self, mail_user):
        self.mail_user = mail_user
    
    def set_mail_type(self, mail_type):
        self.mail_type = mail_type






