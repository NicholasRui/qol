# Object which automatically writes bash scripts for submitting jobs
from qol.bash.BashScript import BashScript

class SlurmBashScript(BashScript):
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
        super().__init__()

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

        # Add tasks to top of file for slurm purposes
        # TODO: rewrite this so comments are aligned with each other
        if self.time is not None:
            self.add_task(f'#SBATCH --time={self.time}     # walltime')
        if self.ntasks is not None:
            self.add_task(f'#SBATCH --ntasks={self.ntasks}     # number of processor cores (i.e. tasks)')
        if self.nodes is not None:
            self.add_task(f'#SBATCH --nodes={self.nodes}     # number of nodes')
        if self.ntasks_per_node is not None:
            self.add_task(f'#SBATCH --ntasks-per-node={self.ntasks_per_node}     # number of tasks per node')
        if self.mem_per_cpu is not None:
            self.add_task(f'#SBATCH --mem-per-cpu={self.mem_per_cpu}     # memory per CPU core')
        if self.job_name is not None:
            self.add_task(f'#SBATCH -J "{self.job_name}"     # job name')
        if self.output is not None:
            self.add_task(f'#SBATCH --output={self.output}')
        if self.error is not None:
            self.add_task(f'#SBATCH --error={self.error}')
        if self.mail_user is not None:
            self.add_task(f'#SBATCH --mail-user={self.mail_user}     # email address')
        if self.mail_type is not None and self.mail_user is not None:
            self.add_task(f'#SBATCH --mail-type={self.mail_type}     # when to email')
        self.add_task('')
