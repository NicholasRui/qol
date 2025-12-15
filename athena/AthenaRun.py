"""
Defines AthenaRun class for scripting AthenaRun files
"""
from qol.athena.AthInput import AthInput
import qol.config as config
from qol.bash.BashScript import BashScript
from qol.bash.SlurmBashScript import SlurmBashScript

import os, shutil

class AthenaRun:
    """
    """
    def __init__(self, run_path, pgen_path, athinput, compile_flags, compile_log_fname=None):
        """
        Docstring for __init__
        
        :param self: Description
        :param run_path: Path to save Athena run
        :param pgen_path: Path from which to copy pgen file (suffix: .cpp)
        :param athinput: AthInput object storing athinput info
        :param compile_flags: list of strings indicating compilation flags to include
        :param compile_log_fname: if specified, copy this configure.log file to directory
        """
        # TODO: allow a configure.log to get passed to avoid multiple runs recompiling on first run each time

        if os.path.exists(run_path):
            raise ValueError(f'Path already exists: {run_path}')
        assert os.path.exists(pgen_path)
        assert pgen_path[-4:] == '.cpp'
        assert isinstance(athinput, AthInput)
        assert isinstance(compile_flags, list)
        for compile_flag in compile_flags:
            assert len(compile_flag) > 1
            assert compile_flag[0] == '-' # assert first character of each compile flag is '-'

        self.run_path = run_path
        self.pgen_path = pgen_path
        self.athinput = athinput
        self.compile_flags = compile_flags
        self.compile_log_fname = compile_log_fname

    def save_directory(self, grant_perms=False,
                       slurm_job_name=None, slurm_job_time=config.slurm_job_time_default,
                       slurm_job_ntasks=config.slurm_job_ntasks_default, slurm_job_nodes=config.slurm_job_nodes_default,
                       slurm_job_ntasks_per_node=config.slurm_job_ntasks_per_node_default,
                       slurm_job_mem_per_cpu=config.slurm_job_mem_per_cpu_default,
                       slurm_job_email_user=True, OMP_NUM_THREADS=config.mesa_OMP_NUM_THREADS,
                       data_path='data/'):
        run_path = self.run_path
        pgen_path = self.pgen_path
        athinput = self.athinput
        compile_flags = self.compile_flags
        compile_log_fname = self.compile_log_fname

        # Create directories
        if os.path.exists(run_path):
            raise ValueError(f'Path already exists: {run_path}')
        os.mkdir(run_path) # make run path
        os.mkdir(os.path.join(run_path, 'output')) # make output path

        # Copy over pgen and write athinput file
        pgen_path_basename = os.path.basename(pgen_path)
        shutil.copy(pgen_path, os.path.join(run_path, pgen_path_basename))
        athinput.save(os.path.join(run_path, 'athinput')) # write athinput

        # Write bash shell for automatically launching job
        bash_script = BashScript()

        bash_script.add_task('# clean and compile (only do if need to compile with new settings)')
        bash_script.add_task('if ! cmp -s configure.log $ATHENA_DIR/configure_hydro.log || ! [ -f $ATHENA_DIR/bin/athena_hydro ] || ! cmp -s full_mesa_blast.cpp $ATHENA_DIR/src/pgen/full_mesa_blast.cpp; then')
        bash_script.add_task('  # copy pgen file')
        bash_script.add_task(f'  cp {pgen_path_basename} $ATHENA_DIR/src/pgen')
        bash_script.add_task('  ')
        bash_script.add_task('  # configure and compile')
        bash_script.add_task('  cd $ATHENA_DIR')
        pgen_name = pgen_path_basename.replace('.cpp', '')
        bash_script.add_task(f"  python configure.py --prob={pgen_name} {' '.join(compile_flags)} # configure") # TODO FLAGS
        bash_script.add_task('  make clean')
        bash_script.add_task('  make -j')
        bash_script.add_task(f"  cp configure.log {os.path.join(run_path, 'configure.log')}")
        bash_script.add_task('  cd -')
        bash_script.add_task('fi')
        bash_script.add_task('')
        bash_script.add_task('# run')
        bash_script.add_task(f"cp bin/athena {os.path.join(run_path, 'athena')}")
        bash_script.add_task('./athena -i athinput')

        bash_script.save(os.path.join(run_path, 'run_athena.sh'))

        # Write slurm job shell
        if slurm_job_name is not None:
            # only email user if say so; get email from config.py file
            if slurm_job_email_user:
                mail_user = config.slurm_job_mail_user
                mail_type = 'BEGIN,FAIL,END'
            else:
                mail_user = mail_type = None
            
            # script to start job
            slurm_bash_script = SlurmBashScript(
                 job_name=slurm_job_name,
                 time=slurm_job_time,
                 ntasks=slurm_job_ntasks, nodes=slurm_job_nodes,
                 ntasks_per_node=slurm_job_ntasks_per_node,
                 mem_per_cpu=slurm_job_mem_per_cpu,
                 output=os.path.join(run_path, 'output.out'),
                 error=os.path.join(run_path, 'error.out'), # absolute paths
                 mail_user=mail_user, # email address
                 mail_type=mail_type # conditions for emailing
                 )
            
            slurm_bash_script.add_task(os.path.join(run_path, 'run_athena.sh'))
            slurm_bash_script.save(os.path.join(run_path, 'slurm_run_athena.sh'))
        
        # if specified, copy given configure.log file over
        # This is to avoid multiple recompilations for multiple runs which have the same
        # compilation settings as each other.
        shutil.copy(compile_log_fname, os.path.join(run_path, 'configure.log'))

        # grant permissions if needed
        if grant_perms:
            os.chmod(os.path.join(run_path, 'run_athena.sh'), 0o755)
            if slurm_job_name is not None:
                os.chmod(os.path.join(run_path, 'slurm_run_athena.sh'), 0o755)
