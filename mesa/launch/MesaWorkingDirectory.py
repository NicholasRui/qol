import qol.config as config
import qol.info as info
import qol.tools.formatter as formatter
from qol.slurm.SlurmBashScript import SlurmBashScript

import numpy as np

import matplotlib.pyplot as plt
import networkx as nx

import os
import shutil

from itertools import chain

class MesaWorkingDirectory:
    """
    Stores information for creating a custom MESA work directory

    Records tasks, which must have the following attributes / methods:
    - rn_string : returns string to put into rn to run this task
    - save : returns command which saves the needed files
    - prereqs : list of prereq files inputed to task
    - products : list of product files outputed by task
    """
    def __init__(self, run_path, mesa_version=config.mesa_version):
        if os.path.exists(run_path):
            raise ValueError(f'Path already exists: {run_path}')

        self.run_path = run_path
        self.mesa_version = mesa_version

        self.history_columns_path = None
        self.profile_columns_path = None

        self.copy_from_path_root_prereqs = []
        self.rel_paths_root_prereq = []
        self.tasks = []

        self.submit_job_path = f'{run_path}/submit_job.sh'
        self.restart_job_path = f'{run_path}/restart_job.sh'

    def add_root_prereq(self, copy_from_abs_path, rel_path):
        """
        Add a model file which serves as a "root prereq"
        (something used by an inlist which is there from the beginning,
                    not outputted as an intermediate product by a task)
        
        Copy from copy_from_abs_path (absolute path) to rel_path (relative path in work/data/ subdirectory)
        """
        if not os.path.exists(copy_from_abs_path):
            raise ValueError(f'Path not found: {copy_from_abs_path}')

        self.copy_from_path_root_prereqs.append(copy_from_abs_path)
        self.rel_paths_root_prereq.append(rel_path)
    
    def add_task(self, task):
        """
        Add a task to the task list, after determining the prereqs already exist
        """
        # Check no conflicts
        self.check_needed_prereqs(task)
        self.check_fname_conflicts(task)

        # Store tasks
        self.tasks.append(task)

    def load_qol_pgstar(self):
        """
        load default qol inlist_pgstar, need to do this before doing use_qol_pgstar
        """
        self.add_root_prereq(copy_from_abs_path=f'{info.qol_path}mesa/resources/r24.08.1/inlist_pgstar', rel_path='inlist_pgstar')

    def check_needed_prereqs(self, task):
        """
        Check if, with addition of new task, prereqs exist, and throw error if not
        Otherwise, do nothing
        """
        task_prereqs = task.data_prereqs

        existing_products = []
        for task in self.tasks:
            existing_products += task.data_products
        existing_products += self.rel_paths_root_prereq

        if not np.in1d(task_prereqs, existing_products).all():
            raise ValueError('task requests prereqs which do not exist')
        
        return
    
    def check_fname_conflicts(self, task):
        """
        make sure the new task doesn't present an fname conflict
        """
        existing_task_names = []
        for task in self.tasks:
            existing_task_names += task.name

        if task.name in existing_task_names:
            raise ValueError('name of task conflicts with existing name')
        
    def copy_history_columns_list(self, abs_path):
        """
        Copy history_columns.list file from abs_path

        None means try to find a preset one defined in qol.
        """
        if os.path.exists(abs_path):
            self.history_columns_path = abs_path
        else:
            raise ValueError(f'No path found: {abs_path}')
        
    def copy_profile_columns_list(self, abs_path=None):
        """
        Copy history_columns.list file from abs_path

        None means try to find a preset one defined in qol.
        """
        if os.path.exists(abs_path):
            self.profile_columns_path = abs_path
        else:
            raise ValueError(f'No path found: {abs_path}')

    def save_directory(self, grant_perms=False, make_flowchart=True, source_sdk=False,
                       slurm_job_name=None, slurm_job_time=config.slurm_job_time_default,
                       slurm_job_ntasks=config.slurm_job_ntasks_default, slurm_job_nodes=config.slurm_job_nodes_default,
                       slurm_job_mem_per_cpu=config.slurm_job_mem_per_cpu_default,
                       slurm_job_email_user=True, OMP_NUM_THREADS=config.mesa_OMP_NUM_THREADS):
        """
        Create MESA directory

        grant_perms: if True, grants permissions for bash files that need to be run
        slurm_job_name: if not None, saves a bash script for sending this job with this arg as job_name
        source_sdk: if True, have submit/restart files do 'source $MESASDK_ROOT/bin/mesasdk_init.sh'
        """
        # CREATE NEW MESA DIRECTORY
        run_path = self.run_path

        # copy important files from $MESA_DIR/star/work, remove files which we want to write explicitly
        if run_path[-1] != '/':
            run_path += '/'
        
        if os.path.exists(run_path):
            raise ValueError(f'Path already exists: {run_path}')
        os.mkdir(run_path)
        os.mkdir(f'{run_path}make')
        os.mkdir(f'{run_path}src')
        os.mkdir(f'{run_path}tasks') # store inlists and python scripts here
        os.mkdir(f'{run_path}data') # store prereqs and products here

        mesadir = config.mesa_paths[self.mesa_version]
        workdir = f'{mesadir}/star/work/'

        shutil.copy(f'{workdir}clean', f'{run_path}clean')
        shutil.copy(f'{workdir}make/makefile', f'{run_path}make/makefile')
        shutil.copy(f'{workdir}mk', f'{run_path}mk')

        # copy run_star_extras, but explicitly substitute standard_run_star_extras.inc for convenience
        if os.path.exists(f'{workdir}src/run.f90'):
            shutil.copy(f'{workdir}src/run.f90', f'{run_path}src/run.f90')
        elif os.path.exists(f'{workdir}src/run.f'):
            shutil.copy(f'{workdir}src/run.f', f'{run_path}src/run.f')
        else:
            raise ValueError('Neither src/run.f90 and src/run.f found')
        
        if os.path.exists(f'{workdir}src/run_star_extras.f90'):
            RSE_in_path = f'{workdir}src/run_star_extras.f90'
            RSE_out_path = f'{run_path}src/run_star_extras.f90'
        elif os.path.exists(f'{workdir}src/run_star_extras.f'):
            RSE_in_path = f'{workdir}src/run_star_extras.f'
            RSE_out_path = f'{run_path}src/run_star_extras.f'
        else:
            raise ValueError('Neither src/run_star_extras.f90 and src/run_star_extras.f found')
        
        newlines = []
        with open(RSE_in_path, 'r') as f:
            lines = f.readlines()
        with open(f'{mesadir}/include/standard_run_star_extras.inc', 'r') as f:
            standard_RSE_text = f.read()

        for line in lines:
            if "include 'standard_run_star_extras.inc'" in line:
                newlines.append(standard_RSE_text)
            else:
                newlines.append(line)
        
        RSE_out_text = ''.join(newlines)

        with open(RSE_out_path,'w') as f:
            f.write(RSE_out_text)
        
        # copy over history and profile columns files, if specified
        if self.history_columns_path is not None:
            shutil.copy(self.history_columns_path, f'{run_path}history_columns.list')
        if self.profile_columns_path is not None:
            shutil.copy(self.profile_columns_path, f'{run_path}profile_columns.list')

        # copy over helper files
        shutil.copy(f'{info.qol_path}mesa/resources/bash/do_one', f'{run_path}do_one')
        shutil.copy(f'{info.qol_path}mesa/resources/bash/re_one', f'{run_path}re_one')
        shutil.copy(f'{info.qol_path}mesa/resources/bash/any_missing', f'{run_path}any_missing')

        # create rn and re
        # rn will call any_missing and do things if False, and terminate if not
        # re will call any_missing and do things if True, and skip if true

        # helper lambda which writes an if statement based on check_exists
        check_missing_template = ''
        check_missing_template += './any_missing <<FILES>>\n'
        check_missing_template += 'missing=$?\n'
        check_missing_template += 'if [ $missing -eq 0 ]; then\n'
        check_missing_template += '    <<NONE_MISSING>>\n'
        check_missing_template += 'else\n'
        check_missing_template += '    <<SOME_MISSING>>\n'
        check_missing_template += 'fi\n\n'

        check_if_missing = lambda fname_list, if_none_missing, if_some_missing: \
          check_missing_template.replace('<<FILES>>', ' '.join([formatter.to_fortran(f'data/{fname}') for fname in fname_list])) \
                                .replace('<<NONE_MISSING>>', if_none_missing) \
                                .replace('<<SOME_MISSING>>', if_some_missing)

        ### Start rn text
        rn_text = ''
        rn_text += '#!/bin/bash\n\n'
        rn_text += '# rn script\n\n'

        rn_text += 'if [ ! -f "star" ]; then\n    echo "QoL Error: Did you forget to do ./mk?"\n    exit 1\n    fi\n\n' # check if ./mk was run
        
        # root prereqs
        rn_text += '# Check root prereqs\n'
        rn_text += check_if_missing(fname_list=self.rel_paths_root_prereq, \
                if_none_missing="echo 'QoL: All root prereqs found, continue!'", \
                if_some_missing="echo 'QoL: SOME ROOT PREREQS MISSING, EXIT'\n    exit 1")

        ### Start re text
        re_text = ''
        re_text += '#!/bin/bash\n\n'
        re_text += '# re script\n\n'
        re_text += 'started=0\n\n'

        re_text += 'if [ ! -f "star" ]; then\n    echo "QoL Error: Did you forget to do ./mk?"\n    exit 1\n    fi\n\n' # check if ./mk was run

        re_text += '# Check root prereqs\n'
        re_text += check_if_missing(fname_list=self.rel_paths_root_prereq, \
                if_none_missing="echo 'QoL: All root prereqs found, continue!'", \
                if_some_missing="echo 'QoL: SOME ROOT PREREQS MISSING, EXIT'\n    exit 1")

        # copy all root prereqs
        for ii, rel_path in enumerate(self.rel_paths_root_prereq):
            copy_from_abs_path = self.copy_from_path_root_prereqs[ii]
            shutil.copy(copy_from_abs_path, f'{self.run_path}/data/{rel_path}')

        # LOOP OVER TASKS AND RUN THEM IN ORDER
        vlevels = [] # vertical level in chart
        hlevels = [] # horizontal level in chart
        sorted_tasks = []

        existing_products = []
        existing_products += self.rel_paths_root_prereq

        level_products = []
        task_ids = [] # dummy index to skip tasks which have been handled
        vlevel, hlevel = 0, 0

        while True:
            new_tasks = 0

            for ii, task in enumerate(self.tasks):
                # if task visited before, skip it
                if ii in task_ids:
                    continue

                # if all prereqs are in existing_products exist,
                # run task, then add task to sorted_tasks and increment
                if np.in1d(task.data_prereqs, existing_products).all():
                    # Add strings to rn and re
                    full_rn_string = f'# Try to run tasks/{task.rel_path}\n'
                    full_re_string = f'# Try to run tasks/{task.rel_path}\n'

                    full_rn_string += check_if_missing(fname_list=task.data_prereqs, \
                            if_none_missing=task.rn_string(), \
                            if_some_missing=f"echo 'QoL: SOME PREREQS MISSING FOR tasks/{task.rel_path}, EXIT'\n    exit 1")
                    
                    # for re_string, add an if statement depending on 'started' variable
                    # if the re-run has started already but needed data is missing, that means something crashed and we shouldn't continue
                    re_string_with_check_started = 'if [ $started -eq 0 ]; then\n'
                    re_string_with_check_started += f'        {task.re_string()}\n'
                    re_string_with_check_started += '        started=1\n'
                    re_string_with_check_started += '    else\n'
                    re_string_with_check_started += f"        echo 'QoL: SOME PREREQS MISSING FOR tasks/{task.rel_path}, EXIT'\n        exit 1\n"
                    re_string_with_check_started += '    fi'
                    full_re_string += check_if_missing(fname_list=task.data_products, \
                            if_none_missing=f"echo 'QoL: All products found, skipping tasks/{task.rel_path}'", \
                            if_some_missing=re_string_with_check_started)
                    
                    rn_text += full_rn_string
                    re_text += full_re_string

                    # Save task
                    task.save(run_path=run_path)

                    task_ids.append(ii)
                    sorted_tasks.append(task)
                    level_products += task.data_products
                    hlevels.append(hlevel)
                    vlevels.append(vlevel)

                    hlevel += 1
                    new_tasks += 1
            
            # if go a full loop without getting any new tasks,
            # go to next vertical level
            if new_tasks == 0:
                existing_products += level_products
                level_products = []

                vlevel += 1
                hlevel = 0
            
            # if all tasks accounted for, break
            if len(sorted_tasks) == len(self.tasks):
                break

        ### End rn text and save rn file
        rn_text += "echo 'QoL: All tasks successfully completed!'\n"
        rn_text += 'exit 0\n'

        with open(f'{run_path}rn', 'w') as f:
            f.write(rn_text)
        
        ### End re text and save re file
        all_products = list(set(chain.from_iterable([task.data_products for task in self.tasks])))

        re_text += check_if_missing(fname_list=all_products, \
                            if_none_missing="echo 'QoL: All tasks successfully completed!'\n    exit 0", \
                            if_some_missing=f"echo 'QoL: Some products missing -- at least one task failed'\n    exit 1")

        with open(f'{run_path}re', 'w') as f:
            f.write(re_text)
        
        if grant_perms:
            os.chmod(f'{run_path}do_one', 0o755)
            os.chmod(f'{run_path}re_one', 0o755)
            os.chmod(f'{run_path}any_missing', 0o755)
            os.chmod(f'{run_path}rn', 0o755)
            os.chmod(f'{run_path}re', 0o755)

        # if desired, make flowchart which plots all of the tasks
        # and shows dependencies of different inputs and outputs
        if make_flowchart:
            # settings
            task_shape = 's'
            task_color = 'xkcd:carolina blue'
            product_shape = 'o'
            product_color = 'xkcd:golden yellow'

            node_size = 50000
            font_size = 20
            pad = 0.7 # how much to bad axis edges by
            scale = 3. # how much to scale up axis by

            # create directed graph where nodes are tasks or prereqs / products
            G = nx.Graph()

            # add nodes and edges
            node_names = []
            node_positions = {}
            node_shapes = {}
            node_colors = {}
            edges = []

            h, v = 0, -0.5
            h_shift = 1 # shift apart tasks and prereqs for visibility

            for rel_path in self.rel_paths_root_prereq: # add nodes for root prereqs
                node_names.append(rel_path)
                node_positions[rel_path] = (h+h_shift, -2*v) # minus sign to make earlier levels higher on plot, 2* to make room for intermediate products
                node_shapes[rel_path] = product_shape
                node_colors[rel_path] = product_color

                h += 1

            v_current, h_product = 0, 0
            for ii, task in enumerate(sorted_tasks): # add nodes and edges for tasks and their products
                h, v = hlevels[ii], vlevels[ii]

                task_node_name = f'({ii+1}) {task.rel_path}'

                node_names.append(task_node_name)
                node_positions[task_node_name] = (h-h_shift, -2*v)
                node_shapes[task_node_name] = task_shape
                node_colors[task_node_name] = task_color

                # for edge purposes, check if next vlevel achieved and update accordingly
                if v != v_current:
                    v_current = v
                    h_product = 0

                # add node for each product
                for product in task.data_products:
                    node_names.append(product)
                    node_positions[product] = (h_product+h_shift, -2*v-1) # minus sign to make earlier levels higher on plot
                    node_shapes[product] = product_shape
                    node_colors[product] = product_color

                    h_product += 1
                
                # add edges
                for prereq in task.data_prereqs:
                    edges.append((task_node_name, prereq))
                for product in task.data_products:
                    edges.append((task_node_name, product))
            
            G.add_nodes_from(node_names)
            G.add_edges_from(edges)

            plt.close()
            fig = plt.figure()

            # draw edges
            nx.draw_networkx_edges(G, node_positions, edge_color='k', width=2)

            # draw nodes
            for node_name in node_names:
                nx.draw_networkx_nodes(G, node_positions, nodelist=[node_name],
                                      node_color=node_colors[node_name],
                                      node_shape=node_shapes[node_name],
                                      node_size=node_size,
                                      edgecolors='k')
            
            # draw labels
            nx.draw_networkx_labels(G, node_positions, font_size=font_size, font_color='k')

            # format figure
            hs = [position[0] for position in node_positions.values()]
            vs = [position[1] for position in node_positions.values()]

            min_h, max_h = np.min(hs), np.max(hs)
            min_v, max_v = np.min(vs), np.max(vs)

            # tune size of figure
            plt.xlim(min_h - pad, max_h + pad)
            plt.ylim(min_v - pad, max_v + pad)
            width = scale * (max_h - min_h + 2 * pad)
            height = scale * (max_v - min_v + 2 * pad)
            fig.set_size_inches(width, height)
            plt.axis('off')
            plt.tight_layout()

            plt.savefig(f'{run_path}/flowchart.pdf')
            plt.close()

        # if desired, save bash scripts
        if slurm_job_name is not None:
            # script to start job
            slurm_bash_script = SlurmBashScript(
                 job_name=slurm_job_name,
                 time=slurm_job_time,
                 ntasks=slurm_job_ntasks, nodes=slurm_job_nodes,
                 mem_per_cpu=slurm_job_mem_per_cpu,
                 output=f'{run_path}/output.out', error=f'{run_path}/error.out', # absolute paths
                 mail_user=config.slurm_job_mail_user, # email address
                 mail_type='BEGIN,FAIL,END' # conditions for emailing
                 )
            
            slurm_bash_script.add_task(f'export OMP_NUM_THREADS={OMP_NUM_THREADS}')

            slurm_bash_script.add_task(f'cd {run_path}')
            if source_sdk:
                slurm_bash_script.add_task('source $MESASDK_ROOT/bin/mesasdk_init.sh')
            slurm_bash_script.add_task(f'./mk')
            slurm_bash_script.add_task(f'./rn')

            slurm_bash_script.save(self.submit_job_path)

            # script to restart job
            if slurm_job_email_user: # only email user if say so; get email from config.py file
                mail_user = config.slurm_job_mail_user
                mail_type = 'BEGIN,FAIL,END'
            else:
                mail_user = mail_type = None

            slurm_bash_script = SlurmBashScript(
                 job_name=slurm_job_name,
                 time=slurm_job_time,
                 ntasks=slurm_job_ntasks, nodes=slurm_job_nodes,
                 mem_per_cpu=slurm_job_mem_per_cpu,
                 output=f'{run_path}/output.out', error=f'{run_path}/error.out', # absolute paths
                 mail_user=mail_user, # email address
                 mail_type=mail_type # conditions for emailing
                 )
            
            slurm_bash_script.add_task(f'cd {run_path}')
            if source_sdk:
                slurm_bash_script.add_task('source $MESASDK_ROOT/bin/mesasdk_init.sh')
            slurm_bash_script.add_task(f'./mk')
            slurm_bash_script.add_task(f'./re')

            slurm_bash_script.save(self.restart_job_path)

