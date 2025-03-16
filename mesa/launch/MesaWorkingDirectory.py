import qol.config as config
import qol.paths as paths

import numpy as np

import matplotlib.pyplot as plt
import networkx as nx

import os
import shutil

class MesaWorkingDirectory:
    """
    Stores information for creating a custom MESA work directory

    Records tasks, which must have the following attributes / methods:
    - rn_string : returns string to put into rn to run this task
    - save : returns command which saves the needed files
    - prereqs : list of prereq files inputed to task
    - products : list of product files outputed by task

    TODO
    need inlist (straight-up inlist file... header?)
    save shell for running job
    """
    def __init__(self, run_path, mesa_version=config.mesa_version):
        if os.path.exists(run_path):
            raise ValueError(f'Path already exists: {run_path}')

        self.run_path = run_path
        self.mesa_version = mesa_version

        self.history_columns_path = None
        self.profile_columns_path = None
        # self.inlist_pgstar = None

        self.copy_from_path_root_prereqs = []
        self.rel_paths_root_prereq = []
        self.tasks = []

    def add_root_prereq(self, copy_from_abs_path, rel_path):
        """
        Add a model file which serves as a "root prereq"
        (something used by an inlist which is there from the beginning,
                    not outputted as an intermediate product by a task)
        
        Copy from copy_from_abs_path (absolute path) to rel_path (relative path in work directory)
        """
        if not os.path.exists(copy_from_abs_path):
            raise ValueError(f'Path not found: {copy_from_abs_path}')

        self.copy_from_path_root_prereqs.append(copy_from_abs_path)
        self.rel_paths_root_prereq.append(rel_path)
    
    def add_task(self, task):
        """
        Add a task to the task list, after determining the prereqs already exist
        """
        self.check_needed_prereqs(task)
        self.check_fname_conflicts(task)

        self.tasks.append(task)

    def check_needed_prereqs(self, task):
        """
        Check if, with addition of new task, prereqs exist, and throw error if not
        Otherwise, do nothing
        """
        task_prereqs = task.prereqs

        existing_products = []
        for task in self.tasks:
            existing_products += task.products
        existing_products += self.rel_paths_root_prereq

        if not np.in1d(task_prereqs, existing_products).all():
            raise ValueError('task requests prereqs which do not exist')
        
        return
    
    def check_fname_conflicts(self, task):
        """
        make sure the new task doesn't present an fname conflict
        """
        existing_rel_paths = []
        for task in self.tasks:
            existing_rel_paths += task.rel_path
        existing_rel_paths += self.rel_paths_root_prereq

        if task.rel_path in existing_rel_paths:
            raise ValueError('rel_path of task conflicts with existing rel_path')
        
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

    # def copy_inlist_pgstar(self, abs_path=None):
    #     """
    #     Copy inlist_pgstar file from abs_path

    #     None means try to find a preset one defined in qol.
    #     """
    #     if os.path.exists(abs_path):
    #         self.inlist_pgstar_path = abs_path
    #     else:
    #         raise ValueError(f'No path found: {abs_path}')

    # def check_and_sort(self):
    #     """
    #     Make sure output and LOGS names do not overlap, and
    #     that each inlist will eventually have its prereqs
    #     """
    #     # loop through until either everything is counted or the number of things doesn't change

    #     ...

    # def save_flowchart(self):
    #     """
    #     Make flowchart which shows what is run in which order,
    #     and which outputs are being used as inputs
    #     """
    #     ...

    def save_directory(self, make_flowchart=True):
        """
        Create MESA directory
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

        mesadir = paths.mesa_paths[self.mesa_version]
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
        
        RSE_out_text = '\n'.join(newlines)

        with open(RSE_out_path,'w') as f:
            f.write(RSE_out_text)
        
        # copy over history and profile columns files, if specified
        if self.history_columns_path is not None:
            shutil.copy(self.history_columns_path, f'{run_path}history_columns.list')
        if self.profile_columns_path is not None:
            shutil.copy(self.profile_columns_path, f'{run_path}profile_columns.list')
        # and inlist_pgstar
        # if self.inlist_pgstar_path is not None:
        #     shutil.copy(self.inlist_pgstar_path, f'{run_path}inlist_pgstar')

        # copy over other files
        shutil.copy(f'{paths.qol_path}/mesa/resources/do_one', f'{run_path}do_one')
        shutil.copy(f'{paths.qol_path}/mesa/resources/re', f'{run_path}re')

        # save all root prereqs
        for ii, rel_path in enumerate(self.rel_paths_root_prereq):
            copy_from_abs_path = self.copy_from_path_root_prereqs[ii]
            shutil.copy(copy_from_abs_path, f'{self.run_path}/{rel_path}')

        # LOOP OVER TASKS AND RUN THEM IN ORDER
        vlevels = [] # vertical level in chart
        hlevels = [] # horizontal level in chart
        sorted_tasks = []

        existing_products = []
        existing_products += self.rel_paths_root_prereq

        level_products = []
        task_ids = [] # dummy index to skip tasks which have been handled
        vlevel, hlevel = 0, 0

        rn_text = ''

        while True:
            new_tasks = 0

            for ii, task in enumerate(self.tasks):
                # if task visited before, skip it
                if ii in task_ids:
                    continue

                # if all prereqs are in existing_products exist,
                # run task, then add task to sorted_tasks and increment
                if np.in1d(task.prereqs, existing_products).all():
                    rn_text += f'{task.rn_string()}\n' # add to rn_text
                    task.save(run_path=run_path)

                    task_ids.append(ii)
                    sorted_tasks.append(task)
                    level_products += task.products
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

        # save rn file
        with open(f'{run_path}rn', 'w') as f:
            f.write(rn_text)
        

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
                for product in task.products:
                    node_names.append(product)
                    node_positions[product] = (h_product+h_shift, -2*v-1) # minus sign to make earlier levels higher on plot
                    node_shapes[product] = product_shape
                    node_colors[product] = product_color

                    h_product += 1
                
                # add edges
                for prereq in task.prereqs:
                    edges.append((task_node_name, prereq))
                for product in task.products:
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

            plt.savefig(f'{run_path}/flowchart.png')
            plt.close()
