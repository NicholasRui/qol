Last updated: 17 Dec 2025 by Nicholas Rui

QoL is a small package to make some of my research tasks easier. I named it QoL for "quality of life," because it makes my life easier. Since it is currently my personal research package, I develop it on-branch at the moment, and make no guarantees for all features working and producing correct results at a given time.

To use it, just clone this repository and then add the parent directory to your ```$PYTHONPATH```. Most of its prerequisites are in ```Anaconda```.

At present, it consists of the following submodules:

#### athena
For automatically creating Athena++ working directories and reading their output. This submodule assumes that you have an environment variable called ```$ATHENA_DIR```.

#### bash
For allowing for the scripted generation of bash files, particularly for those meant to launch jobs to Slurm.

#### mesa
For automatically creating MESA working directories and reading their output. This module implements ```MesaWorkingDirectory``` and ```MesaInlist``` which allow for the scripting and automatic creation of MESA working directories consisting of a potentially large number of sequential tasks (arbitrary combinations of MESA inlists and Python scripts).

This submodule also implements functions such as ```read_data``` and ```read_mod```, which read in MESA ```.data``` and ```.mod``` files and store them as a custom type of ```astropy.table.Table``` object called a ```MesaTable```.

This submodule automatically reads the environment variable ```$MESA_DIR```.

The files ```make_merger_MS_HeWD.py```, ```make_merger_RG_HeWD.py```, and ```call_create_shell_burning_remnant``` are for constructing merger remnant stellar models by generating a shell-burning-supported envelope and manually stitching it to a ``core'' stellar model. If you use these features, please cite Rui & Fuller 2025 (subm.), "Forging neon-distilling white dwarfs in the stellar engulfments of helium white dwarfs".

#### seismology
For automatically computing certain seismology properties of stellar models.

Many of these functions relate to seismic magnetometry. If you use these functions, depending on context, you might consider citing Rui, Fuller, & Hermes 2025.

#### tests
Currently only contains a few tests that the basic MESA preset working directories "compile."

#### tools
Helper functions for computing integrals, formatting numbers as strings, etc.

#### vis
I like pretty plots, but not spaghetti code. This module wraps ```matplotlib``` by creating a custom subtype of ```matplotlib.axes.Axes``` called ```FancyAxes```. It also wraps a lot of ```matplotlib``` functions to make them do change certain formatting automatically.

To use this, instead of calling ```import matplotlib.pyplot as plt```, instead write ```import qol.vis.pyplot as plt```.
