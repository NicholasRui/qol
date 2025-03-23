# Redefine some functions in pyplot to make handling FancyFigure and FancyAxes object as simple as possible

from qol.vis.FancyAxes import FancyAxes
import qol.defaults as defaults

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm, rcParams

# from matplotlib.figure import Figure

# Add font
fe = fm.FontEntry(fname=defaults.vis_fontpath, name=defaults.vis_fontname)
fm.fontManager.ttflist.insert(0, fe)
matplotlib.rcParams['font.family'] = fe.name
params = {'mathtext.default': 'regular' }
plt.rcParams.update(params)

# Import all pyplot functions and only overwrite those which we want to change
# Important: when referencing original pyplot functions in defined ones, always include plt
from matplotlib.pyplot import *

#################################
###### REDEFINED FUNCTIONS ######
#################################

def axes(arg: None | tuple[float, float, float, float] = None,
         **kwargs):
    """
    Redefined to make the axes FancyAxes, and format accordingly
    """
    ax = plt.axes(arg, **kwargs)
    ax.__class__ = FancyAxes
    ax.fancy_formatting()

    return ax

# @_copy_docstring_and_deprecators(Figure.savefig)
def savefig(*args, **kwargs) -> None:
    """
    Set default: bbox_inches='tight'
    """
    kwargs.setdefault('bbox_inches', 'tight')
    res = plt.savefig(*args, **kwargs)

    return res

# @_docstring.interpd
def subplot(*args, **kwargs) -> FancyAxes:
    ax = plt.subplot(*args, **kwargs)
    ax.__class__ = FancyAxes
    ax.fancy_formatting()

    return ax

###########################
###### NEW FUNCTIONS ######
###########################

def cividis() -> None:
    """
    Set the colormap to 'cividis'.

    This changes the default colormap as well as the colormap of the current
    image if there is one. See ``help(colormaps)`` for more information.
    """
    plt.set_cmap("cividis")
