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

def axes(arg: None | tuple[float, float, float, float] = None, blank=False,
         **kwargs):
    """
    Redefined to make the axes FancyAxes, and format accordingly
    """
    ax = plt.axes(arg, **kwargs)
    ax.__class__ = FancyAxes
    ax.blank = blank
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
def subplot(blank=False, *args, **kwargs) -> FancyAxes:
    ax = plt.subplot(*args, **kwargs)
    ax.__class__ = FancyAxes
    ax.blank = blank
    ax.fancy_formatting()

    return ax

def subplots(nrows=1, ncols=1, blank=False, *, sharex=False, sharey=False, squeeze=True, width_ratios=None, height_ratios=None, subplot_kw=None, gridspec_kw=None, **fig_kw):
    fig, axs = plt.subplots(nrows=nrows, ncols=ncols, sharex=sharex, sharey=sharey,
                            squeeze=squeeze, width_ratios=width_ratios, height_ratios=height_ratios,
                            subplot_kw=subplot_kw, gridspec_kw=gridspec_kw, **fig_kw)
    
    if nrows == ncols == 1:
        ax = axs
        ax.__class__ = FancyAxes
        ax.blank = blank
        ax.fancy_formatting()
    elif nrows == 1 or ncols == 1:
        for ax in axs:
            ax.__class__ = FancyAxes
            ax.blank = blank
            ax.fancy_formatting()
    else:
        for _, axs_row in enumerate(axs):
            for ax in axs_row:
                ax.__class__ = FancyAxes
                ax.blank = blank
                ax.fancy_formatting()
    
    return fig, axs

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
