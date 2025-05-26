# Redefine some functions in pyplot to make handling FancyFigure and FancyAxes object as simple as possible

from qol.vis.FancyAxes import FancyAxes
import qol.defaults as defaults

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm, rcParams

# Boilerplate copied imports
from numpy.typing import ArrayLike
from collections.abc import Callable, Hashable, Iterable, Sequence
from matplotlib.typing import (
        ColorType,
        CoordsType,
        HashableList,
        LineStyleType,
        MarkerType,
    )
from typing import Any, BinaryIO, Literal, TypeVar
from matplotlib.colorizer import _ColorizerInterface, ColorizingArtist, Colorizer
from matplotlib.colors import _color_sequences, Colormap
from matplotlib.colors import Normalize
from matplotlib.collections import (
    Collection,
    FillBetweenPolyCollection,
    LineCollection,
    PolyCollection,
    PathCollection,
    EventCollection,
    QuadMesh,
)

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

# Redefining plot and scatter in order to make use of new gca function

#@_copy_docstring_and_deprecators(Figure.gca)
def gca() -> Axes:
    ax = plt.gcf().gca()
    if not isinstance(ax, FancyAxes):
        ax.__class__ = FancyAxes
        ax.blank = False
        ax.fancy_formatting()

    return ax

def sci(im: ColorizingArtist) -> None:
    gca()._sci(im)

def plot(
    *args: float | ArrayLike | str,
    scalex: bool = True,
    scaley: bool = True,
    data=None,
    **kwargs,
) -> list[plt.Line2D]:
    return gca().plot(
        *args,
        scalex=scalex,
        scaley=scaley,
        **({"data": data} if data is not None else {}),
        **kwargs,
    )

# @_copy_docstring_and_deprecators(Axes.scatter)
def scatter(
    x: float | ArrayLike,
    y: float | ArrayLike,
    s: float | ArrayLike | None = None,
    c: ArrayLike | Sequence[ColorType] | ColorType | None = None,
    marker: MarkerType | None = None,
    cmap: str | Colormap | None = None,
    norm: str | Normalize | None = None,
    vmin: float | None = None,
    vmax: float | None = None,
    alpha: float | None = None,
    linewidths: float | Sequence[float] | None = None,
    *,
    edgecolors: Literal["face", "none"] | ColorType | Sequence[ColorType] | None = None,
    colorizer: Colorizer | None = None,
    plotnonfinite: bool = False,
    data=None,
    **kwargs,
) -> PathCollection:
    __ret = gca().scatter(
        x,
        y,
        s=s,
        c=c,
        marker=marker,
        cmap=cmap,
        norm=norm,
        vmin=vmin,
        vmax=vmax,
        alpha=alpha,
        linewidths=linewidths,
        edgecolors=edgecolors,
        colorizer=colorizer,
        plotnonfinite=plotnonfinite,
        **({"data": data} if data is not None else {}),
        **kwargs,
    )
    sci(__ret)
    return __ret

