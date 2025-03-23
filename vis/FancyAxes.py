# Defines base subclass of matplotlib.axes.Axes object
# which automatically makes nice-looking plots
# add helvetica neue font
"""
TODO
- easily change format to put 1, 3, 5, 2, 3, etc., on axes
- easily put a secondary axis on x or y (twinx, twiny)
"""

import qol.defaults as defaults

import matplotlib.pyplot as plt
import numpy as np



class FancyAxes(plt.Axes):
    """
    serves as a base class for other qol-specific axes types
    """
    def __init__(self, fig, rect, *args, **kwargs):
        """
        """
        super().__init__(fig, rect, *args, **kwargs)
        self.fancy_formatting()

    def fancy_formatting(self):
        """
        Formats ticks, grid, etc., to look nicer
        """
        # Put ticks inside axis and on both sides, and configure spines
        self.tick_params(which='both', direction='in',
                         length=defaults.vis_tick_length,
                         width=defaults.vis_axis_border_width,
                         labelsize=defaults.vis_tick_size,
                         zorder=np.inf)
        self.xaxis.set_ticks_position('both')
        self.yaxis.set_ticks_position('both')

        # Format spines
        for spine in self.spines.values():
            spine.set_color(defaults.vis_border_color)
            spine.set_linewidth(defaults.vis_axis_border_width)
            spine.set_zorder(np.inf)
        
        # Format face and grid
        self.set_facecolor(defaults.vis_bg_color)
        self.grid(c='gray', which='both', alpha=defaults.vis_grid_alpha, linewidth=defaults.vis_grid_width)
        self.set_axisbelow(True) # put grid in background : inf zorder for spines/ticks ensures that they're still on top


    # Re-define some matplotlib methods to have better defaults
    def set_title(self, label, fontdict=None, loc='center', pad=None, **kwargs):
        kwargs.setdefault('fontsize', defaults.vis_title_size)
        super().set_title(label, fontdict=fontdict, loc=loc, pad=pad, **kwargs)

    def set_xlabel(self, xlabel, fontdict=None, labelpad=None, *args, loc=None, **kwargs):
        print(xlabel)
        kwargs.setdefault('fontsize', defaults.vis_label_size)
        super().set_xlabel(xlabel, fontdict=fontdict, labelpad=labelpad, *args, loc=loc, **kwargs)

    def set_ylabel(self, ylabel, fontdict=None, labelpad=None, *args, loc=None, **kwargs):
        kwargs.setdefault('fontsize', defaults.vis_label_size)
        super().set_ylabel(ylabel, fontdict=fontdict, labelpad=labelpad, *args, loc=loc, **kwargs)
    







