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

    blank: if True, doing fancy_formatting removes all axes, ticks, etc.
    """
    def __init__(self, fig, rect, blank=False, *args, **kwargs):
        """
        """
        self.blank = blank

        super().__init__(fig, rect, *args, **kwargs)
        self.fancy_formatting()

    def fancy_formatting(self):
        """
        Formats ticks, grid, etc., to look nicer
        """
        if self.blank: # if blank, remove EVERYTHING
            self.axis('off')
            self.set_facecolor('none')
        else:
            # Put ticks inside axis and on both sides, and configure spines
            self.tick_params(which='both', direction='in',
                            length=defaults.vis_tick_length,
                            width=defaults.vis_axis_border_width,
                            labelsize=defaults.vis_tick_size,
                            zorder=np.inf)
            self.xaxis.set_ticks_position('both')
            self.yaxis.set_ticks_position('both')

            # format spines
            for spine in self.spines.values():
                spine.set_color(defaults.vis_border_color)
                spine.set_linewidth(defaults.vis_axis_border_width)
                # spine.set_zorder(np.inf)
            
            # format face and grid
            self.set_facecolor(defaults.vis_bg_color)
            self.grid(c='gray', which='both', alpha=defaults.vis_grid_alpha, linewidth=defaults.vis_grid_width)
            self.set_axisbelow(True) # put grid in background : inf zorder for spines/ticks ensures that they're still on top

            # set zorder of ticks again...
            # self.yaxis.set_zorder(np.inf)
            # self.xaxis.set_zorder(np.inf)

            # grid zorder
            for gridline in self.get_xgridlines() + self.get_ygridlines():
                gridline.set_zorder(-np.inf)

        # set color cycler
        self.set_prop_cycle(defaults.vis_cycler)

    # New convenience methods
    def corner_text(self, text, loc, dx=0.03, dy=0.02, **kwargs):
        """
        show text on a corner, displaced from it by dx, dy in relative axis units (i.e., from 0 to 1)
        loc = 'upper left', 'lower left', 'center center', etc.
        """
        match loc:
            case 'upper left':
                x, y = 0+dx, 1-dy
                ha, va = 'left', 'top'
            case 'center left':
                x, y = 0+dx, 0.5
                ha, va = 'left', 'center'
            case 'lower left':
                x, y = 0+dx, 0+dy
                ha, va = 'left', 'bottom'
            case 'upper center':
                x, y = 0.5, 1-dy
                ha, va = 'center', 'top'
            case 'center center':
                x, y = 0.5, 0.5
                ha, va = 'center', 'center'
            case 'lower center':
                x, y = 0.5, 0+dy
                ha, va = 'center', 'bottom'
            case 'upper right':
                x, y = 1-dx, 1-dy
                ha, va = 'right', 'top'
            case 'center right':
                x, y = 1-dx, 0.5
                ha, va = 'right', 'center'
            case 'lower right':
                x, y = 1-dx, 0+dy
                ha, va = 'right', 'bottom'
    
        self.annotate(text, xy=(x, y), xycoords='axes fraction',
                ha=ha, va=va, **kwargs)

    # Re-define some matplotlib methods to have better defaults
    # TODO: put this in a different file
    def set_title(self, label, fontdict=None, loc='center', pad=None, **kwargs):
        kwargs.setdefault('fontsize', defaults.vis_title_size)
        super().set_title(label, fontdict=fontdict, loc=loc, pad=pad, **kwargs)

    def set_xlabel(self, xlabel, fontdict=None, labelpad=None, *args, loc=None, **kwargs):
        kwargs.setdefault('fontsize', defaults.vis_label_size)
        super().set_xlabel(xlabel, fontdict=fontdict, labelpad=labelpad, *args, loc=loc, **kwargs)

    def set_ylabel(self, ylabel, fontdict=None, labelpad=None, *args, loc=None, **kwargs):
        kwargs.setdefault('fontsize', defaults.vis_label_size)
        super().set_ylabel(ylabel, fontdict=fontdict, labelpad=labelpad, *args, loc=loc, **kwargs)
    
    def legend(self, *args, **kwargs):
        kwargs.setdefault('fontsize', defaults.vis_legend_size)
        super().legend(*args, **kwargs)
    
    def text(self, x, y, s, fontdict=None, **kwargs):
        kwargs.setdefault('fontsize', defaults.vis_text_size)
        super().text(x, y, s, fontdict=fontdict, **kwargs)

    def plot(self, *args, scalex=True, scaley=True, data=None, **kwargs):
        # need to do this to handle aliases
        if 'linewidth' not in kwargs.keys() and \
                'lw' not in kwargs.keys():
            kwargs['linewidth'] = defaults.vis_plot_linewidth
        if 'markersize' not in kwargs.keys() and \
                'ms' not in kwargs.keys():
            kwargs['markersize'] = defaults.vis_scatter_markersize

        super().plot(*args, scalex=scalex, scaley=scaley, data=data, **kwargs)

    def annotate(self, text, xy, xytext=None, xycoords='data', textcoords=None, arrowprops=None, annotation_clip=None, **kwargs):
        kwargs.setdefault('fontsize', defaults.vis_text_size)

        super().annotate(text, xy,
                xytext=xytext,
                xycoords=xycoords,
                textcoords=textcoords,
                arrowprops=arrowprops,
                annotation_clip=annotation_clip,
                **kwargs)

    def scatter(self, x, y, s=None, c=None, *, marker=None, cmap=None, norm=None,
                vmin=None, vmax=None, alpha=None, linewidths=None, edgecolors=None,
                colorizer=None, plotnonfinite=False, data=None, **kwargs):
        marker = marker if marker is not None else defaults.vis_scatter_marker
        if linewidths is None and \
                'linewidth' not in kwargs.keys() and \
                'lw' not in kwargs.keys():
            linewidths = defaults.vis_scatter_linewidth
        s = s if s is not None else defaults.vis_scatter_markersize

        super().scatter(x, y, s=s, c=c, marker=marker, cmap=cmap, norm=norm,
                vmin=vmin, vmax=vmax, alpha=alpha, linewidths=linewidths, edgecolors=edgecolors,
                colorizer=colorizer, plotnonfinite=plotnonfinite, data=data, **kwargs)
