import qol.info as info
from matplotlib import cycler

import os

# default colors
vis_bg_color = '#edf4ff'
vis_border_color = 'k'

vis_cycler_colors = []
vis_cycler_colors.append('xkcd:deep sky blue')
vis_cycler_colors.append('xkcd:reddish pink')
vis_cycler_colors.append('xkcd:marigold')
vis_cycler_colors.append('xkcd:green teal')
vis_cycler_colors.append('xkcd:magenta')

vis_cycler = cycler(color=vis_cycler_colors)

# default line info
vis_tick_length = 5
vis_axis_border_width = 1

vis_grid_width = 0.2
vis_grid_alpha = 0.3

vis_plot_linewidth = 3

# default scatter info
vis_scatter_markersize = 50
vis_scatter_linewidth = 1.5
vis_scatter_marker = 'o'

# default font info
vis_fontpath = os.path.join(info.qol_path, 'vis/resources/fonts/HelveticaNeue-01.ttf')
vis_fontname = 'Helvetica Neue'

vis_title_size = 20
vis_label_size = 18
vis_tick_size = 16
vis_text_size = 14
vis_legend_size = 14



