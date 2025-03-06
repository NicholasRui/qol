


# import matplotlib.pyplot as plt
# import numpy as np

# # TODO Put some of this in a config file...
# import matplotlib
# matplotlib.use('Agg')

# from matplotlib import font_manager as fm, rcParams
# fe = fm.FontEntry(fname='/home/nrui/common/HelveticaNeue-01.ttf', name='Helvetica Neue')
# fm.fontManager.ttflist.insert(0, fe)
# matplotlib.rcParams['font.family'] = fe.name
# params = {'mathtext.default': 'regular' }
# plt.rcParams.update(params)

# # Presets: TODO put in config file

# label_size = 20
# tick_size = 16

# def fancy_background(ax):
#     """
#     Give the ax a fancy background and custom formatted ticks.
#     """
#     ax.tick_params(which='both', direction='in', length=5, width=1, labelsize=tick_size)
#     ax.xaxis.set_ticks_position('both')
#     ax.yaxis.set_ticks_position('both')
#     ax.spines['top'].set_linewidth(1)
#     ax.spines['right'].set_linewidth(1)
#     ax.spines['bottom'].set_linewidth(1)
#     ax.spines['left'].set_linewidth(1)
#     ax.set_facecolor('#F3FAFA')
#     plt.sca(ax)
#     plt.grid(c='gray', zorder=-10000, which='both', alpha=0.3, linewidth=0.2)
#     plt.xticks(fontsize=tick_size); plt.yticks(fontsize=tick_size)






# class PropDiagram:
#     """
#     Make an axis a propagation diagram
#     """
#     def __init__(self, ax=plt.gca(), xvar='r', ynorm='angular_uHz'):

#         self.ax = ax
#         self.xvar = xvar
#         self.ynorm = ynorm
        
#         fancy_background(ax)

#         # Handle x axis information
#         match xvar:
#             case 'r': # radius in Rsun
#                 self.lambda_x = lambda p: 10 ** p['logR']
#                 ax.set_xlabel(r'$r$ ($R_{\!\odot}\!$)', fontsize=label_size)
#             case 'r/R': # fractional radius
#                 self.lambda_x = lambda p: 10 ** (p['logR'] - p['logR'][0])
#                 ax.set_xlabel(r'$r/R$', fontsize=label_size)
#             case 'm': # mass coordinate in Msun
#                 self.lambda_x = lambda p: p['mass']
#                 ax.set_xlabel(r'$m$ ($M_{\!\odot}\!$)', fontsize=label_size)
#             case 'm/M': # fractional mass coordinate
#                 self.lambda_x = lambda p: p['mass'] / p['mass'][0]
#                 ax.set_xlabel(r'$m/M$', fontsize=label_size)
#             case '1-r/R': # fractional radius from surface
#                 self.lambda_x = lambda p: 1 - 10 ** (p['logR'] - p['logR'][0])
#                 ax.set_xlabel(r'$1-r/R$', fontsize=label_size)
#             case _:
#                 raise ValueError(f'Invalid xvar: {xvar}')
        
#         # Handle y axis information
#         match ynorm:
#             case 'angular_uHz':
#                 yfactor = 1e6
#                 ax.set_ylabel(r'$\omega$ ($\mu$Hz)', fontsize=label_size)
#             case 'linear_uHz':
#                 yfactor = 1e6 / 2 / np.pi
#                 ax.set_ylabel(r'$\nu$ ($\mu$Hz)', fontsize=label_size)

#         # # Plot Bruntx
#         # N2thresh = 1e-30
#         # bruntN2 = p['brunt_N2']
#         # bruntN2[bruntN2 < N2thresh] = N2thresh
#         # bruntN = np.sqrt(bruntN2)


#         # method to plot stuff
#         # work out colors
#         # logx, logy?




















