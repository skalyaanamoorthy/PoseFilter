import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.axes_grid.inset_locator as mp

import mpl_toolkits.axes_grid1.inset_locator


fig, ax1 = plt.subplots()
#T_E = np.arange(1,max(T)+1,1)
# The data.
#ax1.plot(T, Cp, 'x', c='b', mew=2, alpha=0.8, label='Experiment')
# The Einstein fit.
#ax1.plot(T_E, CV_E, c='m', lw=2, alpha=0.5, label='Einstein model')
#ax1.set_xlabel(r'$T\,/\mathrm{K}$')
#ax1.set_ylabel(r'$C_p\,/\mathrm{J\,K^{-1}\,mol^{-1}}$')
#ax1.legend(loc=0)

# Create a set of inset Axes: these should fill the bounding box allocated to
# them.
ax2 = plt.axes([0,0,1,1])
# Manually set the position and relative size of the inset axes within ax1
#ip = InsetPosition(ax1, [0.4,0.2,0.5,0.5])

#ip = InsetPosition(ax2, [0.5, 0.1, 0.4, 0.8])
#ax_ins.set_axes_locator(ip)
#ax2.set_axes_locator(ip)

ax_ins = plt.axes([0, 0, 1, 1])
ip = mp.InsetPosition(ax2, [0.5, 0.1, 0.4, 0.2])
#mp.ax_ins.set_axes_locator(ip)
# Mark the region corresponding to the inset axes on ax1 and draw lines


# Some ad hoc tweaks.
#ax1.set_ylim(0,26)
#ax2.set_yticks(np.arange(0,2,0.4))
#ax2.set_xticklabels(ax2.get_xticks(), backgroundcolor='w')
#ax2.tick_params(axis='x', which='major', pad=8)

plt.show()