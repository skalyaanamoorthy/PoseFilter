import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid.inset_locator import (inset_axes, InsetPosition,
                                                  mark_inset)

T, Cp = np.loadtxt('Ta-Cp.txt', unpack=True)
T_E, CV_E = np.loadtxt('Ta-CV_Einstein.txt', unpack=True)
T_D, CV_D = np.loadtxt('Ta-CV_Debye.txt', unpack=True)

fig, ax1 = plt.subplots()
T_E = np.arange(1,max(T)+1,1)
# The data.
ax1.plot(T, Cp, 'x', c='b', mew=2, alpha=0.8, label='Experiment')
# The Einstein fit.
ax1.plot(T_E, CV_E, c='m', lw=2, alpha=0.5, label='Einstein model')
ax1.set_xlabel(r'$T\,/\mathrm{K}$')
ax1.set_ylabel(r'$C_p\,/\mathrm{J\,K^{-1}\,mol^{-1}}$')
ax1.legend(loc=0)

# Create a set of inset Axes: these should fill the bounding box allocated to
# them.
ax2 = plt.axes([0,0,1,1])
# Manually set the position and relative size of the inset axes within ax1
ip = InsetPosition(ax1, [0.4,0.2,0.5,0.5])
ax2.set_axes_locator(ip)
# Mark the region corresponding to the inset axes on ax1 and draw lines
# in grey linking the two axes.
mark_inset(ax1, ax2, loc1=2, loc2=4, fc="none", ec='0.5')

# The data: only display for low temperature in the inset figure.
Tmax = max(T_D)
ax2.plot(T[T<=Tmax], Cp[T<=Tmax], 'x', c='b', mew=2, alpha=0.8,
         label='Experiment')
# The Einstein fit (not very good at low T).
ax2.plot(T_E[T_E<=Tmax], CV_E[T_E<=Tmax], c='m', lw=2, alpha=0.5,
         label='Einstein model')
# The Debye fit.
ax2.plot(T_D, CV_D, c='r', lw=2, alpha=0.5, label='Debye model')
ax2.legend(loc=0)

# Some ad hoc tweaks.
ax1.set_ylim(0,26)
ax2.set_yticks(np.arange(0,2,0.4))
ax2.set_xticklabels(ax2.get_xticks(), backgroundcolor='w')
ax2.tick_params(axis='x', which='major', pad=8)

plt.show()




def Main():
    mol = Chem.MolFromSmiles('c1nccc2n1ccc2')

    print("get res coordinates")
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())

    # Whole position array
    pos = mol.GetConformer(0).GetPositions()

    # Gets the maximum coordinate number (adjusted)
    max_num = GetMaxCoords(pos)

    # Gets the x and y values of a residue with
    x1, y1 = GetResCoords(pos, 3, max_num)

    #GetResCoords(mol, 2, max_val)

    #AtomsWithin(ligand, protein, cutoff)

    #def MainFunction():
    m = Chem.MolFromSmiles('c1nccc2n1ccc2')

    # AllChem.Compute2DCoords(m)
    AllChem.EmbedMolecule(m, AllChem.ETKDG())
    fig = Draw.MolToMPL(m)
    pos = m.GetConformer(0).GetPositions()

    print(pos)

    # Turn off after
    plt.grid(True)

    bbox_props = dict(boxstyle="circle,pad=0.3", ec="b", lw=2)
    # max_c = ((1/2.05) + 0.25)*0.98
    import numpy as np

    # Add for each of the residues that needs to happen
    x1_graph = (complex(x1[0]).real)
    y1_graph = (complex(y1[0]).imag)

    text_add = plt.text(x1_graph , y1_graph, "Met", size=15, bbox=bbox_props)
    coord_min = 0
    coord_max = 0

    fig.savefig('mol.jpeg', bbox_inches='tight', va='center', ha='center')

    # Prints the block
    print(Chem.MolToMolBlock(m))