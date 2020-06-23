

def CreateHeatMap(fingerprint_data, Ligands, type):
  #  figure.set = ()
    sns.set()
    normal_data = np.random.randn(10, 12)

    ax = sns.heatmap(fingerprint_data, cmap="Greens", xticklabels= Ligands, yticklabels= Ligands)
   # ax = sns.heatmap(fingerprint_data)
    ax.set(title = type)
    figure = ax.get_figure()
    figure.savefig(type + ".png")