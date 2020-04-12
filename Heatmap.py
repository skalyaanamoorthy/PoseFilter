import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np


def CreateHeatMap(fingerprint_data, Ligands):
    sns.set()
    LigNum = len(Ligands)

    data = np.empty((LigNum, LigNum), dtype=float)

    l_count = 0
    for line in fingerprint_data:
        col_count = 0
        for item in line.split(','):
          #  print("is num: " + str(item))
            col = col_count-1
         #   print(str(l_count) + " " + str(col))
            data[l_count-1, col] = item
         #   print(item)
            col_count += 1
        l_count += 1

    print(data)

    #ax = sns.heatmap(uniform_data)

    #flights_long = sns.load_dataset("Fingerprint.csv")
    #flights = flights_long.pivot("ligand", "ligand", "passengers")

    # Draw a heatmap with the numeric values in each cell
    #f, ax = plt.subplots(figsize=(9, 6))
    #sns.heatmap(uniform_data, annot=True, fmt="d", linewidths=.5, ax=ax)
    #uniform_data = uniform_data.pivot([], "year", "passengers")
    ax = sns.heatmap(data, cmap="Greens", xticklabels= Ligands, yticklabels= Ligands)
    ax.set(title = 'Heatmap')

    #ax.savefig("Fingerprint.png")
    #sns_plot = sns.pairplot(ax, hue='species', size=2.5)
    ax.figure.savefig("output.png")


# Load the example flights dataset and convert to long-form
LigNum = 5  # 5 total
fingerprint_data = open("Fingerprint.csv", 'r')

CreateHeatMap(fingerprint_data, ["Lig1", "Lig2", "Lig3", "Lig4", "Lig5"])

