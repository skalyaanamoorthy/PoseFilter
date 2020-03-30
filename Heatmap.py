import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
np.random.seed(0)
sns.set()

# Load the example flights dataset and convert to long-form
LigNum = 5  # 5 total
fingerprint_data = open("Fingerprint.csv", 'r')

uniform_data = np.empty((LigNum, LigNum), dtype=float)

l_count = 0
for line in fingerprint_data:
   # print(line)
    if l_count == 0:
        pass

    col_count = 0
    for item in line.split(','):
        if col_count == 0:
            col_count +=1
        else:
            print("is num: " + str(item))
            col = col_count-1
            print(str(l_count) + " " + str(col))
            uniform_data[l_count-1, col] = item
            print(item)
           # print(uniform_data[l_count-1, col])
            col_count += 1

    l_count += 1


print(uniform_data)

#ax = sns.heatmap(uniform_data)

#flights_long = sns.load_dataset("Fingerprint.csv")
#flights = flights_long.pivot("ligand", "ligand", "passengers")

# Draw a heatmap with the numeric values in each cell
#f, ax = plt.subplots(figsize=(9, 6))
#sns.heatmap(uniform_data, annot=True, fmt="d", linewidths=.5, ax=ax)
#uniform_data = uniform_data.pivot([], "year", "passengers")
dp = sns.heatmap(uniform_data, cmap="Greens")
#ax.savefig("Fingerprint.png")
#sns_plot = sns.pairplot(ax, hue='species', size=2.5)
dp.figure.savefig("output.png")


