"""
Plot up location results
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

csvfile = './data/SW_picks_v4.csv'

df = pd.read_csv(csvfile).set_index('Dataset', drop=False)

U_datasets = df.loc["JPL":"BKE", "U"].to_numpy()
U_sum = df.loc["Summary", "U"]
U_std_datasets = df.loc["JPL":"BKE", "U_std"].to_numpy()
U_std_sum = df.loc["Summary", "U_std"]
labels_datasets = df.loc["JPL":"BKE", "Dataset"].to_numpy()

Uarray = np.append(U_datasets, U_sum)
Ustd = np.append(U_std_datasets, U_std_sum)
# labels = np.append(labels_datasets, ['Summary'])
labels = ["Method 1", "Method 2", "Method 3", "Method 4", "Summary"]
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

fig, axs = plt.subplots(3, 1, sharex=True, figsize=(6,9))
labelxoff = 0.95
labelyoff = 0.85

plt.sca(axs[0])
for i in np.arange(len(Uarray)):
    plt.plot(i, Uarray[i], marker='o', color=colors[i],
             label=labels[i])
    plt.errorbar(i, Uarray[i], yerr=Ustd[i], color=colors[i])
# ax = plt.gca()
axs[0].set_xticks(np.arange(len(Uarray)))
axs[0].set_xticklabels(labels)
axs[0].set_ylabel("Group velocity (km/s)")
axs[0].text(labelxoff, labelyoff, 'A', transform=axs[0].transAxes, fontsize=16)
# plt.legend()

Delta_datasets = df.loc["JPL":"BKE", "Delta"].to_numpy()
Delta_sum = df.loc["Summary", "Delta"]
Delta_std_datasets = df.loc["JPL":"BKE", "Delta_std"].to_numpy()
Delta_std_sum = df.loc["Summary", "Delta_std"]
# labels_datasets = df.loc["JPL":"BKE", "Dataset"].to_numpy()

Deltaarray = np.append(Delta_datasets, Delta_sum)
Deltastd = np.append(Delta_std_datasets, Delta_std_sum)
# labels = np.append(labels_datasets, ['Summary'])
# colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

MQS_Delta = 37.01
MQS_Delta_std = 0.94

# fig, axs = plt.subplots(1, 3, sharey=True)
plt.sca(axs[1])
for i in np.arange(len(Deltaarray)):
    plt.plot(i, Deltaarray[i], marker='o', color=colors[i],
             label=labels[i])
    plt.errorbar(i, Deltaarray[i], yerr=Deltastd[i], color=colors[i])
# ax = plt.gca()
# axs[0].set_yticks(np.arange(len(Uarray)))
# axs[0].set_yticklabels(labels)
axs[1].set_ylabel("Distance (degrees)")
xmin, xmax = axs[1].get_xlim()
plt.plot([xmin, xmax], [MQS_Delta, MQS_Delta], 'k--')

axs[1].add_patch(Rectangle((xmin, MQS_Delta - 2.*MQS_Delta_std), xmax-xmin,
                           4.*MQS_Delta_std, color='grey', alpha=0.5))
axs[1].set_xlim((xmin, xmax))
axs[1].text(labelxoff, labelyoff, 'B', transform=axs[1].transAxes, fontsize=16)

t0_datasets = df.loc["JPL":"BKE", "t0"].to_numpy()
t0_sum = df.loc["Summary", "t0"]
t0_std_datasets = df.loc["JPL":"BKE", "t0_std"].to_numpy()
t0_std_sum = df.loc["Summary", "t0_std"]
# labels_datasets = df.loc["JPL":"BKE", "Dataset"].to_numpy()

t0array = np.append(t0_datasets, t0_sum)
t0std = np.append(t0_std_datasets, t0_std_sum)
# labels = np.append(labels_datasets, ['Summary'])
# colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

MQS_t0 = -279.1
MQS_t0_std = 4.6

# fig, axs = plt.subplots(1, 3, sharey=True)
plt.sca(axs[2])
for i in np.arange(len(t0array)):
    plt.plot(i, t0array[i], marker='o', color=colors[i],
             label=labels[i])
    plt.errorbar(i, t0array[i], yerr=t0std[i], color=colors[i])
# ax = plt.gca()
# axs[0].set_yticks(np.arange(len(Uarray)))
# axs[0].set_yticklabels(labels)
axs[2].set_ylabel("t$_0$ (seconds relative to P)")
xmin, xmax = axs[2].get_xlim()
plt.plot([xmin, xmax], [MQS_t0, MQS_t0], 'k--')
axs[2].add_patch(Rectangle((xmin, MQS_t0 - 2.*MQS_t0_std), xmax-xmin,
                           4.*MQS_t0_std, color='grey', alpha=0.5))
axs[2].set_xlim((xmin, xmax))
axs[2].text(labelxoff, labelyoff, 'C', transform=axs[2].transAxes, fontsize=16)

plt.savefig("Figure3.png")
