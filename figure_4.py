import seaborn as sns
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd

fig_folder = "figure_data/fig5/"
old_SEI_folder = "figure_data/fig4_SEI/"
old_Cu_folder = "figure_data/fig2_Cu111/"
lb_fs = 12
tl_fs = 10

sns.set_palette([sns.color_palette("Set1", n_colors=8, desat=.5)[i] for i in [1,4]])

# get the data for part a
old_data = pd.read_csv(old_SEI_folder+"k_vs_V.txt", sep='\t', skiprows=1, names=["V", "MHC_k", "MHC_DOS_k"])
new_data = pd.read_csv(fig_folder+"SEI_k_vs_V_Ef-2.txt", sep='\t', skiprows=1, names=["V", "MHC_k", "MHC_DOS_k"])
# reformat
old_data.drop("MHC_k", axis=1, inplace=True)
old_data["fit"] = "old"
new_data.drop("MHC_k", axis=1, inplace=True)
new_data["fit"] = "new"
SEI_data = pd.concat([old_data, new_data])

# and get data for part b
old_data = pd.read_csv(old_Cu_folder+"k_vs_V.txt", sep='\t', skiprows=1, names=["V", "MHC_k", "MHC_DOS_k"])
new_data = pd.read_csv(fig_folder+"Cu111_k_vs_V_doublelambda.txt", sep='\t', skiprows=1, names=["V", "MHC_k", "MHC_DOS_k"])
old_data.drop("MHC_k", axis=1, inplace=True)
old_data["fit"] = "old"
new_data.drop("MHC_k", axis=1, inplace=True)
new_data["fit"] = "new"
Cu_data = pd.concat([old_data, new_data])

# plot it all
f, (a0, a1) = plt.subplots(2, 1, figsize=(3.7, 5))
# # figure a: two cases of SEI
sns.lineplot(x="V", y="MHC_DOS_k", data=SEI_data, hue='fit', legend=False, linewidth=2, ax=a0)
a0.set_yscale('log')
a0.set_ylabel("Rate constant $k$", fontsize=lb_fs)
a0.set_xlabel("Overpotential $\eta$ [V]", fontsize=lb_fs)
a0.tick_params(labelsize=tl_fs)
a0.set_xlim([-5, 5])
a0.legend(["$E_f^{eq}$ (from fig. 3)", "$E_f^{eq}-2$ eV"], fontsize=tl_fs-1, loc='lower left')

# figure b: reorg energy in Cu111
sns.lineplot(x="V", y="MHC_DOS_k", data=Cu_data, hue='fit', legend=False, linewidth=2, ax=a1)
a1.set_yscale('log')
a1.set_ylabel("Rate constant $k$", fontsize=lb_fs)
a1.set_xlabel("Overpotential $\eta$ [V]", fontsize=lb_fs)
a1.tick_params(labelsize=tl_fs)
a1.set_xlim([-3, 3])
a1.legend(["$\lambda_{orig}$ (from fig. 2)", "$2\lambda_{orig}$ eV"], fontsize=tl_fs-1, loc='lower left')

plt.tight_layout()
plt.savefig("figs/for_paper/fig5.png", dpi=300)
