import seaborn as sns
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd

fig_folder = "figure_data/fig3_SEI/"
lb_fs = 12
tl_fs = 10

#e_color =
sns.set_palette("Set1", n_colors=8, desat=.5)

# get the data
dos = pd.read_csv(fig_folder+"dos_interp.txt", sep='\t', skiprows=1, names=["E", "DOS"])

fit_data = pd.read_csv(fig_folder+"k_vs_V.txt", sep='\t', skiprows=1, names=["V", "MHC_k", "MHC_DOS_k"])
#reformat k data for seaborn
MHC = fit_data[['V', "MHC_k"]]
MHC['fit'] = 'MHC'
MHC.rename(columns={'MHC_k':'k'}, inplace=True)
MHC_DOS = fit_data[['V', "MHC_DOS_k"]]
MHC_DOS['fit'] = 'MHC_DOS'
MHC_DOS.rename(columns={'MHC_DOS_k':'k'}, inplace=True)
fit_data = pd.concat([MHC, MHC_DOS])

# plot it all
f, (a0, a1) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [1,2]}, figsize=(3.7, 4.8))
# figure a
sns.lineplot(x="E", y="DOS", data=dos, legend=False, linewidth=2, ax=a0, color='grey')
a0.set_xlim(-5, 5)
a0.set_xlabel("$E-E_f$ [eV]", fontsize=lb_fs)
a0.set_ylabel("DOS [arb.]", fontsize=lb_fs, labelpad=11)
a0.tick_params(labelsize=tl_fs)
leg_str = "LiF/Li$_2$CO$_3$ DOS"
a0.legend([leg_str])
# figure b
sns.lineplot(x="V", y="k", data=fit_data, hue='fit', legend=False, linewidth=2, ax=a1)
a1.set_yscale('log')
a1.set_ylabel("Rate constant $k$", fontsize=lb_fs)
a1.set_xlabel("Overpotential $\eta$ [V]", fontsize=lb_fs)
a1.tick_params(labelsize=tl_fs)
a1.set_xlim([-5, 5])
#a1.set_ylim([5e-9, 6])
plt.legend(["MHC", "MHC+DOS"], fontsize=tl_fs-1, loc='lower left')
plt.tight_layout()
plt.savefig("figs/for_paper/fig3.png", dpi=300)
