import seaborn as sns
from matplotlib import pyplot as plt
#from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import numpy as np
import pandas as pd

fig_folder = "figure_data/fig1_Li100/"
lb_fs = 12
tl_fs = 10

#e_color =
sns.set_palette("Set1", n_colors=8, desat=.5)

## FIGURE 1
with open(fig_folder+"ecdec_fitparams.txt","r") as f:
     params = f.readlines()[1:]
     MHC_lambda, MHC_A, MHC_DOS_lambda, MHC_DOS_A = [float(p[:-1]) for p in params]

#MHC_lambda, MHC_DOS_lambda

# data for figure a, DOS
dos = pd.read_csv(fig_folder+"dos_interp.txt", sep='\t', skiprows=1, names=["E", "DOS"])

# data for figure b, current / rate constant
exp_data = pd.read_csv(fig_folder+"ecdec.txt", names=["V", "I"])
fit_data = pd.read_csv(fig_folder+"ecdec_fitdata.txt", sep='\t', skiprows=1, names=["V", "MHC_k", "MHC_DOS_k"])
# reformat fit_data for seaborn
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
#a0.set_xlim(-3.5, 3.5)
a0.set_xlim([-1,1])
a0.set_xlabel("$E-E_f$ [eV]", fontsize=lb_fs)
a0.set_ylabel("DOS [arb.]", fontsize=lb_fs, labelpad=8)
a0.tick_params(labelsize=tl_fs)
a0.legend(["Li (100) DOS"], loc='lower right')
# figure b
sns.lineplot(x="V", y="k", data=fit_data, hue='fit', legend=False, linewidth=2, ax=a1)
a1.set_yscale('log')
a1.set_ylabel("Current [mA/cm$^2$]", fontsize=lb_fs)
a1.set_xlabel("Overpotential $\eta$ [V]", fontsize=lb_fs)
a1.tick_params(labelsize=tl_fs)
#a1.set_xlim([-3, 3])
a1.set_xlim([-1, 1])
a1.set_ylim([1, 1e3])
a1.scatter(exp_data.V, exp_data.I, c=[[0,0,0,0]], edgecolors='k', linewidths=1, zorder=10)
a1.legend(["MHC", "MHC+DOS", "expt."], fontsize=tl_fs-1)
plt.tight_layout()
# and the inset
# axins = a1.inset_axes(bounds=[0.62, 0.15, 0.33, 0.5])
# g = sns.lineplot(x="V", y="k", data=fit_data, hue='fit', legend=False, linewidth=1, ax=axins)
# axins.scatter(exp_data.V, exp_data.I, 15, c=[[0,0,0,0]], edgecolors='k', linewidths=0.5, zorder=10)
# axins.set_yscale('log')
# axins.set_ylabel("")
# axins.set_xlabel("")
# axins.set_xlim([-.5, .5])
# axins.yaxis.set_visible(False)
# subfigure labels
plt.text(0.02, 0.7, "a)", fontsize=14, transform=plt.gcf().transFigure, weight='bold')
plt.text(0.02, 0.1, "b)", fontsize=14, transform=plt.gcf().transFigure, weight='bold')
plt.savefig("figs/for_paper/fig1.png", dpi=300)
