# MHC_DOS

Herein you can find the code used for analysis and plotting for [this paper](link_eventually) (that will eventually be a link to the paper).

## Here for the data viz?
Either launch the `dataviz.ipynb` notebook yourself if you're down with that, or for an even easier (but probably slower-performing) way, use Binder!

## Usage
If you just want to reproduce the analysis done in the paper, head over to `export_figure_data.jl` and then run the `figure_*.py` scripts to plot.

If you want to play, I'd suggest you start with the interactive data visualization (see above) to get a sense for what's going on. Otherwise, here's a description of the other contents of this repository:

### Directories
* `DOSes/`: Contains DFT-computed DOSes for each material considered in our manuscript. First column is energy relative to Fermi level in eV, second is number of states
* `exp_data/`: experimental data from [this paper](https://pubs.acs.org/doi/abs/10.1021/acsenergylett.0c00031) for three electrolyte, first column is overpotential in V and second is current in mA/cm^2
* `figs/`: empty on Github, you'll have to populate it yourself by running the scripts! (You may have to manually add in some subdirectories)
* `figure_data/`: Contains all the data used to plot the figures from our manuscript (so you can run the Python scripts directly)

### Julia files
* `check_cutoff.jl`: The file we used to verify that truncating the DOS had no noticeable effect on the rate constant plots
* `compare_materials.jl`: This will generate a bunch of plots for different interfaces at different equilibrium Fermi levels, and also plot the material DOSes
* `export_figure_data.jl`: Runs the analysis done for our manuscript and saves out the data to be plotted
* `fit_exp_data.jl`: Perform fits for different electrolytes
* `functions.jl`: The workhorse file. All the functions for analysis, and a few plotting shortcuts, are defined here.
* `interactive_plots.jl`: Code for generating the interactive plots. Run by the Binder notebook referenced above, but you can run in your local IDE too!

## Dependencies
Analysis was done in Julia 1.4.2 with the dependencies specified in `Project.toml`. Plotting was done in Python 3.6 requires the following packages:
* numpy
* pandas
* matplotlib
* seaborn

## Future Plans
We hope to make the analytical tools developed for this into a proper Julia package eventually!
