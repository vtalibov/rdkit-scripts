# RDKit scripts

Collection of basic RDKit-dependent Python scripts.

* pySMCM - calculates synthetic molecular complexity value
* pyPMIscore - calculates average sum of normalized principal moments of inertia for a set of low energy conformers
* deltaE - IPython Notebook with some analysis of difference in minimized energies between a protein-bound conformation of a fragment-alike small molecule and its conformer with minimal energy. Was used to approximate energy cut-off for metrics whenever multiple conformations are considered (*e.g. pyPMIscore). Has a snippet to plot normalized PMIs for a set of compounds with multiple conformers.
