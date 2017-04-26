
[![ADAM](http://www.fahrenfort.com/images/ADAM_logo.png)](https://github.com/fahrenfort/ADAM)

# Features
- Use the same data for training and testing through a k-fold leave-one-out procedure
- Use different datasets for testing and training
- Gives the option to either do time-frequency (TFR) decomposition first (using Fieldtrip) performing the analysis on a range of frequency bands, or just work with the raw EEG
- The decoding can either be performed on induced TFR power or total TFR power
- Do trial binning if required
- Compute Generalization Across Time (GAT) matrices (King et al, TiCS 2014) or time by frequency matrices, both for classification and for cortical tuning functions (forward modeling, e.g. see Foster et al 2016 or Samaha et al. 2016)
- Compute GAT matrices for raw EEG or for any frequency (if using TFR data)
- Compute, plot and do statistics on spatial topo-maps for any time point based on the Haufe method (NeuroImage, 2014)
- Compute, plot and do statistics on spatial top-maps of weights obtained from the forward model
- Statistical testing, including cluster based permutation testing and FDR on the resulting GAT matrices or on the 2D vectors, either per condition or between conditions
- Visualization of the results (2D graphs, 3D color scale maps, topoplots etc)
- The computation-intensive part of the toolbox has been optimized for use on UNIX computing clusters, to enable fast computation of large datasets (many subjects in parallel), but can also be ran locally on any computer with reasonable specs

