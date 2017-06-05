
[![ADAM](/install/ADAM_header.png)](https://github.com/fahrenfort/ADAM)

# What is it?
ADAM is an open source Matlab Toolbox. It allows you to perform multivariate analyses on your EEG and/or MEG data using backward decoding (BDM) and forward encoding models (FEM).

# Features
- Perform multivariate classification analyses (backward decoding) using an arbitrary number of conditions
- Compute forward encoding models, including Channel Tuning Functions (CTFs), either across time or averaged over time
- Use the same data for training and testing through a k-fold leave-one-out procedure
- Use different datasets for testing and training
- Gives the option to either do time-frequency (TFR) decomposition first (using Fieldtrip) performing the analysis on a range of frequency bands, or just work with the raw EEG
- The decoding can either be performed on induced TFR power or total TFR power
- Do trial binning if required
- Strictly enforces balanced designs
- Compute Generalization Across Time (GAT) matrices (King et al, TiCS 2014) or time by frequency matrices, both for classification and for CTFs (forward modeling, e.g. see Foster et al 2016 or Samaha et al. 2016)
- Compute GAT matrices for raw EEG or for any frequency (if using TFR data)
- Average over training windows, average over testing windows, average over frequency windows
- Compute, plot and do statistics on spatial topomaps for any time point based on the Haufe method (NeuroImage, 2014)
- Compute, plot and do statistics on spatial top-maps of weights obtained from the forward model
- Statistical testing, including cluster based permutation testing and FDR on the resulting GAT matrices or on the 2D vectors, either per condition or between conditions
- Visualization of the results (2D graphs, 3D color scale maps, topoplots etc)
- The computation-intensive part of the toolbox has been optimized for use on UNIX computing clusters, to enable fast computation of large datasets (many subjects in parallel), but can also be ran locally on any computer with reasonable specs

# Why should I use this toolbox?
One of the big advantages of this toolbox is that it takes generic input formats for which many import functions are available (EEGLAB or FieldTrip), allowing researchers to do their own pre-processing any which way they like. The toolbox takes care of the intricacies of multivariate analyses (data handling), allowing a wealth of possibilities as specified above, and always has a group analysis as its endpoint. Although everything is scripted, the scripts are easy to use, doable also for novices.

# Requirements
- A recent version of Matlab including the Statistics Toolbox (>2012b, lower versions might or might not work). No special care was taken to exclude functions from other Matlab Toolboxes. Let me know if you find dependencies on other Toolboxes, I might be able to provide a workaround.
- A recent version of [EEGLAB](https://sccn.ucsd.edu/eeglab/downloadtoolbox.php) (>13, lower versions might or might not work)
- A recent install of [FieldTrip](http://www.fieldtriptoolbox.org/download) (>2015, lower versions might or might not work)
- A reasonably modern computer (>=8GB memory, enough HD space, modern processor, more is better)

# Version
The software is currently in pre-release, meaning that it is available for download but important features are still under development. We are still streamlining the basic analysis pipeline, to make it more consistent and more usable. **If you download the software now, expect no support, but do expect unexpected changes during updates.** We foresee that a first 'official' release will become available within the next few months.

# Cite
When you use the decoding (BDM) feature, please cite:<br>
Fahrenfort, J. J., van Leeuwen, J., Olivers, C. N. L., & Hogendoorn, H. (2017). Perceptual integration without conscious access. *Proceedings of the National Academy of Sciences*, 114(14), 3744–3749. 

When you use the forward modeling (FEM) feature, please cite:<br>
Fahrenfort, J. J., Grubert, A., Olivers, C. N. L., & Eimer, M. (2017). Multivariate EEG analyses support high-resolution tracking of feature-based attentional selection. *Scientific Reports*, 7(1), 1886.


# Manual
We are currently working to supply an up-to-date manual of the toolbox, including some easy to use wrapper functions, which will be presented at the [ICON conference](http://www.icon2017.org) 5-8 August 2017 in Amsterdam
