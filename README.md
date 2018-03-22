
[![ADAM](/install/ADAM_header.png)](http://www.fahrenfort.com/ADAM.htm)

# Download
If you download this toolbox, please go through [http://www.fahrenfort.com/ADAM.htm](http://www.fahrenfort.com/ADAM.htm) where you can leave your e-mail. This way I can keep track of the user base of the toolbox and inform users of serious bugs if they happen.

# What is it?
ADAM is an open source Matlab Toolbox. It allows you to perform multivariate analyses on your EEG and/or MEG data using backward decoding (BDM) and forward encoding models (FEM).

# Features
- Perform multivariate classification analyses (backward decoding) using an arbitrary number of conditions
- Easily compute, plot and compare MVPA results side by side with ERP results
- Compute forward encoding models, including Channel Tuning Functions (CTFs), either across time or averaged over time
- Use the same data for training and testing through a k-fold leave-one-out procedure
- Use different files for training and testing
- Use different event codes for training and testing
- Gives the option to either do time-frequency (TFR) decomposition first (using Fieldtrip) performing the analysis on a range of frequency bands, or just work with the raw EEG
- The decoding can either be performed on induced TFR power or total TFR power
- Do trial binning if required
- Strictly enforces balanced designs
- Computes temporal generalization matrices (King et al, TiCS 2014) or time by frequency matrices, both for classification and for CTFs (forward modeling, e.g. see Foster et al 2016 or Samaha et al. 2016)
- Compute temporal generalization matrices for raw EEG or for any frequency (if using TFR data)
- Average over training windows, average over testing windows, average over frequency windows
- Compute, plot and do statistics on spatial topomaps for any time point based on the Haufe method (NeuroImage, 2014)
- Compute, plot and do statistics on spatial top-maps of weights obtained from the forward model
- Statistical testing, including cluster based permutation testing and FDR on the resulting temporal generalization matrices or on the 2D graphs, either per condition or between conditions
- Visualization of the results (2D graphs, 3D color scale maps, topoplots etc)
- The computation-intensive part of the toolbox has been optimized for use on UNIX computing clusters, to enable fast computation of large datasets (many subjects in parallel), but can also be ran locally on any computer with reasonable specs

# Why should I use this toolbox?
One of the big advantages of this toolbox is that it takes generic input formats for which many import functions are available (EEGLAB or FieldTrip), allowing researchers to do their own pre-processing any which way they like. The toolbox takes care of the intricacies of multivariate analyses (data handling), allowing a wealth of possibilities as specified above, and always has a group analysis as its endpoint. Although everything is scripted, the scripts are easy to use, doable also for novices.

# Requirements
- A recent version of Matlab including the Statistics Toolbox (>=2012b, lower versions might or might not work). No special care was taken to exclude functions from other Matlab Toolboxes. Let me know if you find dependencies on other Toolboxes, I might be able to provide a workaround.
- A recent version of [EEGLAB](https://sccn.ucsd.edu/eeglab/downloadtoolbox.php) (>=13, lower versions might or might not work)
- A recent install of [FieldTrip](http://www.fieldtriptoolbox.org/download) (>=2015, lower versions might or might not work)
- A reasonably modern computer (>=8GB memory, enough HD space, modern processor, more is better)

# Version
The toolbox is currently in version 1.0.0

# Cite
When you use the decoding (BDM) feature, please cite:<br>
Fahrenfort, J. J., van Leeuwen, J., Olivers, C. N. L., & Hogendoorn, H. (2017). Perceptual integration without conscious access. *Proceedings of the National Academy of Sciences*, 114(14), 3744–3749. 

When you use the forward modeling (FEM) feature, please cite:<br>
Fahrenfort, J. J., Grubert, A., Olivers, C. N. L., & Eimer, M. (2017). Multivariate EEG analyses support high-resolution tracking of feature-based attentional selection. *Scientific Reports*, 7(1), 1886.


# Manual
A citable tutorial article covering how to use the decoding apect of the toolbox will become availabel shortly.
