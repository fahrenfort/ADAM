
[![ADAM](/install/ADAM_header.png)](http://www.fahrenfort.com/ADAM.htm)

# Download
If you download this toolbox, please go through [http://www.fahrenfort.com/ADAM.htm](http://www.fahrenfort.com/ADAM.htm) where you can leave your e-mail. This way I can keep track of the user base of the toolbox and inform users of serious bugs if they happen.

# What is it?
ADAM is an open source Matlab Toolbox. It allows you to perform multivariate analyses on your EEG and/or MEG data using backward decoding (BDM) and forward encoding models (FEM).

# Features
- Perform multivariate classification analyses (backward decoding) using any number of conditions
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
- A recent version of Matlab (>=2012b, lower versions might or might not work) including the "Image processing toolbox", "Signal processing toolbox" and "Statistics toolbox". No special care was taken to exclude functions from other Matlab Toolboxes other than EEGLAB and FieldTrip. Let me know if you find dependencies on other Toolboxes, I might be able to provide a workaround.
- [EEGLAB](https://sccn.ucsd.edu/eeglab/downloadtoolbox.php) (>=13, lower versions might or might not work, very recent installs may or may not work. Currently using eeglab14_1_1b)
- [FieldTrip](http://www.fieldtriptoolbox.org/download) (>=2015, lower versions might or might not work, very recent installs may are may not work. Currently using fieldtrip-20220925)
- A reasonably modern computer (>=8GB memory, enough HD space, modern processor, more is better)

# Version
The toolbox is currently in version 1.x.x continuous beta

# Cite
When you use the decoding (BDM) feature, please cite:<br>
Fahrenfort, J. J., van Driel, J., van Gaal, S., & Olivers, C. N. L. (2018). From ERPs to MVPA Using the Amsterdam Decoding and Modeling Toolbox (ADAM). *Frontiers in Neuroscience*, 12. http://doi.org/10.3389/fnins.2018.00368

When you use the forward modeling (FEM) feature, please cite:<br>
Fahrenfort, J. J. (2020) Multivariate methods to track the spatiotemporal profile of feature-based attentional selection using EEG.  *Pollmann (Ed.), Spatial learning and attention guidance.* Neuromethods. New York, Springer. https://osf.io/srmt2/


# Manuals
A citable tutorial article covering how to use the decoding features of the toolbox can be found here: http://doi.org/10.3389/fnins.2018.00368.<br>
A citable tutorial article covering how to use the forward encoding features of the toolbox can be found here: https://osf.io/srmt2/
