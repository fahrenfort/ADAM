[![ADAM](/install/ADAM_header.png)](https://github.com/fahrenfort/ADAM)

## Roadmap
###Planned features

* update core FEM functions to have a consistent API
* improved topomap plotting
* FEM tutorial
* improve statistical testing using prevalence inference, addressing problems of t-testing on MVPA data (Allefeld, C., Görgen, K., & Haynes, J.-D. (2016). Valid population inference for information-based imaging: From the second-level t-test to prevalence inference. NeuroImage, 141, 378–392.)
* implement cluster based permutation for topomaps over time (not just space)
* allow oversampling of trigger values within a class (under consideration) 

##Prior releases

###V0.0.0

* better implementation of class balancing using ADASYN/SMOTE (Chawla, N. V., Bowyer, K. W., Hall, L. O., & Kegelmeyer, W. P. (n.d.). SMOTE: Synthetic Minority Over-sampling Technique. Journal of Artificial Intelligence Research, 16, 321–357.)
* implement AUC (Area Under the Curve) as a dependent measure (right now we have mean accuracy, d' and hr-far, but AUC might be more sensitive).
* updated all core BDM functions to have a (more or less) consistent API
* all core functions are now pre-pended with adam_
* BDM tutorial (1st level, 2nd level, plotting) with example data   
