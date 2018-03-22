[![ADAM](/install/ADAM_header.png)](https://github.com/fahrenfort/ADAM)

## Roadmap
###Planned features

* update core FEM functions to have a consistent API
* improved topomap plotting
* FEM tutorial
* improve statistical testing using prevalence inference, addressing problems of t-testing on MVPA data (Allefeld, C., Görgen, K., & Haynes, J.-D. (2016). Valid population inference for information-based imaging: From the second-level t-test to prevalence inference. NeuroImage, 141, 378–392.)
* implement cluster based permutation for topomaps over time (not just space)
* allow oversampling of trigger values within a class (under consideration) 

##V1.0.0 is the current version

###Implemented prior to V1.0.0 
* whitening (by default on for FEMs, not for BDMs as these use LDA)
* better implementation of class balancing using ADASYN to oversample (Haibo He, Yang Bai, Garcia, E. A., & Shutao Li. (2008). ADASYN: Adaptive synthetic sampling approach for imbalanced learning (pp. 1322–1328). Presented at the 2008 IEEE International Joint Conference on Neural Networks (IJCNN 2008 - Hong Kong), IEEE.)
* implemented AUC (Area Under the Curve) as the default performance measure (on top of balanced accuracy, d' and hr-far).
* updated all core BDM functions to have a (more or less) consistent API
* all core functions are now pre-pended with adam_
* BDM tutorial (1st level, 2nd level, plotting) with example data   
