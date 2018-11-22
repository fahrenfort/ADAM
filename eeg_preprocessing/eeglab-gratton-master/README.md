# Ocular Correction EEGlab Plugin

*WARNING:* This is an old plugin I wrote in 2007 and never updated to reflect changes in EEGlab or MATLAB. PR's to fix this are welcome.

## About

A plugin for the EEGlab Toolbox for MATLAB, implementing the
regression-based algorithm for ocular correction of EEG-data
proposed by Gratton et al. 1982.


    @ARTICLE{Gratton1983,
       author = {Gratton, Gabriele and Coles, Michael G. H. and Donchin, Emanuel},
       title = {A new method for off-line removal of ocular artifact},
       journal = {Electroencephalography and Clinical Neurophysiology},
       year = {1983},
       volume = {55},
       pages = {468--484},
       number = {4}
    }

The algorithm estimates different regression parameters for
blink-segments and non-blink-segments. Therefore, two parameters have
to be specified to distinguish blink vs. non-blink data:
* `blinkcritwin` time window for criterion (in ms)  
* `blinkcritvolt` amplitude criterion for blink detection:
       (EOG(t) - EOG(t-win)) + (EOG(t) - EOG(t+win)) >= crit


## Installation

Copy all m-files into a subdirectory under the 'plugins' path of your
EEGlab-installation. Then restart EEGlab and a menu-item in the
'Tools' menu should appear. For more details about EEGlab-Plugins, refer to
http://www.sccn.ucsd.edu/eeglab/.
