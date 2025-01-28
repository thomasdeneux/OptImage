
# OptImage

OptImage is a Matlab program to load, display, and analyze movies from optical
imaging (wide-field imaging, two-photon imaging, etc.). In particular, it 
allows fast navigation throughout the multi-condition, multi-trial,
multi-session data, and provides tools for standards analyses such as baseline
subtraction, dF/F calculation, ROI analysis, artefact removal (drift, 
heartbeat, bleaching, etc.), standard filters, etc.


## Dependency

OptImage requires the [io](https://github.com/thomasdeneux/io) library to read 
data files of some specific formats (abf, cfd, elphy, blk, mesc, neuroplex, 
scanimage).


## Installation

Add Optimage's "matlab" folder to the Matlab path, as well as the dependent 
"io" folder.

Then you can run the "optimage.m" file to start the program.

