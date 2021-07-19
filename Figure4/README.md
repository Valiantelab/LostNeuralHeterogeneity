Compile and run Figure4_100trials.cpp with inputs of the desired sigma_e and sigma_i values. Number of trials can be changed on line 41 (#define Trials).

Matlab file Plot100TrialAvgs.m will generate figures akin the top row in each panel of Figure 4. Set target=1 (this is an artifact from previous code). Filenames should refer to output from cpp code for excitatory and inhibitory spikes.

Compile and run Figure4_FindFixedPoints.c and choose sigma_e and sigma_i on Lines 31 and 32.

Matlab file PlotEigenvaluesAndFixedPoints.m will generate figures akin to the bottom rows in each panel of Figure 4. Input 'filename' should be the output from Figure4_FindFixedPoints.c, and remaining variables are used in naming the output and by convention should be analogous to what's in the 'filename'.

Remaining Matlab files are various plotting functions that may or may not be used in these sets of code (they are the "library" I keep for all my Matlab plotting code, just in case).

