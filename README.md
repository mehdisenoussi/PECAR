This is a set of MATLAB functions to perform the analysis on the PECAR experiment dataset. (Senoussi, Moreland, Busch & Dugué - 2018 (Submitted))

To perform the analysis you will need to download the MGL toolbox which we used to design the experiment and from which some codes are needed to "unpack" the data.
[MGL-toolbox](http://gru.stanford.edu/doku.php/mgl/download)

You should then download or clone the content of this repository as well as the raw data from [OpenScienceFramework](https://osf.io/2d9sc/?view_only=658a434a48c04ba590cdf1a540cf30dd)

You should place the mgl-master folder as well as the uncompressed "data_pecar" folder in the same directory as the scripts.

Finally you can run the scripts in that order to reproduce the figures from the paper:
- pecar_2afc_task.m
- pecar_p_probe_analysis_allsubjs.m
- pecar_plot_P1andP2_P1minusP2_freqAmp.m
- pecar_plot_6hz_pdiff_discr_padded.m


TODO:
- Clean and add the scripts to generate the Power analysis
