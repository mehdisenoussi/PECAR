## Welcome to the PECAR repository.

This is a set of MATLAB scripts to perform the analysis on the PECAR experiment dataset. (Senoussi, Moreland, Busch & Dugué - 2018 (Submitted))
Here is a link to the preprint: [https://www.biorxiv.org/content/10.1101/443341v2](https://www.biorxiv.org/content/10.1101/443341v2)

To perform the analyses you need to download the MGL toolbox which we used to design the experiment and from which some codes are needed to "unpack" the data:
[MGL-toolbox](http://gru.stanford.edu/doku.php/mgl/download)

You should then download or clone the scripts from this github repository as well as the data from this OpenScienceFramework repository: [data](https://osf.io/2d9sc/?view_only=6ef3f85d9f944d27b23fc7af5a26f087)

You should place the mgl-master folder as well as the uncompressed "data_pecar" folder in the same directory as the scripts.

Finally you can run the scripts in that order to reproduce the figures from the paper:
1) pecar_2afc_task.m
2) pecar_p_probe_analysis_allsubjs.m
3) pecar_plot_probeDiffQuad_P1andP2_P1minusP2_freqSpectra.m
4) pecar_plot_6hz_pdiff_discr_padded.m
5) pecar_probe_report_acc.m
6) pecar_plot_probeSameQuad_P1minusP2_freqSpectra.m
7) pecar_poweranalysis_resampshuffletest.m

Todo:
	- add more comments

Thank you for your interest!

