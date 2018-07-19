## Welcome to the PECAR repository.

This is a set of MATLAB functions to perform the analysis on the PECAR experiment dataset. (Senoussi, Moreland, Busch & Dugué - 2018 (Submitted))

To perform the analysis you need to download the MGL toolbox which we used to design the experiment and from which some codes are needed to "unpack" the data:
[MGL-toolbox](http://gru.stanford.edu/doku.php/mgl/download)

You should then download or clone the scripts from this github repository as well as the data from this OpenScienceFramework repository: [data](https://osf.io/2d9sc/?view_only=658a434a48c04ba590cdf1a540cf30dd)

You should place the mgl-master folder as well as the uncompressed "data_pecar" folder in the same directory as the scripts.

Finally you can run the scripts in that order to reproduce the figures from the paper:
1. **pecar_2afc_task.m** (to plot the results of the 2AFC orientation discrimination task (primary task))
2. **pecar_p_probe_analysis_allsubjs.m** (to compute and save the results of the probe task (secondary task) using the P1-P2 estimation technique and compute surrogates)
3. **pecar_plot_P1andP2_P1minusP2_freqAmp.m** (to plot the P1-P2 results (P1 and P2 (in each condition), P1 minus P2 and FFT of P1 minus P2 + surrogate distribution))
4. **pecar_plot_6hz_pdiff_discr_padded.m** (to plot the 6hz effect using both measure: P1 minus P2 and Discriminant)
5. **pecar_datareformat_for_power_analysis.m** (to reformat the data in order to compute the power of P1-P2 and Discriminant analyses in each condition)
6. **pecar_poweranalysis_resampshuffletest.m** (to perform the bootstrap of the power analysis and plot its results)


## Design of the experiment.
<img src="/pecar_exp_design.png" width="400">

## TODO:
- Add more comments
- Add figures of the experimental design and the results of each script

