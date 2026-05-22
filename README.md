# Track-before-detect-GPP-MB-filter
This repository contains the implementation of the generalised parallel partition (labelled) multi-Bernoulli (GPP-MB) particle filter for track-before-detect described in [[1]](https://ieeexplore.ieee.org/document/7812865). This filter can work with a general track-before-detect measurement model. The shared implementation does not include the label switching algorithm to improve the multi-Bernoulli approximation at each time. This filter can form target trajectories via sequential track building (linking target estimates with the same label, though this has not been implemented in this version).

The code is written in Matlab with C-MEX files for speeding up the computation with particles. The repository contains compiled C-MEX files for Windows. For other platforms, one should compile the code following Matlab C-MEX files compilation procedures. Minor modifications to the C-MEX files may be required for new compilations.

In this demo, the measurement model is not the one in [[1]](https://ieeexplore.ieee.org/document/7812865), but the one [[2]](https://ieeexplore.ieee.org/document/10399382).

The measurement model functions have been taken from the repository from [[2]](https://ieeexplore.ieee.org/document/10399382), which is [3].

Performance is evaluated via the GOSPA metric [[4]](https://ieeexplore.ieee.org/document/8009645) with code originally provided in [5].

To the run the GPP-MB filter with informative birth run

GPP_MB_PF_informative_birth.m

To the run the GPP-MB filter with uninformative birth run

GPP_MB_PF_uninformative_birth.m


Code of the C-files is written in Spanish. Files statutils.c and matrixutils.c were written by Mark Morelande, who has given permission to include them in this repository. These files include some statistical and matrix functions implemented in C, used also in the particle filter for multiple target tracking in [[6]](https://ieeexplore.ieee.org/document/6558005). 


[[1] A. F. García-Fernández, "Track-before-detect labeled multi-Bernoulli particle filter with label switching," in IEEE Transactions on Aerospace and Electronic Systems, vol. 52, no. 5, pp. 2123-2138, October 2016, doi: 10.1109/TAES.2016.150343.](https://ieeexplore.ieee.org/document/7812865)

[[2] E. S. Davies and Á. F. García-Fernández, "Information Exchange Track-Before-Detect Multi-Bernoulli Filter for Superpositional Sensors," in IEEE Transactions on Signal Processing, vol. 72, pp. 607-621, 2024, doi: 10.1109/TSP.2024.3349769](https://ieeexplore.ieee.org/document/10399382)

[3] https://github.com/ESDavies/IPLF-IEMB-filter

[[4] A. S. Rahmathullah, A. F. Garcia-Fernandez and L. Svensson, "Generalized optimal sub-pattern assignment metric," 2017 20th International Conference on Information Fusion (Fusion), Xi'an, China, 2017, pp. 1-8, doi: 10.23919/ICIF.2017.8009645.](https://ieeexplore.ieee.org/document/8009645)

[5] https://github.com/abusajana/GOSPA

[[6] A. F. Garcia-Fernandez, J. Grajal and M. R. Morelande, "Two-Layer Particle Filter for Multiple Target Detection and Tracking," in IEEE Transactions on Aerospace and Electronic Systems, vol. 49, no. 3, pp. 1569-1588, July 2013, doi: 10.1109/TAES.2013.6558005](https://ieeexplore.ieee.org/document/6558005)
