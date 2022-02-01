# Introduction
This directory contains data described in 
*Primary visual cortex straightens natural video trajectories* by Olivier J. Hénaff\*, Yoon Bai\*, Julie Charlton, Ian Nauhaus, Eero P. Simoncelli, and Robbe L. T. Goris (2021).  

# Contents
- curvatureEstimations/: this folder contains the results from curvature estimations (recovery analysis). Results are saved in two formats: (a) CSV spreadsheet and (b) MatLab table. Each row or entry is designated to a single V1 response trajectory elicited from a set of movie frames.
- modelResults/: this folder contains model responses from stimulus images used in V1 recordings. We built an image computable model using a cascade of two linear-nonlinear operations (LN-LN) and saved model responses in a MatLab file, 'aggregated_LNLN_nat_sigG2_sqrt.mat'.
- V1_responses/: this folder contains spike count data recorded from the primary visual cortex (V1) in anesthetized macaques. Details of the data structure are described in the example MatLab script, 'summary.m'.
- stimulus/: this folder contains meta information (stim_info.mat) and pixel values of each image used in the experiment (stim_matrix.mat).
- library/: this folder contains auxiliary code used in the example MatLab script, 'summary.m'.
- figures/: example figures generated from 'summary.m' are saved in this folder.
- summary.m: this MatLab script aggregates all results (curvature estimations & model responses) into a single MatLab table.
- plot_summary.m: MatLab script for plotting summary results. Copy the contents of this folder on your computer. Run 'plot_summary.m' to generate summary figures reported in this project.

