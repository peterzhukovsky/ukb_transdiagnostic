The following code was used to analyse tabular UKB data (ran in MATLAB R2016a, on windows 8.1):

# tabular_ukb_RMatch1to1.m
This scripts includes 
a) data cleaning and creating new variables representing clinical and cognitive variables, including four non-overlapping case groups. Matching controls to cases in terms of age and sex.

b) permutation testing general linear models testing for the effects of clinical groups (F-test) and running case-control post hoc tests (T-statistic) for 

- healthy controls: ('ahc')
- major depression: MDD- ('dep')
- non phobic anxiety disorders: ANX- ('anx')
- comorbid major depression with anxiety: MDD+ANX ('depanx') 
- stress-related disorders: STR- ('str')

Further, effects of MDD, ANX and PTSD polygenic risk scores are tested. 
The resulting brain maps are correlated and the correlations are also permutation tested

c) general linear models (uncorrected p-values) for cognitive data from the UKB cognitive battery

d) demographic table

e) visualizations and tables based on the above results 

We repeat the main analyses several times. First, we repeat the main analysis in a subgroup of White British participants while a) not covarying for polygenic risk scores and b) while covarying for polygenic risk scores. > **tabular_ukb_RMatch1to1_PRS.m**

Second, we repeat the main analysis in unmedicated participants only. > **tabular_ukb_RNoMed.m**

Third, we repeat the main analysis in participants with active MDD (PHQ2>=2) only. > **tabular_ukb_RActiveD.m**

Fourth, we also repeat the main analysis by combining all case groups together to generate case-control maps across all case participants. > **tabular_ukb_R2Group.m**

An analysis of all cases vs controls can be found in **tabular_ukb.m**

# tabular_pls.m
This script includes the following steps:

a) using the same code as in # tabular_ukb.m # we create groups and clean cognitive data, excluding missing data

b) runs a partial least squares regression to predict performance on four cognitive tests (trailmaking, digit-symbol substitution, paired associates learning and fluid intelligence) from partial resting state FC
We run the PLS across all cases together, in each case group separately and in healthy controls separately
Permutation testing is used to test for significance of the % of variance explained by the PLS model in the cognitive outcome variables (Y)

c) generates bootstrapped loadings following Dr Sarah Morgan's scripts:

https://www.pnas.org/content/116/19/9604

https://github.com/SarahMorgan/Morphometric_Similarity_SZ/blob/master/Gene_analyses.md

d) plots the results (specifically the pairwise connectivity loadings contributing to the latent PLS variables) using matrices and chord plots ("wheel-like" graphs). Note that we rearranged the chort plots ordering the ICs according to the Yeo network they best match with.

e) compares PLS results between MDD- and ANX- using Pearson's correlations. 

# network_viz.m
This script was used to generate the colormaps for cortical thickness brain maps. Those colormaps were used with freesurfer to generate a freesurfer parcellation in  Glaser 360 ROI space, similar to https://github.com/peterzhukovsky/brain_ageing

# UKB_ICAd25 

This folder contains the 21 signal components from the 25-dimensional ICA analysis, provided by the UKB: https://www.fmrib.ox.ac.uk/datasets/ukbiobank/index.html

*d25_to_Yeo7_mapping.sh* - this script does the following:

a) maps the components from the volumetric MNI space to the fsaverage surface space (mgh surface maps of the ICA components are also shared in this folder)

b) and checks (in the volumetric MNI 1mm space) what percentage of voxels in each of the components fall under each of the Yeo 7 networks. Yeo 7 networks in volume space are available from https://surfer.nmr.mgh.harvard.edu/fswiki/CorticalParcellation_Yeo2011

# Sample data

We provide random data (n=1000), based on real data from the UKB, which can be used to try out the scripts. The data is included in */random_data/randomised_data.mat*. The randomisation script is also included as *shuffling_script.m* and the *shuffle.m* function. Note that we don't provide medication data as it contains tables of medication regimens. Random data is saved as a *.mat* file that contains among other variables the names of the 360 Glasser regions (and their mapping on to the Harvard-Oxford volumetric atlas, done using fsleyes). We also provide the *ica2yeo7.csv* table used for creating the circular graph and for labelling the connectivity matrices.

# Additional dependencies and small functions:

The code relies on:

MATLAB Statistics toolbox

MATLAB Bioinformatics toolbox

MATLAB Parallization toobox

Visualizations rely on the *circulargraph* function:

https://uk.mathworks.com/matlabcentral/fileexchange/48576-circulargraph

and the *allcomb* function

https://uk.mathworks.com/matlabcentral/fileexchange/10064-allcomb-varargin

We also used plot_arc in creating the circular graph with the networks:

https://uk.mathworks.com/matlabcentral/answers/6322-drawing-a-segment-of-a-circle

This function was used in *drawingcircle.m*

*drawingcircle.m* - a script to generate the colorful outside rim of the circular graphs

*mktbl.m* - a function to generate a table for use with *fitlm* (needed for parallelization in R2016a)

*clean.m* - a function for removing outliers (takes the number of standard deviations from the mean that are considered bounds for outlier detection)

# PRS calculations

PRS calculations followed Wainberg et al 2021: https://pubmed.ncbi.nlm.nih.gov/34008483/
A script used by Dr Michael Wainberg to create the polygenic risk scores can be found at */PRS/PRS_MW.py*

# Runtime

Permutations are notoriously slow as they require running the analysis of interest many time over and over. On a quad code intel i7 CPU with 16 GB RAM, running GLMs with apppx 27,000 participants a permutation block (of n=1000 permutations) should take 2.5-3.5h when using parfor parallelization. Using only one CPU this time would approximately quadruple. 
Other code is relatively quick to run even with larger samples.

A preprint with more information can be found at https://www.researchsquare.com/article/rs-711822/v1
