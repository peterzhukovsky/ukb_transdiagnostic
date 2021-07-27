The following code below was used to analyse tabular UKB data (ran in MATLAB R2016a):

# tabular_ukb.m
This scripts includes 
a) data cleaning and creating new variables representing clinical and cognitive variables, including four non-overlapping case groups
b) permutation testing general linear models testing for the effects of clinical groups (F-test) and running case-control post hoc tests (T-statistic) for 
- major depression: MDD- ('dep')
- non phobic anxiety disorders: ANX- ('anx')
- comorbid major depression with anxiety: MDD+ANX ('depanx') 
- stress-related disorders: STR ('str')

Further, effects of MDD, ANX and PTSD polygenic risk scores are tested. 
The resulting brain maps are correlated 
c) general linear models (uncorrected p-values) for cognitive data from the UKB cognitive battery
d) demographic table
e) visualizations and tables based on the above results 

# tabular_pls.m
This script includes the following steps:
a) using the same code as in # tabular_ukb.m # we create groups and clean cognitive data, excluding missing data
b) runs a partial least squares regression to predict performance on four cognitive tests (trailmaking, digit-symbol substitution, paired associates learning and fluid intelligence) from partial resting state FC
We run the PLS across all cases together, in each case group separately and in healthy controls separately
Permutation testing is used to test for significance of the % of variance explained by the PLS model in the cognitive outcome variables (Y)
c) generates bootstrapped loadings following Dr Sarah Morgan's scripts:

https://www.pnas.org/content/116/19/9604

https://github.com/SarahMorgan/Morphometric_Similarity_SZ/blob/master/Gene_analyses.md

d) plots the results in several ways

# network_viz.m
This script was used to generate the colormaps for cortical thickness brain maps. Those colormaps were used with freesurfer to generate a parcellation (aparc) similar to https://github.com/peterzhukovsky/brain_ageing

The code requires:
MATLAB Statistics toolbox
MATLAB Bioinformatics toolbox
Visualizations rely on the circulargraph function:
https://uk.mathworks.com/matlabcentral/fileexchange/48576-circulargraph
allcomb function
https://uk.mathworks.com/matlabcentral/fileexchange/10064-allcomb-varargin
We also used plot_arc in creating the circular graph with the networks:
https://uk.mathworks.com/matlabcentral/answers/6322-drawing-a-segment-of-a-circle

