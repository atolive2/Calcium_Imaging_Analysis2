%% Correlations Version 2: calculate correlation coefficients based on random data about ROIs using corrcoef
% R = corrcoef(A) returns the matrix of correlation coefficients for A, 
% where the columns of A represent random variables and the rows represent observations.
% https://www.mathworks.com/help/matlab/ref/corrcoef.html

% all of this is copied from correlations_bymodality and
% cluster_analysis_v1 but using random variable data instead of time series

% this code starts with tadcluster_data_20170814, with all fields after
% resp_ROIs deleted.