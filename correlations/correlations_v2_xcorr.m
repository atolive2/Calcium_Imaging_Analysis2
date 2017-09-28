%% Correlations Version 2: cross correlation of time sequences using xcorr
% r = xcorr(x,y) returns the cross-correlation of two discrete-time sequences,
% x and y. Cross-correlation measures the similarity between x and shifted
% (lagged) copies of y as a function of the lag.
% https://www.mathworks.com/help/signal/ref/xcorr.html

% all of this is copied from correlations_bymodality and
% cluster_analysis_v1 but using xcorr function instead of corrcoef

% this code starts with tadcluster_data_20170814, with all fields after
% resp_ROIs deleted.
% tadcluster{1,t}.alldff0 contains time series data for all trials
% tadcluster{1,t}.dff0_bystimtype contains time series data split by trial
% type for high/high combos only

% graphs are saved in C:\Users\Torrey\Desktop\correlations_v2_xcorr_figs

% this is responding ROIs only. Not all experiments have enough ROIs to be
% included (min(ROI count) = 18 for inclusion).
% response = peak is after stim onset (20 frames is a generous buffer), greater than 0.2 df/f0
    %Criteria for elimination (e.g. these represent "bad" trials):
    % -	Peak location is earlier than the stimulus onset
    % -	Area is negative
    % -	Peak is negative
    % - peak is greater than 10
% included in resp_ROI = at least 1 response as defined above.

%% If beginning with basic tadpole data

% import data:
myFolder = 'D:\\Torrey_calcium_imaging\newData_20170924' %'F:/Calcium_Imaging_Analysis/analyzed_compiled/Smoothed_analysis/'; % May need to correct this.
%mkdir([myFolder 'figures']);
if ~isdir(myFolder)
	errorMessage = sprintf('Error: The following folder does not exist:\n%s', myFolder);
	uiwait(warndlg(errorMessage));
	return;
end
filePattern = fullfile(myFolder, 'tadpole*.mat');
matFiles = dir(filePattern)

for k = 1:length(matFiles)
	matFilename = fullfile(myFolder, matFiles(k).name)
	matData = load(matFilename); % Retrieves a structure.
    
	% See if tadpole actually exists in the data structure.
	hasField = isfield(matData, 'tadpole');
	if ~hasField
		% Alert user in popup window.
		warningMessage = sprintf('tadpole is not in %s\n', matFilename);
		uiwait(warndlg(warningMessage));
		% It's not there, so skip the rest of the loop.
		continue; % Go to the next iteration.
	end
	tadcluster{1,k} = matData.tadpole % If you get to here, tadpole existed in the file.
    tadcluster{1,k}.matFilename = matFilename
end

% recombine trial by trial df/f0 into 1 long vector for each ROI
for t = 1:length(tadcluster)
    for i = 1:size(tadcluster{1,t}.smoothed, 1)
        tmp = [];
        for j = 1:size(tadcluster{1,t}.smoothed, 2)
            tmp = [tmp; tadcluster{1,t}.smoothed{i,j}];
        end
        tadcluster{1,t}.alldff0(i,:) = tmp';
    end
end

% recombine each modality separately
%(only use 1-4 here to eliminate variation from stim strength)
tmp = [];
for t = 1:length(tadcluster)
    for s = 1:4
        trials_touse = find(tadcluster{1,t}.stimorder == s)
       %tmp = [];
        
            for j = 1:size(tadcluster{1,t}.df_f0,1) % over each ROI
                for i = 1:length(trials_touse) %over all trials of 1 stim type
                    tmp = [tmp; tadcluster{1,t}.df_f0{j, trials_touse(i)}];
                end
                tadcluster{1,t}.dff0_bystimtype{j, s} = tmp;
                tmp = [];
            end
    end
end

% define responding ROIs
for t = 1:length(tadcluster)
    tadcluster{1,t}.resp_ROIs = find(tadcluster{1,t}.sum_responses_sm > 0);
end

% make ROIcenters, which puts the ROI center info into a roi_num x 2 matrix
for t = 1:length(tadcluster)
    for r = 1:length(tadcluster{1,t}.somaticROICenters)
    tadcluster{1,t}.ROIcenters(r,:) = tadcluster{1,t}.somaticROICenters{1,r}(1).Centroid;
    end
end

%% calculate correlation coefficient for responding ROIs using all trials
% maxlag defaults to 2N-1, with N = greater of the lengths of x and y
% therefore, to prevent crashing Matlab, set maxlag. I will set max lag to
% the length of one trial (so the maximum tested lag is that response 1
% starts at the beginning of the trial and response y starts at the end of
% the trial. 
% chose 'coeff' as the normalization option because it normalizes the
% sequence so that the autocorrelations at 0 lag equal 1.

%respROIdff0_R has dims (lag*2) - 1 by length(resp_ROIs)^2
for t = 1:length(tadcluster)
    set_lag = length(tadcluster{1,t}.df_f0{1,1});
    [tadcluster{1,t}.respROIdff0_R, tadcluster{1,t}.respROIdff0_lag] = xcorr(tadcluster{1,t}.alldff0(tadcluster{1,t}.resp_ROIs,:)', set_lag, 'coeff');
end

% what is the max correlation between 2 ROIs?
% this is the max of each column, and all columns that are an
% autocorrelation will have max = 1
for t = 1:length(tadcluster)
    [tadcluster{1,t}.respROIdff0_maxR(1,:), tadcluster{1,t}.respROIdff0_maxR(2,:)] = max(tadcluster{1,t}.respROIdff0_R);
end

%reshape the array to be size(length(resp_ROIs)) e.g. make it a square
%again for plotting purposes. 
%for reference, xcorr arranges the cols as 1-1, 1-2, 1-3 ... 1-n, 2-1, 2-2,
%etc
for t = 1:length(tadcluster)
    len = length(tadcluster{1,t}.resp_ROIs);
    for r = 1:len
        for c = 1:len
            idx = (r-1)*len + c;
            tadcluster{1,t}.respROIdff0_maxR_sq(r,c) = tadcluster{1,t}.respROIdff0_maxR(1,idx);
            tadcluster{1,t}.respROIdff0_maxRlag_sq(r,c) = tadcluster{1,t}.respROIdff0_lag(1,tadcluster{1,t}.respROIdff0_maxR(2,idx));
        end
    end
end

% Plot correlation coefficients (R) as a matrix
for t = 1:length(tadcluster)
    if length(tadcluster{1,t}.resp_ROIs) > 0
        figure;
        colormap('hot')
        colormap(flipud(colormap))
        imagesc(tadcluster{1,t}.respROIdff0_maxR_sq)
        colorbar
        title(sprintf('tad %d responding cells xcorr (maxR) using all df/f0', t))
        fig_filename = sprintf('tad %d Responding ROIs xcorr maxR alldff0', t);
        saveas(gcf,fig_filename,'png');
        close;
    end
end

% Plot lag of max R as a matrix
% this version plots the absolute value of the lag so I get mirror images
% around the diagonal. 
for t = 1:length(tadcluster)
    if length(tadcluster{1,t}.resp_ROIs) > 0
        %to_plot = abs
        figure;
        colormap('jet')
        %colormap(flipud(colormap))
        imagesc(abs(tadcluster{1,t}.respROIdff0_maxRlag_sq))
        colorbar
        title(sprintf('tad %d responding cells xcorr lag of MaxR using all df/f0', t))
        fig_filename = sprintf('tad %d Responding ROIs xcorr maxR lags', t);
        saveas(gcf,fig_filename,'png');
        close;
    end
end

%% Reapply above analysis to each stimulus modality separately. 

% make dff0_multi, vis, mech, ns into right shape with only respROIs
for t = 1:length(tadcluster)
    for r = 1:length(tadcluster{1,t}.resp_ROIs)
        tadcluster{1,t}.dff0_multi(r, :) = reshape(tadcluster{1,t}.dff0_bystimtype{tadcluster{1,t}.resp_ROIs(r),1},[], 1);
        tadcluster{1,t}.dff0_vis(r, :) = reshape(tadcluster{1,t}.dff0_bystimtype{tadcluster{1,t}.resp_ROIs(r),2},[], 1);
        tadcluster{1,t}.dff0_mech(r, :) = reshape(tadcluster{1,t}.dff0_bystimtype{tadcluster{1,t}.resp_ROIs(r),3},[], 1);
        tadcluster{1,t}.dff0_none(r, :) = reshape(tadcluster{1,t}.dff0_bystimtype{tadcluster{1,t}.resp_ROIs(r),4},[], 1);

    end
end


for t = 1:length(tadcluster)
    if sum(size(tadcluster{1,t}.dff0_bystimtype{1,1})) > 0
        set_lag = length(tadcluster{1,t}.smoothed{1,1});
        [tadcluster{1,t}.respROIdff0_R_MS, tadcluster{1,t}.respROIdff0_lag_MS] = xcorr(tadcluster{1,t}.dff0_multi', set_lag, 'coeff');
        [tadcluster{1,t}.respROIdff0_R_V, tadcluster{1,t}.respROIdff0_lag_V] = xcorr(tadcluster{1,t}.dff0_vis', set_lag, 'coeff');
        [tadcluster{1,t}.respROIdff0_R_M, tadcluster{1,t}.respROIdff0_lag_M] = xcorr(tadcluster{1,t}.dff0_mech', set_lag, 'coeff');
        [tadcluster{1,t}.respROIdff0_R_N, tadcluster{1,t}.respROIdff0_lag_N] = xcorr(tadcluster{1,t}.dff0_none', set_lag, 'coeff');
    end
end

% what is the max correlation between 2 ROIs?
% this is the max of each column, and all columns that are an
% autocorrelation will have max = 1
for t = 1:length(tadcluster)
    if sum(size(tadcluster{1,t}.dff0_bystimtype{1,1})) > 0
    [tadcluster{1,t}.respROIdff0_maxR_MS(1,:), tadcluster{1,t}.respROIdff0_maxR_MS(2,:)] = max(tadcluster{1,t}.respROIdff0_R_MS);
    [tadcluster{1,t}.respROIdff0_maxR_V(1,:), tadcluster{1,t}.respROIdff0_maxR_V(2,:)] = max(tadcluster{1,t}.respROIdff0_R_V);
    [tadcluster{1,t}.respROIdff0_maxR_M(1,:), tadcluster{1,t}.respROIdff0_maxR_M(2,:)] = max(tadcluster{1,t}.respROIdff0_R_M);
    [tadcluster{1,t}.respROIdff0_maxR_N(1,:), tadcluster{1,t}.respROIdff0_maxR_N(2,:)] = max(tadcluster{1,t}.respROIdff0_R_N);
    end
end

%reshape the array to be size(length(resp_ROIs)) e.g. make it a square
%again for plotting purposes. 
%for reference, xcorr arranges the cols as 1-1, 1-2, 1-3 ... 1-n, 2-1, 2-2,
%etc
for t = 1:length(tadcluster)
    if sum(size(tadcluster{1,t}.dff0_bystimtype{1,1})) > 0
        len = length(tadcluster{1,t}.resp_ROIs);
        for r = 1:len
            for c = 1:len
                idx = (r-1)*len + c;
                % max R value
                tadcluster{1,t}.respROIdff0_maxR_sq_MS(r,c) = tadcluster{1,t}.respROIdff0_maxR_MS(1,idx);
                tadcluster{1,t}.respROIdff0_maxR_sq_V(r,c) = tadcluster{1,t}.respROIdff0_maxR_V(1,idx);
                tadcluster{1,t}.respROIdff0_maxR_sq_M(r,c) = tadcluster{1,t}.respROIdff0_maxR_M(1,idx);
                tadcluster{1,t}.respROIdff0_maxR_sq_N(r,c) = tadcluster{1,t}.respROIdff0_maxR_N(1,idx);
                % lag time of max R val
                tadcluster{1,t}.respROIdff0_maxRlag_sq_MS(r,c) = tadcluster{1,t}.respROIdff0_lag_MS(1,tadcluster{1,t}.respROIdff0_maxR_MS(2,idx));
                tadcluster{1,t}.respROIdff0_maxRlag_sq_V(r,c) = tadcluster{1,t}.respROIdff0_lag_V(1,tadcluster{1,t}.respROIdff0_maxR_V(2,idx));
                tadcluster{1,t}.respROIdff0_maxRlag_sq_M(r,c) = tadcluster{1,t}.respROIdff0_lag_M(1,tadcluster{1,t}.respROIdff0_maxR_M(2,idx));
                tadcluster{1,t}.respROIdff0_maxRlag_sq_N(r,c) = tadcluster{1,t}.respROIdff0_lag_N(1,tadcluster{1,t}.respROIdff0_maxR_N(2,idx));
            end
        end
    end
end


% Plot correlation coefficients (R) as a matrix
%multisensory
for t = 1:length(tadcluster)
    if sum(size(tadcluster{1,t}.dff0_bystimtype{1,1})) > 0
    if length(tadcluster{1,t}.resp_ROIs) > 0
        figure;
        colormap('hot')
        colormap(flipud(colormap))
        imagesc(tadcluster{1,t}.respROIdff0_maxR_sq_MS)
        colorbar
        title(sprintf('tad %d responding cells xcorr (maxR) using multi df/f0', t))
        fig_filename = sprintf('tad %d Responding ROIs xcorr maxR multi dff0', t);
        saveas(gcf,fig_filename,'png');
        close;
    end
    end
end

%visual
for t = 1:length(tadcluster)
    if sum(size(tadcluster{1,t}.dff0_bystimtype{1,1})) > 0
    if length(tadcluster{1,t}.resp_ROIs) > 0
        figure;
        colormap('hot')
        colormap(flipud(colormap))
        imagesc(tadcluster{1,t}.respROIdff0_maxR_sq_V)
        colorbar
        title(sprintf('tad %d responding cells xcorr (maxR) using vis df/f0', t))
        fig_filename = sprintf('tad %d Responding ROIs xcorr maxR vis dff0', t);
        saveas(gcf,fig_filename,'png');
        close;
    end
    end
end

%mechanosensory
for t = 1:length(tadcluster)
    if sum(size(tadcluster{1,t}.dff0_bystimtype{1,1})) > 0
    if length(tadcluster{1,t}.resp_ROIs) > 0
        figure;
        colormap('hot')
        colormap(flipud(colormap))
        imagesc(tadcluster{1,t}.respROIdff0_maxR_sq_M)
        colorbar
        title(sprintf('tad %d responding cells xcorr (maxR) using mech df/f0', t))
        fig_filename = sprintf('tad %d Responding ROIs xcorr maxR mech dff0', t);
        saveas(gcf,fig_filename,'png');
        close;
    end
    end
end

%no stimulus
for t = 1:length(tadcluster)
    if sum(size(tadcluster{1,t}.dff0_bystimtype{1,1})) > 0
    if length(tadcluster{1,t}.resp_ROIs) > 0
        figure;
        colormap('hot')
        colormap(flipud(colormap))
        imagesc(tadcluster{1,t}.respROIdff0_maxR_sq_N)
        colorbar
        title(sprintf('tad %d responding cells xcorr (maxR) using no stim df/f0', t))
        fig_filename = sprintf('tad %d Responding ROIs xcorr maxR no stim dff0', t);
        saveas(gcf,fig_filename,'png');
        close;
    end
    end
end

% Plot lag of max R as a matrix
% this version plots the absolute value of the lag so I get mirror images
% around the diagonal. 

%multisensory
for t = 1:length(tadcluster)
    if sum(size(tadcluster{1,t}.dff0_bystimtype{1,1})) > 0
    if length(tadcluster{1,t}.resp_ROIs) > 0
        %to_plot = abs
        figure;
        colormap('jet')
        %colormap(flipud(colormap))
        imagesc(abs(tadcluster{1,t}.respROIdff0_maxRlag_sq_MS))
        colorbar
        title(sprintf('tad %d responding cells xcorr lag of MaxR using multi df/f0', t))
        fig_filename = sprintf('tad %d Responding ROIs xcorr maxR lags multi dff0', t);
        saveas(gcf,fig_filename,'png');
        close;
    end
    end
end

%visual
for t = 1:length(tadcluster)
    if sum(size(tadcluster{1,t}.dff0_bystimtype{1,1})) > 0
    if length(tadcluster{1,t}.resp_ROIs) > 0
        %to_plot = abs
        figure;
        colormap('jet')
        %colormap(flipud(colormap))
        imagesc(abs(tadcluster{1,t}.respROIdff0_maxRlag_sq_V))
        colorbar
        title(sprintf('tad %d responding cells xcorr lag of MaxR using vis df/f0', t))
        fig_filename = sprintf('tad %d Responding ROIs xcorr maxR lags vis dff0', t);
        saveas(gcf,fig_filename,'png');
        close;
    end
    end
end

%mechanosensory
for t = 1:length(tadcluster)
    if sum(size(tadcluster{1,t}.dff0_bystimtype{1,1})) > 0
    if length(tadcluster{1,t}.resp_ROIs) > 0
        %to_plot = abs
        figure;
        colormap('jet')
        %colormap(flipud(colormap))
        imagesc(abs(tadcluster{1,t}.respROIdff0_maxRlag_sq_M))
        colorbar
        title(sprintf('tad %d responding cells xcorr lag of MaxR using mech df/f0', t))
        fig_filename = sprintf('tad %d Responding ROIs xcorr maxR lags mech dff0', t);
        saveas(gcf,fig_filename,'png');
        close;
    end
    end
end

%no stim
for t = 1:length(tadcluster)
    if sum(size(tadcluster{1,t}.dff0_bystimtype{1,1})) > 0
    if length(tadcluster{1,t}.resp_ROIs) > 0
        %to_plot = abs
        figure;
        colormap('jet')
        %colormap(flipud(colormap))
        imagesc(abs(tadcluster{1,t}.respROIdff0_maxRlag_sq_N))
        colorbar
        title(sprintf('tad %d responding cells xcorr lag of MaxR using no stim df/f0', t))
        fig_filename = sprintf('tad %d Responding ROIs xcorr maxR lags no stim dff0', t);
        saveas(gcf,fig_filename,'png');
        close;
    end
    end
end

% distribution of R values for each stim type (1 figure per tadpole)
for t = 1:length(tadcluster)
    if sum(size(tadcluster{1,t}.dff0_bystimtype{1,1})) > 0
        if length(tadcluster{1,t}.resp_ROIs) > 0
        figure;
        subplot(2,2,1)
        hist(tadcluster{1,t}.respROIdff0_maxR_sq_MS,40)
        title('Multi')
        subplot(2,2,2)
        hist(tadcluster{1,t}.respROIdff0_maxR_sq_V,40)
        title('Vis')
        subplot(2,2,3)
        hist(tadcluster{1,t}.respROIdff0_maxR_sq_M,40)
        title('Mech')
        subplot(2,2,4)
        hist(tadcluster{1,t}.respROIdff0_maxR_sq_N,40)
        title('None')
        suptitle(sprintf('tad %d responding cells hist of xcorr R vals by stimtype', t))
        fig_filename = sprintf('tad %d responding cells hist of xcorr Rvals by stimtype', t)
        saveas(gcf,fig_filename,'png');
        close;
        end
    end
end

% distribution of lag of maxR for each stimtype (1 figure per tadpole
for t = 1:length(tadcluster)
    if sum(size(tadcluster{1,t}.dff0_bystimtype{1,1})) > 0
        if length(tadcluster{1,t}.resp_ROIs) > 0
        lim_min = min(tadcluster{1,t}.respROIdff0_lag);
        lim_max = max(tadcluster{1,t}.respROIdff0_lag);
        figure;
        subplot(2,2,1)
        hist(tadcluster{1,t}.respROIdff0_maxRlag_sq_MS,40)
        xlim([lim_min, lim_max])
        title('Multi')
        subplot(2,2,2)
        hist(tadcluster{1,t}.respROIdff0_maxRlag_sq_V,40)
        xlim([lim_min, lim_max])
        title('Vis')
        subplot(2,2,3)
        hist(tadcluster{1,t}.respROIdff0_maxRlag_sq_M,40)
        xlim([lim_min, lim_max])
        title('Mech')
        subplot(2,2,4)
        hist(tadcluster{1,t}.respROIdff0_maxRlag_sq_N,40)
        xlim([lim_min, lim_max])
        title('None')
        suptitle(sprintf('tad %d responding cells hist of xcorr lag of maxR by stimtype', t))
        fig_filename = sprintf('tad %d responding cells hist of xcorr lag of maxR by stimtype', t)
        saveas(gcf,fig_filename,'png');
        close;
        end
    end
end


%% find "highly correlated" ROIs by locating all maxR > 0.5

% based on all df/f0
for t=1:length(tadcluster)
    if sum(size(tadcluster{1,t}.dff0_bystimtype{1,1})) > 0
        if length(tadcluster{1,t}.resp_ROIs) > 0
            for r = 1:size(tadcluster{1,t}.respROIdff0_maxR_sq, 1)
                for c = 1:size(tadcluster{1,t}.respROIdff0_maxR_sq, 2)
                    if tadcluster{1,t}.respROIdff0_maxR_sq(r,c) > 0.5
                        highcorr(r, c) = 1;
                    else
                        highcorr(r, c) = 0;
                    end
                end
            end
            tadcluster{1,t}.respROIdff0_highcorr = logical(highcorr)
            clear('highcorr')
        end
    end
end

% based on multi df/f0
for t=1:length(tadcluster)
    if sum(size(tadcluster{1,t}.dff0_bystimtype{1,1})) > 0
        if length(tadcluster{1,t}.resp_ROIs) > 0
            for r = 1:size(tadcluster{1,t}.respROIdff0_maxR_sq_MS, 1)
                for c = 1:size(tadcluster{1,t}.respROIdff0_maxR_sq_MS, 2)
                    if tadcluster{1,t}.respROIdff0_maxR_sq_MS(r,c) > 0.5
                        highcorr(r, c) = 1;
                    else
                        highcorr(r, c) = 0;
                    end
                end
            end
            tadcluster{1,t}.respROIdff0_highcorr_MS = logical(highcorr)
            clear('highcorr')
        end
    end
end            
            
% based on vis df/f0
for t=1:length(tadcluster)
    if sum(size(tadcluster{1,t}.dff0_bystimtype{1,1})) > 0
        if length(tadcluster{1,t}.resp_ROIs) > 0
            for r = 1:size(tadcluster{1,t}.respROIdff0_maxR_sq_V, 1)
                for c = 1:size(tadcluster{1,t}.respROIdff0_maxR_sq_V, 2)
                    if tadcluster{1,t}.respROIdff0_maxR_sq_V(r,c) > 0.5
                        highcorr(r, c) = 1;
                    else
                        highcorr(r, c) = 0;
                    end
                end
            end
            tadcluster{1,t}.respROIdff0_highcorr_V = logical(highcorr)
            clear('highcorr')
        end
    end
end 
                
% based on mech df/f0
for t=1:length(tadcluster)
    if sum(size(tadcluster{1,t}.dff0_bystimtype{1,1})) > 0
        if length(tadcluster{1,t}.resp_ROIs) > 0
            for r = 1:size(tadcluster{1,t}.respROIdff0_maxR_sq_M, 1)
                for c = 1:size(tadcluster{1,t}.respROIdff0_maxR_sq_M, 2)
                    if tadcluster{1,t}.respROIdff0_maxR_sq_M(r,c) > 0.5
                        highcorr(r, c) = 1;
                    else
                        highcorr(r, c) = 0;
                    end
                end
            end
            tadcluster{1,t}.respROIdff0_highcorr_M = logical(highcorr)
            clear('highcorr')
        end
    end
end 

% based on no stim df/f0
for t=1:length(tadcluster)
    if sum(size(tadcluster{1,t}.dff0_bystimtype{1,1})) > 0
        if length(tadcluster{1,t}.resp_ROIs) > 0
            for r = 1:size(tadcluster{1,t}.respROIdff0_maxR_sq_N, 1)
                for c = 1:size(tadcluster{1,t}.respROIdff0_maxR_sq_N, 2)
                    if tadcluster{1,t}.respROIdff0_maxR_sq_N(r,c) > 0.5
                        highcorr(r, c) = 1;
                    else
                        highcorr(r, c) = 0;
                    end
                end
            end
            tadcluster{1,t}.respROIdff0_highcorr_N = logical(highcorr)
            clear('highcorr')
        end
    end
end 

%% Now assign ROIs to groups based on high correlations

% get indexes of ROIs that are highly correlated, by ROI
for t = 1:length(tadcluster)
   if sum(size(tadcluster{1,t}.dff0_bystimtype{1,1})) > 0
        if length(tadcluster{1,t}.resp_ROIs) > 0
            for r = 1:size(tadcluster{1,t}.respROIdff0_highcorr,1)
                tadcluster{1,t}.correlated_ROIs_alldff0{1,r} = find(tadcluster{1,t}.respROIdff0_highcorr(r,:));
                tadcluster{1,t}.correlated_ROIs_dff0_MS{1,r} = find(tadcluster{1,t}.respROIdff0_highcorr_MS(r,:));
                tadcluster{1,t}.correlated_ROIs_dff0_V{1,r} = find(tadcluster{1,t}.respROIdff0_highcorr_V(r,:));
                tadcluster{1,t}.correlated_ROIs_dff0_M{1,r} = find(tadcluster{1,t}.respROIdff0_highcorr_M(r,:));
                tadcluster{1,t}.correlated_ROIs_dff0_N{1,r} = find(tadcluster{1,t}.respROIdff0_highcorr_N(r,:));
            end
        end
   end
end

% Determine overlap in which other ROIs are correlated with that ROI over all other ROIs
for t = 1:length(tadcluster)
   if sum(size(tadcluster{1,t}.dff0_bystimtype{1,1})) > 0
        if length(tadcluster{1,t}.resp_ROIs) > 0
            for r = 1:size(tadcluster{1,t}.correlated_ROIs_alldff0,2)
                for c = 1:size(tadcluster{1,t}.correlated_ROIs_alldff0,2)
                    tadcluster{1,t}.correlated_ROIs_alldff0_int{r,c} = intersect(tadcluster{1,t}.correlated_ROIs_alldff0{1,r}, tadcluster{1,t}.correlated_ROIs_alldff0{1,c});
                    tadcluster{1,t}.correlated_ROIs_dff0_MS_int{r,c} = intersect(tadcluster{1,t}.correlated_ROIs_dff0_MS{1,r}, tadcluster{1,t}.correlated_ROIs_dff0_MS{1,c});
                    tadcluster{1,t}.correlated_ROIs_dff0_V_int{r,c} = intersect(tadcluster{1,t}.correlated_ROIs_dff0_V{1,r}, tadcluster{1,t}.correlated_ROIs_dff0_V{1,c});
                    tadcluster{1,t}.correlated_ROIs_dff0_M_int{r,c} = intersect(tadcluster{1,t}.correlated_ROIs_dff0_M{1,r}, tadcluster{1,t}.correlated_ROIs_dff0_M{1,c});
                    tadcluster{1,t}.correlated_ROIs_dff0_MN_int{r,c} = intersect(tadcluster{1,t}.correlated_ROIs_dff0_N{1,r}, tadcluster{1,t}.correlated_ROIs_dff0_N{1,c});
                end
            end
        end
   end
end

% Determine the ROIs that are the same across all cells with significant
% overlap with a given ROI
% significan overlap = at least 1/6*total num ROIS in correlated_ROIs_alldff0_int
clear('lens', 'int', 'roi_count')
for t = 1:length(tadcluster)
    if sum(size(tadcluster{1,t}.dff0_bystimtype{1,1})) > 0
        if length(tadcluster{1,t}.resp_ROIs) > 0
            roi_count = (1/6)*length(tadcluster{1,t}.resp_ROIs);
            for row = 1:size(tadcluster{1,t}.correlated_ROIs_alldff0_int,1)
                for ct = 1:size(tadcluster{1,t}.correlated_ROIs_alldff0_int,2)
                    lens(ct) = length(tadcluster{1,t}.correlated_ROIs_alldff0_int{ct,row});
                end
                first_roi = find((lens > roi_count), 1)
                if isempty(first_roi) 
                    continue
                else
                    int = tadcluster{1,t}.correlated_ROIs_alldff0_int{row, first_roi}

                    for col = first_roi:size(tadcluster{1,t}.correlated_ROIs_alldff0_int,2)
                        if lens(col) > roi_count
                            int = intersect(int, tadcluster{1,t}.correlated_ROIs_alldff0_int{row, col});
                        else
                            continue
                        end
                    end
                tadcluster{1,t}.correlated_ROIs_alldff0_common{row} = int;
                end
                clear('lens', 'int')
            end
            clear('roi_count')
        end
    end
end

%% plot tectum-shaped scatter plot of the ROIs remaining in an intersection

% first take the ROI numbers and index them to the actual ROIs 
for t = 1:length(tadcluster)
    if sum(size(tadcluster{1,t}.dff0_bystimtype{1,1})) > 0
        if isfield(tadcluster{1,t}, 'correlated_ROIs_alldff0_common')
            roi_list = tadcluster{1,t}.resp_ROIs;
            for i = 1:length(tadcluster{1,t}.correlated_ROIs_alldff0_common)
                if ~isempty(tadcluster{1,t}.correlated_ROIs_alldff0_common(i))
                    tadcluster{1,t}.correlated_ROIs_alldff0_common_AROI{i} = roi_list(tadcluster{1,t}.correlated_ROIs_alldff0_common{i});
                end
            end
        end
        clear('roi_list')
    end
end


% now generate input data for scatter(x,y) and plot

for t = 1:length(tadcluster)
    if sum(size(tadcluster{1,t}.dff0_bystimtype{1,1})) > 0
        if isfield(tadcluster{1,t}, 'correlated_ROIs_alldff0_common_AROI')
            rois = [];
            list = [];
            for i = 1:length(tadcluster{1,t}.correlated_ROIs_alldff0_common_AROI)
                rois = [rois tadcluster{1,t}.correlated_ROIs_alldff0_common_AROI{i}'] %all ROIs in all groups
                if ~isempty(tadcluster{1,t}.correlated_ROIs_alldff0_common_AROI{i})
                    list = [list tadcluster{1,t}.resp_ROIs(i)]; %the ROIs we're using for the groups
                else
                    continue
                end
            end
            xy_data = tadcluster{1,t}.ROIcenters(rois,:); %get ROI centers for all ROIs being included
            % change size for each ROI-based cluster
            sizes = [];
            N = 1;
            for j = 1:length(tadcluster{1,t}.correlated_ROIs_alldff0_common_AROI)
                if ~isempty(tadcluster{1,t}.correlated_ROIs_alldff0_common{j})
                    len=length(tadcluster{1,t}.correlated_ROIs_alldff0_common{j});
                    sizes = [sizes; N*ones(len,1)] ;
                    N = N+1;
                end
            end
            plot_sizes = (1 ./ sizes) * 1000;
            %assign colors to each cluster
            cmap = colormap(jet(N-1));
            cluster_colors = cmap(sizes, :);
            figure;
            hold on
            scatter(tadcluster{1,t}.ROIcenters(:,1), tadcluster{1,t}.ROIcenters(:,2), '+', 'MarkerEdgeColor', 'k')
            scatter(tadcluster{1,t}.ROIcenters(tadcluster{1,t}.resp_ROIs,1), tadcluster{1,t}.ROIcenters(tadcluster{1,t}.resp_ROIs,2), 'o', 'MarkerEdgeColor', 'b')
            scatter(xy_data(:,1), xy_data(:,2), plot_sizes, cluster_colors, 'filled')
            hold off
            % add text box with list of ROIs
            dim = [0.15 0.1 .2 .2];
            L = length(list);
            str = sprintf('%s\n%s\n%s\n%s', num2str(list(1:floor(L/4))), num2str(list((ceil(L/4)+1):floor(2*L/4))), num2str(list((ceil((2*L/4)+1):floor(3*L/4)))), num2str(list((ceil((3*L/4)+1):end)))) 
            annotation('textbox',dim,'String',str,'LineStyle', 'none', 'FitBoxToText', 'on');
            title(sprintf('tad %d ROIs in a cluster', t))
            fig_filename = sprintf('tad %d ROIs in a cluster', t)
            saveas(gcf,fig_filename,'png');
            close;
            clear('cmap', 'cluster_colors', 'dim', 'str', 'plot_sizes', 'xy_data', 'len', 'L');
        end
    end
end


%% What features are similar among groups of cells that cluster?




%% Are there more correlated cells with multi vs uni?

% Get the number of cells correlated with a given cell using each stim
% modality

for t = 1:length(tadcluster)
   if sum(size(tadcluster{1,t}.dff0_bystimtype{1,1})) > 0
        if length(tadcluster{1,t}.resp_ROIs) > 0
            tadcluster{1,t}.highcorr_numROIs_all = sum(tadcluster{1,t}.respROIdff0_highcorr);
            tadcluster{1,t}.highcorr_numROIs_MS = sum(tadcluster{1,t}.respROIdff0_highcorr_MS);
            tadcluster{1,t}.highcorr_numROIs_V = sum(tadcluster{1,t}.respROIdff0_highcorr_V);
            tadcluster{1,t}.highcorr_numROIs_M = sum(tadcluster{1,t}.respROIdff0_highcorr_M);
            tadcluster{1,t}.highcorr_numROIs_N = sum(tadcluster{1,t}.respROIdff0_highcorr_N);
        end
   end
end

% plot multi vs others with trendline
for t = 1:length(tadcluster)
   if sum(size(tadcluster{1,t}.dff0_bystimtype{1,1})) > 0
        if length(tadcluster{1,t}.resp_ROIs) > 0
            myfit1 = polyfit(tadcluster{1,t}.highcorr_numROIs_MS, tadcluster{1,t}.highcorr_numROIs_all, 1)
            myfit2 = polyfit(tadcluster{1,t}.highcorr_numROIs_MS, tadcluster{1,t}.highcorr_numROIs_V, 1)
            myfit3 = polyfit(tadcluster{1,t}.highcorr_numROIs_MS, tadcluster{1,t}.highcorr_numROIs_M, 1)
            myfit4 = polyfit(tadcluster{1,t}.highcorr_numROIs_MS, tadcluster{1,t}.highcorr_numROIs_N, 1)
            linspace1 = linspace(0, max(tadcluster{1,t}.highcorr_numROIs_MS))
            figure;
            subplot(2,2,1)
            scatter(tadcluster{1,t}.highcorr_numROIs_MS, tadcluster{1,t}.highcorr_numROIs_all)
            hold on
            plot(linspace1, polyval(myfit1, linspace1))
            hold off
            title('multi vs all')
            %xlabel('multi counts')
            ylabel('all counts')

            subplot(2,2,2)
            scatter(tadcluster{1,t}.highcorr_numROIs_MS, tadcluster{1,t}.highcorr_numROIs_V)
            hold on
            plot(linspace1, polyval(myfit2, linspace1))
            hold off
            title('multi vs vis')
            %xlabel('multi counts')
            ylabel('vis counts')

            subplot(2,2,3)
            scatter(tadcluster{1,t}.highcorr_numROIs_MS, tadcluster{1,t}.highcorr_numROIs_M)
            hold on
            plot(linspace1, polyval(myfit3, linspace1))
            hold off
            title('multi vs mech')
            xlabel('multi counts')
            ylabel('mech counts')

            subplot(2,2,4)
            scatter(tadcluster{1,t}.highcorr_numROIs_MS, tadcluster{1,t}.highcorr_numROIs_N)
            hold on
            plot(linspace1, polyval(myfit4, linspace1))
            hold off
            title('multi vs no stim')
            xlabel('multi counts')
            ylabel('no stim counts')

            suptitle(sprintf('tad %d responding cells numROIs high xcorr by cell', t))
            fig_filename = sprintf('tad %d responding cells numROIs high xcorr by cell', t)
            saveas(gcf,fig_filename,'png');
            close;
        end
   end
end

