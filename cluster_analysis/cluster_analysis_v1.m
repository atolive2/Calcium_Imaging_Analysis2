%% Cluster Analysis
% last updated 20170710
% This code takes data by tadpole and hierarchacally clusters it.

%% First, make a new variable with just the useful stuff from tadpole in it. 
% this only uses the smoothed df/f0 data (moving avg, 8).
for t = 1:length(tadpole)
    %basics
    tadcluster{1,t}.stimorder = tadpole{1,t}.stimorder;
    tadcluster{1,t}.expnum = tadpole{1,t}.expnum;
    tadcluster{1,t}.ROIcenters = tadpole{1,t}.ROIlocs_dbl;
    tadcluster{1,t}.df_f0 = tadpole{1,t}.smoothed;
    % by cell, by trial data
    tadcluster{1,t}.area_bytrial = tadpole{1,t}.area_bytrial_sm;
    tadcluster{1,t}.peak_bytrial = tadpole{1,t}.peak_bytrial_sm;
    tadcluster{1,t}.peakloc_bytrial = tadpole{1,t}.peakloc_bytrial_sm;
    % avg trials by stimtype each cell
    tadcluster{1,t}.area_avg = tadpole{1,t}.area_avg_sm;
    tadcluster{1,t}.peak_avg = tadpole{1,t}.peak_avg_sm;
    tadcluster{1,t}.peakloc_avg = tadpole{1,t}.peakloc_avg_sm;
    % calculated info using avg data
    tadcluster{1,t}.MSenh_area = tadpole{1,t}.MSenh_area_sm;
    tadcluster{1,t}.MSenh_peak = tadpole{1,t}.MSenh_peak_sm;
    tadcluster{1,t}.MSenh_peakloc = tadpole{1,t}.MSenh_peakloc_sm;
    tadcluster{1,t}.unimax_peakavg = tadpole{1,t}.unimax_peakavg_sm;
    tadcluster{1,t}.unimax_stimtype = tadpole{1,t}.unimax_stimtype_sm;
    tadcluster{1,t}.multimax_peakavg = tadpole{1,t}.multimax_peakavg_sm;
end

% now save tadcluster as a new variable and open it in a clear workspace

%% Normalize the data by range
% This is a necessary step since Euclidean distance can be dominated by the
% variable with the largest range, skewing the clustering results. 

% I have used range normalization because I cannot assume my data is
% normally distributed (Z scores would work well on normally distributed
% data)
    % RN = (x - xmin)/(xmax - xmin)
    
% area_bytrial was a cell datatype, but everything else is type double so
% we convert area_bytrial to double.
for t = 1:length(tadcluster)
    tadcluster{1,t}.area_bytrial = cell2mat(tadcluster{1,t}.area_bytrial)
end

% convert area to range normalized area
for t = 1:length(tadcluster)
    xmax = max(max(tadcluster{1,t}.area_bytrial))
    xmin = 0 % min(min(tadcluster{1,t}.area_bytrial))
    tadRN{1,t}.area_bytrial = (tadcluster{1,t}.area_bytrial - xmin) / (xmax - xmin);
end

% convert peak to range norm peak
for t = 1:length(tadcluster)
    xmax = max(max(tadcluster{1,t}.peak_bytrial))
    xmin = min(min(tadcluster{1,t}.peak_bytrial))
    tadRN{1,t}.peak_bytrial = (tadcluster{1,t}.peak_bytrial - xmin) / (xmax - xmin);
end

plot(tadRN{1,1}.peak_bytrial)
figure;
plot(tadcluster{1,1}.peak_bytrial)
figure;
plot(tadcluster{1,1}.df_f0{50,:}')
% this is giving too much pull to outliers. 

%% Normalize the data by Z-score
% outliers can be more easily addressed with this method over range
% analysis

% this way normalizes over each trial (each column has mean=0 and SD=1)
for t = 1:length(tadcluster)
    tadRN{1,t}.peak_bytrial_Z = zscore(tadcluster{1,t}.peak_bytrial)
end

figure;
plot(tadRN{1,1}.peak_bytrial_Z')
% Looks like a few trials drive outlier behavior - so maybe we eliminate
% trials with too many large Z-scores from clustering and future analysis?

for t = 1:length(tadcluster)
    tadRN{1,t}.area_bytrial_Z = zscore(tadcluster{1,t}.area_bytrial)
end

for t = 1:length(tadcluster)
    tadRN{1,t}.peakloc_bytrial_Z = zscore(tadcluster{1,t}.peakloc_bytrial)
end


%% Plot Z-scored data in 3D to see what's there by eye

% reshape to get a vector instead of a matrix for Z-scored data
% "The elements in B preserve their columnwise ordering from A." 
    % E.g. 1-124 = first column, 125-248 = second column.
for t = 1:length(tadRN)
    tadRN{1,t}.peak_bytrial_ZV = reshape(tadRN{1,t}.peak_bytrial_Z, 1, []);
end

for t = 1:length(tadRN)
    tadRN{1,t}.peakloc_bytrial_ZV = reshape(tadRN{1,t}.peakloc_bytrial_Z, 1, []);
end

for t = 1:length(tadRN)
    tadRN{1,t}.area_bytrial_ZV = reshape(tadRN{1,t}.area_bytrial_Z, 1, []);
end

% Plot in 3D
for t = 1:length(tadRN)
    figure;
    scatter3(tadRN{1, t}.area_bytrial_ZV,tadRN{1, t}.peakloc_bytrial_ZV,tadRN{1, t}.peak_bytrial_ZV)
    xlabel('area')
    ylabel('time to peak')
    zlabel('peak df/f0')
    title(sprintf('Tadpole %d Z-scored parameters all trials all cells', t))
    fig_filename = sprintf('tadpole%d_Z-score_alltrialsallcells', t)
    savefig(fig_filename)
    close;
end

%% Build the proximity matrix (calculate distance between points) for each tadpole

% pdist function takes an mxn matrix. rows=observations, columns =
% variables. Therefore, collect the Z scored and vectorized peak, peak loc
% and area measurements into 1 matrix
for t = 1:length(tadRN)
    tadRN{1,t}.input_ZV_PPLA = [tadRN{1, t}.peak_bytrial_ZV; tadRN{1, t}.peakloc_bytrial_ZV; tadRN{1,t}.area_bytrial_ZV];
end

% apply pdist to get D, the Euclidean distance between all pairs in the data matrix,  tadRN{1,t}.input_ZV_PPLA.
for t = 1:length(tadRN)
    tadRN{1,t}.Eucliddist_ZV_PPLA = squareform(pdist(tadRN{1,t}.input_ZV_PPLA'));
end

%% Apply heirarchical clustering    
    
% Use linkage function to generate heirarchical tree (defaults =
% single(closest) and Euclidean distance)
for t = 1:length(tadRN)
    tadRN{1,t}.Tree_ZV_PPLA = linkage(tadRN{1,t}.input_ZV_PPLA')
end

% Use dendrogram function to turn Tree_ZV_PPLA into a graph
for t = 1:length(tadRN)
    [tadRN{1,t}.H_ZV_PPLA tadRN{1,t}.T_ZV_PPLA tadRN{1,t}.outperm_ZV_PPLA] = dendrogram(tadRN{1,t}.Tree_ZV_PPLA)
end

%%%%% Whoops, this doesn't work because the recursion limit on dendrogram
%%%%% is 500 and size(tadRN{1,t}.Tree_ZV_PPLA) = [4463 3]

%% Part 2: Repeat analysis with a different set of data that has smaller dimensions. 

% Reduce dimensionality by using average data for each cell (e.g. do not
% include every trial separately)

    
    
%% Idea #2 - calculate cross correlations across df/f0. 

% recombine trial by trial df/f0 into 1 long vector for each ROI
for t = 1:length(tadcluster)
    for i = 1:size(tadcluster{1,t}.df_f0, 1)
        tmp = [];
        for j = 1:size(tadcluster{1,t}.df_f0, 2)
            tmp = [tmp; tadcluster{1,t}.df_f0{i,j}];
        end
        tadcluster{1,t}.alldff0(i,:) = tmp';
    end
end

% % for deleting var when I screwed up - don't run
% for t = 1:length(tadcluster)
% tadcluster{1,t}.alldff0 = [];    
% end

% calculate correlation coefficient for all ROIs using all trials
for t = 1:length(tadcluster)
    [tadcluster{1,t}.alldff0_R, tadcluster{1,t}.alldff0_P] = corrcoef(tadcluster{1,t}.alldff0');
end

% any signficant?
% Is p < 0.05 for any ROI combinations?
for t = 1:length(tadcluster)
    tadcluster{1,t}.sigP = find(tadcluster{1,t}.alldff0_P < 0.05);
    tadcluster{1,t}.sigP_total = length(tadcluster{1,t}.sigP)
end

for t = 1:length(tadcluster)
    sigPcount(t) = tadcluster{1,t}.sigP_total
end

% Plot correlation coefficients as a matrix
for t = 1:length(tadcluster)
    figure;
    colormap('hot')
    colormap(flipud(colormap))
    imagesc(tadcluster{1,t}.alldff0_P)
    colorbar
    title(sprintf('tad %d all cells correlations using all df/f0', t))
    fig_filename = sprintf('tad %d all cells correlation', t)
    saveas(gcf,fig_filename,'png');
    close;
end

%% Make sense of correlation coefficients - do they map onto any obvious
% properties?

% Compare correlation coefficients as a function of proximity

% calculate distance between ROIs
for t = 1:length(tadcluster)
    for i = 1:size(tadcluster{1,t}.ROIcenters,1)
        for j = 1:size(tadcluster{1,t}.ROIcenters)
            tadcluster{1,t}.ROIdist(i,j) = sqrt( (tadcluster{1,t}.ROIcenters(i,1) + tadcluster{1,t}.ROIcenters(j,1))^2 + (tadcluster{1,t}.ROIcenters(i,2) + tadcluster{1,t}.ROIcenters(j,2))^2);
        end
    end
end

% plot distance vs P val for each ROI pair
for t = 1:length(tadcluster)
    figure;
for i = 1:(size(tadcluster{1,t}.alldff0_P,1)-1)
    hold on
    plot(tadcluster{1,t}.ROIdist(i,:), tadcluster{1,t}.alldff0_P(i+1,:), 'o')
end
    hold off
    title(sprintf('tad %d all cells ROI distance vs P value df/f0', t));
    fig_filename = sprintf('tad %d all cells ROI distance vs P value', t);
    xlabel('ROI distance (pixels)')
    ylabel('P value')
    saveas(gcf,fig_filename,'png');
    close;
end

%% Is there more correlation during a specific stimulus type?

% Make df_f0 separated by trial type (only use 1-4 here to eliminate
% variation from stim strength)
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

% % for deleting when I screw up
for t = 1:length(tadcluster)
tadcluster{1,t}.dff0_bystimtype_MS_R = [];    
tadcluster{1,t}.dff0_bystimtype_MS_P = []; 
end

% calculate correlation coefficient for all ROIs, each stimtype separately
% first get a matrix ROI x dff0 for all trials of a given stimtype
for t = 1:length(tadcluster)
    for i = 1:size(tadcluster{1,t}.dff0_bystimtype,1)
        if sum(size(tadcluster{1,t}.dff0_bystimtype{i,1})) > 0 %elimnate tads with no high/high multisensory trials
        tadcluster{1,t}.dff0_multi(i,:) = tadcluster{1,t}.dff0_bystimtype{i,1}';
        tadcluster{1,t}.dff0_vis(i,:) = tadcluster{1,t}.dff0_bystimtype{i,2}';
        tadcluster{1,t}.dff0_mech(i,:) = tadcluster{1,t}.dff0_bystimtype{i,3}';
        tadcluster{1,t}.dff0_none(i,:) = tadcluster{1,t}.dff0_bystimtype{i,4}';
        end
    end
end
    
for t = 1:length(tadcluster)
    if sum(size(tadcluster{1,t}.dff0_bystimtype{1,1})) > 0
    [tadcluster{1,t}.dff0_bystimtype_MS_R, tadcluster{1,t}.dff0_bystimtype_MS_P] = corrcoef(tadcluster{1,t}.dff0_multi');
    [tadcluster{1,t}.dff0_bystimtype_V_R, tadcluster{1,t}.dff0_bystimtype_V_P] = corrcoef(tadcluster{1,t}.dff0_vis');
    [tadcluster{1,t}.dff0_bystimtype_M_R, tadcluster{1,t}.dff0_bystimtype_M_P] = corrcoef(tadcluster{1,t}.dff0_mech');
    [tadcluster{1,t}.dff0_bystimtype_N_R, tadcluster{1,t}.dff0_bystimtype_N_P] = corrcoef(tadcluster{1,t}.dff0_none');
    end
end

% any signficant?
% Is p < 0.05 for any ROI combinations?
for t = 1:length(tadcluster)
    tadcluster{1,t}.sigP = find(tadcluster{1,t}.alldff0_P < 0.05);
    tadcluster{1,t}.sigP_total = length(tadcluster{1,t}.sigP)
end

for t = 1:length(tadcluster)
    sigPcount(t) = tadcluster{1,t}.sigP_total
end

% Plot correlation coefficients as a matrix
% multisensory
for t = 1:length(tadcluster)
    if sum(size(tadcluster{1,t}.dff0_bystimtype{1,1})) > 0
        figure;
        colormap('hot')
        colormap(flipud(colormap))
        imagesc(tadcluster{1,t}.dff0_bystimtype_MS_P)
        colorbar
        title(sprintf('tad %d all cells correlations using Multi df/f0', t))
        fig_filename = sprintf('tad %d all cells Multi correlation', t)
        saveas(gcf,fig_filename,'png');
        close;
    end
end

%visual
for t = 1:length(tadcluster)
    if sum(size(tadcluster{1,t}.dff0_bystimtype{1,1})) > 0
        figure;
        colormap('hot')
        colormap(flipud(colormap))
        imagesc(tadcluster{1,t}.dff0_bystimtype_V_P)
        colorbar
        title(sprintf('tad %d all cells correlations using Vis df/f0', t))
        fig_filename = sprintf('tad %d all cells Vis correlation', t)
        saveas(gcf,fig_filename,'png');
        close;
    end
end

% mech
for t = 1:length(tadcluster)
    if sum(size(tadcluster{1,t}.dff0_bystimtype{1,1})) > 0
        figure;
        colormap('hot')
        colormap(flipud(colormap))
        imagesc(tadcluster{1,t}.dff0_bystimtype_M_P)
        colorbar
        title(sprintf('tad %d all cells correlations using Mech df/f0', t))
        fig_filename = sprintf('tad %d all cells Mech correlation', t)
        saveas(gcf,fig_filename,'png');
        close;
    end
end

% none
for t = 1:length(tadcluster)
    if sum(size(tadcluster{1,t}.dff0_bystimtype{1,1})) > 0
        figure;
        colormap('hot')
        colormap(flipud(colormap))
        imagesc(tadcluster{1,t}.dff0_bystimtype_N_P)
        colorbar
        title(sprintf('tad %d all cells correlations using no stim df/f0', t))
        fig_filename = sprintf('tad %d all cells no stim correlation', t)
        saveas(gcf,fig_filename,'png');
        close;
    end
end

% distribution of P values for each stim type (1 figure per tadpole)
for t = 1:length(tadcluster)
    if sum(size(tadcluster{1,t}.dff0_bystimtype{1,1})) > 0
        figure;
        subplot(2,2,1)
        hist(tadcluster{1,t}.dff0_bystimtype_MS_P,40)
        title('Multi')
        subplot(2,2,2)
        hist(tadcluster{1,t}.dff0_bystimtype_V_P,40)
        title('Vis')
        subplot(2,2,3)
        hist(tadcluster{1,t}.dff0_bystimtype_M_P,40)
        title('Mech')
        subplot(2,2,4)
        hist(tadcluster{1,t}.dff0_bystimtype_N_P,40)
        title('None')
        suptitle(sprintf('tad %d hist of P vals by stimtype', t))
        fig_filename = sprintf('tad %d all cells Pvals of correlation by stimtype', t)
        saveas(gcf,fig_filename,'png');
        close;
    end
end

%% None of the above really gave info that made sense or structured te data meaningfully. 
% Next step: try to "clean up" the data first. 

%% Remove trials with a ton of noise




%% Remove cells that don't respond ever














