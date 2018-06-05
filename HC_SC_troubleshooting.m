%% dissertation High Corr vs Sparse Corr troubleshooting
% histograms of vars from allRespROIs look weird. 
    % want to plot onset time MS, onset time SD MS, MSpeak and MSIndex_peak
    % for figure. 
    % where are my vars in allRespROIs?
        %onsettimeMS = 9
        %onsettimeMSSD = 10
        %peakof MS = 8
        % MSInd by peak = 12

%% Start by re-making allRespROIs with just the variables I need for this

% add a field to allData{1,t} that is all unique highcorr ROIs
for t = 1:length(allData)
    tmp = [];
    if isfield(allData{1,t}, 'uniqueHighCorrROI_MS') & isfield(allData{1,t}, 'uniqueHighCorrROI_V') & isfield(allData{1,t}, 'uniqueHighCorrROI_M')
        allData{1,t}.allhighcorrROI = unique([allData{1,t}.uniqueHighCorrROI_MS; allData{1,t}.uniqueHighCorrROI_V; allData{1,t}.uniqueHighCorrROI_M]);
    else 
        if isfield(allData{1,t}, 'uniqueHighCorrROI_MS')
            tmp = [tmp; allData{1,t}.uniqueHighCorrROI_MS];
        end
        if isfield(allData{1,t}, 'uniqueHighCorrROI_V')
            tmp = [tmp; allData{1,t}.uniqueHighCorrROI_V];
        end
        if isfield(allData{1,t}, 'uniqueHighCorrROI_M')
            tmp = [tmp; allData{1,t}.uniqueHighCorrROI_M];
        end
        tmp1 = unique(tmp)
        allData{1,t}.allhighcorrROI = tmp1;
    end
end
          
idx = 1
for t = 1:length(allData)
    for roi = 1:length(allData{1,t}.resp_ROIs)
        r = allData{1,t}.resp_ROIs(roi);
        allRespROIs(idx, 1) = t;
        allRespROIs(idx, 2) = r;
        allRespROIs(idx, 3) = allData{1,t}.stage;
        allRespROIs(idx, 4) = allData{1,t}.unimean_peak(r);
        allRespROIs(idx, 5) = allData{1,t}.unimean_onsettime(r);
        allRespROIs(idx, 6) = allData{1,t}.unistd_onsettime(r); 
        allRespROIs(idx, 7) = allData{1,t}.resp_reliability(r); 
        allRespROIs(idx, 8) = allData{1,t}.avg_peak2(1,r); %multi peak
        allRespROIs(idx, 9) = allData{1,t}.avg_onsettime2(1,r); %multi avg onsettime
        allRespROIs(idx, 10) = allData{1,t}.std_onsettime2(1,r); %multi std onset time
        allRespROIs(idx, 11) = allData{1,t}.MSInd_peak2(r);
        allRespROIs(idx, 12) = allData{1,t}.MSInd_onsettime2(r);
        allRespROIs(idx, 13) = ismember(r,  allData{1,t}.allhighcorrROIs);
        idx = idx + 1
    end
end

%% index the high corr and sparse corr cells

HC_cells = find(allRespROIs(:,13))
SC_cells = find(allRespROIs(:,13) == 0)
% check
length(HC_cells) % ans = 878
length(SC_cells) % ans = 186

%% Check the whole list of data for each var I am interested in

vars = [8 9 10 11];

for v = 1:length(vars)
    figure;
    histogram(allRespROIs(:, vars(v)),60)
    title(sprintf('var %d', vars(v)))
end
% still has very high proportion of 0s

%% How much NaN, inf?

for v = 1:size(allRespROIs,2)
    %tmp = allRespROIs(:,v)
    %tmp1 = tmp(~isnan(tmp))
    
    diff_all(v,1) = sum(isnan(allRespROIs(:,v)));
    diff_all(v,2) = sum(isinf(allRespROIs(:,v)));
    diff_all(v,3) = sum(isfinite(allRespROIs(:,v)));
end
% 11 and 12 have 50 NaN and 26 inf, everything else is all finite. 
% 20 NaN and 7 inf in the SC cells (159 finite)
% 30 NaN and 19 inf in the HC cells (829 finite)
% seems ok and not source of problem


%% What's the minimum and how many cells have that value?
for v = 1:size(allRespROIs, 2)
    d(v,1) = min(allRespROIs(:,v));
    d(v,2) = length(find(allRespROIs(SC_cells,v) == d(v,1)));
end
% 10-20% of ROIs have the minimum value. 
% 8: 0, 142; 9: 0, 294; 10: 0, 294; 11: -1, 92
% 1 ROI has max val for 8, 9, 10. 11 has 26 inf as the max val. 
% so the problem appears to be with the minimum values. Go back to the
% original data that allRespROIs came from. 

%% Use 1 tadpole to investigate potential issue with the input data (in allData). 
t = 20 %tad 34, all ROIs respond
histogram(allData{1,t}.MSInd_peak2)
figure;
histogram(allData{1,t}.MSInd_peak)
figure;
histogram(MSIn) % after changing calc_MSIndex to only use finite values for max(V,M)

%% Make overlaid histograms as lines to determine if the problem is coming from 1 or 2 tadpoles

edges = 0:0.1:7
figure;
hold on
for t = 1:length(allData)
    [h e] = histcounts(allData{1,t}.avg_onsettime2(1,:))
    plot(e(1:length(h)), h)
end

figure;
hold on
for t = 1:length(allData)
    [h e] = histcounts(allData{1,t}.avg_onsettime(1,:))
    plot(e(1:(end-1)), h)
end
% it may be. 
% went back to avg_data_byRespTrial to investigate calculation. Found that
% if there were no trials with responses, 0 was being inserted. Changed it
% to be NaN if number of trials with response = 0. 

histogram(allData{1,t}.avg_onsettime2(1,:))

% If the MSInd == 0, then the ROI did not respond in a stimulus trial (the only way to get 0 is to put in 0 for MS, M, and V)
