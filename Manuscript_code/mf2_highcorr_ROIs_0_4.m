%% Manuscript figures version 2 -- change threshold for high corr ROI

% Lines 445-568 of manuscript_figs_v2
% changed the code so that the fields get saved to a separate struct, HCData

% for threshold = 0.4, the field names are the same as in allData for
% threshold = 0.5

stims = [ 1 2 3 4];
threshold = 0.4;
%% Get highcorr ROIs 

for k = 1:length(allData)
    if length(allData(k).respROIs_ST) > 1
        for s = 1:length(stims)
            for r = 1:size(allData(k).respROIdff0_maxR_sq, 2)
                for c = 1:size(allData(k).respROIdff0_maxR_sq, 3)
                    if allData(k).respROIdff0_maxR_sq(s, r, c) > threshold
                        highcorr(r, c) = 1;
                    else
                        highcorr(r, c) = 0; 
                    end
                end
            end
            HCData(k).respROIdff0_HC(s, :, :) = logical(highcorr);
            clear('highcorr')
        end
    end
end

%% get indexes of ROIs that are highly correlated, by ROI

for k = 1:length(HCData)
    if length(allData(k).respROIs_ST) > 1
        for s = 1:length(stims)
            for r = 1:size(HCData(k).respROIdff0_HC, 2)
                HCData(k).respROIdff0_correlROIs{s, r} = find(HCData(k).respROIdff0_HC(s, r, :));
            end
        end
    end
end

% determine the overlap in which other ROIs are correlated with that ROI
for k = 1:length(HCData)
    if length(allData(k).respROIs_ST) > 1
        for s = 1:length(stims)
            for r = 1:size(HCData(k).respROIdff0_correlROIs, 2)
                for c = 1:size(HCData(k).respROIdff0_correlROIs, 2)
                    HCData(k).correlROIs_int{s, r, c} = intersect(HCData(k).respROIdff0_correlROIs{s, r}, HCData(k).respROIdff0_correlROIs{s, c});
                end
            end
        end
    end
end
                
%% Determine the ROIs that are the same across all cells with significant
% overlap with a given ROI
% significan overlap = at least 1/6*total num ROIS in correlated_ROIs_alldff0_int

clear('lens', 'int', 'roi_count')
for k = 1:length(allData)
    if length(allData(k).respROIs_ST) > 10
        roi_count = (1/6) * length(allData(k).respROIs_ST);
        for s = 1:length(stims)
            for row = 1:size(HCData(k).correlROIs_int, 2)
                for ct = 1:size(HCData(k).correlROIs_int, 3)                
                    lens(ct) = length(HCData(k).correlROIs_int{s, row, ct});
                end
                first_roi = find((lens > roi_count), 1);
                if isempty(first_roi) 
                    continue
                else
                    int = HCData(k).correlROIs_int{s, row, first_roi};
                    for col = first_roi:size(HCData(k).correlROIs_int, 3)
                        if lens(col) > roi_count
                            int = intersect(int, HCData(k).correlROIs_int{s, row, col});
                        else
                            continue
                        end
                    end
                    HCData(k).correlROIs_common{s,row} = int;
                end
                clear('lens', 'int')
            end
            
        end
    end
    clear('roi_count')
end

% Index ROI numbers to the actual ROIs 
for k = 1:length(allData)
    if length(allData(k).respROIs_ST) > 10
        roi_list = allData(k).respROIs_ST;
        for s = 1:length(stims)
            if size(HCData(k).correlROIs_common, 1) >= s
                for i = 1:size(HCData(k).correlROIs_common,2)
                    if ~isempty(HCData(k).correlROIs_common{s,i})
                        HCData(k).correlROIs_common_AROI{s,i} = roi_list(HCData(k).correlROIs_common{s,i});
                    end
                end
            else
                 HCData(k).correlROIs_common_AROI{s,i} = [];
            end
        end
        clear('roi_list')
    end
end

%% How many HighCorr ROIs do I have? 

% Get unique HC count for each tad
for k = 1:length(allData)
    if ~isempty(HCData(k).correlROIs_common_AROI)
        for s = 1:length(stims)
            tmp = [];
            for q = 1:size(HCData(k).correlROIs_common_AROI, 2)
                tmp = [tmp; HCData(k).correlROIs_common_AROI{s, q}];
            end
            HCData(k).uniqueHC_AROI{s} = unique(tmp);
        end
    end
end
                
% Get HC and NHC count, and percentage in a separate variable
for k = 1:length(allData)
    if ~isempty(HCData(k).correlROIs_common_AROI)
        for s = 1:length(stims)
            highcorr_counts0_4(k, s, 1) = length(HCData(k).correlROIs_common_AROI{s}); %HC ROIs
            highcorr_counts0_4(k, s, 2) = length(allData(k).respROIs_ST); %all ROIs
            highcorr_counts0_4(k, s, 3) = highcorr_counts0_4(k, s, 1) / highcorr_counts0_4(k, s, 2); % proportion HC
        end
    else
        highcorr_counts0_4(k, s, :) = NaN;    
    end
end

%% Get HC and NHC for each stage separately

% HC threshold = 0.4 // highcorr_counts0_4
for k = 1:length(allData)
    stage(k) = allData(k).stage;
end

propHC46 = [];
propHC49 = [];
for k = 1:length(stage)
    if length(allData(k).respROIs_ST) > 10
        if stage(k) == 46
            propHC46 = [propHC46; highcorr_counts0_4(k, :, 3)];
        elseif stage(k) == 49
            propHC49 = [propHC49; highcorr_counts0_4(k, :, 3)];
        else
            continue
        end
    end
end







