%% Determine the ROIs that are the same across all cells with significant
% overlap with a given ROI
% significant overlap = at least threshold*total num ROIS in correlated_ROIs_alldff0_int
%original threshold = 1/6
threshold = 1/2
%% Multisensory
clear('lens', 'int', 'roi_count')
for t = 1:length(allData)
    if isfield(allData{1,t}, 'respROIdff0_R_MS')
        if length(allData{1,t}.resp_ROIs) > 10
            roi_count = threshold*length(allData{1,t}.resp_ROIs);
            for row = 1:size(allData{1,t}.correlated_ROIs_dff0_MS_int,1)
                for ct = 1:size(allData{1,t}.correlated_ROIs_dff0_MS_int,2)
                    lens(ct) = length(allData{1,t}.correlated_ROIs_dff0_MS_int{ct,row});
                end
                first_roi = find((lens > roi_count), 1)
                if isempty(first_roi) 
                    continue
                else
                    int = allData{1,t}.correlated_ROIs_dff0_MS_int{row, first_roi}

                    for col = first_roi:size(allData{1,t}.correlated_ROIs_dff0_MS_int,2)
                        if lens(col) > roi_count
                            int = intersect(int, allData{1,t}.correlated_ROIs_dff0_MS_int{row, col});
                        else
                            continue
                        end
                    end
                allData{1,t}.correlated_ROIs_dff0_MS_common1_2{row} = int;
                end
                clear('lens', 'int')
            end
            clear('roi_count')
        end
    end
end

% Index ROI numbers to the actual ROIs 
for t = 1:length(allData)
   % if sum(size(allData{1,t}.dff0_bystimtype{1,1})) > 0
        if isfield(allData{1,t}, 'correlated_ROIs_dff0_MS_common')
            roi_list = allData{1,t}.resp_ROIs;
            for i = 1:length(allData{1,t}.correlated_ROIs_dff0_MS_common)
                if ~isempty(allData{1,t}.correlated_ROIs_dff0_MS_common(i))
                    allData{1,t}.correlated_ROIs_dff0_MS_common_AROI1_2{i} = roi_list(allData{1,t}.correlated_ROIs_dff0_MS_common{i});
                end
            end
        end
        clear('roi_list')
    %end
end

%% Mechanosensory

clear('lens', 'int', 'roi_count')
for t = 1:length(allData)
    if isfield(allData{1,t}, 'respROIdff0_R_M')
        if length(allData{1,t}.resp_ROIs) > 10
            roi_count = threshold*length(allData{1,t}.resp_ROIs);
            clear('lens')
            for row = 1:size(allData{1,t}.correlated_ROIs_dff0_M_int,1)
                for ct = 1:size(allData{1,t}.correlated_ROIs_dff0_M_int,2)
                    lens(ct) = length(allData{1,t}.correlated_ROIs_dff0_M_int{ct,row});
                end
                first_roi = find((lens > roi_count), 1)
                if isempty(first_roi) 
                    continue
                    clear('lens')
                else
                    int = allData{1,t}.correlated_ROIs_dff0_M_int{row, first_roi}

                    for col = first_roi:size(allData{1,t}.correlated_ROIs_dff0_M_int,2)
                        if lens(col) > roi_count
                            int = intersect(int, allData{1,t}.correlated_ROIs_dff0_M_int{row, col});
                        else
                            continue
                        end
                    end
                allData{1,t}.correlated_ROIs_dff0_M_common1_2{row} = int;
                clear('lens', 'int')
                end
                clear('lens', 'int')
            end
            clear('roi_count')
        end
    end
end

% Index ROI numbers to the actual ROIs 
for t = 1:length(allData)
   % if sum(size(allData{1,t}.dff0_bystimtype{1,1})) > 0
        if isfield(allData{1,t}, 'correlated_ROIs_dff0_M_common')
            roi_list = allData{1,t}.resp_ROIs;
            for i = 1:length(allData{1,t}.correlated_ROIs_dff0_M_common)
                if ~isempty(allData{1,t}.correlated_ROIs_dff0_M_common(i))
                    allData{1,t}.correlated_ROIs_dff0_M_common_AROI1_2{i} = roi_list(allData{1,t}.correlated_ROIs_dff0_M_common{i});
                end
            end
        end
        clear('roi_list')
    %end
end

%% Visual 
%threshold = 1/5
clear('lens', 'int', 'roi_count')
for t = 1:length(allData)
    if isfield(allData{1,t}, 'respROIdff0_R_V')
        if length(allData{1,t}.resp_ROIs) > 10
            roi_count = threshold*length(allData{1,t}.resp_ROIs);
            clear('lens')
            for row = 1:size(allData{1,t}.correlated_ROIs_dff0_V_int,1)
                for ct = 1:size(allData{1,t}.correlated_ROIs_dff0_V_int,2)
                    lens(ct) = length(allData{1,t}.correlated_ROIs_dff0_V_int{ct,row});
                end
                first_roi = find((lens > roi_count), 1)
                if isempty(first_roi) 
                    continue
                else
                    int = allData{1,t}.correlated_ROIs_dff0_V_int{row, first_roi}

                    for col = first_roi:size(allData{1,t}.correlated_ROIs_dff0_V_int,2)
                        if lens(col) > roi_count
                            int = intersect(int, allData{1,t}.correlated_ROIs_dff0_V_int{row, col});
                        else
                            continue
                        end
                    end
                allData{1,t}.correlated_ROIs_dff0_V_common1_2{row} = int;
                end
                clear('lens', 'int')
            end
            clear('roi_count', 'lens')
        end
    end
end

% Index ROI numbers to the actual ROIs 
for t = 1:length(allData)
   % if sum(size(allData{1,t}.dff0_bystimtype{1,1})) > 0
        if isfield(allData{1,t}, 'correlated_ROIs_dff0_V_common')
            roi_list = allData{1,t}.resp_ROIs;
            for i = 1:length(allData{1,t}.correlated_ROIs_dff0_V_common)
                if ~isempty(allData{1,t}.correlated_ROIs_dff0_V_common(i))
                    allData{1,t}.correlated_ROIs_dff0_V_common_AROI1_2{i} = roi_list(allData{1,t}.correlated_ROIs_dff0_V_common{i});
                end
            end
        end
        clear('roi_list')
    %end
end


%% What ROIs?
%%%%%%% Use highcorr_vs_nothighcorr code to create
%%%%%%% allData{1,t}.uniquehighcorrROI

for t = 1:length(allData)
    if isfield(allData{1,t}, 'correlated_ROIs_dff0_MS_common_AROI1_2')
    for r = 1:length(allData{1,t}.correlated_ROIs_dff0_MS_common_AROI1_2)
        tmp = [tmp; cell2mat(allData{1,t}.correlated_ROIs_dff0_MS_common_AROI1_2(r))];
    end
    allData{1,t}.uniqueHighCorrROI1_2 =  unique(tmp)
    tmp = [];
    end
end
tmp = [];
for t = 1:length(allData)
    if isfield(allData{1,t}, 'correlated_ROIs_dff0_M_common_AROI1_2')
    for r = 1:length(allData{1,t}.correlated_ROIs_dff0_M_common_AROI1_2)
        tmp = [tmp; cell2mat(allData{1,t}.correlated_ROIs_dff0_M_common_AROI1_2(r))];
    end
    allData{1,t}.uniqueHighCorrROI_M1_2 =  unique(tmp)
    tmp = [];
    end
end
tmp = [];
for t = 1:length(allData)
    if isfield(allData{1,t}, 'correlated_ROIs_dff0_V_common_AROI1_2')
    for r = 1:length(allData{1,t}.correlated_ROIs_dff0_V_common_AROI1_2)
        tmp = [tmp; cell2mat(allData{1,t}.correlated_ROIs_dff0_V_common_AROI1_2(r))];
    end
    allData{1,t}.uniqueHighCorrROI_V1_2 =  unique(tmp)
    tmp = [];
    end
end
tmp = [];

%% Compare 1/6 (original) to 1/3 

% table of values, by stage
%how many HC ROIs in 1/6 vs 1/5 vs 1/3?
for t = 1:length(allData)
    HC_table(t, 1) = allData{1,t}.stage;
    HC_table(t, 2) = length(allData{1,t}.resp_ROIs)
 % 1/6 threshold
    if isfield(allData{1,t}, 'correlated_ROIs_dff0_MS_common_AROI')
        HC_table(t, 3) = length(allData{1,t}.uniqueHighCorrROI_MS);
    else
        HC_table(t, 3) = 0
    end
    if isfield(allData{1,t}, 'correlated_ROIs_dff0_V_common_AROI')
        HC_table(t, 4) = length(allData{1,t}.uniqueHighCorrROI_V);
    else
        HC_table(t, 4) = 0
    end
    if isfield(allData{1,t}, 'correlated_ROIs_dff0_M_common_AROI')
        HC_table(t, 5) = length(allData{1,t}.uniqueHighCorrROI_M);
    else
        HC_table(t, 5) = 0
    end  
  % 1/5 threshold
    if isfield(allData{1,t}, 'correlated_ROIs_dff0_MS_common_AROI1_5')
        HC_table(t, 6) = length(allData{1,t}.uniqueHighCorrROI1_5);
    else
        HC_table(t, 6) = 0
    end
    if isfield(allData{1,t}, 'correlated_ROIs_dff0_V_common_AROI1_5')
        HC_table(t, 7) = length(allData{1,t}.uniqueHighCorrROI_V1_5);
    else
        HC_table(t, 7) = 0
    end
    if isfield(allData{1,t}, 'correlated_ROIs_dff0_M_common_AROI1_5')
        HC_table(t, 8) = length(allData{1,t}.uniqueHighCorrROI_M1_5);
    else
        HC_table(t, 8) = 0
    end
  % 1/3 threshold
    if isfield(allData{1,t}, 'correlated_ROIs_dff0_MS_common_AROI1_3')
        HC_table(t, 9) = length(allData{1,t}.uniqueHighCorrROI1_3);
    else
        HC_table(t, 9) = 0
    end
    if isfield(allData{1,t}, 'correlated_ROIs_dff0_V_common_AROI1_3')
        HC_table(t, 10) = length(allData{1,t}.uniqueHighCorrROI_V1_3);
    else
        HC_table(t, 10) = 0
    end
    if isfield(allData{1,t}, 'correlated_ROIs_dff0_M_common_AROI1_3')
        HC_table(t, 11) = length(allData{1,t}.uniqueHighCorrROI_M1_3);
    else
        HC_table(t, 11) = 0
    end
     
end

% Change in ROI number 




%%%%%%% Use sfn2017_figures to generate tectum shaped scatter plots by
%%%%%%% tadpole of highcorr ROIs
%% PCA plot by HC status

% calc allhighcorrrROIs for 1/3 
for t = 1:length(allData)
    if isfield(allData{1,t}, 'uniqueHighCorrROI_MS1_3')
        if isfield(allData{1,t}, 'uniqueHighCorrROI_M1_3')
           if isfield(allData{1,t}, 'uniqueHighCorrROI_V1_3') 
                allData{1,t}.allhighcorrROIs1_3 = unique([allData{1,t}.uniqueHighCorrROI_MS1_3; allData{1,t}.uniqueHighCorrROI_M1_3; allData{1,t}.uniqueHighCorrROI_V1_3]);
           else
               allData{1,t}.allhighcorrROIs1_3 = unique([allData{1,t}.uniqueHighCorrROI_MS1_3; allData{1,t}.uniqueHighCorrROI_M1_3]);
           end
        else 
            if isfield(allData{1,t}, 'uniqueHighCorrROI_V1_3')
                allData{1,t}.allhighcorrROIs1_3 = unique([allData{1,t}.uniqueHighCorrROI_MS1_3; allData{1,t}.uniqueHighCorrROI_V1_3]);
            else
                allData{1,t}.allhighcorrROIs1_3 = allData{1,t}.uniqueHighCorrROI_MS1_3;
            end
        end
    else
        if isfield(allData{1,t}, 'uniqueHighCorrROI_M1_3')
           if isfield(allData{1,t}, 'uniqueHighCorrROI_V1_3') 
                allData{1,t}.allhighcorrROIs1_3 = unique([allData{1,t}.uniqueHighCorrROI_M1_3; allData{1,t}.uniqueHighCorrROI_V1_3]);
           else
               allData{1,t}.allhighcorrROIs1_3 = unique([allData{1,t}.uniqueHighCorrROI_M1_3]);
           end
        else 
            if isfield(allData{1,t}, 'uniqueHighCorrROI_V1_3')
                allData{1,t}.allhighcorrROIs1_3 = unique([ allData{1,t}.uniqueHighCorrROI_V1_3]);
            else
                allData{1,t}.allhighcorrROIs1_3 = [];
            end
        end
    end
end


% calc allhighcorrrROIs for 1/5
for t = 1:length(allData)
    if isfield(allData{1,t}, 'uniqueHighCorrROI_MS1_5')
        if isfield(allData{1,t}, 'uniqueHighCorrROI_M1_5')
           if isfield(allData{1,t}, 'uniqueHighCorrROI_V1_5') 
                allData{1,t}.allhighcorrROIs1_5 = unique([allData{1,t}.uniqueHighCorrROI_MS1_5; allData{1,t}.uniqueHighCorrROI_M1_5; allData{1,t}.uniqueHighCorrROI_V1_5]);
           else
               allData{1,t}.allhighcorrROIs1_5 = unique([allData{1,t}.uniqueHighCorrROI_MS1_5; allData{1,t}.uniqueHighCorrROI_M1_5]);
           end
        else 
            if isfield(allData{1,t}, 'uniqueHighCorrROI_V1_5')
                allData{1,t}.allhighcorrROIs1_5 = unique([allData{1,t}.uniqueHighCorrROI_MS1_5; allData{1,t}.uniqueHighCorrROI_V1_5]);
            else
                allData{1,t}.allhighcorrROIs1_5 = allData{1,t}.uniqueHighCorrROI_MS1_5;
            end
        end
    else
        if isfield(allData{1,t}, 'uniqueHighCorrROI_M1_5')
           if isfield(allData{1,t}, 'uniqueHighCorrROI_V1_5') 
                allData{1,t}.allhighcorrROIs1_5 = unique([allData{1,t}.uniqueHighCorrROI_M1_5; allData{1,t}.uniqueHighCorrROI_V1_5]);
           else
               allData{1,t}.allhighcorrROIs1_5 = unique([allData{1,t}.uniqueHighCorrROI_M1_5]);
           end
        else 
            if isfield(allData{1,t}, 'uniqueHighCorrROI_V1_5')
                allData{1,t}.allhighcorrROIs1_5 = unique([ allData{1,t}.uniqueHighCorrROI_V1_5]);
            else
                allData{1,t}.allhighcorrROIs1_5 = [];
            end
        end
    end
end

% calc allhighcorrrROIs for 1/2
for t = 1:length(allData)
    if isfield(allData{1,t}, 'uniqueHighCorrROI_MS1_2')
        if isfield(allData{1,t}, 'uniqueHighCorrROI_M1_2')
           if isfield(allData{1,t}, 'uniqueHighCorrROI_V1_2') 
                allData{1,t}.allhighcorrROIs1_2 = unique([allData{1,t}.uniqueHighCorrROI_MS1_2; allData{1,t}.uniqueHighCorrROI_M1_2; allData{1,t}.uniqueHighCorrROI_V1_2]);
           else
               allData{1,t}.allhighcorrROIs1_2 = unique([allData{1,t}.uniqueHighCorrROI_MS1_2; allData{1,t}.uniqueHighCorrROI_M1_2]);
           end
        else 
            if isfield(allData{1,t}, 'uniqueHighCorrROI_V1_2')
                allData{1,t}.allhighcorrROIs1_2 = unique([allData{1,t}.uniqueHighCorrROI_MS1_2; allData{1,t}.uniqueHighCorrROI_V1_2]);
            else
                allData{1,t}.allhighcorrROIs1_2 = allData{1,t}.uniqueHighCorrROI_MS1_2;
            end
        end
    else
        if isfield(allData{1,t}, 'uniqueHighCorrROI_M1_2')
           if isfield(allData{1,t}, 'uniqueHighCorrROI_V1_2') 
                allData{1,t}.allhighcorrROIs1_2 = unique([allData{1,t}.uniqueHighCorrROI_M1_5; allData{1,t}.uniqueHighCorrROI_V1_2]);
           else
               allData{1,t}.allhighcorrROIs1_2 = unique([allData{1,t}.uniqueHighCorrROI_M1_2]);
           end
        else 
            if isfield(allData{1,t}, 'uniqueHighCorrROI_V1_2')
                allData{1,t}.allhighcorrROIs1_2 = unique([ allData{1,t}.uniqueHighCorrROI_V1_2]);
            else
                allData{1,t}.allhighcorrROIs1_2 = [];
            end
        end
    end
end


%% add cols to allRespROIs
idx = 1
for t = 1:length(allData)
    for roi = 1:length(allData{1,t}.resp_ROIs)
       % allRespROIs(idx, 14) = ismember(r,  allData{1,t}.allhighcorrROIs1_3);
        %allRespROIs(idx, 15) = ismember(r,  allData{1,t}.allhighcorrROIs1_5);
        allRespROIs(idx, 16) = ismember(r,  allData{1,t}.allhighcorrROIs1_2);
        idx = idx + 1
    end
end

%% Make C1 vs C2 plot by HC group 

% Taken from diss_fig4_v2 (PCA version 2)

% C1 and C2 scores colored by HC1/6, HC1/5, HC1/3
% Note that SC cells are not plotted to increase parsability. 

HC1_6 = find(allRespROIs(:,13) == 1)
HC1_5 = find(allRespROIs(:,15) == 1)
HC1_3 = find(allRespROIs(:,14) == 1)
HC1_2 = find(allRespROIs(:,16) == 1)
figure;
hold on 
plot(score_all(HC1_6, 1), score_all(HC1_6, 2), 'o', 'Color', [0 0 0]) %black is original 1/6 threshold
plot(score_all(HC1_5, 1), score_all(HC1_5, 2), 'o', 'Color', [73/255, 0 146/255]) %purple is 1/5 threshold
%plot(score_all(HC1_3, 1), score_all(HC1_3, 2), 'o', 'Color', [146/255 0 0]) %red is 1/3 threshold
hold off
legend({'HC 1/6', 'HC 1/5', 'HC 1/3'})
xlabel('Component 1')
ylabel('Component 2')
set(gca, 'FontSize', 30)
saveas(gcf, 'PCA scores C1 C2 by HC mult thresholds', 'epsc2')
saveas(gcf, 'PCA scores C1 C2 by HC mult thresholds', 'png')


%% Where are the 223 ROIs with HC status > 1/5?

%what stage?
HHC_st = allRespROIs(HC1_5, 3)
length(find(HHC_st == 49))
% 130 from 46, 93 from 49
HHC_st = allRespROIs(HC1_6, 3)
length(find(HHC_st == 46))
% for 1/6, 366 from 49 and 453 from 46

% are they coming from all tads or are a few weirdos driving this?
HHC_t = allRespROIs(HC1_5, 1)
HHC_t = allRespROIs(HC1_6, 1)
length(unique(HHC_t))
%yes a few weirdos are driving this: t=20, 21, 26. 2 46s and 1 49 for HC1_5.
%with 1/6 threshold, 28 tads have HC ROIs. 

