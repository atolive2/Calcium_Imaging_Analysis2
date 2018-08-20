%% Assess high corr ROIs
% start with load('tadcluster_xcorr_assessment_20170907.mat')

% %% Make a struct for high corr ROIs and non high corr ROIs
% 
% for t = 1:length(allData)
%     for r = 1:length(allData{1,t}.resp_ROIs)
%         if ismember(allData{1,t}.resp_ROIs(r), allData{1,t}.correlated_ROIs_dff0_MS_common_AROI)
%             highcorrROIs(c1,1) = t
%             highcorrROIs(c1,2) = allData{1,t}.resp_ROIs(r)
%             highcorrROIs(c1,3) = sum
%         end
%     end
% end

%% define all highcorr ROIs

% all
tmp=[];
for t = 1:length(allData)
    if isfield(allData{1,t}, 'correlated_ROIs_dff0_all_common_AROI')
    for r = 1:length(allData{1,t}.correlated_ROIs_dff0_all_common_AROI)
        tmp = [tmp; cell2mat(allData{1,t}.correlated_ROIs_dff0_all_common_AROI(r))];
    end
    allData{1,t}.uniqueHighCorrROI_all =  unique(tmp)
    tmp = [];
    end
end

%multi
tmp=[];
for t = 1:length(allData)
    if isfield(allData{1,t}, 'correlated_ROIs_dff0_MS_common_AROI')
    for r = 1:length(allData{1,t}.correlated_ROIs_dff0_MS_common_AROI)
        tmp = [tmp; cell2mat(allData{1,t}.correlated_ROIs_dff0_MS_common_AROI(r))];
    end
    allData{1,t}.uniqueHighCorrROI_MS =  unique(tmp)
    tmp = [];
    end
end

% mech
tmp=[];
for t = 1:length(allData)
    if isfield(allData{1,t}, 'correlated_ROIs_dff0_M_common_AROI')
    for r = 1:length(allData{1,t}.correlated_ROIs_dff0_M_common_AROI)
        tmp = [tmp; cell2mat(allData{1,t}.correlated_ROIs_dff0_M_common_AROI(r))];
    end
    allData{1,t}.uniqueHighCorrROI_M =  unique(tmp)
    tmp = [];
    end
end

% visual
tmp=[];
for t = 1:length(allData)
    if isfield(allData{1,t}, 'correlated_ROIs_dff0_V_common_AROI')
    for r = 1:length(allData{1,t}.correlated_ROIs_dff0_V_common_AROI)
        tmp = [tmp; cell2mat(allData{1,t}.correlated_ROIs_dff0_V_common_AROI(r))];
    end
    allData{1,t}.uniqueHighCorrROI_V =  unique(tmp)
    tmp = [];
    end
end


% no stim
tmp=[];
for t = 1:length(allData)
    if isfield(allData{1,t}, 'correlated_ROIs_dff0_N_common_AROI')
    for r = 1:length(allData{1,t}.correlated_ROIs_dff0_N_common_AROI)
        tmp = [tmp; cell2mat(allData{1,t}.correlated_ROIs_dff0_N_common_AROI(r))];
    end
    allData{1,t}.uniqueHighCorrROI_N =  unique(tmp)
    tmp = [];
    end
end


%% Are high corr ROIs more active overall?

highcorrROIs = [];
nothighcorrROIs = [];
for t = 1:length(allData)
    if isfield(allData{1,t}, 'correlated_ROIs_dff0_MS_common_AROI')
        fprintf(num2str(t))
        for r = 1:length(allData{1,t}.resp_ROIs)
            if ismember(allData{1,t}.resp_ROIs(r), allData{1,t}.uniqueHighCorrROI)
                highcorrROIs = [highcorrROIs, allData{1,t}.sum_responses(allData{1,t}.resp_ROIs(r))];
            else
                nothighcorrROIs = [nothighcorrROIs, allData{1,t}.sum_responses(allData{1,t}.resp_ROIs(r))];
            end
        end
    end
end

numHC = length(highcorrROIs)
numNHC = length(nothighcorrROIs)

figure;
hold on 
ecdf(highcorrROIs)
pause
ecdf(nothighcorrROIs)
hold off
title('Total number of responses by high/not high corr st 49')

group = [ones(length(highcorrROIs), 1); 2*ones(length(nothighcorrROIs), 1)];
[p, tbl, stats] = kruskalwallis([highcorrROIs, nothighcorrROIs], group)

%% Is MS index different between high corr and not high corr?
highcorrROIs = [];
nothighcorrROIs = [];
for t = 1:length(allData)
    if isfield(allData{1,t}, 'correlated_ROIs_dff0_MS_common_AROI')
        fprintf(num2str(t))
        for r = 1:length(allData{1,t}.resp_ROIs)
            if ismember(allData{1,t}.resp_ROIs(r), allData{1,t}.uniqueHighCorrROI)
                highcorrROIs = [highcorrROIs, allData{1,t}.MSenh_peak(allData{1,t}.resp_ROIs(r))];
            else
                nothighcorrROIs = [nothighcorrROIs, allData{1,t}.MSenh_peak(allData{1,t}.resp_ROIs(r))];
            end
        end
    end
end

numHC = length(highcorrROIs)
numNHC = length(nothighcorrROIs)

figure;
hold on 
ecdf(highcorrROIs)
pause
ecdf(nothighcorrROIs)
hold off
title('MSIn by high/not high corr st 49')

group = [ones(length(highcorrROIs), 1); 2*ones(length(nothighcorrROIs), 1)];
[p, tbl, stats] = kruskalwallis([highcorrROIs, nothighcorrROIs], group)

%% Is average peak-multi different?

highcorrROIs = [];
nothighcorrROIs = [];
for t = 1:length(allData)
    if isfield(allData{1,t}, 'correlated_ROIs_dff0_MS_common_AROI')
        fprintf(num2str(t))
        for r = 1:length(allData{1,t}.resp_ROIs)
            if ismember(allData{1,t}.resp_ROIs(r), allData{1,t}.uniqueHighCorrROI)
                highcorrROIs = [highcorrROIs, allData{1,t}.peak_avg(1, allData{1,t}.resp_ROIs(r))];
            else
                nothighcorrROIs = [nothighcorrROIs, allData{1,t}.peak_avg(1, allData{1,t}.resp_ROIs(r))];
            end
        end
    end
end

numHC = length(highcorrROIs)
numNHC = length(nothighcorrROIs)

figure;
hold on 
ecdf(highcorrROIs)
pause
ecdf(nothighcorrROIs)
hold off
title('Avg multi peak by high/not high corr st 49')

group = [ones(length(highcorrROIs), 1); 2*ones(length(nothighcorrROIs), 1)];
[p, tbl, stats] = kruskalwallis([highcorrROIs, nothighcorrROIs], group)

%% Is average peak-unimax different?

highcorrROIs = [];
nothighcorrROIs = [];
for t = 1:length(allData)
    if isfield(allData{1,t}, 'correlated_ROIs_dff0_MS_common_AROI')
        fprintf(num2str(t))
        for r = 1:length(allData{1,t}.resp_ROIs)
            if ismember(allData{1,t}.resp_ROIs(r), allData{1,t}.uniqueHighCorrROI)
                highcorrROIs = [highcorrROIs, max([allData{1,t}.peak_avg(2, allData{1,t}.resp_ROIs(r)), allData{1,t}.peak_avg(3, allData{1,t}.resp_ROIs(r))])];
            else
                nothighcorrROIs = [nothighcorrROIs, max([allData{1,t}.peak_avg(2, allData{1,t}.resp_ROIs(r)), allData{1,t}.peak_avg(3, allData{1,t}.resp_ROIs(r))])];
            end
        end
    end
end

numHC = length(highcorrROIs)
numNHC = length(nothighcorrROIs)

figure;
hold on 
ecdf(highcorrROIs)
pause
ecdf(nothighcorrROIs)
hold off
title('Avg uni max peak by high/not high corr st 49')

group = [ones(length(highcorrROIs), 1); 2*ones(length(nothighcorrROIs), 1)];
[p, tbl, stats] = kruskalwallis([highcorrROIs, nothighcorrROIs], group)

%% Are high corr cells spontaneously active together?
% In how many trials are the high corr cells doing the same thing?










