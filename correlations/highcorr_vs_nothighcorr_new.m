%% Assess high corr ROIs
% start with load('tadcluster_xcorr_assessment_20170907.mat')

% %% Make a struct for high corr ROIs and non high corr ROIs
% 
% for t = 1:length(tadcluster)
%     for r = 1:length(tadcluster{1,t}.resp_ROIs)
%         if ismember(tadcluster{1,t}.resp_ROIs(r), tadcluster{1,t}.correlated_ROIs_dff0_MS_common_AROI)
%             highcorrROIs(c1,1) = t
%             highcorrROIs(c1,2) = tadcluster{1,t}.resp_ROIs(r)
%             highcorrROIs(c1,3) = sum
%         end
%     end
% end

%% define all highcorr ROIs

for t = 1:length(tadcluster)
    if isfield(tadcluster{1,t}, 'correlated_ROIs_dff0_MS_common_AROI')
    for r = 1:length(tadcluster{1,t}.correlated_ROIs_dff0_MS_common_AROI)
        tmp = [tmp; cell2mat(tadcluster{1,t}.correlated_ROIs_dff0_MS_common_AROI(r))];
    end
    tadcluster{1,t}.uniqueHighCorrROI =  unique(tmp)
    tmp = [];
    end
end


%% Are high corr ROIs more active overall?

highcorrROIs = [];
nothighcorrROIs = [];
for t = 1:length(tadcluster)
    if isfield(tadcluster{1,t}, 'correlated_ROIs_dff0_MS_common_AROI')
        fprintf(num2str(t))
        for r = 1:length(tadcluster{1,t}.resp_ROIs)
            if ismember(tadcluster{1,t}.resp_ROIs(r), tadcluster{1,t}.uniqueHighCorrROI)
                highcorrROIs = [highcorrROIs, tadcluster{1,t}.sum_responses(tadcluster{1,t}.resp_ROIs(r))];
            else
                nothighcorrROIs = [nothighcorrROIs, tadcluster{1,t}.sum_responses(tadcluster{1,t}.resp_ROIs(r))];
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
for t = 1:length(tadcluster)
    if isfield(tadcluster{1,t}, 'correlated_ROIs_dff0_MS_common_AROI')
        fprintf(num2str(t))
        for r = 1:length(tadcluster{1,t}.resp_ROIs)
            if ismember(tadcluster{1,t}.resp_ROIs(r), tadcluster{1,t}.uniqueHighCorrROI)
                highcorrROIs = [highcorrROIs, tadcluster{1,t}.MSenh_peak(tadcluster{1,t}.resp_ROIs(r))];
            else
                nothighcorrROIs = [nothighcorrROIs, tadcluster{1,t}.MSenh_peak(tadcluster{1,t}.resp_ROIs(r))];
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
for t = 1:length(tadcluster)
    if isfield(tadcluster{1,t}, 'correlated_ROIs_dff0_MS_common_AROI')
        fprintf(num2str(t))
        for r = 1:length(tadcluster{1,t}.resp_ROIs)
            if ismember(tadcluster{1,t}.resp_ROIs(r), tadcluster{1,t}.uniqueHighCorrROI)
                highcorrROIs = [highcorrROIs, tadcluster{1,t}.peak_avg(1, tadcluster{1,t}.resp_ROIs(r))];
            else
                nothighcorrROIs = [nothighcorrROIs, tadcluster{1,t}.peak_avg(1, tadcluster{1,t}.resp_ROIs(r))];
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
for t = 1:length(tadcluster)
    if isfield(tadcluster{1,t}, 'correlated_ROIs_dff0_MS_common_AROI')
        fprintf(num2str(t))
        for r = 1:length(tadcluster{1,t}.resp_ROIs)
            if ismember(tadcluster{1,t}.resp_ROIs(r), tadcluster{1,t}.uniqueHighCorrROI)
                highcorrROIs = [highcorrROIs, max([tadcluster{1,t}.peak_avg(2, tadcluster{1,t}.resp_ROIs(r)), tadcluster{1,t}.peak_avg(3, tadcluster{1,t}.resp_ROIs(r))])];
            else
                nothighcorrROIs = [nothighcorrROIs, max([tadcluster{1,t}.peak_avg(2, tadcluster{1,t}.resp_ROIs(r)), tadcluster{1,t}.peak_avg(3, tadcluster{1,t}.resp_ROIs(r))])];
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










