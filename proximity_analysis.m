%% Analyze "Proximity" for cells by primary modality
% this code is based exclusively on peak analysis. 

myFolder = 'F:/Calcium_Imaging_Analysis/tadpoles_byexp/'; % May need to correct this.
mkdir([myFolder 'figures']);
if ~isdir(myFolder)
	errorMessage = sprintf('Error: The following folder does not exist:\n%s', myFolder);
	uiwait(warndlg(errorMessage));
	return;
end
filePattern = fullfile(myFolder, 'exp*.mat');
matFiles = dir(filePattern)

% collect all data into a single cell array of structs.
[ tadpole ] = get_matFiles( myFolder, matFiles )

% empirically determine the cut of for response/no response
%cutoff = calc_response_threshold(
[ cutoff, std_dev, avg_alltrial ] = calc_response_threshold( tadpole{1}.meanpeak_bytrial, 2 )

% Generate logical array of responses, based on peak
for i = 1:length(tadpole)
    [ tadpole{i}.cat_bytrial, tadpole{i}.sum_byROI ] = cat_response( tadpole{i}.meanpeak_bytrial, 0.2 )
end

% Find only ROIs that respond at least once
for i = 1:length(tadpole)
    tadpole{i}.responsiveROIs = find(tadpole{i}.sum_byROI);
end

% Get the proportion of responses to each stimtype. 
for t = 1:length(tadpole)
    for j = 1:size(tadpole{t}.stimmask,2)
        total_pres = sum(tadpole{t}.stimmask(:,j));
        trials = find(tadpole{t}.stimmask(:,j));
        for i = 1:length(tadpole{t}.responsiveROIs)
            roi = tadpole{t}.responsiveROIs(i);
            tadpole{t}.total_resp(i,j) = sum(tadpole{t}.cat_bytrial(roi,trials));
            tadpole{t}.propResp_bystim(i,j) = tadpole{t}.total_resp(i,j)/total_pres;
        end
        clear('total_pres', 'trials')
    end
end        
%plot to check
for t = 1:length(tadpole)
    figure;
    plot(tadpole{t}.propResp_bystim)
    legend('1', '2', '3', '4', '5', '6', '7', '8', '9')
end

for t = 1:length(tadpole)
    figure;
    hold on
    for i = 1:size(tadpole{t}.df_f0,1)
        for j = 1:size(tadpole{t}.df_f0,2)
            plot(tadpole{t}.df_f0{i,j}(1,:)')
        end
    end
    hold off
    %legend('1', '2', '3', '4', '5', '6', '7', '8', '9')
end

for t = 1:length(tadpole)
    figure;
    hold on
    for i = 1:size(tadpole{t}.signal,1)
        for j = 1:size(tadpole{t}.signal,2)
            if tadpole{t}.signal_good(i,j)
                plot(tadpole{t}.signal{i,j}(1,:)')
            end
        end
    end
    hold off
    %legend('1', '2', '3', '4', '5', '6', '7', '8', '9')
end

% Using only the max stim set (1-4), calculate a "primary" modality.
% primary modality = modality with highest proportion responses 
% code doesn't work. Dimension problem in max(tadpole{1}.propResp_bystim(i,1:4))
for i = 1:size(tadpole{1}.propResp_bystim,1)
    if find(tadpole{1}.propResp_bystim(i,1:4) == max(tadpole{1}.propResp_bystim(i,1:4))) == 1
        [max, tadpole{1}.primary_modality_high(i)] = max(tadpole{1}.propResp_bystim(i,1:4), [], 2)
    else 
        tadpole{1}.primary_modality_high(i) = 0
    end
end


% identify "good" traces, based on having a + signal
for t = 1:length(tadpole)
    for i = 1:size(tadpole{t}.signal,1)
        for j = 1:size(tadpole{t}.signal,2)
            avg_signal(i,j) = mean(tadpole{t}.signal{i,j}(1,:));
            signal_good(i,j) = avg_signal(i,j) > 0.01;
%             if mean(tadpole{t}.signal{i,j}(1,:)) < 0.01
%                 signal_good(i,j) = 0;
%             else
%                 signal_good(i,j) = 1;
%             end
        end
    end
    tadpole{t}.signal_good = logical(signal_good);
end    

% identify traces that should be used for analysis.
% don't use traces based on signal_good, and df_f0 too big or small

% df_f0 ok
for t = 1:length(tadpole)
    for i = 1:size(tadpole{t}.df_f0,1)
        for j = 1:size(tadpole{t}.df_f0,2)
            tadpole{t}.df_f0_good(i,j) = (max(tadpole{t}.df_f0{i,j}(1,:)) < 10) && (min(tadpole{t}.df_f0{i,j}(1,:)) > -2);
        end
    end
end
% combine signal_good and df_f0_good
for t = 1:length(tadpole)
    for i = 1:size(tadpole{t}.df_f0,1)
        for j = 1:size(tadpole{t}.df_f0,2)
            tadpole{t}.use_me(i,j) = tadpole{t}.df_f0_good(i,j) && tadpole{t}.signal_good(i,j)
        end
    end
end















