%% Correct data for diss figs
% started with allData_20180320 in D:\Torrey_calcium_imaging\compare_46-49\analysis_Feb 2018\corrected_for_badtrials


% only include peaks and onset times from trials with a response. 
t = 20
tad =[];
for t = 1:length(allData)
    stims = unique(allData{1,t}.stimorder); %find all types of stims
    for s = 1:length(stims)
        trials = allData{1,t}.stimorder == stims(s); %find each trial of that stimtype
        for r = 1:size(allData{1,t}.peak_bytrial2, 1) %each ROI
            tmp = allData{1,t}.peak_bytrial2(r,trials);% > 0.1
            tmp2 = tmp > 0.1 & tmp < 10;
            if sum(tmp2) > 0
            allData{1,t}.avg_peak2(s,r) = mean(tmp(tmp2));
            tmp3 = allData{1,t}.onsettime(r,trials);
            allData{1,t}.avg_onsettime2(s,r) = mean(tmp3(tmp2));
            allData{1,t}.std_onsettime2(s,r) = std(tmp3(tmp2)); 
            else
                allData{1,t}.avg_peak2(s,r) = NaN;
            
            allData{1,t}.avg_onsettime2(s,r) = NaN;
            allData{1,t}.std_onsettime2(s,r) = NaN;
            tad=[tad, t];
            end
            clear('tmp', 'tmp2', 'tmp3')
        end
    end
end

% for t = 1:length(allData)
%     stims = unique(allData{1,t}.stimorder); %find all types of stims
%     for s = 1:length(stims)
%         trials = allData{1,t}.stimorder == stims(s); %find each trial of that stimtype
%         for r = 1:size(allData{1,t}.boolean_response, 1) %each ROI
%             include = allData{1,t}.boolean_response(r,:) & trials;
%             if length(include) > 0
%                 allData{1,t}.avg_onsettime2(s,r) = mean(allData{1,t}.onsettime(r, include));
%             end
%         end
%     end
% end
% 
% for t = 1:length(allData)
%     stims = unique(allData{1,t}.stimorder); %find all types of stims
%     for s = 1:length(stims)
%         trials = allData{1,t}.stimorder == stims(s); %find each trial of that stimtype
%         for r = 1:size(allData{1,t}.boolean_response, 1) %each ROI
%             include = allData{1,t}.boolean_response(r,:) & trials;
%             if length(include) > 0
%                 allData{1,t}.std_onsettime2(s,r) = std(allData{1,t}.onsettime(r, include));
%             end
%         end
%     end
% end

%% recalculate MSI with this data
for t = 1:length(allData)
    allData{1,t}.MSInd_peak2 = calc_MSIndex(allData{1,t}.avg_peak2);
    allData{1,t}.MSInd_onsettime2 = calc_MSIndex(allData{1,t}.avg_onsettime2);
end


