%% analyse other stimuli

%how many experiments/ROIs have multiple stims?

for t =1: length(tadcluster)
    u = unique(tadcluster{1,t}.stimorder)
    len(t) = length(u)
end

for t = 1:length(tadcluster)
    if len(t) > 4
        respROI(t) = length(tadcluster{1,t}.resp_ROIs)
    end
end
all_respROI = sum(respROI)

% 94 ROIS across 7 tads 
