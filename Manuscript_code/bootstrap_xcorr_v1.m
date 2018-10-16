%% Shuffle and bootstrap xcorr values

% start with allData containing respROIs_ST, smoothed_trunc, and the list
% of trials in order. 

%% Shuffle the trials, calculate xcorr and report maxR, lagmaxR and 0lagR into corresponding fields 

its = 100; % how many bootstrap iterations?
stims = [1 2 3 4] % MS, V, M, NS
% allData(k).incl_trials contains the trial numbers for each ROI pair for
% each stim type
tic
for k = 1:length(allData)
    if length(allData(k).respROIs_ST) > 5
        set_lag = length(allData(k).smoothed_trunc{1,1});
        for s = 1:length(stims) %all 4 stims
            for r1 = 1:length(allData(k).respROIs_ST) %all responding ROIs
                for r2 = 1:length(allData(k).respROIs_ST) %all responding ROIs 
                    all_tr = allData(k).incl_trials{s, r1, r2};   
                    if ~isempty(all_tr)
                        for b = 1:its % number of iterations
                            r1_tr = all_tr(randperm(length(all_tr)));
                            r2_tr = all_tr(randperm(length(all_tr)));
                            tmp1 = [];
                            for tr = 1:length(r1_tr)
                                tmp1 = [tmp1; allData(k).smoothed_trunc{allData(k).respROIs_ST(r1), r1_tr(tr)}];
                            end
                            tmp2 = [];
                            for t2 = 1:length(r2_tr)
                                tmp2 = [tmp2; allData(k).smoothed_trunc{allData(k).respROIs_ST(r2), r2_tr(t2)}];
                            end
                            [boot_xcorr(k).R{s, r1, r2, b}, boot_xcorr(k).lag{s, r1, r2, b}] = xcorr(tmp1, tmp2, set_lag, 'coeff');
                            [boot_xcorr(k).maxR{s, r1, r2}(b), boot_xcorr(k).lagmaxR{s, r1, r2}(b)] = max(boot_xcorr(k).R{s, r1, r2, b}(:,:));
                            boot_xcorr(k).lag0R{s, r1, r2}(b) = boot_xcorr(k).R{s, r1, r2, b}(set_lag+1);
                        end
                        clear('r1_tr', 'r2_tr')
                    end 
                    clear('all_tr')
                end
            end
        end
        toc
    end
end

                    
%% Generate 95% confidence intervals from the bootstrap values for each ROI pair
tic
for k = 1:length(allData)
    if length(allData(k).respROIs_ST) > 5
        for s = 1:length(stims) %all 4 stims
            for r1 = 1:length(allData(t).respROIs_ST) %all responding ROIs
                for r2 = 1:length(allData(t).respROIs_ST) %all responding ROIs 
                    % maxR calculations
                    boot_xcorr(k).maxR_info{s,r1,r2}(1) = mean(boot_xcorr(k).maxR{s,r1,r2}(:)); %mean
                    boot_xcorr(k).maxR_info{s,r1,r2}(2) = std(boot_xcorr(k).maxR{s,r1,r2}(:)); %std
                    boot_xcorr(k).maxR_info{s,r1,r2}(3) = boot_xcorr(k).maxR_info{s,r1,r2}(1) + (2*boot_xcorr(k).maxR_info{s,r1,r2}(2)); %upper 95% CI
                    boot_xcorr(k).maxR_info{s,r1,r2}(4) = boot_xcorr(k).maxR_info{s,r1,r2}(1) - (2*boot_xcorr(k).maxR_info{s,r1,r2}(2)); %lower 95% CI
                    % lag of maxR calculations
                    boot_xcorr(k).lagmaxR_info{s,r1,r2}(1) = mean(boot_xcorr(k).lagmaxR{s,r1,r2}(:)); %mean
                    boot_xcorr(k).lagmaxR_info{s,r1,r2}(2) = std(boot_xcorr(k).lagmaxR{s,r1,r2}(:)); %std
                    boot_xcorr(k).lagmaxR_info{s,r1,r2}(3) = boot_xcorr(k).lagmaxR_info{s,r1,r2}(1) + (2*boot_xcorr(k).lagmaxR_info{s,r1,r2}(2)); %upper 95% CI
                    boot_xcorr(k).lagmaxR_info{s,r1,r2}(4) = boot_xcorr(k).lagmaxR_info{s,r1,r2}(1) - (2*boot_xcorr(k).lagmaxR_info{s,r1,r2}(2)); %lower 95% CI
                    % lag 0 R calculations
                    boot_xcorr(k).lag0R_info{s,r1,r2}(1) = mean(boot_xcorr(k).lag0R{s,r1,r2}(:)); %mean
                    boot_xcorr(k).lag0R_info{s,r1,r2}(2) = std(boot_xcorr(k).lag0R{s,r1,r2}(:)); %std
                    boot_xcorr(k).lag0R_info{s,r1,r2}(3) = boot_xcorr(k).lag0R_info{s,r1,r2}(1) + (2*boot_xcorr(k).lag0R_info{s,r1,r2}(2)); %upper 95% CI
                    boot_xcorr(k).lag0R_info{s,r1,r2}(4) = boot_xcorr(k).lag0R_info{s,r1,r2}(1) - (2*boot_xcorr(k).lag0R_info{s,r1,r2}(2)); %lower 95% CI
                end
            end
        end
        toc
    end
end


%% Compare actual data point to 95% confidence interval for each ROI pair 

for k = 1:length(allData)
    if length(allData(k).respROIs_ST) > 5
        for s = 1:length(stims) %all 4 stims
            for r1 = 1:length(allData(t).respROIs_ST) %all responding ROIs
                for r2 = 1:length(allData(t).respROIs_ST) %all responding ROIs 
                    % 0 if actual is outside of CI
                    if (allData(k).respROIdff0_maxR_sq(s, r1, r2) > boot_xcorr(k).maxR_info{s, r1, r2}(4)) && (allData(k).respROIdff0_maxR_sq(s, r1, r2) < boot_xcorr(k).maxR_info{s, r1, r2}(3))
                        boot_xcorr(k).comparetoCI{s, r1, r2}(1) = 1; % actual data is outside CI
                        if allData(k).respROIdff0_maxR_sq(s, r1, r2) < boot_xcorr(k).maxR_info{s, r1, r2}(4) % actual less than lower CI bound
                        boot_xcorr(k).comparetoCI{s, r1, r2}(2) = boot_xcorr(k).maxR_info{s, r1, r2}(4); - allData(k).respROIdff0_maxR_sq(s, r1, r2); %positive value for difference
                        boot_xcorr(k).comparetoCI{s, r1, r2}(3) = -1; %ID those smaller than lower CI bound with -1
                        elseif allData(k).respROIdff0_maxR_sq(s, r1, r2) > boot_xcorr(k).maxR_info{s, r1, r2}(3) % actual greater than upper CI bound
                            boot_xcorr(k).comparetoCI{s, r1, r2}(2) = allData(k).respROIdff0_maxR_sq(s, r1, r2) - boot_xcorr(k).maxR_info{s, r1, r2}(3); % positive value for difference
                            boot_xcorr(k).comparetoCI{s, r1, r2}(3) = 1; %ID those larger than CI with +1
                        end
                    else % actual is inside CI
                        tmp(1) = boot_xcorr(k).maxR_info{s, r1, r2}(4) - allData(k).respROIdff0_maxR_sq(s, r1, r2); % negative value for difference to lower CI
                        tmp(2) = allData(k).respROIdff0_maxR_sq(s, r1, r2) - boot_xcorr(k).maxR_info{s, r1, r2}(3); % negative value for difference to upper  CI
                        boot_xcorr(k).comparetoCI{s, r1, r2}(2) = min(abs(tmp)); %get the negative value closer to 0 by abs val first
                        boot_xcorr(k).comparetoCI{s, r1, r2}(3) = 0; %ID those inside CI with 0
                    end
                end
            end
        end
    end
end
toc

%% Plot difference between edge of 95% confidence interval and the actual datapoint for each tad
% subplot each stim type
% histogram
% postive = outside (significantly diff from noise)
% negative = inside distribution (not sig diff from noise)

% for k = 1:length(allData)
%     if length(allData(k).respROIs_ST) > 5
%         for s = 1:length(stims) %all 4 stims


%% Plot proportion of pairs outside 95% confidence interval


