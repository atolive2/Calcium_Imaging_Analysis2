%% Measure onset time variation by modality type

% this assumes you are starting with tadcluster_data_20170814.mat which
% contains only tadcluster, with responding ROIs already defined. 

% get the onset frame for trials and rois with a response

% if responding ROI && boolean response = 1 then get onset time using
% stepinfo (http://www.mathworks.com/help/control/ref/stepinfo.html)
% insert into a new matrix with tad num, cell num, trial type, peak and peak time
counter = 1;
for t = 1:length(tadcluster)
    for r = 1:length(tadcluster{1,t}.resp_ROIs) %over each responding ROI
        roi = tadcluster{1,t}.resp_ROIs(r);
        for tt = 1:size(tadcluster{1,t}.df_f0,2)
            if tadcluster{1,t}.boolean_response(roi,tt)
                temp_var(counter).stepinfo = stepinfo(tadcluster{1,t}.df_f0{roi, tt});
                temp_var(counter).tad_num = t;
                temp_var(counter).roi_num = roi;
                temp_var(counter).trial_num = tt;
                temp_var(counter).trial_type = tadcluster{1,t}.stimorder(tt);
                temp_var(counter).abs_peak = tadcluster{1,t}.peak_bytrial(roi, tt);

                if length(tadcluster{1,t}.df_f0{roi, tt}) <=160
                    temp_var(counter).trial_length = 160;
                elseif length(tadcluster{1,t}.df_f0{roi, tt}) > 160
                    temp_var(counter).trial_length = length(tadcluster{1,t}.df_f0{roi, tt});
                end
                temp_var(counter).peakloc = tadcluster{1,t}.peakloc_bytrial(roi, tt);
                temp_var(counter).abs_peakloc = tadcluster{1,t}.peakloc_bytrial(roi, tt) / (trial_length/7); %earlier stage analysis converted peak_loc into 160 frames per trial for 22.9fps
                temp_var(counter).stimonset = tadcluster{1,t}.stim_onset;
                counter = counter+1;
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%add stim onset time to tadcluster%%%%
% exps 1-9 have onset at 0.5s 
% exps 10+ have onset at 2s
for t = 1:length(tadcluster)
    if tadcluster{1,t}.expnum <=9
        tadcluster{1,t}.stim_onset = 0.5;
    elseif tadcluster{1,t}.expnum >9
        tadcluster{1,t}.stim_onset = 2;
    else
        fprintf('error exp %d, t')
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


tads = [temp_var.tad_num]

[occurrences, entries] = hist(tads, unique(tads))

% get onset time as time to 50% of peak in a frame num
counter = 1;
for t = 1:length(tadcluster)
    for r = 1:length(tadcluster{1,t}.resp_ROIs) %over each responding ROI
        roi = tadcluster{1,t}.resp_ROIs(r);
        for tt = 1:size(tadcluster{1,t}.df_f0,2)
            if tadcluster{1,t}.boolean_response(roi,tt)
                half_peak = tadcluster{1,t}.peak_bytrial(roi, tt) / 2;
                temp_var(counter).timeto50peak = find((tadcluster{1,t}.df_f0{roi, tt} > half_peak) , 3);
                counter = counter+1;
            end
        end
    end
end

% now clean up data to remove any 1s or 2s as the 50% peak (these are
% errors from registration. sometimes the first frame is abnormally bright)
for t = 1:length(temp_var)
    if temp_var(t).timeto50peak(1) == (1 || 2)
        if temp_var(t).timeto50peak(2) == 2 
            temp_var(t).timeto50peak_actual = temp_var(t).timeto50peak(3) / (temp_var(t).trial_length/7);
        else
            temp_var(t).timeto50peak_actual = temp_var(t).timeto50peak(2) / (temp_var(t).trial_length/7);
        end
    else
        temp_var(t).timeto50peak_actual = temp_var(t).timeto50peak(1) / (temp_var(t).trial_length/7);
    end
end

% figure out why there are time to 50% peak vals > 7s
bad_trials = find([temp_var.timeto50peak_actual] > 7)
exps_badtrials = [temp_var(bad_trials).tad_num]
%unique(exps_badtrials) = [10, 12]
%tadcluster{1,10} and {1,12} both have 220 frames per trial
% ok, so trials with 136 frames got multiplied to be out of 160 but trials
% with 219/220 did not. %FIXED ABOVE

% Then compare onset times by modality - first high(crash)/high(90%peak) only
% these are trial types 1(MS), 2(V), 3(Mech), 4(no stim)
for s = 1:4
    trials_touse = find([temp_var.trial_type] == s);
    trial_count = length(trials_touse);
    for t = 1:length(trials_touse)
        peaks{1,s} = [temp_var(trials_touse).timeto50peak_actual];
    end
end

% plot these onset times by modality %%%%%%%%%%%%%DOESN'T WORK
figure;
hold on
boxplot([peaks{1,1}, peaks{1,2}])
boxplot(peaks{1,2})
boxplot(peaks{1,3})
boxplot(peaks{1,4})
hold off

% histogram each modality
figure;
subplot(1,3,1)
hist(peaks{1,1})
title('multi')
subplot(1,3,2)
hist(peaks{1,2})
title('vis')
subplot(1,3,3)
hist(peaks{1,3})
title('mech')

%% More useful histogram - time from stim onset
% get difference between stim onset time and time to 50% peak
for t = 1:length(temp_var)
    temp_var(t).diffStimOnset_RespOnset = temp_var(t).timeto50peak_actual - temp_var(t).stimonset;
end

%collect diffStimOnset_RespOnset by stim modality
for s = 1:4
    trials_touse = find([temp_var.trial_type] == s);
    trial_count = length(trials_touse);
    for t = 1:length(trials_touse)
        peaks_diffSO_RO{1,s} = [temp_var(trials_touse).diffStimOnset_RespOnset];
    end
end

% histogram each modality
bin_edges = [-2 -1.5 -1 -.5 0 .5 1 1.5 2 2.5 3 3.5 4 4.5 5 5.5 6]
figure;
subplot(1,3,1)
hist(peaks_diffSO_RO{1,1}, bin_edges)
axis([-2 6 0 200])
title('multi')
ylabel('counts')
xlabel('time after stim onset (s)')
subplot(1,3,2)
hist(peaks_diffSO_RO{1,2}, bin_edges)
axis([-2 6 0 200])
title('vis')
ylabel('counts')
xlabel('time after stim onset (s)')
subplot(1,3,3)
hist(peaks_diffSO_RO{1,3}, bin_edges)
axis([-2 6 0 200])
title('mech')
ylabel('counts')
xlabel('time after stim onset (s)')
suptitle('Response Onset Time by Stimulus Modality')
saveas(gcf, 'Response Onset Time by Stimulus Modality', 'png')

% Remove trials that are negative (e.g. start before stim onset)
for i = 1:length(peaks_diffSO_RO)
    idx = find(peaks_diffSO_RO{1,i} >= 0)
    peaks_diffSO_RO_pos{1,i} = peaks_diffSO_RO{1,i}(idx)
    clear('idx')
end

% stat test if all positive responses come from the same distribution 
% Kruskal-Wallis test. Tests if multiple samples are all drawn from the same 
% populations (or equivalently, from different populations with the same 
% distribution), against the alternative that they are not all drawn from the same population.
Data = cat(2, peaks_diffSO_RO_pos{1,1}, peaks_diffSO_RO_pos{1,2}, peaks_diffSO_RO_pos{1,3})
Groups = [ones(length(peaks_diffSO_RO_pos{1,1}),1); 2*ones(length(peaks_diffSO_RO_pos{1,2}),1); 3*ones(length(peaks_diffSO_RO_pos{1,3}),1)]
[p,tbl,stats] = kruskalwallis(Data, Groups)
c = multcompare(stats)