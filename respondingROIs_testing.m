%% responding ROIs testing

% threshold based on difference between no stim and all stims
%based on meanpeak_bytrial
ind_nostim = find(tadpole.stimorder == 4)
ind_stim = find(tadpole.stimorder ~= 4)

for i = 1:size(tadpole.meanpeak_bytrial,1) %over each ROI
    vals_nostim(i,:)= cell2mat(tadpole.meanpeak_bytrial(i, ind_nostim));
    vals_stim(i,:)= cell2mat(tadpole.meanpeak_bytrial(i, ind_stim));
end

for i = 1:size(vals_stim,1)
    [h(i), p(i)] = ttest2(vals_nostim(i,:), vals_stim(i,:))
end

% oops. h = 1 in 6 cases. Too strigent.

% threshold based on sum of boolean response threshold (so # of times it
% responded
m = mean(tadpole.sum_responses)
s = std(tadpole.sum_responses)
threshold = ceil(m - s) % one stdev below mean
countROIs = length(find(tadpole.sum_responses >= threshold))
respROIs = find(tadpole.sum_responses >= threshold)

%test plot
[unisensory, uniloc] = max(tadpole.peak_avg(2:3,respROIs),[],1); %get val and loc of max unisensory
% area_V = horzcat(area_V, [unisensory(uniloc==1); tadpole.area_avg(uniloc==1)]); % where max unisensory = V (1)
% area_M = horzcat(area_M, [unisensory(uniloc==2); tadpole.area_avg(uniloc==2)]); % where max unisensory = M (2)
% 
% uni=tadpole.peak_avg(max
multi=tadpole.peak_avg(1,respROIs);
figure;
plot(unisensory, multi, 'o')

% from get_respondingROIs
% [ boolean_response, sum_responses ] = get_respondingROIs( area, peak, peakloc )
 [ boolean_response, sum_responses ] = get_respondingROIs(tadpole.area_bytrial, tadpole.meanpeak_bytrial, tadpole.peakloc_bytrial)
m = mean(sum_responses)
s = std(sum_responses)
threshold = ceil(m - s) % one stdev below mean
%threshold = m
countROIs = length(find(sum_responses >= threshold))
respROIs_bool = find(sum_responses >= threshold)

%test plot
[unisensory2, uniloc2] = max(tadpole.peak_avg(2:3,respROIs_bool),[],1); %get val and loc of max unisensory
% area_V = horzcat(area_V, [unisensory(uniloc==1); tadpole.area_avg(uniloc==1)]); % where max unisensory = V (1)
% area_M = horzcat(area_M, [unisensory(uniloc==2); tadpole.area_avg(uniloc==2)]); % where max unisensory = M (2)
% 
% uni=tadpole.peak_avg(max
multi2=tadpole.peak_avg(1,respROIs_bool);
figure;
plot(unisensory2, multi2, 'o')

%% new idea - by proportion responding to each stim
% calculate boolean response sum for each stimulus type
 
order = unique(tadpole.stimorder)
for i = 1:length(order)
    stim = order(i)
    cols = find(tadpole.stimorder == stim)
    sum_resp_bystim(:,i) = sum(tadpole.boolean_response(:,cols),2)
    num_of_stim_pres(i) = length(cols)
end
for s = 1:length(num_of_stim_pres)
    vector = sum_resp_bystim(:,s)
    scalar = num_of_stim_pres(s)
    prop_resp_bystim(:,s) = vector ./ scalar
    clear('vector', 'scalar')
end
single_resp = [1 1 1 1] ./ num_of_stim_pres

% eliminate ROIs with:
    % 1 or fewer responses to all stimuli
    % number of responses to 4 (no stim) is greater than num responses to
        % all other stimuli
included_rois = zeros(size(prop_resp_bystim,1),1);
for i = 1:size(prop_resp_bystim,1)
    if (prop_resp_bystim(i,1) > single_resp(1,1)) || (prop_resp_bystim(i,2) > single_resp(1,2)) || (prop_resp_bystim(i,3) > single_resp(1,3))
        if (prop_resp_bystim(i,1) > prop_resp_bystim(i,4)) || (prop_resp_bystim(i,1) > prop_resp_bystim(i,2)) || (prop_resp_bystim(i,3) > prop_resp_bystim(i,4))
            included_rois(i) = 1
        end
    end
end
roi_count = sum(included_rois)

included_rois_log = logical(included_rois)
peak_included_rois = [];
for i = 1:length(included_rois_log)
    if included_rois_log(i)
        peak_included_rois = [peak_included_rois, tadpole.peak_avg(:,i)]
    end
end

uni_max = max(peak_included_rois(2:3,:))
multi_toplot = peak_included_rois(1,:)
figure;
plot(uni_max, peak_included_rois(1,:), 'o')
% the above analysis does not eliminate the very low peak vals. 

% Next step: use max peak, not average peak to plot. 
% get max response for each modality, using meanpeak_bytrial
order = unique(tadpole.stimorder)
for i = 1:length(order)
    trials = find(tadpole.stimorder == order(i))
    largest_vals(:,i) = max(cell2mat(tadpole.meanpeak_bytrial(:, trials)),[],2)
end

[uni_largestval modality] = max(largest_vals(:,2:3),[],2)
multi_largestval_toplot = largest_vals(:,1)
primary_vis =[];
primary_mech = [];
for i = 1:length(modality)
    if modality(i)==1
        primary_vis = [primary_vis; uni_largestval(i) multi_largestval_toplot(i)]
    elseif modality(i)==2
        primary_mech = [primary_mech; uni_largestval(i) multi_largestval_toplot(i)]
    else
        sprintf('error')
    end
end
figure;
hold on 
%plot(uni_largestval, multi_largestval_toplot, 'o')
plot(primary_vis(:,1), primary_vis(:,2), 'or')
plot(primary_mech(:,1), primary_mech(:,2), 'og')
x = [0 0.6]
y = [0 0.6]
plot(x, y)
hold off
title('Exp 5 All ROIs max of max uni vs max multi')
xlabel('max unisensory max meanpeak, \DeltaF/F_{0}')
ylabel('multisensory max meanpeak, \DeltaF/F_{0}')

% how many ROIs have max multi > max uni?
multi_larger = length(find(uni_largestval > multi_largestval_toplot))
uni_larger = length(find(uni_largestval < multi_largestval_toplot))

% Above uses all ROIs, so what's the difference when only responding ROIs
% are used (included_rois_log)
primary_vis_respROI = [];
primary_mech_respROI = [];
for i = 1:length(modality)
    if included_rois_log(i)
        if modality(i)==1
            primary_vis_respROI = [primary_vis_respROI; uni_largestval(i) multi_largestval_toplot(i)]
        elseif modality(i)==2
            primary_mech_respROI = [primary_mech_respROI; uni_largestval(i) multi_largestval_toplot(i)]
        else
            sprintf('error')
        end
    end
end

figure;
hold on 
%plot(uni_largestval, multi_largestval_toplot, 'o')
plot(primary_vis_respROI(:,1), primary_vis_respROI(:,2), 'or')
plot(primary_mech_respROI(:,1), primary_mech_respROI(:,2), 'og')
x = [0 0.6]
y = [0 0.6]
plot(x, y)
hold off
title('Exp 5 All ROIs max of max uni vs max multi')
xlabel('max unisensory max meanpeak, \DeltaF/F_{0}')
ylabel('multisensory max meanpeak, \DeltaF/F_{0}')

% how many ROIs have max multi > max uni?
multi_largerV = length(find(primary_vis_respROI(:,1) < primary_vis_respROI(:,2)))
multi_largerM = length(find(primary_mech_respROI(:,1) < primary_mech_respROI(:,2)))
%uni_larger = length(find(uni_largestval < multi_largestval_toplot))

% above elimination still includes a bunch of really small ROIs.
% next step is to eliminate all ROIs with a max of max uni < 0.1 (the
% response threshold)
respROIs_bymaxuniV = find((primary_vis_respROI(:,1) > 0.1) | (primary_vis_respROI(:,2) > 0.1))
respROIs_bymaxuniM = find((primary_mech_respROI(:,1) > 0.1) | (primary_mech_respROI(:,1) > 0.1))

figure;
hold on 
%plot(uni_largestval, multi_largestval_toplot, 'o')
plot(primary_vis_respROI(respROIs_bymaxuniV,1), primary_vis_respROI(respROIs_bymaxuniV,2), 'or')
plot(primary_mech_respROI(respROIs_bymaxuniM,1), primary_mech_respROI(respROIs_bymaxuniM,2), 'ob')
x = [0 0.6]
y = [0 0.6]
plot(x, y, 'k')
hold off
title('Exp 5 All ROIs max of max uni vs max multi, either > 0.1')
xlabel('max unisensory max meanpeak, \DeltaF/F_{0}')
ylabel('multisensory max meanpeak, \DeltaF/F_{0}')
legend('visual', 'mechanosensory', 'Y = X')

%% elimination criteria 5: peak average for max uni and multi are both < 0.1

[unisensory5, uniloc5] = max(tadpole.peak_avg(2:3,:),[],1);
multisensory5 = tadpole.peak_avg(1,:)
respROIs_5V = find(unisensory5 > 0.1 | multisensory5 > 0.1)

figure;
plot(max(tadpole.peak_avg(2:3,respROIs_5V)), tadpole.peak_avg(1,respROIs_5V), 'o')
uni_respROI5 = max(tadpole.peak_avg(2:3,respROIs_5V))
multi_respROI5 = tadpole.peak_avg(1,respROIs_5V)



