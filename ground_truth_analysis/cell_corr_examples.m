%% Generate Examples to show that correlations are not an artifact of my recording/analysis methods

%% Open data file
\\files.brown.edu\research\BM_AizenmanLab\Torrey_calciumimaging\_Data Analysis\compare_46-49\analysis
load('allData_workspace_20171109.mat')

%% Get cells that are close together and their correlations (at 0 offset) in a table
% only use respROIs, multisensory responses (since this is for examples)

% put ROI distance and correlation at 0 into a matrix
%t = 1
% some exps don't have ROIdist, so add that field (using original pixel
% ids)
for t = 1:length(allData)
    if ~isfield(allData{1,t}, 'ROIdist')
        if isfield(allData{1,t}, 'somaticROICenters')
            id = t
            for r1 = 1:length(allData{1,t}.somaticROICenters)
                for r2 = 1:length(allData{1,t}.somaticROICenters)
                    r1_loc = allData{1,t}.somaticROICenters{1,r1}.Centroid;
                    r2_loc = allData{1,t}.somaticROICenters{1,r2}.Centroid;
                    allData{1,t}.ROIdist(r1, r2) = sqrt( (r1_loc(1) + r2_loc(1) ) ^2 + (r1_loc(2) + r2_loc(2) ) ^2 );
                end
            end
        end
    end
end

% and some have the centers in a matrix already
for t = 1:length(allData)
    if ~isfield(allData{1,t}, 'ROIdist')
        if isfield(allData{1,t}, 'ROIcenters')
            id = t
            for r1 = 1:length(allData{1,t}.ROIcenters)
                for r2 = 1:length(allData{1,t}.ROIcenters)
                    r1_loc = allData{1,t}.ROIcenters(r1,:);
                    r2_loc = allData{1,t}.ROIcenters(r2,:);
                    allData{1,t}.ROIdist(r1, r2) = sqrt( (r1_loc(1) + r2_loc(1) ) ^2 + (r1_loc(2) + r2_loc(2) ) ^2 );
                end
            end
        end
    end
end

% now we can use the distances and put them in a nice matrix here
for t = 1:length(allData)
    if isfield(allData{1,t}, 'respROIdff0_R_MS')
    
    % Find the 0 offset R value and make a square matrix
    %for reference, xcorr arranges the cols as 1-1, 1-2, 1-3 ... 1-n, 2-1, 2-2,
    row_0offset = ceil(size(allData{1,t}.respROIdff0_R_MS,1) / 2)
    len = length(allData{1,t}.resp_ROIs);
    for r = 1:len
        for c = 1:len
            idx = (r-1)*len + c;
            allData{1,t}.cell_by_corr(r,c,1) = allData{1,t}.respROIdff0_R_MS(row_0offset,idx);
        end
    end

    % import the correct distance (in pixels)
    for r1 = 1:length(allData{1,t}.resp_ROIs)
        for r2 = 1:length(allData{1,t}.resp_ROIs)
            allData{1,t}.cell_by_corr(r1, r2, 2) = allData{1,t}.ROIdist(r1, r2);
        end
    end
    end
end

% Matrix cell_by_corr now contains the R value at 0 offset from xcorr in
% dim3, 1 and the distance between the cells in dim3, 2. 
% Only responding ROIs are included, so this will match up with the xcorr
% calculations used in ppt. 

%% now let's find cells that are nearby (within 600 pixels) AND low
% correlation (less than 0.3) and put logical into cell_by_corr (r1, r2, 3)

for t = 1:length(allData)
    if isfield(allData{1,t}, 'cell_by_corr')
    for r1 = 1:length(allData{1,t}.resp_ROIs)
        for r2 = 1:length(allData{1,t}.resp_ROIs)
            if (allData{1,t}.cell_by_corr(r1, r2, 1) < 0.3) && (allData{1,t}.cell_by_corr(r1, r2, 2) < 600)
                allData{1,t}.cell_by_corr(r1, r2, 3) = 1
            else
                allData{1,t}.cell_by_corr(r1, r2, 3) = 0
            end
        end
    end
    end
end   

% How many cell pairs do we get?
for t = 1:length(allData)
    if isfield(allData{1,t}, 'cell_by_corr')
    num_cells(t) = sum(sum(allData{1,t}.cell_by_corr(:, :, 3)));
    else
        num_cells(t) = nan
    end
end

% whoops, that's 0 pairs across all exps. Need to reassess criteria. 
% with 0.3 and 600, get some pairs

%% What are my ranges?

for t = 1:length(allData)
    if isfield(allData{1,t}, 'cell_by_corr')
% for distances between cells in pixels
ranges(t,1) = min(min(allData{1,t}.cell_by_corr(:, :, 1)));
ranges(t,2) = max(max(allData{1,t}.cell_by_corr(:, :, 1)));
% for correlations (R val)
ranges(t,3) = min(min(allData{1,t}.cell_by_corr(:, :, 2)));
ranges(t,4) = max(max(allData{1,t}.cell_by_corr(:, :, 2)));
    else
        ranges(t, :) = nan
    end
end

% minimums are 200-800 pixels and negative R vals

%% Select cell pairs with low correlation values and plot them
% use cell_by_corr(:,:,3) to select

for t = 1:length(allData)
    if num_cells(t) > 0
        for r1 = 1:length(allData{1,t}.resp_ROIs)
            for r2 = 1:length(allData{1,t}.resp_ROIs)
                if allData{1,t}.cell_by_corr(r1, r2, 3) 
                    figure;
                    hold on
                    plot(allData{1,t}.dff0_multi(r1,:))
                    plot(allData{1,t}.dff0_multi(r2,:))
                    hold off
                    xlabel('time (trials)')
                    ax=gca;
                    xsize = length(allData{1,t}.df_f0{1,1});
                    ax.XTick = [0, xsize, xsize*2, xsize*3, xsize*4, xsize*5, xsize*6, xsize*7, xsize*8, xsize*9, xsize*10, xsize*11, xsize*12, xsize*13, xsize*14, xsize*15, xsize*16, xsize*17, xsize*18, xsize*19, xsize*20];
                    ax.XTickLabel = {'0','1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20'};
                    title(sprintf('Exp %d respROIs %d vs %d, corr=%1.2f, dist=%1.0f', allData{1,t}.expnum, r1, r2, round(allData{1,t}.cell_by_corr(r1, r2, 1), 2), round(allData{1,t}.cell_by_corr(r1, r2, 2))));
                    ylabel('\Delta F / F_{0}');
                    set(gcf,'units','inches','position',[2, 2, 8, 3])
                    fig_filename = sprintf('Exp%d(t%d)respROIs%dvs%d,dist%1.0f', allData{1,t}.expnum, t, r1, r2, round(allData{1,t}.cell_by_corr(r1, r2, 2)));
                    saveas(gcf, fig_filename, 'png')
                    close;
                end
            end
        end
    end
end

                    

