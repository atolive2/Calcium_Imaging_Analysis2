%% maxR vs distance
figure;
hold on
for t = 1:length(allData)
    if length(allData{1,t}.resp_ROIs) > 10
        %figure;
        
         if allData{1,t}.stage == 46
           plot(allData{1,t}.respROIdff0_maxR_sq_all, allData{1,t}.topo_dist, 'go')
        elseif allData{1,t}.stage == 49
            plot(allData{1,t}.respROIdff0_maxR_sq_all, allData{1,t}.topo_dist, 'mo')
        end
    end
end

figure;
hold on
for t = 1:length(allData)
    if length(allData{1,t}.resp_ROIs) > 10
        k = find(triu(allData{1,t}.respROIdff0_maxR_sq_all, 1));
        tmpY = allData{1,t}.topo_dist(k)
        tmpX = allData{1,t}.respROIdff0_maxR_sq_all(k)
        %plot(tmpX, tmpY, 'o')
        [p, S] = polyfit(tmpX, tmpY, 1)
        x1 = 0:0.1:1;
        y1 = polyval(p, x1)
        %hold on
        if allData{1,t}.stage == 46
            plot(x1, y1, 'g')
        elseif allData{1,t}.stage == 49
            plot(x1, y1, 'm')
        end
        R = corrcoef(tmpX, tmpY);
        dist_maxR_corr(t) = R(1,2);
        clear('row', 'col', 'tmpX', 'tmpY', 'x1', 'y1', 'R')
    end
end

figure;
hold on
for t = 1:length(allData)
    if length(allData{1,t}.resp_ROIs) > 10
        k = find(triu(allData{1,t}.respROIdff0_maxR_sq_MS, 1));
        tmpY = allData{1,t}.topo_dist(k)
        tmpX = allData{1,t}.respROIdff0_maxR_sq_MS(k)
        %plot(tmpX, tmpY, 'o')
        [p, S] = polyfit(tmpX, tmpY, 1)
        x1 = 0:0.1:1;
        y1 = polyval(p, x1)
        %hold on
        if allData{1,t}.stage == 46
            plot(x1, y1, 'g')
        elseif allData{1,t}.stage == 49
            plot(x1, y1, 'm')
        end
        R = corrcoef(tmpX, tmpY);
        dist_maxR_corrMS(t) = R(1,2);
        clear('row', 'col', 'tmpX', 'tmpY', 'x1', 'y1', 'R')
    end
end

figure;
hold on
for t = 1:length(allData)
    if length(allData{1,t}.resp_ROIs) > 10
        k = find(triu(allData{1,t}.respROIdff0_maxR_sq_V, 1));
        tmpY = allData{1,t}.topo_dist(k)
        tmpX = allData{1,t}.respROIdff0_maxR_sq_V(k)
        %plot(tmpX, tmpY, 'o')
        [p, S] = polyfit(tmpX, tmpY, 1)
        x1 = 0:0.1:1;
        y1 = polyval(p, x1)
        %hold on
        if allData{1,t}.stage == 46
            plot(x1, y1, 'g')
        elseif allData{1,t}.stage == 49
            plot(x1, y1, 'm')
        end
        R = corrcoef(tmpX, tmpY);
        dist_maxR_corrV(t) = R(1,2);
        clear('row', 'col', 'tmpX', 'tmpY', 'x1', 'y1', 'R')
    end
end


figure;
hold on
for t = 1:length(allData)
    if length(allData{1,t}.resp_ROIs) > 10
        k = find(triu(allData{1,t}.respROIdff0_maxR_sq_M, 1));
        tmpY = allData{1,t}.topo_dist(k)
        tmpX = allData{1,t}.respROIdff0_maxR_sq_M(k)
        %plot(tmpX, tmpY, 'o')
        [p, S] = polyfit(tmpX, tmpY, 1)
        x1 = 0:0.1:1;
        y1 = polyval(p, x1)
        %hold on
        if allData{1,t}.stage == 46
            plot(x1, y1, 'g')
        elseif allData{1,t}.stage == 49
            plot(x1, y1, 'm')
        end
        R = corrcoef(tmpX, tmpY);
        dist_maxR_corrM(t) = R(1,2);
        clear('row', 'col', 'tmpX', 'tmpY', 'x1', 'y1', 'R')
    end
end

% consolidate R vals into a table
for t = 1:length(allData)
    all_corrcoeff(t,1) = allData{1,t}.stage
end
all_corrcoeff(:,2) = dist_maxR_corrMS;
all_corrcoeff(:,3) = dist_maxR_corrV;
all_corrcoeff(:,4) = dist_maxR_corrM;

s49_corrcoeff = all_corrcoeff(st49, :)
s46_corrcoeff = all_corrcoeff(st46, :)
avg_corrcoeff(:,1) = mean(s46_corrcoeff(:,2:4))
avg_corrcoeff(:,2) = mean(s49_corrcoeff(:,2:4))
