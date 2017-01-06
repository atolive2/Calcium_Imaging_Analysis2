% get euclidian distance (in pixels) for Torrey's Tadpoles.

pixelToDistanceScale=2.26; % this is the number of microns per pixel. I assume 1 for this code. N microns per pixel based on a calibration.

for n=1:numel(tadpole)
    pairNum=1;
    tadpole_totalPairs(:,n)=nchoosek(numel(tadpole{n}.somaticROICenters),2);
    for k=1:numel(tadpole{n}.somaticROICenters)
        for j=k+1:numel(tadpole{n}.somaticROICenters)

            roiDXs_byTadpole{n}(:,pairNum)=fix(abs(tadpole{n}.somaticROICenters{j}(1).Centroid(2)-tadpole{n}.somaticROICenters{k}(1).Centroid(2)));
            % faster but not perfect
            roiDYs_byTadpole{n}(:,pairNum)=fix(abs(tadpole{n}.somaticROICenters{j}(1).Centroid(1)-tadpole{n}.somaticROICenters{k}(1).Centroid(1)));
            % we round because only integers should be allowed because we
            % have not transformed from pixels yet; the interpolation done
            % by the roi centroid alogorithm can lead to decimals.
            pairMap{n}(:,:,pairNum)=[k,j];
            pairNum=pairNum+1;
        end
    end
    % assign distance as dy/dx, but handle cases where they are either on
    % the same x or y pixel.
    roiEucDist_byTadpole{n}=zeros(size(roiDXs_byTadpole{n}));
    dNzInds=find(roiDYs_byTadpole{n}~=0 & roiDXs_byTadpole{n}~=0);
    xZInds=find(roiDXs_byTadpole{n}==0);
    yZInds=find(roiDYs_byTadpole{n}==0);
    roiEucDist_byTadpole{n}(dNzInds)=roiDYs_byTadpole{n}(dNzInds)./roiDXs_byTadpole{n}(dNzInds);
    roiEucDist_byTadpole{n}(yZInds)=roiDXs_byTadpole{n}(yZInds);
    roiEucDist_byTadpole{n}(xZInds)=roiDYs_byTadpole{n}(xZInds);
    
end

figure,nhist(roiEucDist_byTadpole{7},'box','minbins',100,'proportion')