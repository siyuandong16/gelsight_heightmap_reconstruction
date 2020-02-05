function center=gray2center(grayim, thresh)
% y,x

grayim(grayim<thresh)=0;
L = logical(grayim);
s = regionprops(L,grayim, 'WeightedCentroid','Area');
areathresh=15;
ind=find([s(:).Area]>areathresh);
center_=[s(ind).WeightedCentroid];
center(:,1)=center_(2:2:end);
center(:,2)=center_(1:2:end);
center(:,3)=[s(ind).Area];