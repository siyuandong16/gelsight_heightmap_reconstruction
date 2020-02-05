function [DiffScore Strongness]=evalueCircle(center, radius, width, I, mask)
% I is the difference image, with marker area marked as 0; mask indicates
% the available pixels

size1=size(I,1);size2=size(I,2);
[xq yq]=meshgrid(1:size2, 1:size1);
xq=xq-center(1);yq=yq-center(2);
rq=xq.^2+yq.^2;

r2=radius^2;r3=(radius+width)^2;
% r0=(radius-width*0.1)^2;
% r1=(radius-width*2)^2;
r1=r2*0.57;
r0=r2;
% r2=(radius+width*0.5)^2;r3=(radius+width*1.5)^2;


m1=((rq<r0) & (rq>r1));
m2=((rq>r2) & (rq<r3));

% t1=I(logical((m2.*mask)));t1=t1(:);

m10=I(m1 & mask);m20=I(m2 & mask);
% Strongness=sum(sum(m1.*I))/sum(sum(m1.*mask));
% DiffScore=sum(sum(m2.*I))/sum(sum(m2.*mask))-Strongness*1.5;
Strongness=mean(m10);
DiffScore=mean(m20)-Strongness;
% Strongness=sum(m10);
% DiffScore=mean(m20)*radius*width*3-Strongness;

% DiffScore=DiffScore+min(m20);
% DiffScore=DiffScore-var(m20)/radius*100;

% DiffScore=DiffScore-var(t1)/5;
% DiffScore=DiffScore-radius/150;