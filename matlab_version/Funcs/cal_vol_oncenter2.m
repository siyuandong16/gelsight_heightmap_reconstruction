function [ut, vt, center_back, AreaChange] = cal_vol_oncenter2(center0, center1, center2)

% center0 is the original position of the center, center 1 is the last
% frame position, in the order of the corresponding for center0; center2 is
% the center of the current frame

n=size(center0,1);
n1=size(center2,1);

no_seq2=[];
for i=1:n1
    dif=abs(center2(i,1)-center1(:,1))+abs(center2(i,2)-center1(:,2));
%     [a b]=min(dif);
    [a b]=min(dif.*(abs(center2(i,3)-center0(:,3))+100));
    no_seq2(i)=b;
end

for i=1:n
    dif=abs(center1(i,1)-center2(:,1))+abs(center1(i,2)-center2(:,2));
    [a b]=min(dif.*(abs(center0(i,3)-center2(:,3))+100));
    a=a/100;
    if center0(i,3)<a    %for small area
       ut(i)=0;vt(i)=0;
        center_back(i,:)=center0(i,:);    
        AreaChange(i)=100;
    elseif i==no_seq2(b)
        vt(i)=center2(b,1)-center0(i,1);
        ut(i)=center2(b,2)-center0(i,2);
        center_back(i,:)=center2(b,:);
        AreaChange(i)=center2(b,3)-center0(i,3);
    else
        ut(i)=0;vt(i)=0;
        center_back(i,:)=center0(i,:);
        AreaChange(i)=100;
    end
end

end