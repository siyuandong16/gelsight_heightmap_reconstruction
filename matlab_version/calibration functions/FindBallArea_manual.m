function [ContactMask, ValidMap, center, Radius]= FindBallArea_manual(I, frame)
% find the circle using manual method 
% The input matrix I is single-channelled, the result of the current image
% sustract the background

IlluThresh=-5; %should be -0.5 for the very soft cases
MarkerAreaThresh=15;
SE=ones(7);
RadiusRange=[40 500]; % the preset range of the circle radius
RadiusTestWidth=50;

size1=size(I,1);size2=size(I,2);
mask1=(I>MarkerAreaThresh); % get the marker area
mask2=~imdilate(mask1,SE);
I_=I.*mask2;

%ini Output
Radius=0;
ValidMap=0;
ContactMask=0; center=0;

%estimate the circle position
stepwidth=20;
sumcol=sum(I_,1);
i0=0;
for i=1:stepwidth:size2-stepwidth
    i0=i0+1;
    t=sum(sum(mask2(:,i:i+stepwidth)));
    sumcol2(i0)=sum(sumcol(i:i+stepwidth))/t;
end
sumbase=(sumcol2(1)+sumcol2(end))/2;
sumcol2=sumcol2-sumbase;
cStart=2;cStart0=1;
cEnd0=length(sumcol2);cEnd=cEnd0-1;
while cStart<cEnd
    cStart=cStart+1;
    if (sumcol2(cStart)<IlluThresh) && (sumcol2(cStart)<sumcol2(cStart-1)) && (sumcol2(cStart)<sumcol2(cStart+1)+0.05)
        if (sumcol2(cStart)<sumcol2(cStart-2)-0.4)               % is the start
            cStart0=0;break;
        end
    elseif cStart>4 && (sumcol2(cStart)<sumcol2(cStart-3)-0.4) &&(sumcol2(cStart)<-0.3)
            cStart0=0;break;
    end
end
while cStart<cEnd
    cEnd=cEnd-1;
    if (sumcol2(cEnd)<IlluThresh) && (sumcol2(cEnd)<sumcol2(cEnd-1)+0.05) && (sumcol2(cEnd)<sumcol2(cEnd+1))
        if (sumcol2(cEnd)<sumcol2(cEnd+2)-0.4)
            cEnd0=0;
            break;
        end  
   elseif cEnd<cEnd0-3 && (sumcol2(cEnd)<sumcol2(cEnd+3)-0.4) && (sumcol2(cEnd)<-0.3)
            cEnd0=0;
            break;
    end
end
if ~(cEnd0 || cStart0)
    center1X=round((cEnd+cStart)/2*stepwidth);
    Radius1=(cEnd-cStart)/2*stepwidth;
else
    display('no circle');
    return;
end


sumrow=sum(I_,2);
i0=0;
for i=1:stepwidth:size1-stepwidth
    i0=i0+1;
    t=sum(sum(mask2(:,i:i+stepwidth)));
    sumrow2(i0)=sum(sumrow(i:i+stepwidth))/t;
end
sumrow2=sumrow2-sumbase;
cStart=1;cStart0=1;
cEnd0=length(sumrow2);cEnd=cEnd0;
if sumrow2(cStart0)>-1
while cStart<cEnd0-1
    cStart=cStart+1;
    if sumrow2(cStart)<-5
        cStart0=0;break;
    elseif (sumrow2(cStart)<sumrow2(cStart+1)+0.04) && (sumrow2(cStart)<IlluThresh) && (sumrow2(cStart)<sumrow2(cStart-1)-0.6) && (sumrow2(cStart-1)>-1)
        cStart0=0;break;
    elseif (sumrow2(cStart)<IlluThresh) && (sumrow2(cStart)<sumrow2(cStart-1)) && (sumrow2(cStart)<sumrow2(cStart+1)+0.04)
        if cStart==2
            continue;
        end
        if (sumrow2(cStart)<sumrow2(cStart-2)-0.4)               % is the start
            cStart0=0;break;
        end
    end
end
end
if sumrow2(cEnd0)>-1
while cStart<=cEnd && cEnd>3
    cEnd=cEnd-1;
    if cEnd<-5
        cEnd0=0; break;
    elseif (sumrow2(cEnd)<IlluThresh) && (sumrow2(cEnd)<sumrow2(cEnd-1)+0.04) && (sumrow2(cEnd)<sumrow2(cEnd+1)-0.6) 
        cEnd0=0; break;
    elseif (sumrow2(cEnd)<IlluThresh) && (sumrow2(cEnd)<sumrow2(cEnd-1)+0.04) && (sumrow2(cEnd)<sumrow2(cEnd+1))
        if(cEnd0-cEnd==1)
            continue;
        end
        if (sumrow2(cEnd)<sumrow2(cEnd+2)-0.4) 
            cEnd0=0;
            break;
        end
    end
end
end

if ~(cStart0 || cEnd0) && cStart~=cEnd
    center1Y=round((cEnd+cStart)/2*stepwidth);
elseif cStart==cEnd
    if sumrow2(cEnd+3)>sumrow2(cEnd-3)
        center1Y=cEnd*stepwidth-Radius1;
    else
        center1Y=cEnd*stepwidth+Radius1;
    end
elseif ~cStart0
    center1Y=cStart*stepwidth+Radius1;
elseif ~cEnd0
    center1Y=cEnd*stepwidth-Radius1;
else
    center1Y=size(I,1)/2;
end

display([Radius1 center1X center1Y]);

%%
kstep=3;

hf=figure;
disIm=frame;
[xq yq]=meshgrid(1:size(frame,2),1:size(frame,1));
center=[center1X center1Y];
BallBord(1)=center(1)-Radius1;BallBord(2)=center(1)+Radius1;
BallBord(3)=center(2)-Radius1;BallBord(4)=center(2)+Radius1;
rq=(xq-center(1)).^2+(yq-center(2)).^2;
r2=Radius1^2;
a=frame(:,:,1);
a(rq<r2)=a(rq<r2)-100;
disIm(:,:,1)=a;
himage=imshow(disIm,[]);


while ishandle(hf)
    k=waitforbuttonpress;
    if k   %keyboard
       c=get(hf,'CurrentCharacter');
%        display(c);
       if c==27     %ESC
           break;
       elseif c==30 || c=='w'
           center(2)=center(2)-kstep;
       elseif c==31 || c=='s'
           center(2)=center(2)+kstep;
       elseif c==28 || c=='a'
           center(1)=center(1)-kstep;
       elseif c==29 || c=='d'
           center(1)=center(1)+kstep;
       elseif c=='+'
           Radius1=Radius1+kstep*3;
           r2=Radius1^2;
       elseif c=='-'
           Radius1=Radius1-kstep*3;
           r2=Radius1^2;
       end
       BallBord(1)=center(1)-Radius1;BallBord(2)=center(1)+Radius1;
       BallBord(3)=center(2)-Radius1;BallBord(4)=center(2)+Radius1;
    else
        p=get(gca,'CurrentPoint');
        xt=p(1,1);yt=p(1,2);
        t=abs([BallBord(1:2)-xt BallBord(3:4)-yt]);
        [t1 t2]=min(t);
        if t2<3
            BallBord(t2)=xt;
            center(1)=(BallBord(1)+BallBord(2))*0.5;
            Radius1=(BallBord(2)-BallBord(1))*0.5;
            r2=Radius1^2;
            BallBord(3)=center(2)-Radius1;
            BallBord(4)=center(2)+Radius1;
        else
            BallBord(t2)=yt;
            center(2)=(BallBord(4)+BallBord(3))*0.5;
            Radius1=(BallBord(4)-BallBord(3))*0.5;
            r2=Radius1^2;
            BallBord(1)=center(1)-Radius1;
            BallBord(2)=center(1)+Radius1;
        end
    end
    % redraw the circle
    rq=(xq-center(1)).^2+(yq-center(2)).^2;
    a=frame(:,:,1);
    a(rq<r2)=a(rq<r2)-130;
    disIm(:,:,1)=a;
    set(himage,'CData',disIm);  
end
close(hf);
Radius=Radius1;
ContactMask=(rq<r2);
ValidMap=mask2;
