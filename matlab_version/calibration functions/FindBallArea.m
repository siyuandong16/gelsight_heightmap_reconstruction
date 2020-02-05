function [ContactMask, ValidMap, center, Radius]= FindBallArea(I)
% find the circle in the pressing sphere sample case. 
% The input matrix I is single-channelled, the result of the current image
% sustract the background

TEST_DESPLAY=1;     %whether display the ball area


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

%% loop to find the exact center and the radius
stepradius=2;  
stepCenter=2;

[center, Radius]=loopCenter(center1X, center1Y, stepCenter*2, I_, mask2, RadiusTestWidth, Radius1, stepradius*2, RadiusRange);
[center, Radius]=loopCenter(center(1), center(2), stepCenter, I_, mask2, RadiusTestWidth, Radius, stepradius, RadiusRange);

Radius=Radius-3;

[xq yq]=meshgrid(1:size2, 1:size1);
xq=xq-center(1);yq=yq-center(2);
rq=xq.^2+yq.^2;
ContactMask=(rq<(Radius^2));
ValidMap=mask2;
%% display
if TEST_DESPLAY
 
displayIm=zeros(size1, size2, 3);
displayIm(:,:,1)=-I_/50;
displayIm(:,:,2)=(rq<(Radius^2))/2.5;
figure,imshow(displayIm);

end

function [ResCenter, ResRadius]=loopCenter(centerx, centery, cstep, I, mask, rwidth, r0, rstep, rrange)
LoopFlag=1;
centerx0=centerx;centery0=centery;
[R0, CurrentScore]=loopRadius(centerx, centery, I, mask, rwidth, r0, rstep, rrange);        

while LoopFlag
    LoopFlag=0;
    [R1, score1]=loopRadius(centerx0+cstep, centery0, I, mask, rwidth, r0, rstep, rrange);        
    [R2, score2]=loopRadius(centerx0-cstep, centery0, I, mask, rwidth, r0, rstep, rrange);        
    if score1>CurrentScore && score1>score2
        LoopFlag=1;
        while score1>CurrentScore
            CurrentScore=score1;
            centerx0=centerx0+cstep;R0=R1;
            [R1, score1]=loopRadius(centerx0+cstep, centery0, I, mask, rwidth, r0, rstep, rrange);
        end
    elseif score2>CurrentScore && score1<score2
        LoopFlag=1;
        while score2>CurrentScore
            CurrentScore=score2;
            centerx0=centerx0-cstep;R0=R2;
            [R2, score2]=loopRadius(centerx0-cstep, centery0, I, mask, rwidth, r0, rstep, rrange);
        end
    end
    [R1, score1]=loopRadius(centerx0, centery0+cstep, I, mask, rwidth, r0, rstep, rrange);        
    [R2, score2]=loopRadius(centerx0, centery0-cstep, I, mask, rwidth, r0, rstep, rrange);        
    if score1>CurrentScore && score1>score2
        LoopFlag=1;
        while score1>CurrentScore
            CurrentScore=score1;
            centery0=centery0+cstep;R0=R1;
            [R1, score1]=loopRadius(centerx0, centery0+cstep, I, mask, rwidth, r0, rstep, rrange);
        end
    elseif score2>CurrentScore && score1<score2
        LoopFlag=1;
        while score2>CurrentScore
            CurrentScore=score2;
            centery0=centery0-cstep;R0=R2;
            [R2, score2]=loopRadius(centerx0, centery0-cstep, I, mask, rwidth, r0, rstep, rrange);
        end
    end    
end
ResCenter=[centerx0 centery0];ResRadius=R0;


function [R, maxScore]=loopRadius(centerx, centery, I, mask, width, r0, rstep, rrange)
r=r0;
[currentScore t]=evalueCircle([centerx centery], r0, width, I, mask);
[s1 t]=evalueCircle([centerx centery], r-rstep, width, I, mask);
[s2 t]=evalueCircle([centerx centery], r+rstep, width, I, mask);
isFound=0;
if s1>currentScore && s1>s2 %decresae
    r=r-rstep;
    while r>=rrange(1)
        [s2 t]=evalueCircle([centerx centery], r-rstep, width, I, mask);
        if s2<s1 
           isFound=1;
           R=r;maxScore=s1;
           break;
        else
            s1=s2;r=r-rstep;
        end
    end
elseif s1<s2 && s2>currentScore %increase
    r=r+rstep;
    while r<=rrange(2)
        [s1 t]=evalueCircle([centerx centery], r+rstep, width, I, mask);
        if s2>s1 
           isFound=1;
           R=r;maxScore=s2;
           break;
        else
            s2=s1;r=r+rstep;
        end
    end    
else  %equal, or both large 
    isFound=1;R=r0;maxScore=currentScore;
%     R=R*1.02;
end
    
if ~isFound
    R=0;maxScore=0;
end

