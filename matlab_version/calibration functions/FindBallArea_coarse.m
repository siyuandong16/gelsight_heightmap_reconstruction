function [ContactMask, ValidMap, center, Radius]= FindBallArea_coarse(I,frame, MANUAL)
% find the circle in the pressing sphere sample case. 
% The input matrix I is single-channelled, the result of the current image
% sustract the background
% MANUAL: whether to adjust the center and radius manually

TEST_DESPLAY=0;     %whether display the ball area
if ~exist('MANUAL')
    MANUAL=false;
end
if ~exist('frame')
    frame(:,:,1)=(I+50)*1.5;
    frame(:,:,2)=(I+50)*1.5;
    frame(:,:,3)=(I+50)*1.5;
end
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

%% estimate the circle position
stepwidth=20;
thresh=0.35;
sumcol=sum(I_,1);
sumcol0=sumcol-max(sumcol);
t=find(sumcol0<(min(sumcol0)*thresh));
center1X=(t(1)+t(end))/2;
Radius1=(t(end)-t(1))/2;

sumrow=sum(I_,2);
sumrow0=sumrow-max(sumrow);
t=find(sumrow0<(min(sumrow0)*thresh));
center1Y=(t(1)+t(end))/2;


% display([Radius1 center1X center1Y]);

%% loop to find the exact center and the radius
stepradius=2;  
stepCenter=2;

[center, Radius]=loopCenter(center1X, center1Y, stepCenter*2, I_, mask2, RadiusTestWidth, Radius1, stepradius*2, RadiusRange);
[center, Radius]=loopCenter(center(1), center(2), stepCenter, I_, mask2, RadiusTestWidth, Radius, stepradius, RadiusRange);

Radius=Radius-3
%% display
if TEST_DESPLAY
    [xq,yq]=meshgrid(1:size2, 1:size1);
    xq=xq-center(1);yq=yq-center(2);
    rq=xq.^2+yq.^2;
    displayIm=zeros(size1, size2, 3);
    displayIm(:,:,1)=-I_/150;
    displayIm(:,:,2)=(rq<(Radius^2))/2.5;
    figure;imshow(displayIm);
end

if MANUAL
%     center(1) = 240;
%     center(1) = 320;
%     Radius1 = 30;
    kstep = 2;
    Radius1=Radius;
    hf=figure;
    disIm=uint8(frame);
    [xq,yq]=meshgrid(1:size(frame,2),1:size(frame,1));
    % center=[center1X center1Y];
    BallBord(1)=center(1)-Radius1;BallBord(2)=center(1)+Radius1;
    BallBord(3)=center(2)-Radius1;BallBord(4)=center(2)+Radius1;
    rq=(xq-center(1)).^2+(yq-center(2)).^2;
    r2=Radius1^2;
    a=frame(:,:,1);
    a(rq<r2)=a(rq<r2)-100;
    disIm(:,:,1)=a;
    himage=imshow(disIm);


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
               Radius1=Radius1+kstep;
               r2=Radius1^2;
           elseif c=='-'
               Radius1=Radius1-kstep;
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
        a=uint8(frame(:,:,1));
        a(rq<r2)=a(rq<r2)-130;
        disIm(:,:,1)=a;
        set(himage,'CData',disIm);  
    end
    close(hf);
    Radius=Radius1;
    ContactMask=(rq<r2);
end

%% 
[xq,yq]=meshgrid(1:size2, 1:size1);
xq=xq-center(1);yq=yq-center(2);
rq=xq.^2+yq.^2;
ContactMask=(rq<(Radius^2));
ValidMap=mask2;


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

