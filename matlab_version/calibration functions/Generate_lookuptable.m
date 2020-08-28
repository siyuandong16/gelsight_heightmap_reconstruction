
% Generate lookup table from the video
% To make the calibration video: press a sphere onto the GelSight sensor
% surface, roll it on different positions. Maker sure the 1st frame (or at
% least one frame) is blank with nothing contacting the Gel
% Pick up one or more frames in the videos to calibrate the sensor
% Modify the Ball's radius and the sensor's display ratio for each test
%
% Mod by Wenzhen Yuan (yuanwenzhen@gmail.com), Feb 2017
clear;close all;
BallRad=4.76/2;% Ball's radius, in mm   
border=20;
BALL_MANUAL=1;      % whether to find the ball manually

type='new';  % choose whether it's the new sensor or the old one
name1=[type '_0303'];
%%
filename = [3,6,7,8,9,10,12,13,14,15,16,17,18,19,26,33,34,35,36,37,38,39,40,41,42,43,45,46,47,48];
ref = imread('C:\Users\siyua\Desktop\CalibrationData60\ref.jpg');
f0 = double(ref(border+1:end-border,border+1:end-border,:));
for i = 1:length(filename) 
    img(:,:,:,i) = imread(['C:\Users\siyua\Desktop\CalibrationData60\sample_' num2str(filename(i)) '.jpg']);
end 
imshow(int(f0),[]);
%% generate lookup table 
READ_RADIUS=0;

bins=80;
gradmag=[];gradir=[];countmap=[];
if strcmp(type,'new')
    %load('evaluation_data_new_sensor_0303.mat');
    Pixmm=0.106;
    zeropoint=-90;
    lookscale=180;
elseif strcmp(type,'old')
    %load('evaluation_data_old_sensor_v23.mat');
    zeropoint=-1e-3;
    lookscale=5e-3;
    Pixmm=0.0254;
end
BallRad_pix=BallRad/Pixmm;
             
for Frn=1:size(img,4)
    frame = img(:,:,:,Frn);
    display(['Calibration on Frame' num2str(Frn)]);
    frame_=frame(border+1:end-border,border+1:end-border,:);
    I=double(frame_)-double(f0);
    
    dI=(min(I,[],3)-max(I,[],3))/2;
    if ~READ_RADIUS
        [ContactMask, validMask, touchCenter, Radius]= FindBallArea_coarse(dI,frame_,BALL_MANUAL);
        validMask=validMask & ContactMask;
    else
        touchCenter=touchCenter_record(:,:,Frn);
        Radius=radius_record(Frn);
        if Radius>BallRad_pix
            Radius=BallRad_pix;
        end
        Radius2=Radius^2;
        size1=size(I,1); size2=size(I,2);
        validMask=zeros(size1, size2);
        [xq,yq]=meshgrid(1:size2, 1:size1);xq=xq-touchCenter(1);yq=yq-touchCenter(2);
        rq=xq.^2+yq.^2;
        validMask(rq<Radius2)=1;
    end
    if strcmp(type, 'old')
        nomarkermask=max(-I,[],3)<50;  % for old sensor
        nomarkermask=imerode(nomarkermask,strel('disk',3)); 
        validMask=validMask & nomarkermask;
        [gradmag, gradir,countmap]=LookuptableFromBall_RL(I,f0, bins , touchCenter, BallRad, Pixmm, validMask, gradmag, gradir, countmap, zeropoint, lookscale);
    elseif strcmp(type, 'new')
        nomarkermask=min(-I,[],3)<30;   % for new sensor
        nomarkermask=imerode(nomarkermask,strel('disk',3)); 
        validMask=validMask & nomarkermask;      
        [gradmag, gradir,countmap]=LookuptableFromBall_Bnz(I,f0, bins , touchCenter, BallRad, Pixmm, validMask, gradmag, gradir, countmap, zeropoint, lookscale);
    end   
end
[GradMag, GradDir]=LookuptableSmooth(bins, gradmag, gradir, countmap);
LookupTable.bins=bins;
LookupTable.GradMag=GradMag;
LookupTable.GradDir=GradDir;
LookupTable.Zeropoint=zeropoint;
LookupTable.Scale=lookscale;
LookupTable.Pixmm=Pixmm;
LookupTable.FrameSize=size(frame);
savename=[name1 '.mat'];
save(savename,'LookupTable');    %Save the loopup table file
display('Calibration Done');
