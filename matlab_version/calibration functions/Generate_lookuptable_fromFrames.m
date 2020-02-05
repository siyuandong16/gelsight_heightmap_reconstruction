
% Generate lookup table from images in a single folder
% To make the calibration video: press a sphere onto the GelSight sensor
% surface, roll it on different positions. Maker sure the 1st frame (or at
% least one frame) is blank with nothing contacting the Gel
% Pick up one or more frames in the videos to calibrate the sensor
% Modify the Ball's radius and the sensor's display ratio for each test
%
% Mod by Wenzhen Yuan (yuanwenzhen@gmail.com), July 2017 
% clear;close all;
BallRad=6.35/2;% Ball's radius, in mm   
BALL_MANUAL=1;      % whether to find the ball manually
type='new';  % choose whether it's the new sensor or the old one
name1=[type '_gelsight_5'];
%%
v = VideoReader('calibration_video.mp4');
t = 1;
while hasFrame(v)
    video = readFrame(v);
    figure(1);imshow(video);
    w = waitforbuttonpress;
    if w == 0
        disp(['frame' num2str(t)]);
        img(:,:,:,t) = video;
        t = t+1; 
    end
end
%%   
num = 1;
for i = [230, 240, 320, 340, 410, 420,  540, 620, 630, 640, 750, 840, 950, 1010, 1110, 1130, 1190, 1290,  1300, 1380, 1460, 1540, 1630]
    img(:,:,:,num) = imread(['C:/Users/Siyuan/Desktop/Heightmap_reconstruction_with_GelSlim-master/Heightmap_reconstruction_with_GelSlim-master/new_data/new_data/sample2_' num2str(i) '.jpg']); 
    num = num + 1;
end 
%%
frame0 = imread(['C:/Users/Siyuan/Desktop/Heightmap_reconstruction_with_GelSlim-master/Heightmap_reconstruction_with_GelSlim-master/new_data/new_data/ref.jpg']); 
% frame0_ref = imread(['C:/Users/Siyuan/Desktop/gelsight_image/raw_image_' num2str(1) '.jpg']);
% frame0_ref2 = imread(['C:/Users/Siyuan/Desktop/gelsight_image/raw_image_' num2str(1) '.jpg']);
pad = 20;
% frame0(133-pad:133+pad, 298-pad:298+pad,:) = frame0_ref(133-pad:133+pad, 298-pad:298+pad,:);
% frame0(196-pad:196+pad, 160-pad:160+pad,:) = frame0_ref(196-pad:196+pad, 160-pad:160+pad,:);
figure;imshow(frame0,[]);
%%
v = VideoReader('calibration_video.mp4');
t = 1;
while hasFrame(v)
    video = readFrame(v);
    figure(1);imshow(video);
    w = waitforbuttonpress;
    if w == 0
        disp(['frame' num2str(t)]);
        frame0 = video;
        t = t+1; 
    end
end
%% generate lookup table 
% READ_RADIUS=1;
border=20;
bins=80;
gradmag=[];gradir=[];countmap=[];
gradx=[];grady=[];
if strcmp(type,'new')
    Pixmm=0.0283; %for 960x720
    Pixmm = 0.0671; 
%     Pixmm=0.0324;  % for 640x480

%     zeropoint=-90;
%     lookscale=180;
    zeropoint=-90;
    lookscale=180;
elseif strcmp(type,'old')
    zeropoint=-1e-3;
    lookscale=5e-3;
    Pixmm=0.0254;
end
BallRad_pix=BallRad/Pixmm;


for i=1:size(img,4)
    f0 = iniFrame(double(frame0), border);
    imshow(f0);
    frame_=img(border+1:end-border,border+1:end-border,:,i);
    I=double(frame_)-f0;
    dI=min(-I,[],3);
    [ContactMask, validMask, touchCenter, Radius]= FindBallArea_coarse(dI,frame_,BALL_MANUAL);
    validMask=validMask & ContactMask;  
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
% LookupTable.GradX=-cos(GradDir).*GradMag;
% LookupTable.GradY=sin(GradDir).*GradMag;
LookupTable.Zeropoint=zeropoint;
LookupTable.Scale=lookscale;
LookupTable.Pixmm=Pixmm;
LookupTable.FrameSize=size(img(:,:,:,1));
savename=[name1 '.mat'];
save(savename,'LookupTable');    %Save the loopup table file
display('Calibration Done');

%% Test
clear LookupTable;
lookupfile='new_gelsight_5.mat';
% inputpath='./CaliIm_0702/';
i=5;
border=0;

load(lookupfile);
% load([inputpath 'anno.mat']);
% frame=imread(Anno(i).FramePath);
% frame0=imread(Anno(i).F0Path);
num = 201;
pad = 20;
% frame = imread(['C:/Users/Siyuan/Desktop/gelsight_image/trial_3/raw_image_' num2str(num) '.jpg']); 
% frame = img(20:320-20,20:427-20,:,1);
% frame0 = frame0(20:320-20,20:427-20,:);
frame = img(:,:,:,1);
f0 = iniFrame(frame0, border);
disp(size(f0));
disp(size(frame0));
border=20;
frame_=frame(border:end-border,border:end-border,:);
disp(size(frame_));
I=double(frame_)-f0;

[ImGradX, ImGradY, ImGradMag, ImGradDir]=matchGrad_Bnz(LookupTable, I, f0);

hm=fast_poisson2(ImGradX, ImGradY);
disp(size(hm))
% hm(hm<0) = 0; 
figure,subplot(1,2,1);imshow(frame);
subplot(1,2,2);mesh(hm*1.5);axis equal
figure;imshow(ImGradX,[]);
figure;imshow(ImGradY,[]);