 %% Important Instructions:
%       hm -  3D height map created from GelSight, where brightness
%       represents the height. Please adjust the display range of the
%       imshow function in line 130 for proper display contrast (Default is
%       [0 30]), i.e.,
%
%           subplot 122; h12 = imshow(hm,[0 30]);title('3D Height Map','FontSize',15);
%
%       Check readme.txt for details about how to use the software.
%       Note that the current version of the software is experimental only and NOT for distribution.
%
%   Author: Rui Li (rui@mit.edu). 11/02/2014.
%   Copyright 2014 MIT.

if exist('flagVid','var')
    if flagVid == 1
        stop(vid);
        delete(vid);
    end
end

clear; close all;clc;

%% addpath and load lookuptable
addpath(genpath('calibration functions'));
addpath(genpath('Funcs'));
% load('new_0717Exp2.mat')
load('new_shanluo.mat')
% load('new_0217.mat')
% load('new_0717Exp2_short.mat')
%% Create video input object
camera = 'logitech';
% camera = 'webcam';

if strcmp(camera,'logitech')
    % FOR LOGITECH CAMERA
    
    vid = videoinput('winvideo',1,'RGB24_640x480'); 
    set(vid,'TriggerRepeat',Inf);
    vid.FrameGrabInterval = 1;
    src = getselectedsource(vid);
    set(src,'Tag','GelSight realtime 3D');
    src.ExposureMode = 'manual';
    src.Exposure = -0;
    atobj=getselectedsource(vid);
    
    vid.ReturnedColorspace = 'rgb';
    src.BacklightCompensation = 'off';
    imSize=vid.ROIPosition;
    col = imSize(3);row = imSize(4);
else
    disp('Please use correct camera type! logitech or webcam.');
    return;
end

flagVid = 0;
border=10;
% Set trigger type to manual, to prevent automatic triggering (and automatic stopping, unless 'triggerrepeat' is set to 'inf') of the video object
triggerconfig(vid, 'manual'); % Default is 'immediate', which triggers directly after start()...

if ~isrunning(vid)
    start(vid);
    flagVid = 1;
end

% Wait for peekdata to return a frame
while isempty(peekdata(vid,1))
end

i=0;
showscale=4;
MarkerAreaThresh=30;
rawim = peekdata(vid,1);
imdown = impyramid(rawim,'reduce');  % Downsample image
[f0, f00]=iniFrame(imdown, border);%% cut the border and do some filter 
f01=sum(f0,3);
frame0=imdown(border+1:end-border,border+1:end-border,:,1);
I=f0-double(frame0);
dI=min(I,[],3);
flowcenter=gray2center(dI, MarkerAreaThresh);  % initialize the center of the markers in the initialization frame
center_last=flowcenter; 


while true
    loop = tic;
    i=i+1
    
    %% Read images from camera
    rawim = peekdata(vid,1);
    imdown = impyramid(rawim,'reduce');  % Downsample image
    imwobd= imdown(border+1:end-border,border+1:end-border,:);
    I = double(imwobd)-f0;
    dI=min(-I,[],3);
%     MarkerCenter=gray2center(dI, MarkerAreaThresh); 
%     [ut,vt, center_last, AreaChange]=cal_vol_oncenter2(flowcenter,center_last,MarkerCenter); 
    
    %% 3D reconstruction
    [ImGradX, ImGradY, ImGradMag, ImGradDir]=matchGrad_Bnz(LookupTable, I, f0);
    hm=fast_poisson2(ImGradX, ImGradY);
    
    % Smooth the image and enhance the edges. Optional.
    %         hm1 = hm - imfilter(hm, h);
    %         hm2 = hm + hm1;
    %         hm = (hm2/2+hmPrev/2);
    
    maxHm = max(hm(:));
    hm(hm<0) = 0;
    if maxHm < 2
        hm = zeros(size(hm));
    end
    
    %% Display results
    if i < 2
        figure1 = figure(1);
        subplot 121; h11 = imshow(imwobd);hold on;
%         h13 = quiver(flowcenter(:,2),flowcenter(:,1), ut'*showscale, vt'*showscale,'LineWidth',2,'Color','y','AutoScale','off');
        subplot 122; h12 = mesh(hm);title('3D Height Map','FontSize',15);view([-53,29])
        axis equal;
        drawnow;        
    else
        set(h11,'CData',imwobd);
%         set(h13,'UData',ut'*showscale);
%         set(h13,'VData',vt'*showscale);
%         set(h13,'XData',flowcenter(:,2));
%         set(h13,'YData',flowcenter(:,1));
        set(h12,'ZData',hm);set(gca,'zlim',[-5 30]);caxis([-5 30])
        drawnow;
    end
%     hmrecord(:,:,i) = hm; 
    display(['Loop time is: ', num2str(toc(loop)), 's']);
end

