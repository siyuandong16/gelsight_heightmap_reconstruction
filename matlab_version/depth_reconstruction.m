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

clear; close all;clc;
addpath(genpath('calibration functions'));
addpath(genpath('Funcs'));
% load('new_0717Exp2.mat')
load('gelslim.mat')
% load('new_0217.mat')
% load('new_0717Exp2_short.mat')
%% Create video input object
flagVid = 0;
border=20;

i=0;
showscale=4;
MarkerAreaThresh=30;
ref = imread('C:\Users\siyua\Desktop\CalibrationData60\ref.jpg');
[f0, f00]=iniFrame(ref, border);%% cut the border and do some filter 
f01=sum(f0,3);
frame0=ref(border+1:end-border,border+1:end-border,:,1);
I=f0-double(frame0);
dI=min(I,[],3);
flowcenter=gray2center(dI, MarkerAreaThresh);  % initialize the center of the markers in the initialization frame
center_last=flowcenter; 

% Read images from camera
rawim = imread('C:\Users\siyua\Desktop\CalibrationData60\reconstruct_6.jpg');
imwobd = rawim(border+1:end-border,border+1:end-border,:);
I = double(imwobd)-f0;
dI=min(-I,[],3);
figure;imshow(imwobd,[]);
%     MarkerCenter=gray2center(dI, MarkerAreaThresh); 
%     [ut,vt, center_last, AreaChange]=cal_vol_oncenter2(flowcenter,center_last,MarkerCenter); 

% 3D reconstruction
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

figure1 = figure(1);
subplot 121; h11 = imshow(imwobd);hold on;
subplot 122; h12 = mesh(hm);title('3D Height Map','FontSize',15);view([-53,29])
axis equal;
drawnow;        
%% show gradient 
figure1 = figure(1);
subplot 121; imshow(ImGradX,[]);
subplot 122; imshow(ImGradY,[]);
figure;mesh(hm);axis equal;
%%

