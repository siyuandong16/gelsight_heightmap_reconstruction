% main_generate Calibration file

inputpath='./data/';  % the folder of the input video files
outputpath='./CaliIm_0702/';

if ~exist(inputpath)
    error('input folder not existed')
end
if ~exist(outputpath)
    mkdir(outputpath);
end

border=30;

addpath('Func');
videoIndex=dir(inputpath);
videoLen=length(videoIndex);
% Anno=[];
Frn=0;
for i=1:5
% for i=5:5
    videoname=videoIndex(i).name;
    if videoname(1)=='.'
        continue;
    end
    Frn=Frn+1;
    
    %% Read each video and get the frames
    v = VideoReader([inputpath videoname]);
    frame0=readFrame(v);%first frame
    MaxDiff=1e6;
    frame1=frame0;
%     ii=0;
    while hasFrame(v)
%         ii=ii+1;
        frame=readFrame(v);
        d=sum(abs(frame(:)-frame0(:)));
        if d>MaxDiff
            MaxDiff=d;
            frame1=frame;
%             display(ii);
        end                
    end    
    %find the center and radius
    f0= iniFrame(frame0, border);
    frame_=frame1(border+1:end-border,border+1:end-border,:);
    I=double(frame_)-f0;
%     [ContactMask, ValidMap, center, Radius]= FindBallArea_manual(I, frame_);
    [ContactMask, validMask, touchCenter, Radius]= FindBallArea_coarse(dI,frame_,1);
    Anno(Frn).center=center+border;
    Anno(Frn).Radius=Radius;
    Anno(Frn).FramePath=[outputpath videoname(1:end-4) '_1.jpg'];
    Anno(Frn).F0Path=[outputpath videoname(1:end-4) '_0.jpg'];
    save([outputpath 'anno.mat'], 'Anno');
    imwrite(frame0, [outputpath videoname(1:end-4) '_0.jpg']);
    imwrite(frame1, [outputpath videoname(1:end-4) '_1.jpg']);  
%     display(ii);
    
end
save([outputpath 'anno.mat'], 'Anno');