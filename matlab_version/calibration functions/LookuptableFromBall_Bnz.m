function [gradmag, graddir, countmap]=LookuptableFromBall_Bnz(dI,f0, bins, center, BallRad, Pixmm, validmask,gradmag, graddir, countmap, zeropoint, lookscale, f01)
% Calibration function for GelSight-Bnz sensor
% dI is the differencew image, 3 channel;f0 is the initialized image
% dI shall be positive mostly
% bins could be 60;
% need to be fixed after calculated at last: call LookuptableSmooth

% Wenzhen Yuan (yuanwenzhen@gmail.com)  Feb 2017


BallRad=BallRad/Pixmm;  % in pix

%% parameter ini
if ~exist('zeropoint')
    zeropoint=-90;
end
if ~exist('lookscale')
    lookscale=180;
end
if ~exist('f01')
    t=mean(f0(:));
    f01=1+((t./f0)-1)*2;
end

%% geometry
sizey=size(dI,1);sizex=size(dI,2);
[xq, yq]=meshgrid(1:sizex, 1:sizey);
xq=(xq-center(1));    % in mm
yq=(yq-center(2));    % in mm
disp(xq)
disp(yq) 
validid=find(validmask);
xvalid=xq(validid);yvalid=yq(validid);
rvalid=sqrt(xvalid.^2+yvalid.^2);
if max(rvalid-BallRad)>0
    display('Contact Radius is too large. Ignoring the exceeding area');
    rvalid(rvalid>BallRad)=BallRad-0.001;
end

gradxseq=asin(rvalid/BallRad);gradyseq=atan2(-yvalid, -xvalid);

%% colorvalid
binm=bins-1;
sizet=sizex*sizey;sizet2=sizet*2;
r1=dI(validid).*f01(validid);g1=dI(validid+sizet).*f01(validid+sizet);b1=dI(validid+sizet2).*f01(validid+sizet2);
disp("mean");
disp(min(r1));disp(max(r1));
disp(min(g1));disp(max(g1));
disp(min(b1));disp(max(b1));
r2=(r1-zeropoint)/lookscale;r2(r2<0)=0;r2(r2>1)=1;
g2=(g1-zeropoint)/lookscale;g2(g2<0)=0;g2(g2>1)=1;
b2=(b1-zeropoint)/lookscale;b2(b2<0)=0;b2(b2>1)=1;
r3=floor(r2*binm)+1;
g3=floor(g2*binm)+1;
b3=floor(b2*binm)+1;
%% save lookuptable
if ~(exist('gradmag') && size(gradmag,1))    % no previous table
    gradmag=zeros(bins, bins, bins);
    countmap=gradmag;graddir=gradmag;
end

for i=1:length(r1)
    t1=countmap(r3(i), g3(i), b3(i));
    countmap(r3(i), g3(i), b3(i))=countmap(r3(i), g3(i), b3(i))+1;
    if t1    
        gradmag(r3(i), g3(i), b3(i))=(gradmag(r3(i), g3(i), b3(i))*t1+gradxseq(i))/(t1+1);
        a1=graddir(r3(i), g3(i), b3(i));a2=gradyseq(i);
        if a2-a1>pi 
            a2=a2-pi*2;
        elseif a1-a2>pi
            a2=a2+pi*2;
        end
        graddir(r3(i), g3(i), b3(i))=(a1*t1+a2)/(t1+1);

    else
        gradmag(r3(i), g3(i), b3(i))=gradxseq(i);
        graddir(r3(i), g3(i), b3(i))=gradyseq(i);
    end    
end
i=0;

