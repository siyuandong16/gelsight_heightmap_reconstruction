function [gradxout, gradyout]=LookuptableSmooth(bins, gradx, grady, countmap)
    if ~countmap(1,1,1)
        countmap(1,1,1)=1;
        gradx(1,1,1)=0;grady(1,1,1)=0;
    end
    
    validid=find(countmap);
    if length(validid)==bins^3  %no interpolation needed
        gradxout0=gradx;gradyout0=grady;
        return;
    end
    [xout, yout,zout]=meshgrid(1:bins, 1:bins, 1:bins);
%     gradxout = interp3(xout(validid),yout(validid),zout(validid),gradx(validid),...
%         xout,yout,zout, 'Cubic','spline');
%     gradyout = interp3(xout(validid),yout(validid),zout(validid),gradx(validid),...
%         xout,yout,zout, 'Cubic', 'spline');
    
    % not know what happened if the boundry is NaN
%     i=0;


%%%%%%%%%%%%%% NEED TO CONSIDER THE 
%% Method: get the current value to the nearest
gradxout0=gradx;gradyout0=grady;
invalidid=find(~countmap);

% invalidid=find(isnan(gradxout) | isnan(gradyout));
% validid=find(~isnan(gradxout) & ~isnan(gradyout));
xvalid=xout(validid);yvalid=yout(validid);zvalid=zout(validid);
for i=1:length(invalidid)
    t3=invalidid(i);
    t1=(xvalid-xout(t3)).^2+(yvalid-yout(t3)).^2+(zvalid-zout(t3)).^2;
    [t1 t2]=min(t1);
    gradxout0(t3)=gradx(validid(t2));
    gradyout0(t3)=grady(validid(t2));
end
gradxout=gradxout0;
gradyout=gradyout0;
% gradxout=tan(gradxout0).*cos(gradyout0);
% gradyout=tan(gradxout0).*sin(gradyout0);
