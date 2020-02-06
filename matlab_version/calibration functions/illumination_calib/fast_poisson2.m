% function [img_direct] = poisson_solver_function(gx,gy,boundary_image)
% Inputs; Gx and Gy -> Gradients
% Boundary Image -> Boundary image intensities
% Gx Gy and boundary image should be of same size
% code borrowed from Ramesh Raskar,
% http://www.merl.com/people/raskar/photo/code.pdf 
function [img_direct] = fast_poisson2(gx,gy)

%     ticPoisson = tic;


	[ydim,xdim] = size(gx);
	gxx = zeros(ydim,xdim); 
	gyy = zeros(ydim,xdim);

% 	f = zeros(ydim,xdim); 
	j = 1:ydim-1; 
	k = 1:xdim-1;

	% Laplacian
	gyy(j+1,k) = gy(j+1,k) - gy(j,k); 
	gxx(j,k+1) = gx(j,k+1) - gx(j,k);
    
	f = gxx + gyy; 
    
%     subplot 223;imshow(gxx,[]);title('x Laplacian');
%     subplot 224;imshow(gyy,[]);title('y Laplacian');
%     figure(1);imshow(f,[]);title('sum of Laplacian');
    
    %% comment from here for 2nd derivative ONLY
    
	clear j k gxx gyy gyyd gxxd


	% DST Sine Transform algo starts here
	f2 = f(2:end-1,2:end-1); 
	clear f1

	%compute sine transform
	tt = dst(f2); 
	f2sin = dst(tt')'; 
	clear f2

	%compute Eigen Values
	[x,y] = meshgrid(1:xdim-2,1:ydim-2); 
	denom = (2*cos(pi*x/(xdim-1))-2) + (2*cos(pi*y/(ydim-1)) - 2);
    disp(denom)
	%divide
	f3 = f2sin./denom; 
	clear f2sin x y

	%compute Inverse Sine Transform
	tt = idst(f3); clear f3; img_tt = idst(tt')'; 
	clear tt

	% put solution in inner points; outer points obtained from boundary image
	img_direct = zeros(ydim, xdim);
% 	img_direct(2:end-1,2:end-1) = 0;
	img_direct(2:end-1,2:end-1) = img_tt;
    
%     t=toc(ticPoisson);
%     disp(['fast Poisson time is ' num2str(t)]);

% 	return;

