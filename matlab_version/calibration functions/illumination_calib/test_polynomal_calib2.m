clear; clc; close all;


calib = load('data3_calib_poly_single_02-12-20_13:48:50.mat'); % load training data

%
u = readNPY('test_image_single/x.npy');
v = readNPY('test_image_single/y.npy');

% rgb
r = readNPY('test_image_single/r.npy');
g = readNPY('test_image_single/g.npy');
b = readNPY('test_image_single/b.npy');
rgb = double([r(:)'; g(:)'; b(:)']);
rgb_01 = rgb; %(rgb - calib.nrmlz.rgb.mean_val)./calib.nrmlz.rgb.std_dev;

% mask
mask = double(1 - readNPY('test_image_single/mask.npy'));

[~, N] = size(rgb_01); 

%%%%%%%%%%%%%%%%% Compute  %%%%%%%%%%%%%%%%%%%%
r = rgb_01(1, :); g = rgb_01(2, :); b = rgb_01(3, :); 
xdata = [ones(1,N); r; g; b; r.*g.*b; r.*b.^2; r.*g.^2; b.*r.^2; b.*g.^2; 
    g.*r.^2; g.*b.^2; r.^3; b.^3; g.^3]; 
ghat_01 = calib.A*xdata; 

% denormalize
ghat = ghat_01; %ghat_01.*calib.nrmlz.g.std_dev + calib.nrmlz.g.mean_val;

% true values
gx = readNPY('gradient_gt_single/gx_t.npy');
gy = readNPY('gradient_gt_single/gy_t.npy');


mask2 = (u > 143) & (u < 219) & (v > 148) & (v < 208); 

ghat_img = zeros([size(mask), 2]); 
ghat_img(:, :, 1) = reshape(ghat(1, :), size(mask)); %.*mask2;
ghat_img(:, :, 2) = reshape(ghat(2, :), size(mask)); %.*mask2;


%% Plotting
figure(1); clf;

subplot(2,2,1);
imshow(ghat_img(:, :, 1), []);
% title("G_x")

subplot(2,2,2);
imshow(gx, [])


subplot(2,2,3);
imshow(ghat_img(:, :, 2), []);
% title("G_y")

subplot(2,2,4);
imshow(gy, [])

figure(2); clf; 
subplot(1, 2, 1)
imshow(ghat_img(:, :, 1) - gx, []);
subplot(1,2,2)
imshow(ghat_img(:, :, 2) - gy, []);


depth_img = fast_poisson2(ghat_img(:, :, 1),ghat_img(:, :, 2));
figure(3); 
mesh(depth_img);


depth_img_t = fast_poisson2(gx,gy);
figure(4); 
mesh(depth_img_t);
