clear; clc; close all;


datadir = './data_test_2/'; % test date set

calib = load('data3_calib_poly_02-10-20_14:51:14'); % load training data

r = readNPY([datadir 'r.npy']);
g = readNPY([datadir 'g.npy']);
b = readNPY([datadir 'b.npy']);
rgb = (double([r(:)'; g(:)'; b(:)']) - calib.nrmlz.rgb.min_val) ...
    ./calib.nrmlz.rgb.range; 


%%%%%%%%%%%%%%%%% Compute Error %%%%%%%%%%%%%%%%%%%%
xdata = poly_nd(calib.hp.nfit, rgb(1, :), rgb(2, :), rgb(3, :))';

ghat = calib.A*xdata; 


ghat_img = zeros([size(r), 2]); 
ghat_img(:, :, 1) = reshape(ghat(1, :), size(r));
ghat_img(:, :, 2) = reshape(ghat(2, :), size(r));

figure(1); clf; hold on;
subplot(1,2,1);
imshow(ghat_img(:, :, 1) - min(ghat(1, :)));
title("G_y")
subplot(1,2,2); hold on;
imshow(ghat_img(:, :, 2) - min(ghat(2, :)));
title("G_x")

depth_img = fast_poisson2(ghat_img(:, :, 2),ghat_img(:, :, 1));

figure(2); 
subplot(1,2,1); hold on; view(90,0)
mesh(depth_img);
zlabel('Depth'); ylabel('Y')
subplot(1,2,2); hold on; view(0,0)
mesh(depth_img); xlabel('X')