clear; clc; close all;


datadir = './data_test_2/'; % test date set

calib = load('data3_calib_soft_02-10-20_14:04:08.mat'); % load training data
beta = calib.hp.beta;

u = readNPY([datadir 'x.npy']);
v = readNPY([datadir 'y.npy']);
uv = (double([u(:)'; v(:)']) - calib.nrmlz.uv.min_val) ...
    ./calib.nrmlz.uv.range;

r = readNPY([datadir 'r.npy']);
g = readNPY([datadir 'g.npy']);
b = readNPY([datadir 'b.npy']);
rgb = (double([r(:)'; g(:)'; b(:)']) - calib.nrmlz.rgb.min_val) ...
    ./calib.nrmlz.rgb.range;

[~, nc] = size(calib.c);
[~, N] = size(uv); 


%%%%%%%%%%%%%%%%% Compute Error %%%%%%%%%%%%%%%%%%%%
err_uv = zeros(nc, N);

for i = 1:nc
    err_uv(i, :) = sum(bsxfun(@minus, uv, calib.c(:, i)).^2);
end

% soft assign cluster
z = zeros(nc, N);
for i = 1:nc
    z(i, :) = exp(-beta*err_uv(i, :))./sum(exp(-beta*err_uv));
end

% predict gradient 
ghat = zeros(2, N); 
for i = 1:nc    
    ghat = ghat +  calib.A(:, :, i) * rgb.*z(i,:);
end

ghat_img = zeros([size(u), 2]); 
ghat_img(:, :, 1) = reshape(ghat(1, :), size(u));
ghat_img(:, :, 2) = reshape(ghat(2, :), size(u));

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