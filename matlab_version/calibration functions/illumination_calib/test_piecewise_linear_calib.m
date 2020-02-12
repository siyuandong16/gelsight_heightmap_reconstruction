clear; clc; close all;


calib = load('data3_calib3.mat');
datadir = './test_image/';


u = readNPY([datadir 'x.npy']);
v = readNPY([datadir 'y.npy']);
uv = double([u(:)'; v(:)']);

r = readNPY([datadir 'r.npy']);
g = readNPY([datadir 'g.npy']);
b = readNPY([datadir 'b.npy']);
rgb = double([r(:)'; g(:)'; b(:)']);


mask = double(1 - readNPY([datadir 'mask.npy']));

[~, nc] = size(calib.c);
[~, N] = size(uv); 

gx = readNPY([datadir 'gx_t.npy']);
gy = readNPY([datadir 'gy_t.npy']);

%%%%%%%%%%%%%%%%% Compute Error %%%%%%%%%%%%%%%%%%%%
err_uv = zeros(nc, N);
% err_g = zeros(nc, N);

for i = 1:nc
    %         cluster_mask = cluster_ind == i;
    err_uv(i, :) = sum(bsxfun(@minus, uv, calib.c(:, i)).^2);
%     err_g(i, :) = sum((g - calib.A(:,:,i)*rgb).^2);
end

% assign cluster
[~, cluster_ind] = min(err_uv);

% predict gradient 
ghat = zeros(2, N); 
for i = 1:nc    
    ghat(:, cluster_ind == i) = calib.A(:, :, i) * rgb(:, cluster_ind == i);
end

ghat_img = zeros([size(u), 2]); 
ghat_img(:, :, 1) = reshape(ghat(1, :), size(u)).*mask;
ghat_img(:, :, 2) = reshape(ghat(2, :), size(u)).*mask;


depth_img = fast_poisson2(ghat_img(:, :, 1),ghat_img(:, :, 2));
depth_img_t = fast_poisson2(gx,gy);

figure(1); clf; hold on;
subplot(1,2,1);
imshow(ghat_img(:, :, 1) - min(ghat(1, :)));
title("G_y")
subplot(1,2,2); hold on;
imshow(ghat_img(:, :, 2) - min(ghat(2, :)));
title("G_x")


figure(2); 
% subplot(1,2,1); hold on; 
mesh(depth_img);
% zlabel('Depth'); ylabel('Y')
% subplot(1,2,2); hold on; view(0,0)
% mesh(depth_img); xlabel('X')

figure(3); 
% subplot(1,2,1); hold on; view(90,0)
mesh(depth_img_t);


% data = load('