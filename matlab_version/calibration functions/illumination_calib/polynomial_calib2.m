clear;clc; close all;

red = readNPY('training_data/r.npy');
green = readNPY('training_data/g.npy');
blue = readNPY('training_data/b.npy');
rgb = double([red';green';blue']);

gx = readNPY('training_data/gx.npy');
gy = readNPY('training_data/gy.npy');
grad = double([gx';gy']);
 
% data dimensions
[nrgb, N] = size(rgb);
[ng, ~] = size(grad);


% compute normalization
nrmlz = struct();
% 
% nrmlz.rgb.mean_val = mean(rgb, 2);
% nrmlz.rgb.std_dev = std(rgb, [], 2);
% 
% nrmlz.g.mean_val = mean(grad, 2);
% nrmlz.g.std_dev = std(grad, [], 2);
 
% % normalize data to [0, 1]
rgb_01 = rgb; %(rgb  - nrmlz.rgb.mean_val)./nrmlz.rgb.std_dev;
g_01 = grad; %(grad - nrmlz.g.mean_val)./nrmlz.g.std_dev;

% build xdata
r = rgb_01(1, :); g = rgb_01(2, :); b = rgb_01(3, :); 
xdata = [ones(1,N); r; g; r.*g.*b; r.*b.^2; r.*g.^2; b.*r.^2; b.*g.^2; 
    g.*r.^2; g.*b.^2; r.^3; b.^3; g.^3]; 

A = (xdata'\g_01')';
g_01_hat = A*xdata;


figure(1); clf;
subplot(1,2,1);
plot(g_01(1,:), g_01_hat(1, :), 'b.');
% xlim([-3, 3]); ylim([-3, 3])
axis square;
subplot(1,2,2);
plot(g_01(2,:), g_01_hat(2, :), 'b.');
% xlim([-3, 3]); ylim([-3, 3])
axis square; 

nskip = 1; 
figure(2); 
scatter3(rgb(1, 1:nskip:end), rgb(2, 1:nskip:end), ...
    rgb(3, 1:nskip:end), [], g_01(1, 1:nskip:end))
xlabel('R'); 
ylabel('G');
zlabel('B'); 
caxis([-0.3, 0.3])
colorbar;

figure(3); 
scatter3(rgb(1, 1:nskip:end), rgb(2, 1:nskip:end), ...
    rgb(3, 1:nskip:end), [], g_01_hat(1, 1:nskip:end))
xlabel('R'); 
ylabel('G');
zlabel('B'); 
caxis([-0.3, 0.3])
colorbar;

g_01_err = sum((g_01 - g_01_hat).^2);
disp(mean(g_01_err))
disp(max(g_01_err))

format_out = 'mm-dd-yy_HH:MM:SS';
save(['data3_calib_poly_single_' datestr(datetime('now'), format_out)], ...
    'A', 'nrmlz');

