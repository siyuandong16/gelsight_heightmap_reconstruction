clear;clc; close all;

%% Hyper parameters

hp = struct();
hp.wuv = 0;
hp.wg = 1;
hp.beta = 40;

%% load data

data = load('data3.mat');
uv = double([data.data_struct.x'; data.data_struct.y']);
rgb = [data.data_struct.r'; data.data_struct.g'; data.data_struct.b'];
g = [data.data_struct.gx'; data.data_struct.gy'];

% data dimensions
[nuv, N] = size(uv);
[nrgb, ~] = size(rgb);
[ng, ~] = size(g);

% compute normalization
nrmlz = struct();

nrmlz.uv.min_val = min(uv, [], 2); 
nrmlz.uv.range = max(uv, [], 2) - nrmlz.uv.min_val; 

nrmlz.rgb.min_val = min(rgb, [], 2); 
nrmlz.rgb.range = max(rgb, [], 2) - nrmlz.rgb.min_val; 

nrmlz.g.min_val = min(g, [], 2); 
nrmlz.g.range = max(g, [], 2) - nrmlz.g.min_val; 


% normalize data to [0, 1]
uv_01 = (uv - nrmlz.uv.min_val)./nrmlz.uv.range;
rgb_01 = (rgb - nrmlz.rgb.min_val)./nrmlz.rgb.range;
g_01 = (g - nrmlz.g.min_val)./nrmlz.g.range;


% initial cluster locations
c = double([0.5, 0.25, 0.75, 0.75;
            0.5, 0.5,  0.25, 0.75]);

% number of clusters
[~,nc] = size(c);

% true map from (r, g, b) to (gx, gy)
A = rand(ng, nrgb, nc);


%% soft k means aglorithm

ct = 1; Nmax = 10;
% c = rand(nuv, nc);
cmap = winter(Nmax);

while ct < Nmax
    fprintf("----------- Iteration: %f ---------- \n", ct)
    
    %%%%%%%%%%%%%%%%% Compute Error %%%%%%%%%%%%%%%%%%%%
    err_uv = zeros(nc, N);
    err_g = zeros(nc, N);
    
    for i = 1:nc
        err_uv(i, :) = sum(bsxfun(@minus, uv_01, c(:, i)).^2);
        err_g(i, :) = sum((g_01 - A(:,:,i)*rgb_01).^2);
    end
    
    err_tot = hp.wuv * err_uv + hp.wg * err_g;
    
    %%%%%%%%%%%%%%%% Plotting %%%%%%%%%%%%%%%%%%%%%%
    figure(1);
    subplot(2, 2, 1);  hold on;
    plot(err_uv(:), '*', 'color', cmap(ct, :))
    %     ylim([0, 2])
    title("pixel distance")
    
    subplot(2, 2, 3);  hold on;
    plot(err_g(:), '*', 'color', cmap(ct, :))
    %     ylim([0, 2])
    title('gradient error')
    
        
    subplot(2,2,2);
    plot(c(1, :), c(2, :), '*', 'color', cmap(ct, :))
    xlim([0, 1]); ylim([0, 1])
    
    colors = [1, 0, 0, 0.2;
        0, 1, 0, 0.2;
        0, 0, 1, 0.2;
        0, 0, 0, 0.2];
    
    [uu, vv] = meshgrid(linspace(0, 1, 100), linspace(0, 1, 100));
    error_plot = 0 * uu;
    
    subplot(2,2,4); hold on;
    for i = 1:nc
        dist_err_i = (uu - c(1, i)).^2 + (vv - c(2, i)).^2;
        error_plot = error_plot + exp(-hp.beta*dist_err_i);
    end
    contourf(uu, vv, error_plot)
    xlim([0, 1]); ylim([0, 1])
%     colorbar;
    
    
    fprintf("Average error: %f \r", (mean(err_tot(:))))
    
    %     fprintf("Centroid error (distance): (%f, %f) \r", vecnorm(c - ctrue, 2))
    %     fprintf("Amatrix Error: (%f, %f) \r", ...
    %         norm(A(:, :, 1) - Atrue(:, :, 1), 'fro'), ...
    %         norm(A(:, :, 2) - Atrue(:, :, 2), 'fro'))
    
    %%%%%%%%%%% Soft assign pixels to clusters %%%%%%%%%
    
    z = zeros(nc, N);
    for i = 1:nc
        z(i, :) = exp(-hp.beta*err_tot(i, :))./sum(exp(-hp.beta*err_tot));
    end
    
    %%%%%%%%%%%%%% Compute new new c, A, b %%%%%%%%%%%%%%
    for i = 1:nc
        c(:, i) = sum(bsxfun(@times, uv_01, z(i, :)), 2)/sum(z(i, :));
        A(:, :, i) = lscov(rgb_01', g_01', 1 + 0*z(i, :))';
    end
    
    ct = ct + 1;
    
end


format_out = 'mm-dd-yy_HH:MM:SS';
save(['data3_calib_soft_' datestr(datetime('now'), format_out)], ...
    'A', 'c', 'hp', 'nrmlz');
%% Generate data

% % image size
% umin = 0; umax = 1;
% vmin = 0; vmax = 1;
%
% % data dimensions
% nuv = 2; nrgb = 3; ng = 2;
% n = 25; N = n^2;
%
% % uv data
% [uu, vv] = meshgrid(linspace(umin, umax, n),linspace(vmin, vmax, n));
% uv = [uu(:)'; vv(:)'];
%
% % rgb data
% rgb = rand(nrgb, N);
%
% % true cluster locations
% ctrue = [(umax-umin)/4, 3*(umax-umin)/4;
%     (vmax - vmin)/2, (vmax - vmin)/2];
%
% % number of clusters
% [~,nc] = size(ctrue);
%
% % true map from (r, g, b) to (gx, gy)
% Atrue = rand(ng, nrgb, nc);
%
% % true cluster assignment
% err_uv_true = zeros(nc, N);
% for i = 1:nc
%     err_uv_true(i, :) = sum(bsxfun(@minus, uv, ctrue(:, i)).^2);
% end
%
% [~, cluster_ind_true] = min(err_uv_true);
%
% % gradient data
% g = zeros(ng, N);
% for i = 1:nc
%     rgb_i = rgb(:, cluster_ind_true == i);
%     g(:,  cluster_ind_true == i) = Atrue(:, :, i)*rgb_i;
% end
%
% % initialize linear model
% A = repmat(eye(ng, nrgb), 1, 1, nc);



