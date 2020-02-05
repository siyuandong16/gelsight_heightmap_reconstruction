clear;clc; close all;

%% Hyper parameters

w1 = 0.0001;
w2 = 1;

%% Image information

% image size
umin = 0; umax = 1;
vmin = 0; vmax = 1;

% data dimensions
nuv = 2; nrgb = 3; ng = 2;
n = 25; N = n^2;

% fake data
u = linspace(umin, umax, n);
v = linspace(vmin, vmax, n);
[uu, vv] = meshgrid(u,v);
uv = [uu(:)'; vv(:)'];

g = zeros(ng, N);
rgb = rand(nrgb, N);


%% Initialize


% initial (u, v) positions
% c = [umin, umin, 0.5*(umin + umax), umax, umax;
%     vmin, vmax, 0.5*(vmin + vmax), vmin, vmax ];

ctrue = [(umax-umin)/4, 3*(umax-umin)/4;
    (vmax - vmin)/2, (vmax - vmin)/2];


[~,nc] = size(ctrue); % number of clusters



%% Generate gradient data
Atrue = rand(ng, nrgb, nc);

err_uv_true = zeros(nc, N);
for i = 1:nc
    err_uv_true(i, :) = sum(bsxfun(@minus, uv, ctrue(:, i)).^2);
end

[~, cluster_ind_true] = min(err_uv_true);

for i = 1:nc
    rgb_i = rgb(:, cluster_ind_true == i);
    g(:,  cluster_ind_true == i) = Atrue(:, :, i)*rgb_i;
end

% initialize linear model
A = repmat(eye(ng, nrgb), 1, 1, nc);

%% Initial assignment

err_uv_0 = zeros(nc, N);
err_g_0 = zeros(nc, N);

for i = 1:nc
    
    err_uv_0(i, :) = sum(bsxfun(@minus, uv, ctrue(:, i)).^2);
    err_g_0(i, :) = sum((g - A(:,:,i)*rgb).^2);
    
end

err_tot_0 = w1 * err_uv_0 + w2 * err_g_0;

[~, cluster_ind] = min(err_tot_0);

%% k means aglorithm

ct = 1; Nmax = 10;
c = rand(nuv, nc);

while ct < Nmax
    fprintf("----------- Iteration: %f ---------- \n", ct)
    
    %%%%%%%%%%%%%%%%% Compute Error %%%%%%%%%%%%%%%%%%%%
    err_uv = zeros(nc, N);
    err_g = zeros(nc, N);
    
    for i = 1:nc
        %         cluster_mask = cluster_ind == i;
        err_uv(i, :) = sum(bsxfun(@minus, uv, c(:, i)).^2);
        err_g(i, :) = sum((g - A(:,:,i)*rgb).^2);
    end
    
    figure(ct); clf; hold on;
    plot(err_uv(:), 'r*')
    plot(err_g(:), 'bo')
    ylim([0, 3])
    legend('pixel distance', 'gradient')
    
       
    err_tot = w1 * err_uv + w2 * err_g;
    
    mean(err_tot(:));
    
    disp(c - ctrue)
    disp(A - Atrue)    
    
    %%%%%%%%%%% Assign pixels to clusters %%%%%%%%%%%%%%
    [~, cluster_ind] = min(err_tot);
    
    %%%%%%%%%%%%%% Compute new new c, A, b %%%%%%%%%%%%%%
    for i = 1:nc
        ndata = sum(cluster_ind == i);
        c(:, i) = mean(uv(:, cluster_ind_true == i), 2);
        A(:, :, i) =  ([rgb(:, cluster_ind_true == i)]'\g(:, cluster_ind_true == i)')';
    end
    
    
    ct = ct + 1;
    
    
end



