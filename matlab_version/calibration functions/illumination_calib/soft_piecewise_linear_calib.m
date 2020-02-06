clear;clc; close all;

%% Hyper parameters

w1 = 0.0001;
w2 = 1;
beta = 10;

%% Generate data

% image size
umin = 0; umax = 1;
vmin = 0; vmax = 1;

% data dimensions
nuv = 2; nrgb = 3; ng = 2;
n = 25; N = n^2;

% uv data
[uu, vv] = meshgrid(linspace(umin, umax, n),linspace(vmin, vmax, n));
uv = [uu(:)'; vv(:)'];

% rgb data
rgb = rand(nrgb, N);

% true cluster locations
ctrue = [(umax-umin)/4, 3*(umax-umin)/4;
    (vmax - vmin)/2, (vmax - vmin)/2];

% number of clusters
[~,nc] = size(ctrue);

% true map from (r, g, b) to (gx, gy)
Atrue = rand(ng, nrgb, nc);

% true cluster assignment 
err_uv_true = zeros(nc, N);
for i = 1:nc
    err_uv_true(i, :) = sum(bsxfun(@minus, uv, ctrue(:, i)).^2);
end

[~, cluster_ind_true] = min(err_uv_true);

% gradient data
g = zeros(ng, N);
for i = 1:nc
    rgb_i = rgb(:, cluster_ind_true == i);
    g(:,  cluster_ind_true == i) = Atrue(:, :, i)*rgb_i;
end

% initialize linear model
A = repmat(eye(ng, nrgb), 1, 1, nc);

%% soft k means aglorithm

ct = 1; Nmax = 10;
c = rand(nuv, nc);
cmap = winter(Nmax);

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
  
    %%%%%%%%%%%%%%%% Plotting %%%%%%%%%%%%%%%%%%%%%%
    figure(1);
    subplot(2, 1, 1);  hold on;
    plot(err_uv(:), '*', 'color', cmap(ct, :))
%     ylim([0, 2])
    title("pixel distance")
    
    subplot(2, 1, 2);  hold on;
    plot(err_g(:), '*', 'color', cmap(ct, :))
%     ylim([0, 2])
    title('gradient error')
    
    err_tot = w1 * err_uv + w2 * err_g;
    
    fprintf("Average error: %f \r", (mean(err_tot(:))))
    
    fprintf("Centroid error (distance): (%f, %f) \r", vecnorm(c - ctrue, 2))
    fprintf("Amatrix Error: (%f, %f) \r", ...
        norm(A(:, :, 1) - Atrue(:, :, 1), 'fro'), ... 
        norm(A(:, :, 2) - Atrue(:, :, 2), 'fro'))    
    
    %%%%%%%%%%% Soft assign pixels to clusters %%%%%%%%%
    
    z = zeros(nc, N); 
    for i = 1:nc
        z(i, :) = exp(-beta*err_tot(i, :))./sum(exp(-beta*err_tot));
    end
    
    %%%%%%%%%%%%%% Compute new new c, A, b %%%%%%%%%%%%%%
    for i = 1:nc
        c(:, i) = sum(uv.*z(i, :), 2)/sum(z(i, :));
        
%         lhs = zeros(nrgb, ng); 
%         rhs = zeros(nrgb, nrgb); 
%         for j = 1:N
%             lhs = lhs + z(i,j)*rgb(:, j)*g(:, j)'; 
%             rhs = rhs + z(i,j)*rgb(:, j)*rgb(:, j)'; 
%             
%         end
        A(:, :, i) = lscov(rgb', g', z(i, :))';
%         A(:, :, i) = (rhs' \ lhs)';
    end
    
    ct = ct + 1;    
    
end



