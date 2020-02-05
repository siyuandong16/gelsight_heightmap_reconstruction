clear; clc; close all; 

%% LOAD TRAINING DATA

train_datadir = 'data3'; 
train_npy_data = dir(fullfile(train_datadir, '*.npy')); 

train_data_struct = struct();

for i=1:numel(train_npy_data)
    
    fname = train_npy_data(i).name;
    fpath = fullfile(train_datadir, fname);
    train_data_struct.(strtok(fname, '.')) = double(readNPY(fpath));
    
end

xtrain = [train_data_struct.r, train_data_struct.g, train_data_struct.b, ...
    train_data_struct.x, train_data_struct.y];

%% FITTING 

poly_order = 1; 
p_gx = polyfitn(xtrain, train_data_struct.gx, poly_order);
p_gy = polyfitn(xtrain, train_data_struct.gy, poly_order);

%% LOAD TESTING DATA

test_datadir = 'data_test_2'; 
test_npy_data = dir(fullfile(test_datadir, '*.npy')); 

test_data_struct = struct();

for i=1:numel(test_npy_data)
    
    fname = test_npy_data(i).name;
    fpath = fullfile(test_datadir, fname);
    test_data_struct.(strtok(fname, '.')) = double(readNPY(fpath));
    
end

img_size = size(test_data_struct.x);

testx = [test_data_struct.r(:), test_data_struct.g(:), test_data_struct.b(:), ...
    test_data_struct.x(:), test_data_struct.y(:)];

gx_test = reshape(polyvaln(p_gx, testx), img_size);
gy_test = reshape(polyvaln(p_gy, testx), img_size);


%% PLOTTING 
figure(1); clf; hold on;
subplot(1,2,1);
imshow(gx_test)
subplot(1,2,2);
imshow(gy_test) 
