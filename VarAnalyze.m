close all;
clear all;
clc
%% simulation parameters
dx = 0.01;
dy = 0.5;
dz = 0.5;
Nx = 201;
Ny = 21;
Nz = 21;
x_arr = -dx*(Nx-1)/2:dx:dx*(Nx-1)/2;
y_arr = -dy*(Ny-1)/2:dy:dy*(Ny-1)/2;
z_arr = -dz*(Nz-1)/2:dz:dz*(Nz-1)/2;
mua=0.3;
mus=50;
g = 0.9;
sigma = 0;
gruneisen_coef = 0.5;

%% create grid
[X,Y,Z] = meshgrid(x_arr,y_arr,z_arr);

scattering_matrix = mus + mus.*sigma.*randn(size(X));
negative_idx = scattering_matrix<0;
scattering_matrix(negative_idx)=0;
absorption_matrix = mua + mua.*sigma.*randn(size(X));
negative_idx = absorption_matrix<0;
absorption_matrix(negative_idx)=0;
mu_matrix = sqrt(3.*absorption_matrix.*(absorption_matrix+(1-g).*scattering_matrix));
mu_expected_value = mean(mu_matrix(:));
mu_variance = var(mu_matrix(:));


%% Regression

cutoff = 30;
cutoff_end = 90;
H = load(['VarResults\H_var_', num2str(sigma)]).H;
analyze_H_print(H,X,Y,Z,mua,mus,g,cutoff,cutoff_end, "ValoMC", sigma)
fprintf("expected mu  %d\n", mu_expected_value)
fprintf("std mu  %d\n", sqrt(mu_variance))


%% K-Wave Regression
H_k_wave_recon = load(['VarResults\H_recon_var_', num2str(sigma)]).H_k_wave_recon;
sensor_data = load(['VarResults\Sensor_Data_var_', num2str(sigma)]).sensor_data;
cutoff = 30;
cutoff_end = 70;
analyze_H_print(H_k_wave_recon, X, Y, Z, mua, mus, g, cutoff, cutoff_end, 'K-Wave', sigma);


%% comparing source and reconstruction
H_mid = H(Z==0 & Y==0);
H_k_wave_recon_mid = H_k_wave_recon(Z==0 & Y==0);
figure;
plot(x_arr, H_mid, "DisplayName","H");
hold on
plot(x_arr, H_k_wave_recon_mid, "DisplayName","H k-wave");
legend;
xlabel("Depth [mm]")
legend()
title("Absorbtion Before and After Reconstruction")
filename = join(['BeforeAndAfter_', num2str(sigma), '.png'], '');
saveas(gcf, filename);