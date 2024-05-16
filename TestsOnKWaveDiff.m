%% Initialize
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
cutoff_end = 70;
H = load(['VarResults\H_var_', num2str(sigma)]).H;
analyze_H(H,X,Y,Z,mua,mus,g,cutoff,cutoff_end, "ValoMC")
fprintf("expected mu  %d\n", mu_expected_value)
fprintf("std mu  %d\n", sqrt(mu_variance))


%% K-Wave Regression
H_k_wave_recon = load(['VarResults\H_recon_var_', num2str(sigma)]).H_k_wave_recon;
sensor_data = load(['VarResults\Sensor_Data_var_', num2str(sigma)]).sensor_data;
cutoff = 30;
cutoff_end = 70;
analyze_H(H_k_wave_recon, X, Y, Z, mua, mus, g, cutoff, cutoff_end, 'K-Wave');



%% smoothed regression
H_smooth = smooth3(H);
analyze_H(H_smooth,X,Y,Z,mua,mus,g,cutoff,cutoff_end, "Smooth")
fprintf("expected mu  %d\n", mu_expected_value)
fprintf("std mu  %d\n", sqrt(mu_variance))
%% mean regression
H_mean = squeeze(mean(H, [1,3]));
H_mean = repmat(H_mean,21,1,21);
analyze_H(H_mean,X,Y,Z,mua,mus,g,cutoff,cutoff_end, "Mean")
fprintf("expected mu  %d\n", mu_expected_value)
fprintf("std mu  %d\n", sqrt(mu_variance))

%% mean k-wave
H_k_wave = permute(H_mean, [2,1,3]);

% p0
kgrid = kWaveGrid(Nx, dx*1e-3, Ny, dy*1e-3, Nz, dz*1e-3);
medium.sound_speed = 1500;    % [m/s]
medium.density = 1000;        % [kg/m^3]
source.p0 = gruneisen_coef.*H_k_wave;

% sensor
sensor.mask = zeros(Nx, Ny, Nz);
sensor.mask(1,:,:) = 1;
% input arguments
input_args = {'PlotLayout', true, 'PlotPML', false, ...
    'DataCast', 'gpuArray-single', 'CartInterp', 'nearest', 'PMLInside', false, 'DataRecast', true};


% run the simulation
sensor_data = kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:});
%%
% reshape sensor data to y, z, t
sensor_data_rs = reshape(sensor_data, Ny, Nz, kgrid.Nt);

% reconstruct the initial pressure
p_xyz = kspacePlaneRecon(sensor_data_rs, kgrid.dy, kgrid.dz, kgrid.dt, ...
    medium.sound_speed, 'DataOrder', 'yzt', 'PosCond', true, 'Plot', true);



% define a second k-space grid using the dimensions of p_xyz
[Nx_recon, Ny_recon, Nz_recon] = size(p_xyz);
kgrid_recon = kWaveGrid(Nx_recon, kgrid.dt * medium.sound_speed, Ny_recon, dy, Nz_recon, dz);

% resample p_xyz to be the same size as source.p0
p_xyz_rs = interp3(kgrid_recon.y, kgrid_recon.x - min(kgrid_recon.x(:)),kgrid_recon.z, p_xyz, kgrid.y, kgrid.x - min(kgrid.x(:)), kgrid.z);

%%
H_k_wave_recon = permute(p_xyz_rs, [2,1,3])./gruneisen_coef;
cutoff = 30;
cutoff_end = 70;
analyze_H(H_k_wave_recon,X,Y,Z,mua,mus,g,cutoff,cutoff_end, "K-Wave")




%% comparing source and reconstruction
H_mid = H(Z==0 & Y==0);
H_k_wave_recon_mid = 2.1.*H_k_wave_recon(Z==0 & Y==0);
figure;
plot(x_arr, H_mid, "DisplayName","H");
hold on
plot(x_arr, H_k_wave_recon_mid, "DisplayName","H k-wave");
legend;
%% Making the sensors smaller

H_k_wave = permute(H, [2,1,3]);

% p0
kgrid = kWaveGrid(Nx, dx*1e-3, Ny, dy*1e-3, Nz, dz*1e-3);
medium.sound_speed = 1500;    % [m/s]
medium.density = 1000;        % [kg/m^3]
source.p0 = gruneisen_coef.*H_k_wave;

% sensor
sensor.mask = zeros(Nx, Ny, Nz);
sensor.mask(1,9:13,9:13) = 1;
% input arguments
input_args = {'PlotLayout', true, 'PlotPML', false, ...
    'DataCast', 'gpuArray-single', 'CartInterp', 'nearest', 'PMLInside', false, 'DataRecast', true};

% run the simulation
sensor_data = kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:});
%%
% reshape sensor data to y, z, t
sensor_data_rs = reshape(sensor_data, 5, 5, kgrid.Nt);

% reconstruct the initial pressure
p_xyz = kspacePlaneRecon(sensor_data_rs, kgrid.dy, kgrid.dz, kgrid.dt, ...
    medium.sound_speed, 'DataOrder', 'yzt', 'PosCond', true, 'Plot', true);



% define a second k-space grid using the dimensions of p_xyz
[Nx_recon, Ny_recon, Nz_recon] = size(p_xyz);
kgrid_recon = kWaveGrid(Nx_recon, kgrid.dt * medium.sound_speed, Ny_recon, dy, Nz_recon, dz);

% resample p_xyz to be the same size as source.p0
p_xyz_rs = interp3(kgrid_recon.y, kgrid_recon.x - min(kgrid_recon.x(:)),kgrid_recon.z, p_xyz, kgrid.y, kgrid.x - min(kgrid.x(:)), kgrid.z);

%%
H_k_wave_recon = permute(p_xyz_rs, [2,1,3])./gruneisen_coef;
cutoff = 30;
cutoff_end = 70;
analyze_H_new(H_k_wave_recon,X,Y,Z,mua,mus,g,cutoff,cutoff_end, "K-Wave")


%% Zeroing Frequencies polar symmetry

H_fft = fftshift(fftn(H))

threshold = 2e8;

evanescent_idx = medium.sound_speed.*sqrt((kgrid.kx.^2 + kgrid.ky.^2 + kgrid.kz.^2))>threshold;
H_fft(evanescent_idx) = 0;
H_reconstructed = ifftn(ifftshift(H_fft));


analyze_H(abs(H_reconstructed),X,Y,Z,mua,mus,g,cutoff,cutoff_end, "up_to_w")




%% Zeroing y,z Frequencies

H_fft = fftshift(fftn(H))

threshold = 1e6;

unseen_idx = medium.sound_speed.*sqrt((kgrid.ky.^2 + kgrid.kz.^2))>threshold;

H_fft(unseen_idx) = 0;
H_reconstructed = smooth3(ifftn(ifftshift(H_fft)));


analyze_H(abs(H_reconstructed),X,Y,Z,mua,mus,g,cutoff,cutoff_end, "ky,kz0")

%% Zeroing y,z Frequencies using geometry

H_fft = fftshift(fftn(H))

h = 0.5;
W = 5;
unseen_idx = sqrt((kgrid.ky.^2 + kgrid.kz.^2))>kgrid.kx.*(W/h);
H_fft(unseen_idx) = 0;
H_reconstructed = smooth3(ifftn(ifftshift(H_fft)));


analyze_H(abs(H_reconstructed),X,Y,Z,mua,mus,g,cutoff,cutoff_end, "Geometry")

%% Combining Zeroing

H_fft = fftshift(fftn(H))

evanescent_threshold = 2e8;
unseen_threshold = 1e5;

evanescent_idx = medium.sound_speed.*sqrt((kgrid.kx.^2 + kgrid.ky.^2 + kgrid.kz.^2))>evanescent_threshold;
H_fft(evanescent_idx) = 0;
unseen_idx = medium.sound_speed.*sqrt((kgrid.ky.^2 + kgrid.kz.^2))>unseen_threshold;
H_fft(unseen_idx) = 0;
H_reconstructed = ifftn(ifftshift(H_fft));


analyze_H(abs(H_reconstructed),X,Y,Z,mua,mus,g,cutoff,cutoff_end, "combinedEffects")



