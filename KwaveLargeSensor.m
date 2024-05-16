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
% analyze_H_print(H,X,Y,Z,mua,mus,g,cutoff,cutoff_end, "ValoMC", sigma)
fprintf("expected mu  %d\n", mu_expected_value)
fprintf("std mu  %d\n", sqrt(mu_variance))



%% Making the sensors bigger
jump = 5;
dx_sensor = 0.05;
dy_sensor = 0.5;
dz_sensor = 0.5;
Nx_sensor = 41;
Ny_sensor = (Ny - 1)*jump + 1;
Nz_sensor = (Nz - 1)*jump + 1;
x_arr_sensor = -dx_sensor*(Nx_sensor-1)/2:dx_sensor:dx_sensor*(Nx_sensor-1)/2;
y_arr_sensor = -dy_sensor*(Ny_sensor-1)/2:dy_sensor:dy_sensor*(Ny_sensor-1)/2;
z_arr_sensor = -dz_sensor*(Nz_sensor-1)/2:dz_sensor:dz_sensor*(Nz_sensor-1)/2;
[X_sensor,Y_sensor,Z_sensor] = meshgrid(x_arr_sensor,y_arr_sensor,z_arr_sensor);

H_sensor = interp3(X,Y,Z, H, X_sensor, Y_sensor, Z_sensor);
out_of_bounds_idx = abs(X_sensor)>abs(max(x_arr(:))) | abs(Y_sensor)>abs(max(y_arr(:))) | abs(Z_sensor)>abs(max(z_arr(:)));
H_sensor(out_of_bounds_idx) = 0;
%%
figure;
slice(X_sensor, Y_sensor, Z_sensor, H_sensor, 0, 0, 0);
xlabel('x [mm]');
ylabel('y [mm]');
zlabel('z [mm]');
view(315,25);
hold


snapnow;

%%

H_k_wave = permute(H_sensor, [2,1,3]);

% p0
kgrid = kWaveGrid(Nx_sensor, dx_sensor*1e-3, Ny_sensor, dy_sensor*1e-3, Nz_sensor, dz_sensor*1e-3);
medium.sound_speed = 1500;    % [m/s]
medium.density = 1000;        % [kg/m^3]
source.p0 = gruneisen_coef.*H_k_wave;

% sensor
mask1 = zeros(Nx_sensor, Ny_sensor, Nz_sensor);
mask2 = zeros(Nx_sensor, Ny_sensor, Nz_sensor);
mask1(1,1:jump:Ny_sensor,:) = 1;
mask2(1,:,1:jump:Nz_sensor) = 1;
sensor.mask = mask1.*mask2;

% input arguments
input_args = {'PlotLayout', true, 'PlotPML', false, ...
    'DataCast', 'gpuArray-single', 'CartInterp', 'nearest', 'PMLInside', false, 'DataRecast', true};

% run the simulation
sensor_data = kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:});
%%

figure;
plot(sensor_data(1,:))

%%
% reshape sensor data to y, z, t
sensor_data_rs = reshape(sensor_data, Ny, Nz, kgrid.Nt);

% reconstruct the initial pressure
p_xyz = kspacePlaneRecon(sensor_data_rs, kgrid.dy, kgrid.dz, kgrid.dt, ...
    medium.sound_speed, 'DataOrder', 'yzt', 'PosCond', true, 'Plot', true);



% define a second k-space grid using the dimensions of p_xyz
[Nx_recon, Ny_recon, Nz_recon] = size(p_xyz);
kgrid_recon = kWaveGrid(Nx_recon, kgrid.dt * medium.sound_speed, Ny_recon, dy_sensor, Nz_recon, dz_sensor);

% resample p_xyz to be the same size as source.p0
p_xyz_rs = interp3(kgrid_recon.y, kgrid_recon.x - min(kgrid_recon.x(:)),kgrid_recon.z, p_xyz, kgrid.y, kgrid.x - min(kgrid.x(:)), kgrid.z);
H_k_wave_recon = permute(p_xyz_rs, [2,1,3])./gruneisen_coef;
H_k_wave_recon = interp3(X_sensor,Y,Z, H_k_wave_recon, X_sensor, Y_sensor, Z_sensor);
%%

cutoff = 30;
cutoff_end = 70;
analyze_H_print(H_k_wave_recon,X,Y,Z,mua,mus,g,cutoff,cutoff_end, "K-Wave", sigma)

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