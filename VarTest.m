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
sigma = 0.1;
gruneisen_coef = 0.5;

%% create grid
vmcmesh = createGridMesh(x_arr, y_arr, z_arr); % function provided by ValoMC
vmcmedium = createMedium(vmcmesh);
[X,Y,Z] = meshgrid(x_arr,y_arr,z_arr); % Matlab function
scattering_matrix = mus + mus.*sigma.*randn(size(X));
negative_idx = scattering_matrix<0;
scattering_matrix(negative_idx)=0;
absorption_matrix = mua + mua.*sigma.*randn(size(X));
mu_matrix = sqrt(3.*absorption_matrix.*(absorption_matrix+(1-g).*scattering_matrix));
mu_expected_value = mean(mu_matrix(:));
mu_variance = var(mu_matrix(:));
negative_idx = absorption_matrix<0;
absorption_matrix(negative_idx)=0;

vmcmedium.scattering_coefficient = scattering_matrix;
vmcmedium.absorption_coefficient = absorption_matrix;
vmcmedium.scattering_anisotropy = g;        
vmcmedium.refractive_index = 1;



%% create source
clear vmcboundary;
vmcboundary = createBoundary(vmcmesh, vmcmedium);
lightsource = findBoundaries(vmcmesh, 'direction', [-1.5 0 0], [-0.5 0 0], 4.5);
vmcboundary.lightsource(lightsource) = {'direct'};
options.photon_count = 1e6;
solution = ValoMC(vmcmesh, vmcmedium, vmcboundary, options);

%% Regression

cutoff = 30;
cutoff_end = 90;
H = vmcmedium.absorption_coefficient .* solution.grid_fluence*1e6;
analyze_H(H,X,Y,Z,mua,mus,g,cutoff,cutoff_end, "ValoMC")
fprintf("expected mu  %d\n", mu_expected_value)
fprintf("std mu  %d\n", sqrt(mu_variance))


%% k-wave simulation

H_k_wave = permute(H, [2,1,3]);

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
    'DataCast', 'single', 'CartInterp', 'nearest', 'PMLInside', false};

% run the simulation
sensor_data = kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:});
%%
% reshape sensor data to y, z, t
sensor_data_rs = reshape(sensor_data, Ny, Nz, kgrid.Nt);

% reconstruct the initial pressure
p_xyz = kspacePlaneRecon(sensor_data_rs, kgrid.dy, kgrid.dz, kgrid.dt, ...
    medium.sound_speed, 'DataOrder', 'yzt', 'PosCond', true, 'Plot', true);


%% k-wave reconstruction
% define a second k-space grid using the dimensions of p_xyz
[Nx_recon, Ny_recon, Nz_recon] = size(p_xyz);
kgrid_recon = kWaveGrid(Nx_recon, kgrid.dt * medium.sound_speed, Ny_recon, dy, Nz_recon, dz);

% resample p_xyz to be the same size as source.p0
p_xyz_rs = interp3(kgrid_recon.y, kgrid_recon.x - min(kgrid_recon.x(:)),kgrid_recon.z, p_xyz, kgrid.y, kgrid.x - min(kgrid.x(:)), kgrid.z);


%% K-Wave Regression
H_k_wave_recon = permute(p_xyz_rs, [2,1,3])./gruneisen_coef;
cutoff = 30;
cutoff_end = 90;
analyze_H(H_k_wave_recon,X,Y,Z,mua,mus,g,cutoff,cutoff_end, "K-Wave")


%% Analyzing straight from sensor
cutoff = 30;
cutoff_end = 90;
x_arr = -dx*(Nx-1)/2:dx:dx*(Nx-1)/2;
H_from_sensor = squeeze(sensor_data_rs(11,11,:));
x_arr = X(1,:,1);
lm = fitlm(x_arr(cutoff:cutoff_end), log(H_from_sensor(cutoff:cutoff_end)));
intercept = lm.Coefficients.Estimate(1);
mu_reg = -lm.Coefficients.Estimate(2);
simple_reg_energy = exp(-mu_reg*x_arr+ intercept);
fprintf("mu from simple linear regression %d\n", mu_reg)
figure;
plot(x_arr(cutoff:cutoff_end), H_from_sensor(cutoff:cutoff_end), "DisplayName","Real Data")
hold on
plot(x_arr(cutoff:cutoff_end), simple_reg_energy(cutoff:cutoff_end), "DisplayName","Linear Regression")
title("Absorbtion for Different Regressions")
xlabel("Depth [mm]")
legend()
title("From Sensor")


figure;
plot(x_arr, H_from_sensor, "DisplayName","data")
hold on
plot(x_arr, log(H_from_sensor), "DisplayName","log")
hold on
line([x_arr(1)+(cutoff-1)*dx x_arr(1)+(cutoff-1)*dx], ylim, 'Color', 'r', 'LineStyle', '--', 'LineWidth', 2);
hold on
line([x_arr(1)+(cutoff_end-1)*dx x_arr(1)+(cutoff_end-1)*dx], ylim, 'Color', 'r', 'LineStyle', '--', 'LineWidth', 2);
xlabel("Depth [mm]")
legend()
title(["Absorbtion in middle ", name])

%% smoothed regression
H = load("H_21X21X201_Var0.1.mat").H;
H_smooth = smooth(H);
analyze_H(H_smooth,X,Y,Z,mua,mus,g,cutoff,cutoff_end, "Smooth")
fprintf("expected mu  %d\n", mu_expected_value)
fprintf("std mu  %d\n", sqrt(mu_variance))

H_mean = squeeze(mean(H, [1,3]));
H_mean = repmat(H_mean,21,1,21);
analyze_H(H_mean,X,Y,Z,mua,mus,g,cutoff,cutoff_end, "Mean")
fprintf("expected mu  %d\n", mu_expected_value)
fprintf("std mu  %d\n", sqrt(mu_variance))