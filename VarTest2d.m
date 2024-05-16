
%% simulation parameters
dx = 0.01;
dy = 0.5;
dz = 0.5;
Nx = 201;
Ny = 21;
x_arr = -dx*(Nx-1)/2:dx:dx*(Nx-1)/2;
y_arr = -dy*(Ny-1)/2:dy:dy*(Ny-1)/2;
mua=0.3;
mus=50;
g = 0.9;
sigma = 0;
gruneisen_coef = 0.5;

%% create grid
vmcmesh = createGridMesh(x_arr, y_arr); % function provided by ValoMC
vmcmedium = createMedium(vmcmesh);
[X,Y] = meshgrid(x_arr,y_arr); % Matlab function

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
lightsource = findBoundaries(vmcmesh, 'direction', [-1.5 0 ], [-0.5 0 ], 4.5);
vmcboundary.lightsource(lightsource) = {'direct'};
options.photon_count = 1e9;
solution = ValoMC(vmcmesh, vmcmedium, vmcboundary, options);

%% Regression

cutoff = 30;
cutoff_end = 90;
H = vmcmedium.absorption_coefficient .* solution.grid_fluence*1e6;
% H = load(['VarResults\H_var_', num2str(sigma)]).H;
analyze_H_2d(H,X,Y,mua,mus,g,cutoff,cutoff_end, "ValoMC")
fprintf("expected mu  %d\n", mu_expected_value)
fprintf("std mu  %d\n", sqrt(mu_variance))
% save(['VarResults\H_var_', num2str(sigma)], 'H')

%% k-wave simulation

H_k_wave = permute(H, [2,1]);

% p0
kgrid = kWaveGrid(Nx, dx*1e-3, Ny, dy*1e-3);
medium.sound_speed = 1500;    % [m/s]
medium.density = 1000;        % [kg/m^3]
source.p0 = gruneisen_coef.*H_k_wave;

% sensor
sensor.mask = zeros(Nx, Ny);
sensor.mask(1,:) = 1;
% input arguments
input_args = {'PlotLayout', true, 'PlotPML', false, ...
    'DataCast', 'gpuArray-single', 'CartInterp', 'nearest', 'PMLInside', false, 'DataRecast', true};

% run the simulation
sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});

% reshape sensor data to y, z, t
sensor_data_rs = reshape(sensor_data, Ny, kgrid.Nt);

%%
% reconstruct the initial pressure
p_xy = kspaceLineRecon(sensor_data_rs, kgrid.dy, kgrid.dt, ...
    medium.sound_speed, 'DataOrder', 'yt', 'PosCond', true, 'Plot', true);



% define a second k-space grid using the dimensions of p_xyz
[Nx_recon, Ny_recon, Nz_recon] = size(p_xy);
kgrid_recon = kWaveGrid(Nx_recon, kgrid.dt * medium.sound_speed, Ny_recon, dy, Nz_recon, dz);

% resample p_xyz to be the same size as source.p0
p_xy_rs = interp2(kgrid_recon.y, kgrid_recon.x - min(kgrid_recon.x(:)), p_xy, kgrid.y, kgrid.x - min(kgrid.x(:)));


%% K-Wave Regression
H_k_wave_recon = permute(p_xy_rs, [2,1])./gruneisen_coef;
cutoff = 30;
cutoff_end = 90;
analyze_H_2d(H_k_wave_recon,X,Y,mua,mus,g,cutoff,cutoff_end, "K-Wave")
save(['VarResults\H_recon_var_', num2str(sigma)], 'H_k_wave_recon')
save(['VarResults\Sensor_Data_var_', num2str(sigma)], 'sensor_data')

