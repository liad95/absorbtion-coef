%% simulation parameters
dx = 0.01;
dy = 0.5;
dz = 0.5;
Nx = 200;
Ny = 20;
Nz = 20;
x_arr = -dx*Nx/2:dx:dx*(Nx/2-1);
y_arr = -dy*Ny/2:dy:dy*(Ny/2-1);
z_arr = -dz*Nz/2:dz:dz*(Nz/2-1);
mua=0.3;
mus=50;
g = 0.9;
sigma = 1;
gruneisen_coef = 0.5;

%% create grid
vmcmesh = createGridMesh(x_arr, y_arr, z_arr); % function provided by ValoMC
vmcmedium = createMedium(vmcmesh);
[X,Y,Z] = meshgrid(x_arr,y_arr,z_arr); % Matlab function

vmcmedium.scattering_coefficient = mus + mus.*sigma.*randn(size(X));
vmcmedium.absorption_coefficient = mua + mus.*sigma.*randn(size(X));
vmcmedium.scattering_anisotropy = g;        
vmcmedium.refractive_index = 1;



%% create source
clear vmcboundary;
vmcboundary = createBoundary(vmcmesh, vmcmedium);
lightsource = findBoundaries(vmcmesh, 'direction', [-1.5 0 0], [-0.5 0 0], 4.5);
vmcboundary.lightsource(lightsource) = {'direct'};
options.photon_count = 1e7;
solution = ValoMC(vmcmesh, vmcmedium, vmcboundary, options);

%% Regression
cutoff = 50;
cutoff_end = 110;
H = vmcmedium.absorption_coefficient .* solution.grid_fluence*1e6;
analyze_H(H,X,Y,Z,mua,mus,g,cutoff,cutoff_end)


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
cutoff = 50;
cutoff_end = 110;
analyze_H(H_k_wave_recon,X,Y,Z,mua,mus,g,cutoff,cutoff_end)