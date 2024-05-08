%% simulation parameters
dx = 0.01;
dy = 0.02;
dz = 0.02;
Nx = 200;
Ny = 20;
Nz = 20;
x_arr = -dx*Nx/2:dx:dx*(Nx/2-1);
y_arr = -dy*Ny/2:dy:dy*(Ny/2-1);
z_arr = -dz*Nz/2:dz:dz*(Nz/2-1);
mua=0.3;
mus=50;
g = 0.9;
gruneisen_coef = 0.5;

%% create grid
vmcmesh = createGridMesh(x_arr, y_arr, z_arr); % function provided by ValoMC
vmcmedium = createMedium(vmcmesh);
[X,Y,Z] = meshgrid(x_arr,y_arr,z_arr); % Matlab function

vmcmedium.scattering_coefficient = mus;
vmcmedium.absorption_coefficient = mua*ones(size(X));  %refractive index is now a three dimensional array
vmcmedium.scattering_anisotropy = g;        
vmcmedium.refractive_index = 1;



%% create source
clear vmcboundary;
vmcboundary = createBoundary(vmcmesh, vmcmedium);
lightsource = findBoundaries(vmcmesh, 'direction', [-1.5 0 0], [-0.5 0 0], 4.5);
vmcboundary.lightsource(lightsource) = {'direct'};
options.photon_count = 1e6;
solution = ValoMC(vmcmesh, vmcmedium, vmcboundary, options);

%% show absorbtion map

H = vmcmedium.absorption_coefficient .* solution.grid_fluence*1e6;
figure;
slice(X, Y, Z, H, 0, 0, 0);
xlabel('x [mm]');
ylabel('y [mm]');
zlabel('z [mm]');
view(125,25);
hold

snapnow;



%% Find mu_eff theoretical

center = 0;
H_mid = H(Z==center & Y==center);


mu_expected = sqrt(3*mua*(mua+(1-g)*mus));
cutoff = 25;
cutoff_end = 175;

%% full regression 1d
ft = fittype('fullfitfunc(x, mu_eff, D, a)');
[full1d_parameters,gof] = fit(x_arr(cutoff:cutoff_end)',H_mid(cutoff:cutoff_end),ft);
full1d_reg_energy = fullfitfunc_1d(x_arr, full1d_parameters.mu_eff, full1d_parameters.D, full1d_parameters.a);
%% full regression
ft = fittype('fullfitfunc(x, mu_eff, D, a)');
[full_parameters,gof] = fit(x_arr(cutoff:cutoff_end)',H_mid(cutoff:cutoff_end),ft);
full_reg_energy = fullfitfunc(x_arr, full_parameters.mu_eff, full_parameters.D, full_parameters.a);

%% single exponent regression

ft = fittype('singleexpfitfunc(x, mu_eff, l, a)');
[single_parameters,full_gof] = fit(x_arr(cutoff:cutoff_end)',H_mid(cutoff:cutoff_end),ft);
single_reg_energy = singleexpfitfunc(x_arr, single_parameters.mu_eff, single_parameters.l, single_parameters.a);

% Dont understand why we get the mu in negative. Maybe because of the
% absolute value?

%% simple linear regression 3D
cutoff_x = x_arr(1)+(cutoff-1)*dx;
cutoff_end_x = x_arr(1)+(cutoff_end-1)*dx;
radius = 3;
relevant_idx = (Z<=center+radius & Y<=center+radius & Z>=center-radius & Y>=center-radius ...
    & X<=cutoff_end_x & X>=cutoff_x);
x_3d_arr = X(relevant_idx)';
H_3d = H(relevant_idx);
lm_3d = fitlm(x_3d_arr, log(H_3d));
intercept_3d = lm_3d.Coefficients.Estimate(1);
mu_reg_3d = -lm_3d.Coefficients.Estimate(2);
simple_reg_3d_energy = exp(-mu_reg_3d*x_arr+ intercept_3d);

%% simple linear regression
lm = fitlm(x_arr(cutoff:cutoff_end), log(H_mid(cutoff:cutoff_end)));
intercept = lm.Coefficients.Estimate(1);
mu_reg = -lm.Coefficients.Estimate(2);
simple_reg_energy = exp(-mu_reg*x_arr+ intercept);

%% regression comparisons

fprintf("mu from full regression %d\n", full_parameters.mu_eff)
fprintf("mu from full 1d regression %d\n", full1d_parameters.mu_eff)
fprintf("mu from single exp regression %d\n", single_parameters.mu_eff)
fprintf("mu from ground truth %d\n", mu_expected)
fprintf("mu from simple linear regression %d\n", mu_reg)
fprintf("mu from simple 3d linear regression %d\n", mu_reg_3d)



figure;
plot(x_arr(cutoff:cutoff_end), H_mid(cutoff:cutoff_end), "DisplayName","Real Data")
hold on
plot(x_arr(cutoff:cutoff_end), full_reg_energy(cutoff:cutoff_end), "DisplayName","Full Regression")
hold on
plot(x_arr(cutoff:cutoff_end), single_reg_energy(cutoff:cutoff_end), "DisplayName","Single Exp Regression")
hold on
plot(x_arr(cutoff:cutoff_end), simple_reg_energy(cutoff:cutoff_end), "DisplayName","Linear Regression")
hold on
plot(x_arr(cutoff:cutoff_end), simple_reg_3d_energy(cutoff:cutoff_end), "DisplayName","Linear 3D Regression")
hold on
plot(x_arr(cutoff:cutoff_end), full1d_reg_energy(cutoff:cutoff_end), "DisplayName","Full 1D Regression")
title("Absorbtion for Different Regressions")
xlabel("Depth [mm]")
legend()



figure;
plot(x_arr, H_mid, "DisplayName","data")
hold on
plot(x_arr, log(H_mid), "DisplayName","log")
hold on
line([x_arr(1)+(cutoff-1)*dx x_arr(1)+(cutoff-1)*dx], ylim, 'Color', 'r', 'LineStyle', '--', 'LineWidth', 2);
hold on
line([x_arr(1)+(cutoff_end-1)*dx x_arr(1)+(cutoff_end-1)*dx], ylim, 'Color', 'r', 'LineStyle', '--', 'LineWidth', 2);
xlabel("Depth [mm]")
legend()

H = vmcmedium.absorption_coefficient .* solution.grid_fluence*1e6;
figure;
slice(X, Y, Z, H, 0, 0, 0);
xlabel('x [mm]');
ylabel('y [mm]');
zlabel('z [mm]');
view(315,25);
hold

snapnow;


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


%% Preparing Regression
H_k_wave_recon = permute(p_xyz_rs, [2,1,3])./gruneisen_coef;


center = 0;
H_k_wave_recon_mid = H_k_wave_recon(Z==center & Y==center);
cutoff = 25;
cutoff_end = 50;
%% Find mu_eff theoretical
mu_expected = sqrt(3*mua*(mua+(1-g)*mus));
%% full regression 1d
ft = fittype('fullfitfunc(x, mu_eff, D, a)');
[full1d_parameters,gof] = fit(x_arr(cutoff:cutoff_end)',H_k_wave_recon_mid(cutoff:cutoff_end),ft);
full1d_reg_energy = fullfitfunc_1d(x_arr, full1d_parameters.mu_eff, full1d_parameters.D, full1d_parameters.a);
%% full regression
ft = fittype('fullfitfunc(x, mu_eff, D, a)');
[full_parameters,gof] = fit(x_arr(cutoff:cutoff_end)',H_k_wave_recon_mid(cutoff:cutoff_end),ft);
full_reg_energy = fullfitfunc(x_arr, full_parameters.mu_eff, full_parameters.D, full_parameters.a);

%% single exponent regression

ft = fittype('singleexpfitfunc(x, mu_eff, l, a)');
[single_parameters,full_gof] = fit(x_arr(cutoff:cutoff_end)',H_k_wave_recon_mid(cutoff:cutoff_end),ft);
single_reg_energy = singleexpfitfunc(x_arr, single_parameters.mu_eff, single_parameters.l, single_parameters.a);

% Dont understand why we get the mu in negative. Maybe because of the
% absolute value?

%% simple linear regression 3D
cutoff_x = x_arr(1)+(cutoff-1)*dx;
cutoff_end_x = x_arr(1)+(cutoff_end-1)*dx;
radius = 3;
relevant_idx = (Z<=center+radius & Y<=center+radius & Z>=center-radius & Y>=center-radius ...
    & X<=cutoff_end_x & X>=cutoff_x);
x_3d_arr = X(relevant_idx)';
H_k_wave_recon_3d = H_k_wave_recon(relevant_idx);
lm_3d = fitlm(x_3d_arr, log(H_k_wave_recon_3d));
intercept_3d = lm_3d.Coefficients.Estimate(1);
mu_reg_3d = -lm_3d.Coefficients.Estimate(2);
simple_reg_3d_energy = exp(-mu_reg_3d*x_arr+ intercept_3d);

%% simple linear regression
lm = fitlm(x_arr(cutoff:cutoff_end), log(H_k_wave_recon_mid(cutoff:cutoff_end)));
intercept = lm.Coefficients.Estimate(1);
mu_reg = -lm.Coefficients.Estimate(2);
simple_reg_energy = exp(-mu_reg*x_arr+ intercept);

%% regression comparisons

fprintf("mu from full regression %d\n", full_parameters.mu_eff)
fprintf("mu from full 1d regression %d\n", full1d_parameters.mu_eff)
fprintf("mu from single exp regression %d\n", single_parameters.mu_eff)
fprintf("mu from ground truth %d\n", mu_expected)
fprintf("mu from simple linear regression %d\n", mu_reg)
fprintf("mu from simple 3d linear regression %d\n", mu_reg_3d)



figure;
plot(x_arr(cutoff:cutoff_end), H_k_wave_recon_mid(cutoff:cutoff_end), "DisplayName","Real Data")
hold on
plot(x_arr(cutoff:cutoff_end), full_reg_energy(cutoff:cutoff_end), "DisplayName","Full Regression")
hold on
plot(x_arr(cutoff:cutoff_end), single_reg_energy(cutoff:cutoff_end), "DisplayName","Single Exp Regression")
hold on
plot(x_arr(cutoff:cutoff_end), simple_reg_energy(cutoff:cutoff_end), "DisplayName","Linear Regression")
hold on
plot(x_arr(cutoff:cutoff_end), simple_reg_3d_energy(cutoff:cutoff_end), "DisplayName","Linear 3D Regression")
hold on
plot(x_arr(cutoff:cutoff_end), full1d_reg_energy(cutoff:cutoff_end), "DisplayName","Full 1D Regression")
title("Absorbtion for Different Regressions")
xlabel("Depth [mm]")
legend()



figure;
plot(x_arr, H_k_wave_recon_mid, "DisplayName","data")
hold on
plot(x_arr, log(H_k_wave_recon_mid), "DisplayName","log")
hold on
line([x_arr(1)+(cutoff-1)*dx x_arr(1)+(cutoff-1)*dx], ylim, 'Color', 'r', 'LineStyle', '--', 'LineWidth', 2);
hold on
line([x_arr(1)+(cutoff_end-1)*dx x_arr(1)+(cutoff_end-1)*dx], ylim, 'Color', 'r', 'LineStyle', '--', 'LineWidth', 2);
xlabel("Depth [mm]")
legend()

figure;
slice(X, Y, Z, H_k_wave_recon, 0, 0, 0);
xlabel('x [mm]');
ylabel('y [mm]');
zlabel('z [mm]');
view(315,25);
hold

snapnow;

%%
figure;
plot(x_arr(cutoff:cutoff_end), log(H_k_wave_recon_mid(cutoff:cutoff_end)), "DisplayName","Real Data")
hold on
plot(x_arr(cutoff:cutoff_end), log(simple_reg_energy(cutoff:cutoff_end)), "DisplayName","Linear Regression")
legend;

%%
figure;
plot(x_arr, H_k_wave_recon_mid, 'DisplayName', 'Reconstructed H')
hold on
plot(x_arr, H_mid, 'DisplayName', 'H')
legend;


%%
H = load('H_200_deep.mat');
H = H.H;
analyze_H(H,X,Y,Z,mua,mus,g,cutoff,cutoff_end)