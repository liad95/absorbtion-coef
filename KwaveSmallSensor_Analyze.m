function KwaveSmallSensor2d_Analyze(Radius)

%% simulation parameters
dx = 0.05;
dy = 0.5;
Nx = 201;
Ny = 101;
x_arr = -dx*(Nx-1)/2:dx:dx*(Nx-1)/2;
y_arr = -dy*(Ny-1)/2:dy:dy*(Ny-1)/2;
mua=0.01;
mus=10;
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
negative_idx = absorption_matrix<0;
absorption_matrix(negative_idx)=0;


mu_3d_matrix = sqrt(3.*absorption_matrix.*(absorption_matrix+(1-g).*scattering_matrix));
mu_expected_value_3d = mean(mu_3d_matrix(:));
mu_variance_3d = var(mu_3d_matrix(:));


vmcmedium.scattering_coefficient = scattering_matrix;
vmcmedium.absorption_coefficient = absorption_matrix;
vmcmedium.scattering_anisotropy = g;        
vmcmedium.refractive_index = 1;


%% Regression

cutoff = 30;
cutoff_end = 90;
H = load(['temp\H_rad_', num2str(Radius)]).H;
analyze_H_2d_print(H,X,Y,mua,mus,g,cutoff,cutoff_end, "ValoMC", Radius)

fileID = fopen(join(['H_rad_', "ValoMC", '_', num2str(Radius), '.txt']), 'a');
fprintf(fileID, "expected mu 3d %d\n", mu_expected_value_3d)
fprintf(fileID, "std mu 3d %d\n", sqrt(mu_variance_3d))
fclose(fileID);

% H_k_wave_recon = load(['temp\H_recon_rad_', num2str(sensor_size)]).H_k_wave_recon;
% cutoff = 40;
% cutoff_end = 90;
% analyze_H_2d_print(H_k_wave_recon,X,Y,mua,mus,g,cutoff,cutoff_end, "K-Wave", sensor_size)


end
