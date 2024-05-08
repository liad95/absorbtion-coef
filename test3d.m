function[]=test3d(mus, mua)
%% simulation parameters
dx = 0.1;
dy = 0.1;
dz = 0.1;
Nx = 41;
Ny = 41;
Nz = 41;
x_arr = 0:dx:(dx*Nx);
y_arr = 0:dy:(dy*Ny);
z_arr = 0:dz:(dz*Nz);
% mua=0.2;
% mus=1.0;
g = 0.9;
folder = "3D";

%% create grid
vmcmesh = createGridMesh(x_arr, y_arr, z_arr); % function provided by ValoMC
vmcmedium = createMedium(vmcmesh);
[X,Y,Z] = meshgrid(x_arr,y_arr,z_arr); % Matlab function

vmcmedium.scattering_coefficient = mus;
vmcmedium.absorption_coefficient = mua*ones(size(X));  %refractive index is now a three dimensional array
vmcmedium.scattering_anisotropy = g;        
vmcmedium.refractive_index = 1;

vmcboundary = createBoundary(vmcmesh, vmcmedium);

%% create source

lightsource = findBoundaries(vmcmesh, 'direction', [2 2 0], [2 2 10], 1);
vmcboundary.lightsource(lightsource) = {'direct'};
solution = ValoMC(vmcmesh, vmcmedium, vmcboundary);

%% show absorbtion map

H = vmcmedium.absorption_coefficient .* solution.grid_fluence;
figure("Visible","off")
slice(X, Y, Z, H, 2, 2, 0);
xlabel('x [mm]');
ylabel('y [mm]');
zlabel('z [mm]');
view(125,25);
hold

snapnow;



%% Find mu_eff theoretical


H_mid = flip(H(X==2 & Y==2));


mu_expected = sqrt(3*mua*(mua+0.1*mus));
D_expected = 1/(3*(mua+0.1*mus));
cutoff = 20;

%% full regression
ft = fittype('fullfitfunc(x, mu_eff, D, a)');
[full_parameters,gof] = fit(z_arr(cutoff:end)',H_mid(cutoff:end),ft, "Lower", [0,0,0]);
full_reg_energy = fullfitfunc(z_arr, full_parameters.mu_eff, full_parameters.D, full_parameters.a);

%% single exponent regression

ft = fittype('singleexpfitfunc(x, mu_eff, l, a)');
[single_parameters,full_gof] = fit(z_arr(cutoff:end)',H_mid(cutoff:end),ft,"Lower", [0,0]);
single_reg_energy = singleexpfitfunc(z_arr, single_parameters.mu_eff, single_parameters.l, single_parameters.a);

%% simple linear regression
lm = fitlm(z_arr(cutoff:end), log(H_mid(cutoff:end)));
intercept = lm.Coefficients.Estimate(1);
mu_reg = -lm.Coefficients.Estimate(2);
simple_reg_energy = exp(-mu_reg*z_arr+ intercept);

%% regression comparisons
fprintf("mu from simple linear regression %d\n", mu_reg)
fprintf("mu from full regression %d\n", full_parameters.mu_eff)
fprintf("mu from single exp regression %d\n", single_parameters.mu_eff)
fprintf("mu from ground truth %d\n", mu_expected)


figure_folder = "Figures";
figure_name = "AbsorbReg-" + "mua" + string(mua) + "mus" + string(mus) + ".fig";
full_file_path = fullfile(folder, figure_folder, figure_name);

figure("Visible","off")
plot(z_arr(cutoff:end), H_mid(cutoff:end), "DisplayName","Real Data")
hold on
plot(z_arr(cutoff:end), full_reg_energy(cutoff:end), "DisplayName","Full Regression")
hold on
plot(z_arr(cutoff:end), single_reg_energy(cutoff:end), "DisplayName","Single Exp Regression")
hold on
plot(z_arr(cutoff:end), simple_reg_energy(cutoff:end), "DisplayName","Linear Regression")
title("Absorbtion for Different Regressions")
xlabel("Depth [mm]")
legend()
savefig(full_file_path)


figure_name = "Cutoff-" + "mua" + string(mua) + "mus" + string(mus) + ".fig";
full_file_path = fullfile(folder, figure_folder, figure_name);
figure("Visible","off")
plot(z_arr, H_mid, "DisplayName","data")
hold on
plot(z_arr, log(H_mid), "DisplayName","log")
hold on
line([cutoff*dz cutoff*dz], ylim, 'Color', 'r', 'LineStyle', '--', 'LineWidth', 2);
xlabel("Depth [mm]")
legend()
savefig(full_file_path)

%% Algorithm parameters
alg_cutoff = cutoff;
alg_cutoff_end = 41;
%% Finding mu effective from point
l = D_expected*3;
mu_found_full = mu_recon(H_mid, z_arr', mu_expected, D_expected);
mu_found = -log(abs(z_arr'-l).*H_mid*(4*pi*D_expected))./abs(z_arr'-l); % single exponent

mu_found_sd0 = -log(abs(z_arr').*H_mid*(4*pi*D_expected))./abs(z_arr'); % skin depth = 0




%% loop
mu_eff_guess = 1;
D_alg = mua^2/mu_eff_guess;
for i = 1:500
%     mu_found_no_info = -log(abs(z_arr).*H_mid*(4*pi*D_alg))./abs(z_arr);
    mu_found_no_info = mu_recon(H_mid, z_arr', mu_eff_guess, D_alg);
    mu_eff_guess = mean(mu_found_no_info(alg_cutoff:alg_cutoff_end));
    D_alg = mua^2/mu_eff_guess;
end


% l_alg = D_alg*3;
% mu_found_alg = -log(abs(z_arr-l_alg).*H_mid*(4*pi*D_alg))./abs(z_arr-l_alg);
% figure("Visible","off")
% plot(z_arr, mu_found, "DisplayName", "semi-complete")

%% algorithm comparisons

fprintf("algorithm mu %d\n", mu_eff_guess)
fprintf("full expressiom mu %d\n", mean(mu_found_full(alg_cutoff:alg_cutoff_end)))
fprintf("single exponent mu %d\n", mean(mu_found(alg_cutoff:alg_cutoff_end)))
fprintf("skin depth 0 mu %d\n", mean(mu_found_sd0(alg_cutoff:alg_cutoff_end)))

figure_name = "muReconstruction-" + "mua" + string(mua) + "mus" + string(mus) + ".fig";
full_file_path = fullfile(folder, figure_folder, figure_name);
figure("Visible","off")
plot(z_arr, mu_found_full, "DisplayName", "full")
hold on
plot(z_arr, mu_found, "DisplayName", "single exponent")
hold on
plot(z_arr, mu_found_sd0, "DisplayName", "skin depth = 0")
legend()
title("mu effective for different reconstructions")
xlabel("z [cm]")
savefig(full_file_path)


%% write to file



result_folder = "Results";
file_name = "mua" + string(mua) + "mus" + string(mus) + ".csv";
full_file_path = fullfile(folder, result_folder, file_name);
parameters_matrix = [mu_expected, mu_reg, full_parameters.mu_eff, single_parameters.mu_eff, mu_eff_guess, mean(mu_found_full(alg_cutoff:alg_cutoff_end)), mean(mu_found(alg_cutoff:alg_cutoff_end)), mean(mu_found_sd0(alg_cutoff:alg_cutoff_end))];
csvwrite(full_file_path, parameters_matrix);


