%% simulation parameters
dx = 0.01;
dy = 0.01;
Nx = 400;
Ny = 400;

x_arr = 0:dx:(dx*Nx);
y_arr = 0:dy:(dy*Ny);
mua=0.2;
mus=1.0;
g = 0.9;

%% create grid
vmcmesh = createGridMesh(x_arr, y_arr); % function provided by ValoMC
vmcmedium = createMedium(vmcmesh);
[X,Y] = meshgrid(x_arr,y_arr); % MATLAB function


vmcmedium.scattering_coefficient = mus;
vmcmedium.absorption_coefficient = mua*ones(size(X));  %refractive index is now a three dimensional array
vmcmedium.scattering_anisotropy = g;        
vmcmedium.refractive_index = 1;

vmcboundary = createBoundary(vmcmesh, vmcmedium);

%% create source

lightsource = findBoundaries(vmcmesh, 'direction', [2 0], [2 -5], 0.1);
vmcboundary.lightsource(lightsource) = {'direct'};

solution = ValoMC(vmcmesh, vmcmedium, vmcboundary);

%% show absorbtion map

figure
H = vmcmedium.absorption_coefficient .* solution.grid_fluence;
imagesc(x_arr, y_arr, (H), [min((H(:))) ...
        max((H(:)))]);

xlabel('[mm]');
ylabel('[mm]');
c = colorbar;                             % create a colorbar
c.Label.String = 'Fluence [J/mm^2]';






%% Find mu_eff theoretical


H_mid = H(X==2);


mu_expected = sqrt(3*mua*(mua+0.1*mus));
D_expected = 1/(3*(mua+0.1*mus));
cutoff = 200;

%% full regression 2D
ft = fittype('fullfitfunc2D(x, mu_eff, D, a)');
[full2d_parameters,gof] = fit(y_arr(cutoff:end)',H_mid(cutoff:end),ft, "Lower", [0,0,0]);
full2d_reg_energy = fullfitfunc2D(y_arr, full2d_parameters.mu_eff, full2d_parameters.D, full2d_parameters.a);

%% full regression
ft = fittype('fullfitfunc(x, mu_eff, D, a)');
[full_parameters,gof] = fit(y_arr(cutoff:end)',H_mid(cutoff:end),ft, "Lower", [0,0,0]);
full_reg_energy = fullfitfunc(y_arr, full_parameters.mu_eff, full_parameters.D, full_parameters.a);

%% single exponent regression

ft = fittype('singleexpfitfunc(x, mu_eff, l, a)');
[single_parameters,full_gof] = fit(y_arr(cutoff:end)',H_mid(cutoff:end),ft,"Lower", [0,0]);
single_reg_energy = singleexpfitfunc(y_arr, single_parameters.mu_eff, single_parameters.l, single_parameters.a);

%% simple linear regression
lm = fitlm(y_arr(cutoff:end), log(H_mid(cutoff:end)));
intercept = lm.Coefficients.Estimate(1);
mu_reg = -lm.Coefficients.Estimate(2);
simple_reg_energy = exp(-mu_reg*y_arr+ intercept);

%% regression comparisons
fprintf("mu from simple linear regression %d\n", mu_reg)
fprintf("mu from full regression %d\n", full_parameters.mu_eff)
fprintf("mu from full regression for 2D %d\n", full2d_parameters.mu_eff)
fprintf("mu from single exp regression %d\n", single_parameters.mu_eff)
fprintf("mu from ground truth %d\n", mu_expected)




figure;
plot(y_arr(cutoff:end), H_mid(cutoff:end), "DisplayName","Real Data")
hold on
plot(y_arr(cutoff:end), full_reg_energy(cutoff:end), "DisplayName","Full Regression")
hold on
plot(y_arr(cutoff:end), full2d_reg_energy(cutoff:end), "DisplayName","Full 2D Regression")
hold on
plot(y_arr(cutoff:end), single_reg_energy(cutoff:end), "DisplayName","Single Exp Regression")
hold on
plot(y_arr(cutoff:end), simple_reg_energy(cutoff:end), "DisplayName","Linear Regression")
title("Absorbtion for Different Regressions")
xlabel("Depth [mm]")
legend()



figure;
plot(y_arr, H_mid, "DisplayName","data")
hold on
plot(y_arr, log(H_mid), "DisplayName","log")
hold on
line([cutoff*dy cutoff*dy], ylim, 'Color', 'r', 'LineStyle', '--', 'LineWidth', 2);
xlabel("Depth [mm]")
legend()




