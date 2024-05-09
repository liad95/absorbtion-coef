function analyze_H(H,X,Y,Z,mua,mus,g,cutoff,cutoff_end,name)
%% Find mu_eff theoretical
H_mid = H(Z==0 & Y==0);
x_arr = X(1,:,1);
dx = abs(x_arr(2)-x_arr(1));
mu_expected = sqrt(3*mua*(mua+(1-g)*mus));

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
relevant_idx = (Z<=radius & Y<=radius & Z>=-radius & Y>=-radius ...
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
title(["Comparison of Absorbtion for Different Reconstructions ", name])


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
title(["Absorbtion in middle ", name])

figure;
slice(X, Y, Z, H, 0, 0, 0);
xlabel('x [mm]');
ylabel('y [mm]');
zlabel('z [mm]');
view(315,25);
hold
title(["Absorbtion ", name])

snapnow;


end