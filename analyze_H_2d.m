function analyze_H_2d(H,X,Y,mua,mus,g,cutoff,cutoff_end,name)
%% Find mu_eff theoretical
H_mid = H(Y==0);
x_arr = X(1,:);
y_arr = Y(:,1);
dx = abs(x_arr(2)-x_arr(1));
mu_expected = sqrt(2*mua*(mua+(1-g)*mus));

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

%% simple linear regression 2D
cutoff_x = x_arr(1)+(cutoff-1)*dx;
cutoff_end_x = x_arr(1)+(cutoff_end-1)*dx;
radius = 3;
relevant_idx = (Y<=radius & Y>=-radius ...
    & X<=cutoff_end_x & X>=cutoff_x);
x_2d_arr = X(relevant_idx)';
H_2d = H(relevant_idx);
lm_2d = fitlm(x_2d_arr, log(H_2d));
intercept_2d = lm_2d.Coefficients.Estimate(1);
mu_reg_2d = -lm_2d.Coefficients.Estimate(2);
simple_reg_2d_energy = exp(-mu_reg_2d*x_arr+ intercept_2d);

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
fprintf("mu from simple 2d linear regression %d\n", mu_reg_2d)



figure;
plot(x_arr(cutoff:cutoff_end), H_mid(cutoff:cutoff_end), "DisplayName","Real Data")
hold on
plot(x_arr(cutoff:cutoff_end), full_reg_energy(cutoff:cutoff_end), "DisplayName","Full Regression")
hold on
plot(x_arr(cutoff:cutoff_end), single_reg_energy(cutoff:cutoff_end), "DisplayName","Single Exp Regression")
hold on
plot(x_arr(cutoff:cutoff_end), simple_reg_energy(cutoff:cutoff_end), "DisplayName","Linear Regression")
hold on
plot(x_arr(cutoff:cutoff_end), simple_reg_2d_energy(cutoff:cutoff_end), "DisplayName","Linear 2D Regression")
hold on
plot(x_arr(cutoff:cutoff_end), full1d_reg_energy(cutoff:cutoff_end), "DisplayName","Full 1D Regression")
title("Absorbtion for Different Regressions")
xlabel("Depth [mm]")
legend()
title(["Comparison of Absorbtion for Different Reconstructions ", name])


figure;
[ax, h1, h2] = plotyy(x_arr, H_mid, x_arr, log(H_mid));
set(h1, "DisplayName","data");
set(h2, "DisplayName","log");
hold on
line([x_arr(1)+(cutoff-1)*dx x_arr(1)+(cutoff-1)*dx], ylim, 'Color', 'r', 'LineStyle', '--', 'LineWidth', 2);
hold on
line([x_arr(1)+(cutoff_end-1)*dx x_arr(1)+(cutoff_end-1)*dx], ylim, 'Color', 'r', 'LineStyle', '--', 'LineWidth', 2);
ylabel(ax(1), 'Signal'); % Label for the left y-axis
ylabel(ax(2), 'Log Signal'); % Label for the right y-axis
xlabel("Depth [mm]")
legend;
title(["Absorbtion in middle ", name])

figure;
image(x_arr, y_arr, H);
xlabel('x [mm]');
ylabel('y [mm]');
colorbar;
title(["Absorbtion ", name])


end