function analyze_H_new(H,X,Y,Z,mua,mus,g,cutoff,cutoff_end,name, sigma)
%% Find mu_eff theoretical
H_mid = H(Z==0 & Y==0);
x_arr = X(1,:,1);
dx = abs(x_arr(2)-x_arr(1));
mu_expected = sqrt(3*mua*(mua+(1-g)*mus));


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

fileID = fopen(join(['mus_', name, '_', num2str(sigma), '.txt']), 'w');
fprintf(fileID, "mu from ground truth %d\n", mu_expected)
fprintf(fileID, "mu from simple linear regression %d\n", mu_reg)
fprintf(fileID, "mu from simple 3d linear regression %d\n", mu_reg_3d)
fprintf(fileID, "mu from single exp regression %d\n", single_parameters.mu_eff)

fprintf("mu from ground truth %d\n", mu_expected)
fprintf("mu from simple linear regression %d\n", mu_reg)
fprintf("mu from simple 3d linear regression %d\n", mu_reg_3d)
fprintf("mu from single exp regression %d\n", single_parameters.mu_eff)

fclose(fileID);



figure;
plot(x_arr(cutoff:cutoff_end), H_mid(cutoff:cutoff_end), "DisplayName","Real Data")
hold on
plot(x_arr(cutoff:cutoff_end), simple_reg_energy(cutoff:cutoff_end), "DisplayName","Linear Regression")
hold on
plot(x_arr(cutoff:cutoff_end), simple_reg_3d_energy(cutoff:cutoff_end), "DisplayName","Linear 3D Regression")
hold on
plot(x_arr(cutoff:cutoff_end), single_reg_energy(cutoff:cutoff_end), "DisplayName","Single Exp Regression")
xlabel("Depth [mm]")
legend()
title(join(["Absorbtion for Different Regressions" , name]))
filename = join(['Reconstructions_', name, '_', num2str(sigma), '.png'], '');
saveas(gcf, filename);


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
title(join(["Absorbtion in Center ", name]))
filename = join(['CenterAbsorbtion_', name, '_', num2str(sigma), '.png'], '');
saveas(gcf, filename);

figure;
slice(X, Y, Z, H, 0, 0, 0);
xlabel('x [mm]');
ylabel('y [mm]');
zlabel('z [mm]');
view(315,25);
hold
title(join(["Absorbtion ", name]))
filename = join(['Absorbtion_', name, '_', num2str(sigma), '.png'], '');
saveas(gcf, filename);

snapnow;


end