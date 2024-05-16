function analyze_H_new(H,X,Y,Z,mua,mus,g,cutoff,cutoff_end,name)
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
ft = fittype('exp(-mu_eff*x+a);');
[simple_3d_parameters,full_gof] = fit(x_arr(cutoff:cutoff_end)',H_3d(cutoff:cutoff_end),ft);
simple_reg_3d_energy = exp(-simple_3d_parameters.mu_eff*x_arr+ simple_3d_parameters.a);

%% simple linear regression
ft = fittype('exp(-mu_eff*x+a);');
[simple_parameters,full_gof] = fit(x_arr(cutoff:cutoff_end)',H_mid(cutoff:cutoff_end),ft);
simple_reg_energy = exp(-simple_parameters.mu_eff*x_arr+ simple_parameters.a);

%% regression comparisons



fprintf("mu from ground truth %d\n", mu_expected)
fprintf("mu from simple linear regression %d\n", simple_parameters.mu_eff)
fprintf("mu from simple 3d linear regression %d\n", simple_3d_parameters.mu_eff)
fprintf("mu from single exp regression %d\n", single_parameters.mu_eff)



figure;
plot(x_arr(cutoff:cutoff_end), H_mid(cutoff:cutoff_end), "DisplayName","Real Data")
hold on
plot(x_arr(cutoff:cutoff_end), simple_reg_energy(cutoff:cutoff_end), "DisplayName","Linear Regression")
hold on
plot(x_arr(cutoff:cutoff_end), simple_reg_3d_energy(cutoff:cutoff_end), "DisplayName","Linear 3D Regression")
hold on
plot(x_arr(cutoff:cutoff_end), single_reg_energy(cutoff:cutoff_end), "DisplayName","Single Exp Regression")
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
slice(X, Y, Z, H, 0, 0, 0);
xlabel('x [mm]');
ylabel('y [mm]');
zlabel('z [mm]');
view(315,25);
hold
title(["Absorbtion ", name])

snapnow;


end