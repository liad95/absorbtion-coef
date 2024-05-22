function analyze_H_2d_print(H,X,Y,mua,mus,g,cutoff,cutoff_end,name, sigma)
%% Find mu_eff theoretical
H_mid = H(Y==0);
x_arr = X(1,:);
y_arr = Y(:,1);
dx = abs(x_arr(2)-x_arr(1));
mu_expected_3d = sqrt(3*mua*(mua+(1-g)*mus));
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

fileID = fopen(join(['mus_', name, '_', num2str(sigma), '.txt']), 'w');
fprintf(fileID, "mu from ground truth 3d %d\n", mu_expected_3d)
fprintf(fileID, "mu from simple linear regression %d\n", mu_reg)
fprintf(fileID, "mu from simple 3d linear regression %d\n", mu_expected_3d)
fclose(fileID);


fprintf("mu from ground truth 3d %d\n", mu_expected_3d)
fprintf("mu from simple linear regression %d\n", mu_reg)
fprintf("mu from simple 2d linear regression %d\n", mu_reg_2d)



figure;
plot(x_arr(cutoff:cutoff_end), H_mid(cutoff:cutoff_end), "DisplayName","Real Data", "LineWidth",2, "LineStyle",":")
hold on
plot(x_arr(cutoff:cutoff_end), simple_reg_energy(cutoff:cutoff_end), "DisplayName","Center Line Fit", "LineWidth",2, "LineStyle","-")
hold on
plot(x_arr(cutoff:cutoff_end), simple_reg_2d_energy(cutoff:cutoff_end), "DisplayName","Cylinder Fit", "LineWidth",2, "LineStyle","-")
xlabel("Depth [mm]")
ylabel("Signal")
legend()

filename = join(['Reconstructions_', name, '_', num2str(sigma), '.png'], '');
set(gca, 'FontSize', 16, 'OuterPosition', [0 0 1 1])
set(gcf, 'PaperPositionMode', 'auto')
saveas(gcf, filename);


figure;

[ax, h1, h2] = plotyy(x_arr, H_mid, x_arr, log(H_mid));
set(h1, "DisplayName","Data", "LineWidth",2, "LineStyle",":")
set(h2, "DisplayName","Log", "LineWidth",2, "LineStyle",":")
y_edges = ylim;
hold on
patch([cutoff_x cutoff_end_x cutoff_end_x cutoff_x], [y_edges(1) y_edges(1) y_edges(2) y_edges(2)], 'green', 'FaceAlpha', 0.05, 'EdgeColor', 'green', 'DisplayName', 'Regression Area');
ylabel(ax(1), 'Signal'); % Label for the left y-axis
ylabel(ax(2), 'Log Signal'); % Label for the right y-axis
xlabel("Depth [mm]")
legend;
filename = join(['CenterAbsorbtion_', name, '_', num2str(sigma), '.png'], '');
set(gca, 'FontSize', 16, 'OuterPosition', [0 0 1 1])
set(gcf, 'PaperPositionMode', 'auto')
saveas(gcf, filename);

figure;
image(x_arr, y_arr, H);
xlabel('x [mm]');
ylabel('y [mm]');
colorbar;
filename = join(['Absorbtion_', name, '_', num2str(sigma), '.png'], '');
set(gca, 'FontSize', 16, 'OuterPosition', [0 0 1 1])
set(gcf, 'PaperPositionMode', 'auto')

saveas(gcf, filename);



end