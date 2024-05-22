error = linspace(0,1,1000);
mu_eff = 0.176;
l = 1;
R_1 = sqrt((log(error)/mu_eff).^2 - 2*(log(error)/mu_eff).*(1-l));
R_5 = sqrt((log(error)/mu_eff).^2 - 2*(log(error)/mu_eff).*(5-l));
R_10 = sqrt((log(error)/mu_eff).^2 - 2*(log(error)/mu_eff).*(10-l));
figure;
plot(R_1, error, "DisplayName", "x=1[mm]")
hold on
plot(R_5, error, "DisplayName", "x=5[mm]")
hold on
plot(R_10, error, "DisplayName", "x=10[mm]")
ylabel("\epsilon - normalized error")
xlabel("R [mm]")
legend()
set(gca, 'FontSize', 16, 'OuterPosition', [0 0 1 1])
set(gcf, 'PaperPositionMode', 'auto')
filename = join(['Illumnation Radius vs Normalized Error.png'], '');
saveas(gcf, filename);
