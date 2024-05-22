%% Sensor Size
sensor_sizes = [5,11,21,31,51,71,91];
mu_ideal = 0.174;

mu_linear_R = 0.1.*[1.79, 1.74, 1.90, 1.83, 2.21, 2.00, 1.98];
% mu_cylinder_R = 0.1.*[3.52 1.20 1.86 1.88 1.82];




plot(sensor_sizes,20*log10(abs(mu_linear_R)./mu_ideal), 'm+--', DisplayName="Center Line Regression" , LineWidth=2)
hold on;
legend(Location='best');
xlabel("Sensor Grid Size [pixels]")
ylabel("dB")
grid;
set(gca, 'FontSize', 16, 'OuterPosition', [0 0 1 1])
set(gcf, 'PaperPositionMode', 'auto')
filename = join(['mu plot.png'], '');
saveas(gcf, filename);




%% Illumination Radius
D = [10 20 30 40 49];
mu_ideal = 0.174;

mu_linear_R = 0.1.*[2.30 1.87 1.78 1.767 1.76];
mu_cylinder_R= 0.1.*[2.36 1.88 1.78 1.764 1.76];


z = 6;
l=1;
y_1 = sqrt((z-l).^2+(D/2).^2);
y_2 = exp(-mu_ideal.*y_1)./(1-exp(-mu_ideal.*y_1));
y = 1-y_2.*((z-l)./y_1 - 1);




plot(D,20*log10(abs(mu_linear_R)./mu_ideal), 'm+--', DisplayName="Center Line Regression" , LineWidth=2)
hold on;
plot(D,20*log10(abs(mu_cylinder_R)./mu_ideal), 'bd--', DisplayName="Cylinder Regerssion" , LineWidth=2)
hold on;
plot(D,20*log10(y), 'r-', DisplayName="Theoretical Value for x=6" , LineWidth=2)
hold on;
legend(Location='best');
xlabel("Diameter [mm]")
ylabel("dB")
grid;
set(gca, 'FontSize', 16, 'OuterPosition', [0 0 1 1])
set(gcf, 'PaperPositionMode', 'auto')
filename = join(['mu plot.png'], '');
saveas(gcf, filename);



%% 2 layers
mu_a_2 = [0.01, 0.03, 0.05, 0.1, 0.2, 0.5];
mu_ideal_1 = 0.174;
mu_ideal_2 = sqrt(3.*mu_a_2.*(mu_a_2+1));

% Close
mu_linear_kWave_Close = 0.1.*[1.79 2.57  3.25   4.62 6.73 11.5];
mu_cylinder_kWave_Close = 0.1.*[1.79 2.57 3.24  4.61 6.73 11.5];


% Far
mu_linear_kWave_Far = 0.1.*[1.71 2.15 2.46 2.95 3.52 4.19];
mu_cylinder_kWave_Far = 0.1.*[1.71 2.14 2.45 2.95 3.49 4.15];


plot(mu_a_2,20*log10(mu_linear_kWave_Close./mu_ideal_2), 'm+--', DisplayName="2nd Layer Close" , LineWidth=2)
hold on;
plot(mu_a_2,20*log10(mu_linear_kWave_Far./mu_ideal_1), 'bd--', DisplayName="2nd Layer Far" , LineWidth=2)
% hold on;
% yline(mu_ideal, 'r--', 'LineWidth', 1, DisplayName="mu ideal")
legend(Location='best');
xlabel("\mu_a")
ylabel("dB")
grid;
set(gca, 'FontSize', 16, 'OuterPosition', [0 0 1 1])
set(gcf, 'PaperPositionMode', 'auto')
filename = join(['mu plot.png'], '');
saveas(gcf, filename);


%% Ellipse


mult = [1,2,4];
mu_ideal = 0.174;

% Ellipse D
mu_linear_kWaveD = 0.1.*[1.69 1.66 1.64];
mu_cylinder_kWaveD = 0.1.*[1.68 1.66 1.64];

% Ellipse C
mu_linear_kWaveC = 0.1.*[1.81 1.74 1.74];
mu_cylinder_kWaveC = 0.1.*[1.80 1.71 1.72];

%Ellipse B
mu_linear_kWaveB = 0.1.*[1.71 1.66 1.54];
mu_cylinder_kWaveB = 0.1.*[1.70 1.65 1.56];

%Ellipse A
mu_linear_kWaveA = 0.1.*[1.89 2.07 2.66];
mu_cylinder_kWaveA = 0.1.*[1.80 2.03 2.66];

plot(mult,20*log10(mu_linear_kWaveA/mu_ideal), 'm+--', DisplayName="Test A" , LineWidth=2)
hold on;
plot(mult,20*log10(mu_linear_kWaveB/mu_ideal), 'bd--', DisplayName="Test B" , LineWidth=2)
hold on;
plot(mult,20*log10(mu_linear_kWaveC/mu_ideal), 'ko--', DisplayName="Test C", LineWidth=2)
hold on
plot(mult,20*log10(mu_linear_kWaveD/mu_ideal), 'g^--', DisplayName="Test D", LineWidth=2)
legend(Location='best');
xlabel("q")
ylabel("dB")
grid;
set(gca, 'FontSize', 16, 'OuterPosition', [0 0 1 1])
set(gcf, 'PaperPositionMode', 'auto')
filename = join(['mu plot.png'], '');
saveas(gcf, filename);

%%

sigmas = [0, 0.1, 0.3, 0.5, 1];
% Correlated Variation
% mu_linear_kWave = [0.176 0.178 0.191 0.124 0.211];
% mu_cylinder_kWave = [0.175 0.178 0.189 0.127 0.211];

% 2d variation
% mu_expected = [0.174 0.173 0.169 0.161 0.145];
mu_linear_kWave = [0.177 0.185 0.178 0.248 0.558];
mu_cylinder_kWave = [0.176 0.185 0.175 0.239 0.566];
figure;


plot(sigmas,mu_expected, 'm+--', DisplayName="mu expected value" , LineWidth=2)
hold on;
plot(sigmas,mu_cylinder_kWave, 'bd--', DisplayName="mu cylinder" , LineWidth=2)
hold on;
plot(sigmas,mu_linear_kWave, 'ko--', DisplayName="mu linear", LineWidth=2)
hold on
yline(mu_ideal, 'r--', 'LineWidth', 1, DisplayName="mu ideal")
legend(Location='best');
xlabel("sigma")
set(gca, 'FontSize', 16, 'OuterPosition', [0 0 1 1])
set(gcf, 'PaperPositionMode', 'auto')
filename = join(['mu plot.png'], '');
saveas(gcf, filename);

%% db
% mu_cylinder_kWave_db = abs((mu_cylinder_kWave-mu_ideal)./mu_ideal);
% mu_linear_kWave_db = abs((mu_linear_kWave-mu_ideal)./mu_ideal);
mu_cylinder_kWave_db = 20*log10(mu_cylinder_kWave/mu_ideal);
mu_linear_kWave_db = 20*log10((mu_linear_kWave)./mu_ideal);
% mu_expected_db = abs((mu_expected-mu_ideal)./mu_ideal);
figure;
% plot(sigmas,mu_expected_db, 'ms--', DisplayName="\mu expected value" , LineWidth=2)
% hold on;
plot(sigmas,mu_cylinder_kWave_db, 'bd--', DisplayName="Cylinder Regression" , LineWidth=2)
hold on;
plot(sigmas,mu_linear_kWave_db, 'ko--', DisplayName="Center Line Regression", LineWidth=2)
grid;
hold on
legend(Location='best');
xlabel("sigma")
ylabel("dB")
set(gca, 'FontSize', 16, 'OuterPosition', [0 0 1 1])
set(gcf, 'PaperPositionMode', 'auto')
filename = join(['mu plot.png'], '');
saveas(gcf, filename);
