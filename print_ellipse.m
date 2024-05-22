%% simulation parameters
dx = 0.05;
dy = 0.5;
Nx = 201;
Ny = 101;
x_arr = -dx*(Nx-1)/2:dx:dx*(Nx-1)/2;
y_arr = -dy*(Ny-1)/2:dy:dy*(Ny-1)/2;
mua=0.01;
ellipse_mult = 4;
ellipse_mua = (ellipse_mult)*mua;
ellipse_mus = 0;
x_center = 120;
y_center = 30;
x_radius = 7;
y_radius = 15;
angle = 45;

mus=10;
sigma = 0.1;
g = 0.9;
gruneisen_coef = 0.5;




ellipse = create_ellipse_2d(Ny, Nx, x_center, y_center, x_radius, y_radius, angle);
close all;
figure;
imagesc(x_arr, y_arr, ellipse);
xlabel("x [mm]")
ylabel("y [mm]")
set(gca, 'FontSize', 16, 'OuterPosition', [0 0 1 1])
set(gcf, 'PaperPositionMode', 'auto')
colorbar;
filename = join(['Ellipse.png'], '')
saveas(gcf, filename);