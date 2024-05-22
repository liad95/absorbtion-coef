%% simulation parameters
dx = 0.05;
dy = 0.5;
Nx = 201;
Ny = 101;
x_arr = -dx*(Nx-1)/2:dx:dx*(Nx-1)/2;
y_arr = -dy*(Ny-1)/2:dy:dy*(Ny-1)/2;
mua_layer = 1;

x_layer = 70;
[X,Y] = meshgrid(x_arr,y_arr); % Matlab function
absorption_matrix = zeros(size(X));
layer_idx = X>x_arr(x_layer);
absorption_matrix(layer_idx) = mua_layer;

figure;
imagesc(x_arr, y_arr ,absorption_matrix)

xlabel("x [mm]")
ylabel("y [mm]")
set(gca, 'FontSize', 16, 'OuterPosition', [0 0 1 1])
set(gcf, 'PaperPositionMode', 'auto')
colorbar;
filename = join(['Layers Close.png'], '')
saveas(gcf, filename);