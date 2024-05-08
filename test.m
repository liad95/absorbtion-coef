%% initialization

addpath('ValoMC/')
cd ValoMC
compile_vmc_mex
%% 1d sim
N = 10;
Na = 7;
Ns = 7;
mua  = logspace(log10(0.001), log10(0.2), Na);
mus = logspace(log10(0.01), log10(20), Ns)
mua_rect = mua;
dx = 0.1e-3;
dy = 0.1e-3;
Nx = 150;
Ny = 150;
lightsource_type='direct';
options.photon_count = 1e5;
locations = floor(linspace(2,148,N));
distances = dy*locations;
H_rect = zeros(N, Na, Ns);
x = zeros(N, Na, Ns);
for k=1:Na
    for l=1:Ns
        for i=1:N
            H = oneDsim(locations(i),mua(k),mus(l), mua_rect(k), dx, dy, Nx, Ny, lightsource_type, options);
            H_rect(i, k, l) = mean(H(:, locations(i)));
            fprintf("mua %d, mus %d, N %d", k, l, i)
%             figure;
%             imagesc(H, [min(H(:)) ...
%                     max(H(:))]);
%              
%             colormap default;
%             
%             c = colorbar;  % create a colorbar
%             axis image;
%             title('Absorbed energy [J/m3]');
%         
        end
    end
end

%% analyzing results for complex fit


slopes = zeros(Na, Ns);
estimated_slopes = zeros(Na, Ns)
fo = fitoptions( 'Method', 'NonlinearLeastSquares',...
               'Lower',[0,0.1, 0.1, 0.1],...
               'Upper',[1e-2,2, 10, 20],...
               'StartPoint',[8e-3, 0.5, 2, -2]);
ft = fittype('fitfunc(x,a,b,z0,z1)', 'options', fo)
for k=1:Na
    for l=1:Ns
        estimated_slopes(k, l) = sqrt(3*mua(k)*(mua(k)+0.1*mus(l)));
        results = fit(distances', log(H_rect(:, k, l)), ft);
%         if results.b>0
%             slopes(k, l) = results.b
%         else
%             slopes(k, l) = 0
%         end
        slopes(k, l) = results.b
        
        figure;
        plot(distances, H_rect(: ,k, l), 'DisplayName','Data')
        hold on
        plot(distances, fitfunc(distances, results.a, results.b, results.z0, results.z1), 'DisplayName','Fit')
        legend;

        xlabel("distance [m] ")
        ylabel("Heat")
        title("Heat as function of object depth")
    end
end
figure;
surf(mua, mus, slopes)
hold on
surf(mua, mus, estimated_slopes)

%% analyzing results
slopes = zeros(Na, Ns);
estimated_slopes = zeros(Na, Ns)
for k=1:Na
    for l=1:Ns
        lm = fitlm(distances(1:end-1), log(H_rect(1:end-1, k, l)));
        coefficients_table = lm.Coefficients;
        slopes(k, l)  = -1e-3*coefficients_table.Estimate(2);
        estimated_slopes(k, l) = sqrt(3*mua(k)*(mua(k)+0.1*mus(l)));
        figure;
        plot(distances(1:end-1), log(H_rect(1:end-1 ,k, l)))
        xlabel("distance [m] ")
        ylabel("Heat")
        title("log Heat as function of object depth")
    end
end
figure;
surf(mua, mus, slopes)
hold on
surf(mua, mus, estimated_slopes)

%% analyzing results
% need to take into account factor of 1/r in the case of isotropic input
figure;
plot(distances(3:end), H_rect(3:end))
xlabel("distance [m] ")
ylabel("Heat")
title("Heat as function of object depth")

figure;
plot(distances(3:end), log(H_rect(3:end)))
xlabel("distance [m] ")
ylabel("Heat")
title("log Heat as function of object depth")

% Perform linear regression
lm = fitlm(distances(3:end), log(H_rect(3:end)));

% Get the coefficients table
coefficients_table = lm.Coefficients;

% Extract the slope (coefficient) from the table
slope = coefficients_table.Estimate(2); % Coefficient for the x variable

% Display the slope
disp(['Slope: ', num2str(-1e-3*slope)]);
mu_eff = sqrt(3*mua*(mua+0.1*mus));
disp(['mu effective ', num2str(mu_eff)]);



%% 2d sim
N = 10;
Na = 7;
Ns = 7;
mua  = logspace(log10(0.001), log10(0.2), Na);
mus = logspace(log10(0.01), log10(20), Ns)
mua_rect = 3*mua;
dx = 0.1e-3;
dy = 0.1e-3;
Nx = 150;
Ny = 150;
rad = 3;
lightsource_type='isotropic';
options.photon_count = 2e5;
locations = floor(linspace(2,148,N));
distances = dy*locations;
H_rect = zeros(N, Na, Ns);
x = zeros(N, Na, Ns);
for k=1:Na
    for l=1:Ns
        for i=1:N
            H = circleInMedium(locations(i),mua(k),mus(l), mua_rect(k), dx, dy, Nx, Ny, lightsource_type, options ,rad);
%             figure;
%             imagesc(H, [min(H(:)), max(H(:))]);
%             colormap default;
%             c = colorbar;  % create a colorbar
%             axis image;
%             title('Absorbed energy [J/m3]');

            H = H.*makeDisc(Nx, Ny, floor(Nx/2), locations(i), rad);
%             figure;
%             imagesc(H, [min(H(:)), max(H(:))]);
%             colormap default;
%             c = colorbar;  % create a colorbar
%             axis image;
%             title('Absorbed energy in disc location [J/m3]');
            circle_loc = find(H);
            H_rect(i, k, l) = mean(H(circle_loc));
            fprintf("mua %d, mus %d, N %d", k, l, i)
        end
    end
end

%% 2d K-Wave sim
mua  = 0.05;
mus = 5;
mua_rect = 3*mua;
dx = 0.1e-3;
dy = 0.1e-3;
Nx = 150;
Ny = 150;
rad = 3;
lightsource_type='isotropic';
locations = 25;
distances = dy*locations;
kgrid = makeGrid(Nx, dx, Ny, dy);
medium.sound_speed = 1500;    % [m/s]

H = circleInMedium(locations,mua,mus, mua_rect, dx, dy, Nx, Ny, lightsource_type, rad);


%% testing k-wave
source.p0 = ones(Nx, Ny);
% source.p0 = smooth(source.p0, true);
% Define a circular sensor
sensor_radius = 6e-3;       % [m]
num_sensor_points = 50;     % number of sensor points   
mask = zeros(Nx, Ny)
mask(1, :) = 1;
sensor.mask = mask;
sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, 'PMLInside', false);

%% 2d K-Wave sim
N = 10;
Na = 1;
Ns = 1;
mua  = linspace(0.05, 0.2, Na);
mus = linspace(5,20,Ns);
mua_rect = 10*mua;
dx = 0.1e-3;
dy = 0.1e-3;
Nx = 150;
Ny = 150;
rad = 3;
lightsource_type='isotropic';
locations = floor(linspace(25,125,N));
distances = dy*locations;
H_rect = zeros(N, Na, Ns);
x = zeros(N, Na, Ns);
kgrid = makeGrid(Nx, dx, Ny, dy);
medium.sound_speed = 1500;    % [m/s]

for k=1:Na
    for l=1:Ns
        for i=1:N
            H = circleInMedium(locations(i),mua(k),mus(l), mua_rect(k), dx, dy, Nx, Ny, lightsource_type, rad);
            source.p0 = H;
            source.p0 = smooth(source.p0, true);
            kgrid.makeTime(medium.sound_speed);
            % Define a circular sensor
            sensor_radius = 6e-3;       % [m]
            num_sensor_points = 50;     % number of sensor points
            mask = zeros(Nx, Ny)
            mask(1, :) = 1;
            sensor.mask = mask;
            sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, 'PMLInside', false);
            p_xy = kspaceLineRecon(sensor_data.', dy, kgrid.dt, medium.sound_speed, 'Plot', true);
            % define a second k-space grid using the dimensions of p_xy
            [Nx_recon, Ny_recon] = size(p_xy);
            kgrid_recon = kWaveGrid(Nx_recon, kgrid.dt * medium.sound_speed, Ny_recon, dy);
            
            % resample p_xy to be the same size as source.p0
            p_xy_rs = interp2(kgrid_recon.y, kgrid_recon.x - min(kgrid_recon.x(:)), p_xy, kgrid.y, kgrid.x - min(kgrid.x(:)));




            figure;
            imagesc(source.p0, [min(source.p0(:)), max(source.p0(:))]);
            colormap default;
            c = colorbar;  % create a colorbar
            axis image;
            title('Initial Pressure [J/m3]');


            
            
            imagesc(p_xy_rs, [min(p_xy_rs(:)), max(p_xy_rs(:))]);
            colormap default;
            c = colorbar;  % create a colorbar
            axis image;
            title('Reconstructed Pressure [J/m3]');

            

            p_xy_rs = p_xy_rs.*makeDisc(Nx, Ny, floor(Nx/2), locations(i), rad);
            figure;
            imagesc(p_xy_rs, [min(p_xy_rs(:)), max(p_xy_rs(:))]);
            colormap default;
            c = colorbar;  % create a colorbar
            axis image;
            title('Reconstructed Pressure in disc location [J/m3]');
            circle_loc = find(p_xy_rs);
            H_rect(i, k, l) = mean(p_xy_rs(circle_loc));
        end
    end
end


%% analyzing results
% need to take into account factor of 1/r in the case of isotropic input
figure;
plot(distances, H_rect)
xlabel("distance [m] ")
ylabel("Heat")
title("Heat as function of object depth")

figure;
plot(distances(1:8), log(H_rect(1:8)))
xlabel("distance [m] ")
ylabel("Heat")
title("log Heat as function of object depth")

% Perform linear regression
lm = fitlm(distances(1:8), log(H_rect(1:8)));

% Get the coefficients table
coefficients_table = lm.Coefficients;

% Extract the slope (coefficient) from the table
slope = coefficients_table.Estimate(2); % Coefficient for the x variable

% Display the slope
disp(['Slope: ', num2str(-1e-3*slope)]);
mu_eff = sqrt(3*mua*(mua+0.1*mus));
disp(['mu effective ', num2str(mu_eff)]);

%% Function for fit
function y = func(x,a,b,z0,z1)
y = a*(exp(-b*(x-z0))./(x-z0) - exp(-b*(x-z1))./(x-z1))
end