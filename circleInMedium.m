
function [H] = oneDsim(y_loc,mua,mus, mua_rect, dx, dy, Nx, Ny, lightsource_type, options, rad)

%% Simulating the photoacoustic effect using k-Wave: kwavetest.m
% This example demonstrates simulation of a pressure field generated through 
% the absorption of an externally introduced light pulse.
% The light propagation is simulated using ValoMC and the propagation and
% detection of pressure wavefield is simulated using k-Wave, see 
% http://www.k-wave.org/documentation/k-wave_initial_value_problems.php.
% The example also shows how the computation grid of k-Wave and mesh of
% ValoMC can be made compatible.
% Note that k-Wave uses SI units (e.g. [m]) and ValoMC works in millimetre-scale 
% (e.g. [mm]).
%
% k-Wave is an open source acoustics toolbox for MATLAB and C++ developed
% by Bradley Treeby and Ben Cox (University College London) and Jiri Jaros 
% (Brno University of Technology). The software is designed for time domain 
% acoustic and ultrasound simulations in complex and tissue-realistic media. 
%
% k-Wave homepage: http://www.k-wave.org/.
% B. E. Treeby and B. T. Cox: "k-Wave: MATLAB toolbox for the simulation
% and reconstruction of photoacoustic wave-fields", Journal of Biomedical 
% Optics, 15(2):021314, 2010.
%
%
% <html>
% <font color="red">Note: there was an incorrectly explained unit conversion earlier in this example.
% See the text shown in red. </font>
% </html>
%

%% k-Wave initialisation
% The initialisation is done as normal in k-Wave.
% Care must be taken at the initialization ValoMC, to make a 
% matching computational simulation area for (see ValoMC initialization)
% the photon propagation simulation.


% Create the k-Wave grid
% Nx = 150;           % number of grid points in the x (row) direction
% Ny = 150;           % number of grid points in the y (column) direction
% dx = 0.1e-3;        % grid point spacing in the x direction [m]
% dy = 0.1e-3;        % grid point spacing in the y direction [m]
kgrid = makeGrid(Nx, dx, Ny, dy);
discs = makeDisc(Nx, Ny, floor(Nx/2), y_loc, rad);


% Define the acoustic properties 
disc_indices = find(discs==1);


%% Create a ValoMC mesh
% ValoMC uses triangles and tetrahedrons as the basis shape, whereas
% in k-Wave pixels and voxels are used. createGridMesh can be used to
% create a straightforward mapping between the triangles and pixels.
% *Note that K-Wave uses matrices in the format matrix(X,Y), whereas* 
% *MATLAB t(c.f. meshgrid) and ValoMC uses matrix(Y,X)*
% Therefore x and y should be swapped when moving between ValoMC
% arrays and k-Wave arrays

vmcmesh = createGridMesh(kgrid.y_vec*1e3, kgrid.x_vec*1e3); % [m to mm]
%% Define optical coefficients
% For users accustomed to k-Wave, the optical coefficients can be set 
% in similar fashion as in k-Wave, i.e. using multidimensional arrays.
% If one of the arrays defining the medium is given as multidimensional 
% array to ValoMC, the code will assume that the mesh was created
% using 'createGridMesh' and the output fluence will also given as 
% a two dimensional array in solution.grid_fluence. See the example 
% 'Working with pixel and  voxel data' on how to achieve similar 
% assignments using one dimensional indexing.

vmcmedium.scattering_coefficient = mus*ones(Nx, Ny);
vmcmedium.absorption_coefficient = mua*ones(Nx, Ny);
% Define the acoustic properties    

vmcmedium.absorption_coefficient(disc_indices) = mua_rect;
vmcmedium.scattering_anisotropy = 0.9;           % scattering anisotropy parameter [unitless]
vmcmedium.refractive_index = 1.0*ones(Nx, Ny);
% vmcmedium.refractive_index(:,Ny/2:end) = 1.4;    % refractive index [unitless]

% Define the Gruneisen parameter describing photoacoustic efficiency
vmcmedium.gruneisen_parameter = 0.02*ones(Nx, Ny);      % [unitless]

%% Create a light source 

% Set a light source with a width of 2 mm and cosinic directional profile 
% in -x direction
% boundary_with_lightsource = findBoundaries(vmcmesh, 'direction', ...
%                                            [0 0], ...
%                                            [-1000 0], ...
%                                            250);
                                    
vmcboundary.lightsource(floor(Nx/2)) = {lightsource_type};
% vmcboundary.lightsource_gaussian_sigma(boundary_with_lightsource) = 0.1;


%% Run the Monte Carlo simulation
solution = ValoMC(vmcmesh, vmcmedium, vmcboundary, options);
H = vmcmedium.absorption_coefficient .* solution.grid_fluence;
% imagesc(kgrid.y_vec*1e3, kgrid.x_vec*1e3, H, [min(H(:)) ...
%         max(H(:))]);
%  
%  xlabel('y-position [mm]');
%  ylabel('x-position [mm]');
% colormap default;
% 
% c = colorbar;  % create a colorbar
% axis image;
% title('Absorbed energy [J/m3]');



