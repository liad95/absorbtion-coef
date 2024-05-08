%% Create the geometry
% Create a circular mesh using Toast++ and import it to ValoMC

rad = 10;
nsect = 8;
nring = 8;
nbnd = 2;


[vtx,idx,eltp] = mkcircle(rad,nsect,nring,nbnd);

toastmesh = toastMesh(vtx,idx,eltp);
vmcmesh = importToastMesh(toastmesh);


%% Set the optical coefficients
% Note that the optical coefficients are given for each node in Toast,
% whereas in ValoMC they are uniform values for each element.

% absorption coefficient [1/mm]
mua_bkg = 0.01;  
% scattering coefficient [1/mm]
mus_bkg = 1.0;                   	        
% scattering anisotropy parameter [unitless]
scattering_anisotropy_bkg = 0.0; 	        
% reduced scattering coefficient [1/mm]
mus_reduced = mus_bkg*(1-scattering_anisotropy_bkg); 

% Set optical coefficients for Toast.

nnd = toastmesh.NodeCount;
toast_mua = ones(nnd,1) * mua_bkg;  % absorption coefficient [1/mm]
toast_mus = ones(nnd,1) * mus_reduced; % reduced scattering coefficient [1/mm]

% Set optical coefficients for ValoMC. The refractive index is set but it
% does not affect the solution as there is no mismatch on the boundary.

nne = size(vmcmesh.H,1); % number of elements

% absorption coefficient [1/mm]
vmcmedium.absorption_coefficient = ones(nne,1)*mua_bkg;
 % scattering coefficient [1/mm]
vmcmedium.scattering_coefficient = ones(nne,1)*mus_bkg;
 % refractive index [unitless]
vmcmedium.refractive_index = ones(nne,1)*1;
% scattering anisotropy parameter [unitless]
vmcmedium.scattering_anisotropy = ones(nne,1)*scattering_anisotropy_bkg;

% Create the boundary so that there is no mismatch
vmcboundary = createBoundary(vmcmesh, vmcmedium);

%% Create the source and detector positions
%
% A collimated lightsoure (pencil beam) can be approximated in DA
% by placing an isotropic source at a distance 1/mus' from the surface,
% where mus' is the reduced scattering coefficient.

% Build source/detector locations for Toast
l = 1/(mus_reduced + mua_bkg)

nq = 10;
Q(1,:) = (rad - l) * [0 1]; % source position
for ii=1:nq
  M(ii,:) = rad*((ii-1)/nq) * [0 1]; % detector position
end

toastmesh.SetQM(Q,M);


% Build source/detector locations for ValoMC

% Sources and detectors are placed to the nearest boundary elements

source_boundary_elements = findBoundaries(vmcmesh, 'location', Q);
detector_boundary_elements = findBoundaries(vmcmesh, 'location', M);


%% Plot the source/detector arrangement
% Note that the lightsources have a finite width in ValoMC, which introduces
% a small discretisation error in the comparison.

figure('rend','painters','pos',[10 10 1000 1000])
hold on
patch('Faces', vmcmesh.H, 'Vertices',vmcmesh.r, 'FaceVertexCData', 0, 'FaceColor', 'flat', 'LineWidth', 0.1);


h3 = plot(Q(:, 1),Q(:,2),'ro','MarkerFaceColor','r');
h4 = plot(M(:, 1),M(:,2),'bs','MarkerFaceColor','b');

title('Source/detector setup');
legend([h3 h4], {'Toast source', 'Toast detector'});
hold off



%% Solve the FEM system with Toast
% The system matrix is constructed manually using 2D coefficients (by
% default, Toast uses formulas derived from the radiative transfer equation
% for 3D geometry). For more detailed information about 2D and 3D
% coefficients, see e.g. T. Tarvainen: Computational Methods for Light
% Transport in Optical Tomography, PhD thesis, University of Kuopio, 2006.


% Create isotropic sources 
qvec = toastmesh.Qvec('Isotropic','Point');

% To obtain comparable results, the measurement vectors have a sharp
% Gaussian profile (a narrow detector). Conversion factor between fluence
% and exitance is 2/pi in 2D
mvec = 2/pi*toastmesh.Mvec('Gaussian',0.5,0);

% 2D diffusion coefficient
toast_kap = 1./(2.*(toast_mua+toast_mus));

S1 = toastmesh.SysmatComponent('PFF',toast_mua);
S2 = toastmesh.SysmatComponent('PDD',toast_kap);

% The boundary term is multiplied by 2/pi in 2D
S3 = toastmesh.SysmatComponent ('BndPFF', ones(nnd,1)*2/pi);
K = S1+S2+S3;

Phi = K\(qvec);                % solve the fluence
Y_toast  = mvec.' * Phi;       % compute the exitance on each detector



%% Plot the solutions



figure;
plot(log(flip(Y_toast(2:end))));
fit = fitlm(rad*((1:10 - 1)/10), log(flip(Y_toast(2:end))))
title('Toast result');




