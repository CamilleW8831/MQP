clc; clear; close all;

% This code is configured to run RAW DATA

% USER INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% File name:
% this should be a string referring to a .mat file containing three
% matrices with spatially orientinted magnetic field components labelled
% x, y, z (otherwise names must be changed where the load file is
% referenced below)
fileName = "Pure_MXene.mat"; 

% Pixel size:
xDelta = 5.86e-6; % length of pixels in x-direction (horizontal) [m]
yDelta = 5.86e-6; % length of pixels in y-direction (vertical) [m]

% Measurement height:
zp = 5e-7; % estimated height of the measurement plane [m]

% Hanning window: 
% the Hanning window offers the option of filtering out
% noisier high-frequency components with a tradeoff of resolution.

% "ON" to activate the hanning window
% "OFF" to ignore the hanning window
% To avoid unncessary processing, if the hanning window is not used, set
% hannWindow= "OFF"
hannWindow= "OFF";

% hannRatio defines a radius of spatial frequencies allowed to remain within 
% the reconstruction based on a percentage of the largest magnitude.
% (e.g. hannRatio= 0.8 corresponds to only freqeuncies smaller than 80% of 
% the Nyquist frequency magnitude to remain)
hannRatio= 0.8;

% Magnetic field plotting:
% "ON" to plot magnetic fields
% "OFF" to only plot current density components and magnitude
plotMagField= "OFF";

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% INVERSE PROCEDURE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SETUP: CONSTANTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mu0 = 4 * pi * 1e-7; % permeability of free space [m*kg/s^2*A^2] OR [T*m/A]

% SETUP: LOAD FILES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load the magnetic field data:
loadBFields= load(fileName);

Bx = loadBFields.x * 1e-6; % magnetic field x-component [T] 
By = loadBFields.y * 1e-6; % magnetic field y-component [T]
Bz = loadBFields.z * 1e-6; % magnetic field z-component [T]

% SETUP: DEFINE THE GRID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[M,N] = size(Bx); % dimensions of the magnetic field matrices
x = 0 : xDelta : xDelta*(N-1); % [m]
y = 0 : yDelta : yDelta*(M-1); % [m]
[X,Y] = meshgrid(x,y); % organizes points into a grid

% STEP 1: ENTER THE FOURIER DOMAIN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate spatial frequency indices for Fourier transform.
kxDelt = 2 * pi * 1 / (N * xDelta); % [rad/m]
kyDelt = 2 * pi * 1 / (M * xDelta); % [rad/m]

% kx and ky are the spatial frequency vectors
% -N/2 : (N/2 - 1)
kx = kxDelt * (-N/2) : kxDelt : kxDelt * (N/2-1); % [rad/m]
ky = kyDelt * (-M/2) : kyDelt : kyDelt * (M/2-1); % [rad/m]
[Kx, Ky] = meshgrid(kx, ky); % organizes points into a grid

% precompute values:
Kmag = sqrt(Kx.^2 + Ky.^2); % the magnitude of spatial frequencies [rad/m]
g= mu0/2 * exp(-Kmag * zp); % Fourier green's function [T*m/A]

% apply the 2D fourier transform to Bx, By, and Bz, shift the output so
% that it is centered around kx=0, ky=0. This aligns the correct Fourier 
% magnitudes with the manner in which the spatial frequencies were set up 
% (shifted configuration)
bx = fftshift(fft2(Bx));
by = fftshift(fft2(By));
bz = fftshift(fft2(Bz));

% HANNING WINDOW
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if hannWindow == "YES"

    kx_nyq = pi / xDelta; % x Nyquist frequency [rad/m]
    ky_nyq = pi / yDelta; % y Nyquist frequenct [rad/m]
    k_nyq = sqrt(kx_nyq^2 + ky_nyq^2);  % Nyquist radius [rad/m]

    kmax = hannRatio * k_nyq; % frequency limit [rad/m]
    
    % apply the Hanning window
    HannWindow = 0.5 * (1 + cos(pi * Kmag / kmax));
    HannWindow(Kmag > kmax) = 0; % values above kmax are set to 0
    
    % Apply window to Fourier-transformed magnetic field components
    bx = HannWindow .* bx;
    by = HannWindow .* by;
    bz = HannWindow .* bz;

end

% STEP 2: SOLVE SYSTEM OF EQUATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% iterate every pair of spatial frequencies
for v = 1:length(ky)
    for u = 1:length(kx)

        if Kmag(v,u)==0 % adjusts the system to avoid singularities when Kmag=0

            A = [0, g(v,u); -g(v,u), 0];

            b = [bx(v,u); by(v,u)];

        else % all other cases (when Kmag =/ 0)

            A = [0, g(v,u), 0; ...
                -g(v,u), 0, 0; ...
                -1i * g(v,u) * Ky(v,u)/Kmag(v,u), 1i * g(v,u) * Kx(v, u)/Kmag(v,u), 1];

            b = [bx(v,u); by(v,u); bz(v,u)];

        end

        % Solve for j 
        j = A \ b; 

        % Save jx and jy
        jx_invert(v, u) = j(1);
        jy_invert(v, u) = j(2);

    end
end

% STEP 3: EXIT THE FOURIER DOMAIN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% unshift the magnitudes and apply an inverse 2D Fourier transform to
% obtain the components of current density in the spatial domain (only the
% real part is logical).
Jx_OUT= real(ifft2(ifftshift(jx_invert))); % [A/m^2]
Jy_OUT= real(ifft2(ifftshift(jy_invert))); % [A/m^2] 

% calculate the magnitude of current density
Jnorm_OUT= sqrt(Jx_OUT.^2 + Jy_OUT.^2); % [A/m^2]

% PLOT THE CURRENT DENSITY COMPONENTS AND MAGNITUDE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Jnorm
figure();
imagesc(x, y, Jnorm_OUT);
% labels:
title("Reconstructed Current Density [A/m^2]");
xlabel("x [m]"); 
ylabel("y [m]");
% graph/text settings: 
set(gca, 'Ydir', 'normal', 'FontSize', 12, 'TickLength',[0, 0]); 
set(gcf, 'Color', 'w', 'Position', [100, 350, 400, 400]);
axis equal tight;
% color bar bounds and other settings:
colormap(magma);
colorbar;
clim_min = prctile(Jnorm_OUT(:), 1);
clim_max = prctile(Jnorm_OUT(:), 95);
clim([clim_min,clim_max]);
c= colorbar; c.TickLength= 0;

% Jx
figure();
imagesc(x, y, Jx_OUT); 
% labels:
title("Reconstructed J_x [A/m^2]");
xlabel("x [m]"); 
ylabel("y [m]");
% graph/text settings: 
set(gca, 'Ydir', 'normal', 'FontSize', 12, 'TickLength',[0, 0]); 
set(gcf, 'Color', 'w', 'Position', [550, 350, 400, 400]);
axis equal tight;
% color bar bounds and other settings:
colormap(magma); 
colorbar;
clim_min = prctile(Jx_OUT(:), 1);
clim_max = prctile(Jx_OUT(:), 95);
clim([clim_min,clim_max]);
c= colorbar; c.TickLength= 0;

% Jy
figure();
imagesc(x, y, Jy_OUT); 
% labels:
title("Reconstructed J_y [A/m^2]");
xlabel("x [m]"); 
ylabel("y [m]");
% graph/text settings: 
set(gca, 'Ydir', 'normal', 'FontSize', 12, 'TickLength',[0, 0]); 
set(gcf, 'Color', 'w', 'Position', [1000, 350, 400, 400]);
axis equal tight;
% color bar bounds and other settings:
colormap(magma); 
colorbar;
clim_min = prctile(Jy_OUT(:), 1);
clim_max = prctile(Jy_OUT(:), 95);
clim([clim_min,clim_max]);
c= colorbar; c.TickLength= 0;

% PLOT MAGNETIC FIELDS (OPTIONAL)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if plotMagField == "ON"

    % Bx
    figure();
    imagesc(x, y, Bx); 
    % labels:
    title("B_x [T]");
    xlabel("x [m]"); 
    ylabel("y [m]");
    % graph/text settings: 
    set(gca, 'Ydir', 'normal', 'FontSize', 12, 'TickLength',[0, 0]); 
    set(gcf, 'Color', 'w', 'Position', [100, 100, 400, 400]);
    axis equal tight;
    % color bar bounds and other settings:
    colormap(turbo); 
    colorbar;
    clim_min = prctile(Bx(:), 1);
    clim_max = prctile(Bx(:), 95);
    clim([clim_min,clim_max]);
    c= colorbar; c.TickLength= 0;
    
    % By
    figure();
    imagesc(x, y, By); 
    % labels:
    title("B_y [T]");
    xlabel("x [m]"); 
    ylabel("y [m]");
    % graph/text settings: 
    set(gca, 'Ydir', 'normal', 'FontSize', 12, 'TickLength',[0, 0]); 
    set(gcf, 'Color', 'w', 'Position', [550, 100, 400, 400]);
    axis equal tight;
    % color bar bounds and other settings:
    colormap(turbo);
    colorbar;
    clim_min = prctile(By(:), 1);
    clim_max = prctile(By(:), 95);
    clim([clim_min,clim_max]);
    c= colorbar; c.TickLength= 0;
    
    % Bz
    figure();
    imagesc(x, y, Bz); 
    % labels:
    title("B_z [T]");
    xlabel("x [m]"); 
    ylabel("y [m]");
    % graph/text settings: 
    set(gca, 'Ydir', 'normal', 'FontSize', 12, 'TickLength',[0, 0]); 
    set(gcf, 'Color', 'w', 'Position', [1000, 100, 400, 400]);
    axis equal tight;
    % color bar bounds and other settings:
    colormap(turbo); 
    colorbar;
    clim_min = prctile(Bz(:), 1);
    clim_max = prctile(Bz(:), 95);
    clim([clim_min,clim_max]);
    c= colorbar; c.TickLength= 0;

end
