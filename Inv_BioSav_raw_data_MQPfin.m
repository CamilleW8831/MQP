clc; clear; close all;

% This code is configured to run RAW DATA

% INPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % File name, defining magnetic field inputs
loadMagField= load("Pure_MXene.mat");
Bx = loadMagField.x * 1e-6; % [T] 
By = loadMagField.y * 1e-6; % [T]
Bz = loadMagField.z * 1e-6; % [T]

zp = 5 * 1e-7; % This is a guess for the height of the measurement plane
xDelt = 5.86e-6; % pixel size in x
yDelt = 5.86e-6; % pixel size in y

% DEFINE THE GRID %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[M,N]= size(Bx); % dimensions of the grid
x= 0 : xDelt : xDelt*(N-1); % [m]
y= 0 : yDelt : yDelt*(M-1); % [m]
[X,Y]= meshgrid(x,y);

% CONSTANTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mu0 = 4 * pi * 1e-7; % Permeability of free space [m*kg/s^2*A^2] OR [T*m/A]

% RECONSTRUCT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate spatial frequency indices for Fourier transform.
kxDelt = 2 * pi * 1 / (N * xDelt); % [rad/m]
kyDelt = 2 * pi * 1 / (M * xDelt); % [rad/m]

% kx and ky are the spatial frequency vectors
% -N/2 : (N/2 - 1)
kx = kxDelt * (-N/2) : kxDelt : kxDelt * (N/2-1); % [rad/m]
ky = kyDelt * (-M/2) : kyDelt : kyDelt * (M/2-1); % [rad/m]
[Kx, Ky] = meshgrid(kx, ky);
Kmag = sqrt(Kx.^2 + Ky.^2);

% Fourier Transform
bx = fftshift(fft2(Bx));
by = fftshift(fft2(By));
bz = fftshift(fft2(Bz));

% % % %%% TRYING WINDOW
% % % % Determine kmax to include 80% of radial frequency range
% % % kx_nyq = pi / xDelt;
% % % ky_nyq = pi / yDelt;
% % % k_nyq = sqrt(kx_nyq^2 + ky_nyq^2);  % Radial Nyquist frequency
% % % kmax_ratio = 0.8;  % include 80% of radial frequency content
% % % kmax = kmax_ratio * k_nyq;
% % % 
% % % % Hann taper in radial frequency domain: 0.5 * (1 + cos(pi * |K| / kmax)) for |K| <= kmax
% % % HannFilter = 0.5 * (1 + cos(pi * Kmag / kmax));
% % % HannFilter(Kmag > kmax) = 0;  % zero out values beyond cutoff
% % % 
% % % % Apply window to Fourier-transformed magnetic field components
% % % bx = HannFilter .* bx;
% % % by = HannFilter .* by;
% % % bz = HannFilter .* bz;
% % % %%%

g= mu0/2 * exp(-Kmag * zp); % Fourier green's function

for v = 1:length(ky)
    for u = 1:length(kx)

        if Kmag(v,u)==0
            A = [0, g(v,u); -g(v,u), 0];
            b = [bx(v,u); by(v,u)];
        else
            A = [0, g(v,u), 0; -g(v,u), 0, 0; -1i * g(v,u) * Ky(v,u)/Kmag(v,u), 1i * g(v,u) * Kx(v, u)/Kmag(v,u), 1];
            b = [bx(v,u); by(v,u); bz(v,u)];
        end

        % Solve for j using a least squares solution to handle the overdetermined system
        j = A \ b; 

        % Save jx and jy
        jx_invert(v, u) = j(1);
        jy_invert(v, u) = j(2);
    end
end

Jx_OUT= real(ifft2(ifftshift(jx_invert)));
Jy_OUT= real(ifft2(ifftshift(jy_invert)));

Jnorm_OUT= sqrt(Jx_OUT.^2 + Jy_OUT.^2);

clim_min = 0;
clim_max = 11;
% Plot the reconstructed current norm with arrows indicating direction.
figure();
imagesc(x, y, Jnorm_OUT); 
xlabel("x [m]"); 
ylabel("y [m]");
title("Reconstructed Current Density [A/m]");
set(gca, 'Ydir', 'normal', 'FontSize', 12, 'TickLength',[0, 0]); 
set(gcf, 'Color', 'w');
axis equal tight;
colormap(magma); 
clim([0,clim_max]);
c= colorbar;
c.TickLength= 0;

clim_min = prctile(Jx_OUT(:), 1);
clim_max = prctile(Jx_OUT(:), 95);
figure();
imagesc(x, y, Jx_OUT); 
xlabel("x [m]"); 
ylabel("y [m]");
title("Reconstructed J_x [A/m]");
set(gca, 'Ydir', 'normal', 'FontSize', 12, 'TickLength',[0, 0]); 
set(gcf, 'Color', 'w');
axis equal tight;
colormap(magma); 
clim([clim_min,clim_max]);
colorbar;

clim_min = prctile(Jy_OUT(:), 1);
clim_max = prctile(Jy_OUT(:), 95);
figure();
imagesc(x, y, Jy_OUT); 
xlabel("x [m]"); 
ylabel("y [m]");
title("Reconstructed J_y [A/m]");
set(gca, 'Ydir', 'normal', 'FontSize', 12, 'TickLength',[0, 0]); 
set(gcf, 'Color', 'w');
axis equal tight;
colormap(magma); 
clim([clim_min,clim_max]);
c= colorbar;
c.TickLength= 0;

clim_min = prctile(Bx(:), 1);
clim_max = prctile(Bx(:), 95);
figure();
imagesc(x, y, Bx); 
xlabel("x [m]"); 
ylabel("y [m]");
title("B_x [T]");
set(gca, 'Ydir', 'normal', 'FontSize', 12, 'TickLength',[0, 0]); 
set(gcf, 'Color', 'w');
axis equal tight;
colormap(turbo); 
clim([-100*1e-6,100*1e-6]);
c= colorbar;
c.TickLength= 0;

clim_min = prctile(By(:), 1);
clim_max = prctile(By(:), 95);
figure();
imagesc(x, y, By); 
xlabel("x [m]"); 
ylabel("y [m]");
title("B_y [T]");
set(gca, 'Ydir', 'normal', 'FontSize', 12, 'TickLength',[0, 0]); 
set(gcf, 'Color', 'w');
axis equal tight;
colormap(turbo); 
clim([-100*1e-6,100*1e-6]);
c= colorbar;
c.TickLength= 0;

clim_min = prctile(Bz(:), 1);
clim_max = prctile(Bz(:), 95);
figure();
imagesc(x, y, Bz); 
xlabel("x [m]"); 
ylabel("y [m]");
title("B_z [T]");
set(gca, 'Ydir', 'normal', 'FontSize', 12, 'TickLength',[0, 0]); 
set(gcf, 'Color', 'w');
axis equal tight;
colormap(turbo); 
clim([-100*1e-6,100*1e-6]);
c= colorbar;
c.TickLength= 0;
