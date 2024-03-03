% Define parameters
nx = 200;            % Number of grid points in the x-direction
ny = 200;            % Number of grid points in the y-direction
nt = 1000;           % Number of time steps
dx = 10.0;           % Grid spacing in the x-direction (meters)
dy = 10.0;           % Grid spacing in the y-direction (meters)
dt = 0.001;          % Time step (seconds)
c0 = 1500.0;         % P-wave velocity (meters/second)
rho = 1000.0;        % Density (kg/m^3)
source_x = nx / 2;   % X-coordinate of the seismic source
source_y = ny / 2;   % Y-coordinate of the seismic source
source_frequency = 10; % Source frequency (Hz)

% Initialize velocity and pressure fields
u = zeros(nx, ny);
v = zeros(nx, ny);
p = zeros(nx, ny);

% Create a simple source waveform (Ricker wavelet)
t = (0:nt-1) * dt;
source_waveform = 2 * pi * source_frequency * (1.0 - 0.25 * (2 * pi * source_frequency * t).^2) .* exp(-0.125 * (2 * pi * source_frequency * t).^2);

% Main time-stepping loop
for it = 1:nt
    % Update velocity components using finite difference scheme
    for i = 2:nx-1
        for j = 2:ny-1
            u(i,j) = u(i,j) - (dt / (rho * dx)) * (p(i,j) - p(i-1,j));
            v(i,j) = v(i,j) - (dt / (rho * dy)) * (p(i,j) - p(i,j-1));
        end
    end
    
    % Inject the source waveform into the pressure field
    p(source_x, source_y) = p(source_x, source_y) + source_waveform(it) / (dx * dy);
    
    % Update pressure field using finite difference scheme
    for i = 2:nx-1
        for j = 2:ny-1
            p(i,j) = p(i,j) - rho * c0^2 * dt * ((u(i+1,j) - u(i,j)) / dx + (v(i,j+1) - v(i,j)) / dy);
        end
    end
    
    % Visualize the wavefield (optional)
    imagesc(p);
    colormap(jet);
    colorbar;
    drawnow;
end

% Plot the final wavefield (optional)
figure;
imagesc(p);
colormap(jet);
colorbar;
title('Final Wavefield');
