clear; close all; clc

%% ShockTube_VanLeer.m
% This MATLAB script solves the 1D Euler equations for a shock tube problem
% using a Van Leer–style flux splitting scheme. It then animates the evolution
% of density and pressure profiles over the simulation time.

clear; clc; close all;

%% Problem Setup
gamma = 1.4;        % Ratio of specific heats
L = 2.0;            % Length of the shock tube in meters
nx = 200;           % Number of computational cells
dx = L / nx;        % Spatial resolution (cell width)
CFL = 0.5;          % CFL number for stability
t_final = 0.25;     % Final simulation time (seconds)

% Create a grid (cell centers)
x = linspace(dx/2, L - dx/2, nx)';

%% Initial Conditions
% Left state (x < 1 m): high density and pressure; Right state: low density and pressure.
rho = zeros(nx,1);
p   = zeros(nx,1);
u   = zeros(nx,1);  % Initial velocity is zero everywhere

rho(x < 1) = 1.4;
rho(x >= 1) = 2.8;
p(x < 1)   = 1;
p(x >= 1)  = 2;

% Initial energy per cell: E = p/(gamma-1) + 0.5*rho*u^2.
E = p/(gamma-1) + 0.5*rho.*u.^2;

% Conserved variable vector: U = [density, momentum, energy]
U = [rho, rho.*u, E];

%% Time Integration Setup
t = 0;
U_history = {};     % To store snapshots for animation
time_history = [];  % Time history of the solution

% Main time-stepping loop
while t < t_final
    % Extend U with ghost cells using simple boundary extrapolation
    U_ext = zeros(nx+2, 3);
    U_ext(2:end-1,:) = U;
    U_ext(1,:) = U(1,:);       % Left boundary
    U_ext(end,:) = U(end,:);   % Right boundary
    
    % Compute primitive variables from extended U
    rho_ext = U_ext(:,1);
    u_ext   = U_ext(:,2) ./ rho_ext;
    E_ext   = U_ext(:,3);
    p_ext   = (gamma-1) * (E_ext - 0.5*rho_ext.*u_ext.^2);
    
    % Compute speed of sound in the extended domain and determine dt by CFL condition
    c_ext = sqrt(gamma * p_ext ./ rho_ext);
    dt = CFL * dx / max(abs(u_ext) + c_ext);
    if t + dt > t_final
        dt = t_final - t;
    end
    
    % Compute the flux splitting for the extended domain
    [F_plus, F_minus] = van_leer_flux_splitting(rho_ext, u_ext, p_ext, E_ext, gamma);
    
    % Compute numerical fluxes at cell interfaces:
    % flux_{i+1/2} = F_plus(i) + F_minus(i+1), for i = 1,...,nx+1.
    flux = F_plus(1:end-1,:) + F_minus(2:end,:);
    
    % Finite volume update:
    % U^{n+1}(i) = U^n(i) - (dt/dx)*(flux_{i+1/2} - flux_{i-1/2})
    U_new = U - (dt/dx) * (flux(2:end,:) - flux(1:end-1,:));
    U = U_new;
    
    % Update time and record snapshots
    t = t + dt;
    time_history(end+1,1) = t;
    U_history{end+1} = U;
end

%% Post-Processing: Extract Density and Pressure at Each Time Step
num_frames = length(U_history);
density_frames = zeros(nx, num_frames);
pressure_frames = zeros(nx, num_frames);
for k = 1:num_frames
    U_snap = U_history{k};
    rho_snap = U_snap(:,1);
    u_snap   = U_snap(:,2) ./ rho_snap;
    E_snap   = U_snap(:,3);
    p_snap   = (gamma-1) * (E_snap - 0.5*rho_snap.*u_snap.^2);
    
    density_frames(:, k) = rho_snap;
    pressure_frames(:, k) = p_snap;
end

%% Animation: Plot Density and Pressure Evolution
figure;
subplot(2,1,1);
h1 = plot(x, density_frames(:,1), 'b-', 'LineWidth', 2);
xlabel('x (m)');
ylabel('Density');
title(sprintf('Density at t = %.3f s', time_history(1)));
ylim([0, 1.2*max(density_frames(:,1))]);

subplot(2,1,2);
h2 = plot(x, pressure_frames(:,1), 'r-', 'LineWidth', 2);
xlabel('x (m)');
ylabel('Pressure');
title(sprintf('Pressure at t = %.3f s', time_history(1)));
ylim([0, 1.2*max(pressure_frames(:,1))]);

for k = 1:num_frames
    % Update the plots with the current frame data
    set(h1, 'YData', density_frames(:,k));
    set(h2, 'YData', pressure_frames(:,k));
    
    subplot(2,1,1);
    title(sprintf('Density at t = %.3f s', time_history(k)));
    subplot(2,1,2);
    title(sprintf('Pressure at t = %.3f s', time_history(k)));
    
    drawnow;
    pause(0.05);  % Adjust pause to control animation speed
end

%% Function Definition: Van Leer Flux Splitting
function [F_plus, F_minus] = van_leer_flux_splitting(rho, u, p, E, gamma)
    % Computes the positive and negative flux splitting based on the Van Leer method.
    % Input:
    %   rho  - Density vector
    %   u    - Velocity vector
    %   p    - Pressure vector
    %   E    - Energy vector
    %   gamma- Ratio of specific heats
    % Output:
    %   F_plus, F_minus - Split fluxes with shape (n x 3) for each conserved variable.
    
    % Calculate the local speed of sound and Mach number.
    c = sqrt(gamma * p ./ rho);
    M = u ./ c;
    
    % Initialize splitting factors.
    a_plus = zeros(size(M));
    a_minus = zeros(size(M));
    b_plus = zeros(size(M));
    b_minus = zeros(size(M));
    
    % For |M| <= 1.
    mask = abs(M) <= 1;
    a_plus(mask) = 0.25 * (M(mask) + 1).^2;
    a_minus(mask) = -0.25 * (M(mask) - 1).^2;
    b_plus(mask) = 0.25 * (M(mask) + 1).^2 .* (2 - M(mask));
    b_minus(mask) = 0.25 * (M(mask) - 1).^2 .* (2 + M(mask));
    
    % For M > 1.
    mask = M > 1;
    a_plus(mask) = M(mask);
    a_minus(mask) = 0;
    b_plus(mask) = 0.5 * (M(mask) + 1);
    b_minus(mask) = 0;
    
    % For M < -1.
    mask = M < -1;
    a_plus(mask) = 0;
    a_minus(mask) = M(mask);
    b_plus(mask) = 0;
    b_minus(mask) = 0.5 * (1 - M(mask));
    
    % Compute the positive flux splitting components.
    % These expressions follow from one common formulation for the Euler equations.
    F_plus = zeros(length(rho), 3);
    F_plus(:,1) = rho .* c .* a_plus;
    F_plus(:,2) = rho .* u .* c .* a_plus + p .* b_plus;
    F_plus(:,3) = (E + p) .* c .* a_plus + p .* u .* b_plus;
    
    % Compute the negative flux splitting components.
    F_minus = zeros(length(rho), 3);
    F_minus(:,1) = rho .* c .* a_minus;
    F_minus(:,2) = rho .* u .* c .* a_minus + p .* b_minus;
    F_minus(:,3) = (E + p) .* c .* a_minus + p .* u .* b_minus;
end
