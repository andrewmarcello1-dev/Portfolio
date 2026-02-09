clear; close all; clc

%% PARAMETERS
alp = 0;                                                                   % Angle of attach (deg)
M_inf = 0.5;                                                               % Free-stream Mach #
T_inf = 300;                                                               % Free-stream temperature (K)
p_inf = 1e5;                                                               % Free-stream pressure (Pa)
R = 287;                                                                   % Ideal gas constant of air (J/kg*K)
gamma = 1.4;                                                               % Ratio of specific heats for air
a_inf = sqrt(gamma*R*T_inf);                                               % Speed of sound (m/s)
v_inf = M_inf*a_inf;                                                       % Free-stream velocity (m/s)
rho_inf = p_inf/(R*T_inf);                                                 % Free-stream density (kg/m^3)

%% AIRFOIL GENERATION

% x-grid
Nx1 = 50;                                                                  % Number of points in downstream
Nx2 = 41;                                                                  % Number of points on airfoil
Nx3 = 50;                                                                  % Number of points in downstream
gf = 0.8;                                                                  % Growth factor
P = 1+gf;                                                                  % Stretching exponent
delx = 0.025;                                                              % Initial offset from airfoil
dely = 0.025;

x1_min = -30;                                                              % Starting point of grid section 1
x1_max = -delx;                                                            % Ending Point of grid section 1
x1_initial = linspace(0,1,Nx1);                                            % Initialize 
x1 = x1_min+(x1_max-x1_min)*(1-(1-x1_initial).^P);                         % Use power law formula to create non-uniform grid in upstream

x2 = linspace(0, 1, Nx2);                                                  % Uniform grid on airfoil

x3_min = 1+delx;                                                           % Repeat section 1 process for section 3
x3_max = 31;
x3_initial = linspace(0,1,Nx3);
x3 = x3_min+(x3_max-x3_min)*(x3_initial.^P);                               % Non-uniform grid in downstream

x_points = [x1,x2,x3];                                                     % Combine the three sections

% y-grid
Ny = 50;                                                                   % Number of points perpendicular to airfoil
y_min = dely;                                                              % Repeat procession for y used for section 3 in x
y_max = 30;
y_initial = linspace(0,1,Ny);
y = y_min+(y_max-y_min)*(y_initial.^P);                                    
y_points = [0 y];

x_airfoil = linspace(0,1,41);                                              % x-points along airfoil
y_airfoil = sqrt(4.182^2-(x_airfoil-0.5).^2)-4.152;                         % y-points along airfoil
dydx_airfoil = 0.5*(4.182^2-(x_airfoil-0.5).^2).^(-0.5).*(-2*x_airfoil+1);  % Compute dy/dx for airfoil for boundary conditions later

%% DISCRETEZATION
[X, Y] = meshgrid(x_points, y_points);                                     % Create mesh grid
mask1 = (X(1,:) < 0);                                                      % Find where X is less then 0 for boundary condition on symmetry place
mask2 = (X(1,:) >= 0) & (X(1,:) <= 1);                                     % Find where the airfoil is to input wall condition
mask3 = (X(1,:) > 1);                                                      % Repeat for values greater than 1 meter for symmetry plane
Y(1,mask2) = y_airfoil;                                                    % Replace y-points within airfoil to y-point on airfoil surface  
Y(2,mask2) = y_airfoil+0.01;

%% INITIALIZE PHI SOLUTION
phi = zeros(size(X));                                                      % Initialize solution matrix

y_plus = Y(3:end,:)-Y(2:end-1,:);                                          % Determine delta y constants for every point
y_minus = Y(2:end-1,:)-Y(1:end-2,:);
y0 = (Y(3:end,:)-Y(1:end-2,:))/2;

x_plus = X(:,3:end)-X(:,2:end-1);                                          % Determine delta x constants for every point
x_minus = X(:,2:end-1)-X(:,1:end-2);
x0 = (X(:,3:end)-X(:,1:end-2))/2;

Y1 = Y(2,:)-Y(1,:);
y_plus_airfoil = Y1(mask2);

phi(:,1) = v_inf*X(:,1);                                                   % Initialize upstream boundary condition
phi(:,end) = v_inf*X(:,end);                                               % Initialize downstream boundary condition
phi(end,:) = v_inf*X(end,:);                                               % Initialize upper surface of domain with boundary condition

phi(1,mask1) = phi(2,mask1);                                               % d(phi)/dy = 0
phi(1,mask2) = phi(2,mask2)-y_plus_airfoil*v_inf.*dydx_airfoil;                       % Apply boundary condition at airfoil
phi(1,mask3) = phi(2,mask3);                                               % d(phi)/dy = 0

%% PHI SOLUTION
B = 1-M_inf^2;                                                             % Beta constant for efficieny

d = zeros(size(phi));                                                      % Initialize variables for GS Point Relaxation Method
b = zeros(size(phi));
f = zeros(size(phi));
a = zeros(size(phi));
g = zeros(size(phi));

for i = 1:length(d(1,:))-2                                                 % For loop to determine point relaxation constants at every internal point
    for j = 1:length(d(:,1))-2
        d(j+1,i+1) = B*y0(j,i)/x_plus(j,i)+B*y0(j,i)/x_minus(j,i)+x0(j,i)/y_plus(j,i)+x0(j,i)/y_minus(j,i);
        b(j+1,i+1) = -B*y0(j,i)/x_minus(j,i);
        a(j+1,i+1) = -B*y0(j,i)/x_plus(j,i);
        f(j+1,i+1) = -x0(j,i)/y_minus(j,i);
        g(j+1,i+1) = -x0(j,i)/y_plus(j,i);
    end
end

res = 1;                                                                   % Initialize residual error
n = 0;                                                                     % Initialize iteration count
eps = 1e-8;                                                                % Residual error tolerance

residuals = [];                                                            % Store residuals
iterations = [];                                                           % Store iteration count

while res > eps

    n = n+1;                                                               % Track iterations
    phi_old = phi;                                                         % Store old solution for residuals

    for i = 2:length(d(1,:))-1                                             % For loop to preform point relaxation method at every internal point
        for j = 2:length(d(:,1))-1
            phi(j,i) = (phi(j,i-1)*-b(j,i) + phi_old(j,i+1)*-a(j,i)+phi(j-1,i)*-f(j,i)+phi_old(j+1,i)*-g(j,i))/d(j,i);
        end
    end
    
    phi(:,1) = v_inf*X(:,1);                                               % Initialize upstream boundary condition
    phi(:,end) = v_inf*X(:,end);                                           % Initialize downstream boundary condition
    phi(end,:) = v_inf*X(end,:);                                           % Initialize upper surface of domain with boundary condition
    
    phi(1,mask1) = phi(2,mask1);                                           % d(phi)/dy = 0
    phi(1,mask2) = phi(2,mask2)-y_plus_airfoil*v_inf.*dydx_airfoil;        % Apply boundary condition at airfoil
    phi(1,mask3) = phi(2,mask3);                                           % d(phi)/dy = 0

    residuals = [residuals res];                                           % Store resitual error for plotting
    iterations = [iterations n];                                           % Store iterations for plotting

    res = norm(phi - phi_old,'fro') / norm(phi_old,'fro');                 % Calculate residual error frobenius norm for computational efficiency
    fprintf('Iteration: %d, Residual: %e\n', n, res);                      % Display iteration and residual error
end

%% VELOCITY
u = zeros(size(phi));                                                      % Initialize velocity matrices
v = zeros(size(phi));

u(:,2:end-1) = (phi(:,3:end) - phi(:,1:end-2))./(2*x0);                    % Central differencing for internal u
u(:,1) = (phi(:,2) - phi(:,1)) / (X(1,2)-X(1,1));                          % Forward differencing for upstream
u(:,end) = (phi(:,end) - phi(:,end-1)) / (X(1,end)-X(1,end-1));            % Backward differencing for downstream

v(2:end-1,:) = (phi(3:end,:) - phi(1:end-2,:))./(2*y0);                    % Central differencing for internal v
v(1,:) = (phi(2,:) - phi(1,:)) ./ (Y(2,:)-Y(1,:));                         % Forward differencing for bottom of domain
v(end,:) = (phi(end,:) - phi(end-1,:)) ./ (Y(end,:)-Y(end-1,:));           % Backward differencing for top of domain

U = (u.^2+v.^2).^(1/2);                                                    % Calculate norm of u and v to find magnitude of velocity

%% COEFFICIENT OF PRESSURE
p = p_inf*(1-(gamma-1)/2*M_inf^2*(U.^2/v_inf^2-1)).^(gamma/(gamma-1));     % Calculate pressure
cp = 2*(p-p_inf)/(gamma*p_inf*M_inf^2);                                    % Calculate coefficient of pressure

%% TEMPERATURE & MACH NUMBER
To = T_inf*(1+(gamma-1)/2*M_inf^2);                                        % Stagnation Temperature (K)
Cp = 1000;                                                                 % Specfic heat capacity of air (J/kg*K)
T = To - 1/(2*Cp)*U.^2;                                                    % Calculate Temperature (K)
Mach = U./sqrt(gamma*R*T);                                                 % Calculate Mach number based on speed of sound at each point

%% PLOTS

% Airfoil with Grid
figure;
plot(X, Y, 'k.', 'MarkerSize', 1);
hold on;
plot(x_airfoil, y_airfoil, 'r', 'LineWidth', 2);
hold off;

axis equal;
xlabel('x');
ylabel('y');
title('Uniform Grid Around Airfoil with Solid Wall');
grid on;

% Iteration vs. Residual
figure;
semilogy(iterations, residuals, 'k-');
grid on;
xlabel('Iteration');
ylabel('Residual Error');
title('Convergence History of the Solver');

% Pressure coefficient along airfoil
figure;
plot(x_airfoil,cp(1,mask2),'o-r')
grid on;
xlabel('Chord (m)')
ylabel('Pressure Coefficient on Airfoil')

% Mach Number Contours
figure;
contourf(X, Y, Mach, 20, 'LineStyle', 'none');
colorbar;
xlabel('x (m)');
ylabel('y (m)');
title('Mach Number Contours');

% Temperature Contours
figure;
contourf(X, Y, T, 20, 'LineStyle', 'none');
colorbar;
xlabel('x (m)');
ylabel('y (m)');
title('Temperature Contours (K)');

% Pressure Contours
figure;
contourf(X, Y, abs(p), 20, 'LineStyle', 'none');
colorbar;
xlabel('x (m)');
ylabel('y (m)');
title('Pressure Contours (Pa)');
