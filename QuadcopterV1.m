%% Main Function
function res = QuadcopterV1()
% Constants
g = -9.81;                      % Acceleration (m/s^2)
m = 1.2;                        % Mass (kg)
L = 0.024;                      % Length (m)
k = 3e-6;                       % Thrust coefficient
b = 1e-7;                       % Drag coefficient
kd = 0.25;                      % Global drag coefficient
I = [16.64 0.76 2.16;           % Moment of Inertia Matrix
     0.76 20.43 -0.17; 
     2.16 -0.17 30.18];

% Time variables
tstart = 0;                     % Time begining (s)
tend = 10;                      % Time ending (s)
dt = 0.005;                     % Time step (s)
ts = tstart:dt:tend;            % Create time vector
N = numel(ts);                  % Number of points in the simulation

% Initialize inputs
x = [0; 0; 10];                 % Position in 3 dimensitons x,y,z (m)
xdot = zeros(3,1);              % Velocity of quadcopter in 3 dimensions x,y,z (m/s)
theta = zeros(3,1);             % Empty matrix for angular orientation from origin (rad)


% Here is where we'd need to use the helper functions to compute forces,
% torques, and accelerations


[T, M] = ode45(@derivs, [tstart, tend], [x, xdot, theta])

    function res = derivs(t, V)
        x = V(1);
        theta = V(2);
        

    end

res  = 
end

%% Helper Functions

% Compute thrusts
function T = thrust(w, k)

end

% Compute acceleration in inertial frame
function a = acceleration(inputs, angles, vels, m, g, k, kd)

end

% Compute torques
function tau = torques(w, L, b, k)
 tau = [L*k*(w(1)^2 - w(3)^3); 
        L*k*(w(2)^2 - w(4)^2);
        b*(w(1)^2 - w(2)^2 + w(3)^2 -w(4)^2)];
end

% Compute acceleration in body frame
function omegad = angular_acceleration(w, omega, I, L, b, k)
tau = torques(w, L, b, k);
omegad = I\(tao - cross(omega, (I*omega)));

end

% Convert roll, pitch, yaw derivatives to omega
function omega = thetadot2omega(thetadot, angles)
    phi = angles(1);
    theta = angles(2);
    psi = angles(3);
    W = [
        1, 0, -sin(theta)
        0, cos(phi), cos(theta)*sin(phi)
        0, -sin(phi), cos(theta)*cos(phi)
    ];
    omega = W * thetadot;
end

% Convert omega to roll, pitch, yaw derivatives
function thetadot = omega2thetadot(omega, angles)
    phi = angles(1);
    theta = angles(2);
    psi = angles(3);
    W = [
        1, 0, -sin(theta)
        0, cos(phi), cos(theta)*sin(phi)
        0, -sin(phi), cos(theta)*cos(phi)
    ];
    thetadot = W\omega;
end

% Compute rotation matrix for given angles
function R = rotation(angles)
    phi = angles(3);
    theta = angles(2);
    psi = angles(1);

    R = zeros(3);
    R(:, 1) = [
        cos(phi) * cos(theta)
        cos(theta) * sin(phi)
        - sin(theta)
    ];
    R(:, 2) = [
        cos(phi) * sin(theta) * sin(psi) - cos(psi) * sin(phi)
        cos(phi) * cos(psi) + sin(phi) * sin(theta) * sin(psi)
        cos(theta) * sin(psi)
    ];
    R(:, 3) = [
        sin(phi) * sin(psi) + cos(phi) * cos(psi) * sin(theta)
        cos(psi) * sin(phi) * sin(theta) - cos(phi) * sin(psi)
        cos(theta) * cos(psi)
    ];
end

