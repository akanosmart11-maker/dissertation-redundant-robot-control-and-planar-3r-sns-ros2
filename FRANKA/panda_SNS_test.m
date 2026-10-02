% Panda Milestone 2 - SNS algorithm on 7-DOF Franka
% Direct analogue of planar Milestone 2
% Task: 3D straight-line end-effector translation (position only, m=3)
% Redundancy: n-m = 7-3 = 4 degrees

clear all; clc; close all;

%% Setup
T_total = 3.0;       % trajectory duration (s)
dt      = 0.001;     % sampling time (s)
N       = T_total/dt;
Kp      = 50;        % proportional gain

% Joint velocity limits (tight, to force saturations)
Vmax =  0.5 * ones(7,1);
Vmin = -0.5 * ones(7,1);

% Start configuration - Panda ready pose
q0 = [0; -pi/4; 0; -3*pi/4; 0; pi/2; pi/4];

% Get start and end EE positions
T0   = panda_fk(q0);
p_start = T0(1:3,4);
p_end   = p_start + [0.1; 0.0; -0.1];  % 10cm in x, 10cm down in z

fprintf('Start EE: [%.4f, %.4f, %.4f]\n', p_start(1),p_start(2),p_start(3));
fprintf('End   EE: [%.4f, %.4f, %.4f]\n', p_end(1),p_end(2),p_end(3));

%% Quintic timing
tau_vec = linspace(0, 1, N);
gamma   = 6*tau_vec.^5 - 15*tau_vec.^4 + 10*tau_vec.^3;
dgamma  = (1/T_total) * (30*tau_vec.^4 - 60*tau_vec.^3 + 30*tau_vec.^2);

%% Storage
q_hist      = zeros(7, N);
p_hist      = zeros(3, N);
p_des_hist  = zeros(3, N);
err_hist    = zeros(1, N);
iter_hist   = zeros(1, N);
q_hist(:,1) = q0;

%% Main control loop
q = q0;
for k = 1:N-1
    % Desired position and velocity
    p_des  = p_start + gamma(k)   * (p_end - p_start);
    v_des  = dgamma(k) * (p_end - p_start);

    % Current EE position
    T_curr = panda_fk(q);
    p_curr = T_curr(1:3,4);

    % Task velocity with feedback (position task only, m=3)
    x_dot = v_des + Kp*(p_des - p_curr);

    % Jacobian - use position rows only (top 3 rows)
    J_full = panda_jacobian(q);
    J      = J_full(1:3, :);    % 3x7 position Jacobian

    % SNS algorithm
    [q_dot, s, iter] = SNS_algorithm(J, x_dot, Vmin, Vmax, 7, 3);

    % Integrate
    q = q + q_dot * dt;

    % Store
    q_hist(:,k+1)   = q;
    p_hist(:,k)     = p_curr;
    p_des_hist(:,k) = p_des;
    err_hist(k)     = norm(p_curr - p_des)*1000;  % mm
    iter_hist(k)    = iter;
end

% Final step
T_final = panda_fk(q);
p_hist(:,N)     = T_final(1:3,4);
p_des_hist(:,N) = p_end;
err_hist(N)     = norm(T_final(1:3,4) - p_end)*1000;

%% Results
fprintf('\n=== PANDA SNS RESULTS ===\n');
fprintf('Max EE error:    %.3f mm\n', max(err_hist));
fprintf('Final EE error:  %.4f mm\n', err_hist(end));
fprintf('Max iterations:  %d\n',      max(iter_hist));
fprintf('Mean iterations: %.2f\n',    mean(iter_hist));

%% Plots
t = linspace(0, T_total, N);
figure('Name','Panda SNS - EE Tracking');

subplot(2,2,1);
plot3(p_des_hist(1,:), p_des_hist(2,:), p_des_hist(3,:), 'b--', ...
      p_hist(1,:),     p_hist(2,:),     p_hist(3,:),     'r-', 'LineWidth',1.5);
legend('Desired','Actual'); grid on; xlabel('x(m)'); ylabel('y(m)'); zlabel('z(m)');
title('3D End-Effector Path');

subplot(2,2,2);
plot(t, err_hist, 'r', 'LineWidth',1.5);
xlabel('Time (s)'); ylabel('Error (mm)'); title('EE Tracking Error'); grid on;

subplot(2,2,3);
plot(t, q_hist', 'LineWidth',1.2);
xlabel('Time (s)'); ylabel('Joint angle (rad)'); title('Joint Trajectories'); grid on;
legend('q1','q2','q3','q4','q5','q6','q7');

subplot(2,2,4);
plot(t, iter_hist, 'k', 'LineWidth',1.5);
xlabel('Time (s)'); ylabel('Iterations'); title('SNS Iterations per Step'); grid on;

sgtitle('Franka Panda - SNS Algorithm Validation');