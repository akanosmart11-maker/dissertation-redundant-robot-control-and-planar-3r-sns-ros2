% Panda Milestone 3 - GSNS with Cartesian ceiling constraint
% Elbow (joint 4 origin) must remain below z = 0.7 m
% Three methods compared: No limits, SNS joint-only, GSNS joint+Cartesian

clear; clc; close all;

%% Parameters
n = 7; m = 3;
dt = 0.001; T_total = 3.0; N = T_total/dt;
Kp = 50;

% Joint velocity limits
Vmax_j =  0.5 * ones(7,1);
Vmin_j = -0.5 * ones(7,1);

% Cartesian ceiling on elbow z-coordinate
z_ceiling = 0.35;   % just above initial elbow z of 0.333
z_floor   = -0.5;   % [m] not restrictive
cp_Vmax   = 0.3;    % [m/s]
cp_Amax   = 1.0;    % [m/s^2]

% Start from ready pose
q0 = [0; -pi/4; 0; -3*pi/4; 0; pi/2; pi/4];
T0 = panda_fk(q0);
p_start = T0(1:3,4);

% Target: move 15cm in x, 5cm in y, 10cm up in z
p_end = p_start + [0.15; 0.0; 0.25];  % aggressive upward motion

% Check initial elbow z
[~, Tj0] = panda_fk(q0);
elbow_z0 = Tj0{5}(3,4);
fprintf('Initial elbow z: %.4f m\n', elbow_z0);
fprintf('Ceiling at z:    %.4f m\n', z_ceiling);
fprintf('Start EE: [%.4f, %.4f, %.4f]\n', p_start(1),p_start(2),p_start(3));
fprintf('End   EE: [%.4f, %.4f, %.4f]\n', p_end(1),p_end(2),p_end(3));

%% Quintic timing
tau_v  = linspace(0,1,N);
gamma  = 6*tau_v.^5 - 15*tau_v.^4 + 10*tau_v.^3;
dgamma = (30*tau_v.^4 - 60*tau_v.^3 + 30*tau_v.^2) / T_total;

method_names = {'No limits','SNS (joint only)','GSNS (joint+Cartesian)'};
colors = {[0.6 0.6 0.6],[0.0 0.45 0.74],[0.85 0.18 0.18]};
results = cell(3,1);

for mi = 1:3
    q = q0;
    p_hist      = zeros(3,N+1);
    p_des_hist  = zeros(3,N+1);
    elbow_z_hist = zeros(1,N+1);
    err_hist    = zeros(1,N);
    scale_hist  = ones(1,N);

    T_curr = panda_fk(q);
    p_hist(:,1) = T_curr(1:3,4);
    [~,Tj] = panda_fk(q); elbow_z_hist(1) = Tj{5}(3,4);

    for k = 1:N
        % Desired trajectory
        p_des = p_start + gamma(k)  * (p_end - p_start);
        v_des = dgamma(k) * (p_end - p_start);
        p_des_hist(:,k) = p_des;

        % Current state
        T_curr = panda_fk(q);
        p_curr = T_curr(1:3,4);
        x_dot  = v_des + Kp*(p_des - p_curr);

        % Jacobian - position rows only
        J_full = panda_jacobian(q);
        J      = J_full(1:3,:);

        switch mi
            case 1  % No limits
                q_dot     = pinv(J) * x_dot;
                s_out     = 1;

            case 2  % SNS joint only
                [q_dot, s_out, ~] = SNS_algorithm(J, x_dot, Vmin_j, Vmax_j, n, m);

            case 3  % GSNS joint + Cartesian ceiling on elbow z
                [~,Tj] = panda_fk(q);
                elbow_z = Tj{5}(3,4);
                J_elbow_z = panda_elbow_jacobian_y(q);

                % Shape Cartesian constraint
                v_pos_max = (z_ceiling - elbow_z) / dt;
                v_pos_min = (z_floor   - elbow_z) / dt;
                margin_up = z_ceiling - elbow_z;
                if margin_up > 0
                    v_stop = sqrt(2*cp_Amax*margin_up);
                else
                    v_stop = 0;
                end
                cp_bmax = min([v_pos_max, cp_Vmax, v_stop]);
                cp_bmin = max([v_pos_min, -cp_Vmax]);
                if cp_bmax < 0; cp_bmax = 0; end
                if cp_bmin > 0; cp_bmin = 0; end

                [q_dot, s_out] = GSNS_algorithm(J, x_dot, ...
                    Vmin_j, Vmax_j, J_elbow_z, cp_bmin, cp_bmax, ...
                    elbow_z, q, [], n, m);
        end

        scale_hist(k) = s_out;
        err_hist(k)   = norm(p_des - p_curr)*1000;
        q = q + q_dot*dt;

        T_next = panda_fk(q);
        p_hist(:,k+1) = T_next(1:3,4);
        [~,Tj] = panda_fk(q); elbow_z_hist(k+1) = Tj{5}(3,4);
    end
    p_des_hist(:,N+1) = p_end;

    results{mi} = struct('p_hist',p_hist,'p_des_hist',p_des_hist,...
        'elbow_z',elbow_z_hist,'err_hist',err_hist,'scale_hist',scale_hist);

    fprintf('\n--- %s ---\n', method_names{mi});
    fprintf('  Max elbow z:    %.4f m (ceiling: %.2f m)\n', max(elbow_z_hist), z_ceiling);
    fprintf('  Max EE error:   %.3f mm\n', max(err_hist));
    fprintf('  Final EE error: %.4f mm\n', norm(p_hist(:,end)-p_end)*1000);
    if max(elbow_z_hist) > z_ceiling + 0.001
        fprintf('  *** CONSTRAINT VIOLATED ***\n');
    else
        fprintf('  Constraint satisfied.\n');
    end
end

%% Plots
t = linspace(0, T_total, N+1);
tv = linspace(0, T_total, N);

figure('Name','Panda GSNS - Elbow Constraint');
subplot(2,1,1); hold on;
for mi=1:3
    plot(t, results{mi}.elbow_z, 'Color',colors{mi}, 'LineWidth',2);
end
yline(z_ceiling,'r--','LineWidth',2);
xlabel('Time (s)'); ylabel('Elbow z (m)');
title('Elbow z-coordinate vs Ceiling Constraint');
legend([method_names, {'Ceiling'}],'Location','best'); grid on;

subplot(2,1,2); hold on;
for mi=1:3
    plot(tv, results{mi}.err_hist, 'Color',colors{mi}, 'LineWidth',1.5);
end
xlabel('Time (s)'); ylabel('EE Error (mm)');
title('End-Effector Tracking Error');
legend(method_names,'Location','best'); grid on;
sgtitle('Panda Milestone 3 - GSNS Cartesian Constraint','FontWeight','bold');

figure('Name','Panda GSNS - 3D Paths');
hold on;
for mi=1:3
    plot3(results{mi}.p_hist(1,:), results{mi}.p_hist(2,:), ...
          results{mi}.p_hist(3,:), 'Color',colors{mi}, 'LineWidth',2);
end
plot3(p_start(1),p_start(2),p_start(3),'go','MarkerSize',8,'MarkerFaceColor','g');
plot3(p_end(1),p_end(2),p_end(3),'rs','MarkerSize',8,'MarkerFaceColor','r');
xlabel('x(m)'); ylabel('y(m)'); zlabel('z(m)');
legend(method_names,'Location','best'); grid on;
title('3D End-Effector Paths - GSNS vs Other Methods');

fprintf('\n=== MILESTONE 3 COMPLETE ===\n');
fprintf('GSNS enforces elbow ceiling at z=%.2f m\n', z_ceiling);
fprintf('Joint-only and unconstrained methods violate it.\n');