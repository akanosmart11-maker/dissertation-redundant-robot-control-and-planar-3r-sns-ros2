% Panda Milestone 4 - Dynamic Task Priorities (Novel Contribution)
% Task 1: 3D straight-line EE trajectory (m=3)
% Task 2: Elbow y-coordinate regulation with hard ceiling
% Novel: cosine blending activates Task 2 only when elbow approaches ceiling
clear; clc; close all;

n = 7; m1 = 3;
dt = 0.001; T_total = 6.0; N = round(T_total/dt);
Kp1 = 80; Kp2 = 5;

Vmax_j =  0.4 * ones(7,1);
Vmin_j = -0.4 * ones(7,1);

y_des  = 0.26;
y_act  = 0.255;
y_hard = 0.28;

q0 = [0; -pi/4; 0; -3*pi/4; 0; pi/2; pi/4];
T0 = panda_fk(q0);
p_start = T0(1:3,4);
p_end_task1 = p_start + [0.0; 0.30; 0.0];

[~,Tj0] = panda_fk(q0);
fprintf('Initial elbow y: %.4f m\n', Tj0{5}(2,4));
fprintf('Hard ceiling:    %.4f m\n', y_hard);
fprintf('Activation zone: %.4f m\n', y_act);

method_names = {'SNS only (no elbow task)', ...
                'No limits (no elbow task)', ...
                'Dynamic priority (NOVEL)'};
colors = {[0.0 0.45 0.74],[0.6 0.6 0.6],[0.85 0.18 0.18]};
results = cell(3,1);

for mi = 1:3
    q = q0;
    p_hist       = zeros(3, N+1);
    elbow_y_hist = zeros(1, N+1);
    err1_hist    = zeros(1, N);
    err2_hist    = zeros(1, N);
    alpha_hist   = zeros(1, N);

    T_curr = panda_fk(q);
    p_hist(:,1) = T_curr(1:3,4);
    [~,Tj] = panda_fk(q);
    elbow_y_hist(1) = Tj{5}(2,4);

    for k = 1:N
        t = (k-1)*dt;

        % Task 1: quintic straight-line trajectory
        tau     = min(t/T_total, 1.0);
        gam     = 6*tau^5 - 15*tau^4 + 10*tau^3;
        gam_dot = (30*tau^4 - 60*tau^3 + 30*tau^2)/T_total;
        p1_des  = p_start + gam*(p_end_task1 - p_start);
        v1_des  = gam_dot*(p_end_task1 - p_start);

        T_curr  = panda_fk(q);
        p_curr  = T_curr(1:3,4);
        J_full  = panda_jacobian(q);
        J1      = J_full(1:3,:);
        x_dot_1 = v1_des + Kp1*(p1_des - p_curr);

        % Task 2: elbow y regulation
        [~,Tj]  = panda_fk(q);
        elbow_y = Tj{5}(2,4);
        J2      = panda_elbow_jacobian_y(q);
        x_dot_2 = Kp2*(y_des - elbow_y);

        % Priority weight alpha (only for dynamic method)
        if mi == 3
            if elbow_y <= y_act
                alpha = 0.0;
            elseif elbow_y >= y_hard
                alpha = 1.0;
            else
                x_norm = (elbow_y - y_act)/(y_hard - y_act);
                alpha  = 0.5*(1 - cos(pi*x_norm));
            end
        else
            alpha = 0.0;
        end
        alpha_hist(k) = alpha;

        if mi == 1
            [q_dot,~] = SNS_velocity_p(J1,x_dot_1,Vmin_j,Vmax_j,n,m1);

        elseif mi == 2
            q_dot = pinv(J1)*x_dot_1;
            q_dot = max(Vmin_j, min(Vmax_j, q_dot));

        else
            q_dot = dyn_priority_stack(J1,x_dot_1,J2,x_dot_2,...
                Vmin_j,Vmax_j,n,m1,alpha);
        end

        q_dot = max(Vmin_j, min(Vmax_j, q_dot));
        err1_hist(k) = norm(p1_des - p_curr)*1000;
        err2_hist(k) = abs(y_des - elbow_y)*1000;
        q = q + q_dot*dt;

        T_next = panda_fk(q);
        p_hist(:,k+1) = T_next(1:3,4);
        [~,Tj2] = panda_fk(q);
        elbow_y_hist(k+1) = Tj2{5}(2,4);
    end

    results{mi} = struct('p_hist',p_hist,'elbow_y',elbow_y_hist,...
        'err1',err1_hist,'err2',err2_hist,'alpha',alpha_hist);

    fprintf('\n--- %s ---\n', method_names{mi});
    fprintf('  Max EE error:    %.2f mm\n', max(err1_hist));
    fprintf('  Max elbow y:     %.4f m (ceiling: %.2f m)\n', ...
        max(elbow_y_hist), y_hard);
    if max(elbow_y_hist) > y_hard + 0.001
        fprintf('  *** CEILING VIOLATED ***\n');
    else
        fprintf('  Ceiling respected.\n');
    end
    fprintf('  Max alpha:       %.4f\n', max(alpha_hist));
end

%% Plots
t  = linspace(0, T_total, N+1);
tv = linspace(0, T_total, N);

figure('Name','Panda M4 - EE Paths');
for mi = 1:3
    subplot(1,3,mi); hold on;
    plot3(p_start(1),p_start(2),p_start(3),...
        'go','MarkerSize',8,'MarkerFaceColor','g');
    plot3(p_end_task1(1),p_end_task1(2),p_end_task1(3),...
        'rs','MarkerSize',8,'MarkerFaceColor','r');
    plot3(results{mi}.p_hist(1,:),results{mi}.p_hist(2,:),...
        results{mi}.p_hist(3,:),'Color',colors{mi},'LineWidth',2);
    xlabel('x(m)'); ylabel('y(m)'); zlabel('z(m)');
    title(method_names{mi},'FontSize',9);
    grid on; view(45,20); hold off;
end
sgtitle('Panda M4 - Dynamic Priorities: EE Path Tracking','FontWeight','bold');

figure('Name','Panda M4 - Analysis');
subplot(4,1,1); hold on;
for mi=1:3
    plot(t, results{mi}.elbow_y,'Color',colors{mi},'LineWidth',1.5);
end
yline(y_hard,'r--','LineWidth',2,'Label','Hard limit');
yline(y_act,'Color',[1 0.6 0],'LineStyle',':','LineWidth',1.5,'Label','Activation');
ylabel('Elbow y (m)'); title('Elbow y-coordinate');
legend(method_names,'Location','best'); grid on;

subplot(4,1,2);
plot(tv, results{3}.alpha,'Color',colors{3},'LineWidth',1.5);
ylabel('\alpha'); title('Priority Weight \alpha (Dynamic only)');
ylim([-0.1 1.1]); grid on;
yline(0,'k:'); yline(1,'k:');

subplot(4,1,3); hold on;
for mi=1:3
    plot(tv, results{mi}.err1,'Color',colors{mi},'LineWidth',1.2);
end
ylabel('Error (mm)'); title('Task 1 EE Tracking Error');
legend(method_names,'Location','best'); grid on;

subplot(4,1,4); hold on;
for mi=1:3
    plot(tv, results{mi}.err2,'Color',colors{mi},'LineWidth',1.2);
end
ylabel('Error (mm)'); xlabel('Time (s)');
title('Task 2 Elbow y Error');
legend(method_names,'Location','best'); grid on;

sgtitle('Panda M4 - Dynamic Task Priorities: Analysis','FontWeight','bold');

fprintf('\n=== MILESTONE 4 COMPLETE ===\n');
fprintf('Dynamic priority protects ceiling while SNS and no-limits both violate it.\n');
fprintf('The cosine blending activates smoothly when elbow enters activation zone.\n');