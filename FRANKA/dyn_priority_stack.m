function q_dot = dyn_priority_stack(J1,xd1,J2,xd2,Vmin,Vmax,n,m1,alpha)

% When alpha=0: pure Task 1 only (no Task 2 influence at all)
if alpha < 1e-6
    [q_dot,~] = SNS_velocity_p(J1,xd1,Vmin,Vmax,n,m1);
    q_dot = max(Vmin, min(Vmax, q_dot));
    return;
end

% When alpha=1: pure Task 2 priority
if alpha > 1-1e-6
    [qd_t2,~] = SNS_velocity_p(J2,xd2,Vmin,Vmax,n,1);
    P2  = eye(n) - pinv(J2)*J2;
    J1p = J1*P2;
    if rank(J1p,1e-6) >= m1
        q_dot = qd_t2 + pinv(J1p)*(xd1 - J1*qd_t2);
    else
        q_dot = qd_t2;
    end
    q_dot = max(Vmin, min(Vmax, q_dot));
    return;
end

% Intermediate alpha: compute both orderings and blend
% Ordering 1: Task 1 > Task 2
[qd_t1,~] = SNS_velocity_p(J1,xd1,Vmin,Vmax,n,m1);
P1  = eye(n) - pinv(J1)*J1;
J2p = J2*P1;
if rank(J2p,1e-6) >= 1
    qd_12 = qd_t1 + pinv(J2p)*(xd2 - J2*qd_t1);
else
    qd_12 = qd_t1;
end
qd_12 = max(Vmin, min(Vmax, qd_12));

% Ordering 2: Task 2 > Task 1
[qd_t2,~] = SNS_velocity_p(J2,xd2,Vmin,Vmax,n,1);
P2  = eye(n) - pinv(J2)*J2;
J1p = J1*P2;
if rank(J1p,1e-6) >= m1
    qd_21 = qd_t2 + pinv(J1p)*(xd1 - J1*qd_t2);
else
    qd_21 = qd_t2;
end
qd_21 = max(Vmin, min(Vmax, qd_21));

% Blend
q_dot = (1-alpha)*qd_12 + alpha*qd_21;
q_dot = max(Vmin, min(Vmax, q_dot));
end