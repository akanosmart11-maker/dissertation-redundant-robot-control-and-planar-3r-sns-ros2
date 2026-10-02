function [q_dot, s_star] = GSNS_algorithm(J, x_dot, ...
    Vmin_j, Vmax_j, J_cp, cp_bmin, cp_bmax, cp_pos, q, L, n, m)

has_cartesian = ~isempty(J_cp);

if has_cartesian
    d_cp = size(J_cp, 1);
    n_constraints = n + d_cp;
    A = [eye(n); J_cp];
    bmin = [Vmin_j; cp_bmin];
    bmax = [Vmax_j; cp_bmax];
else
    n_constraints = n;
    A = eye(n);
    bmin = Vmin_j;
    bmax = Vmax_j;
    d_cp = 0;
end

q_dot_N = zeros(n,1);
s_star = 0;
P = eye(n);
Alim = [];
a_dot_N = [];
q_dot_star = zeros(n,1);
max_iter = n + d_cp;

for iter = 1:max_iter
    JP = J * P;
    if rank(JP, 1e-6) < m; break; end
    JP_pinv = pinv(JP);
    q_dot = q_dot_N + JP_pinv * (x_dot - J * q_dot_N);
    a_dot = A * q_dot;
    violations = (a_dot < bmin - 1e-10) | (a_dot > bmax + 1e-10);
    if ~any(violations)
        s_star = 1;
        q_dot_star = q_dot;
        break;
    end
    alpha = A * JP_pinv * x_dot;
    beta  = a_dot - alpha;
    [s_k, k_crit] = getTaskScalingFactor_GSNS(alpha, beta, bmin, bmax, n_constraints);
    if s_k > s_star
        s_star = s_k;
        q_dot_star = max(Vmin_j, min(Vmax_j, ...
            q_dot_N + JP_pinv * (s_k * x_dot - J * q_dot_N)));
    end
    Ak = A(k_crit,:);
    if a_dot(k_crit) > bmax(k_crit)
        sat_val = bmax(k_crit);
    else
        sat_val = bmin(k_crit);
    end
    Alim    = [Alim;    Ak];
    a_dot_N = [a_dot_N; sat_val];
    P = eye(n) - pinv(Alim) * Alim;
    if rank(J*P, 1e-6) < m
        if size(Alim,1) > 1
            Alim    = Alim(1:end-1,:);
            a_dot_N = a_dot_N(1:end-1);
            P = eye(n) - pinv(Alim) * Alim;
        end
        break;
    end
    q_dot_N = pinv(Alim) * a_dot_N;
end

q_dot = max(Vmin_j, min(Vmax_j, q_dot_star));
end

function [s, k_crit] = getTaskScalingFactor_GSNS(alpha, beta, bmin, bmax, nc)
s_values = ones(nc,1);
for h = 1:nc
    L_h = bmin(h) - beta(h);
    U_h = bmax(h) - beta(h);
    if alpha(h) < 0 && L_h < 0
        if alpha(h) < L_h
            s_values(h) = L_h / alpha(h);
        else
            s_values(h) = 1;
        end
    elseif alpha(h) > 0 && U_h > 0
        if alpha(h) > U_h
            s_values(h) = U_h / alpha(h);
        else
            s_values(h) = 1;
        end
    else
        s_values(h) = 0;
    end
end
[s, k_crit] = min(s_values);
end