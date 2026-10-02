function [q_dot, s_star, iterations] = SNS_algorithm(J, x_dot, Vmin, Vmax, n, m)

W = eye(n);
q_dot_N = zeros(n,1);
s = 1;
s_star = 0;
q_dot_star = zeros(n,1);
iterations = 0;
max_iter = n + 1;
limit_exceeded = true;

while limit_exceeded && iterations <= max_iter
    iterations = iterations + 1;
    limit_exceeded = false;
    JW = J * W;
    if rank(JW) < m
        break;
    end
    JW_pinv = pinv(JW);
    q_dot = q_dot_N + JW_pinv * (s * x_dot - J * q_dot_N);
    if all(q_dot >= Vmin - 1e-10) && all(q_dot <= Vmax + 1e-10)
        s_star = s;
        q_dot_star = q_dot;
        break;
    end
    a = JW_pinv * x_dot;
    b = q_dot - s * a;
    s_values = ones(n,1);
    for i = 1:n
        if abs(a(i)) < 1e-10
            if b(i) < Vmin(i)-1e-10 || b(i) > Vmax(i)+1e-10
                s_values(i) = 0;
            end
        else
            s_lo = (Vmin(i)-b(i))/a(i);
            s_hi = (Vmax(i)-b(i))/a(i);
            if s_lo > s_hi; tmp=s_lo; s_lo=s_hi; s_hi=tmp; end
            if s_hi < 0 || s_lo > 1
                s_values(i) = 0;
            else
                s_values(i) = min(1, max(0, s_hi));
            end
        end
    end
    [s_k, j] = min(s_values);
    if s_k > s_star
        s_star = s_k;
        q_dot_star = b + s_k * a;
    end
    W_test = W; W_test(j,j) = 0;
    if rank(J * W_test) < m
        break;
    end
    W(j,j) = 0;
    if q_dot(j) > Vmax(j)
        q_dot_N(j) = Vmax(j);
    else
        q_dot_N(j) = Vmin(j);
    end
    limit_exceeded = true;
end

q_dot = max(Vmin, min(Vmax, q_dot_star));
end