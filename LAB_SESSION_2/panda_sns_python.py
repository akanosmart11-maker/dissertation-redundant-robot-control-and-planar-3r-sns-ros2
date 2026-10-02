import numpy as np
from panda_fk_python import panda_jacobian

def sns_algorithm(J, x_dot, Vmin, Vmax):
    """
    SNS Algorithm - Flacco et al. 2015
    Input:
        J     - m x n Jacobian (use position rows only: J[:3,:])
        x_dot - m x 1 desired task velocity
        Vmin  - n x 1 joint velocity lower bounds
        Vmax  - n x 1 joint velocity upper bounds
    Output:
        q_dot - n x 1 joint velocity command
        s     - task scaling factor
    """
    n = J.shape[1]
    m = J.shape[0]

    W = np.eye(n)
    q_dot_N = np.zeros(n)
    s_star = 0.0
    q_dot_star = np.zeros(n)

    for _ in range(n + 1):
        JW = J @ W
        if np.linalg.matrix_rank(JW) < m:
            break

        JW_pinv = np.linalg.pinv(JW)
        q_dot = q_dot_N + JW_pinv @ (x_dot - J @ q_dot_N)

        # Check if within bounds
        if np.all(q_dot >= Vmin - 1e-10) and np.all(q_dot <= Vmax + 1e-10):
            s_star = 1.0
            q_dot_star = q_dot
            break

        # Find scaling factor
        a = JW_pinv @ x_dot
        b = q_dot - a
        s_vals = np.ones(n)

        for i in range(n):
            if abs(a[i]) < 1e-10:
                if b[i] < Vmin[i] - 1e-10 or b[i] > Vmax[i] + 1e-10:
                    s_vals[i] = 0.0
            else:
                sL = (Vmin[i] - b[i]) / a[i]
                sU = (Vmax[i] - b[i]) / a[i]
                if sL > sU:
                    sL, sU = sU, sL
                if sU < 0 or sL > 1:
                    s_vals[i] = 0.0
                else:
                    s_vals[i] = min(1.0, max(0.0, sU))

        j_crit = np.argmin(s_vals)
        s_k = s_vals[j_crit]

        if s_k > s_star:
            s_star = s_k
            q_dot_star = np.clip(
                q_dot_N + JW_pinv @ (s_k * x_dot - J @ q_dot_N),
                Vmin, Vmax)

        W_test = W.copy()
        W_test[j_crit, j_crit] = 0
        if np.linalg.matrix_rank(J @ W_test) < m:
            break

        W[j_crit, j_crit] = 0
        if q_dot[j_crit] > Vmax[j_crit]:
            q_dot_N[j_crit] = Vmax[j_crit]
        else:
            q_dot_N[j_crit] = Vmin[j_crit]

    return np.clip(q_dot_star, Vmin, Vmax), s_star


if __name__ == "__main__":
    # Quick test at ready pose
    import numpy as np
    q = [0, -np.pi/4, 0, -3*np.pi/4, 0, np.pi/2, np.pi/4]
    J = panda_jacobian(q)
    J_pos = J[:3, :]  # position rows only

    x_dot = np.array([0.05, 0.0, 0.0])  # 5cm/s in x
    Vmin = -0.5 * np.ones(7)
    Vmax =  0.5 * np.ones(7)

    q_dot, s = sns_algorithm(J_pos, x_dot, Vmin, Vmax)
    print("SNS test - q_dot:", q_dot)
    print("Scaling factor s:", s)
    print("All within limits:", np.all(q_dot >= Vmin) and np.all(q_dot <= Vmax))