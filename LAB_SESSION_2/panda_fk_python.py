import numpy as np

def panda_fk(q):
    """
    Franka Panda forward kinematics - standard DH convention
    Input:  q - list or array of 7 joint angles in radians
    Output: T_ee - 4x4 end-effector transform
            T_joints - list of 8 cumulative transforms (T_joints[0]=base)
    """
    # DH parameters [a, d, alpha] for joints 1-7
    DH = [
        [0,       0.333,  0],
        [0,       0,     -np.pi/2],
        [0,       0.316,  np.pi/2],
        [0.0825,  0,      np.pi/2],
        [-0.0825, 0.384, -np.pi/2],
        [0,       0,      np.pi/2],
        [0.088,   0.107,  np.pi/2],
    ]

    T_joints = [np.eye(4)]  # T_joints[0] = base frame
    T_current = np.eye(4)

    for i in range(7):
        a     = DH[i][0]
        d     = DH[i][2]  # wait - standard DH: [a, d, alpha]
        alpha = DH[i][2]
        d     = DH[i][1]
        th    = q[i]

        ct = np.cos(th); st = np.sin(th)
        ca = np.cos(alpha); sa = np.sin(alpha)

        Ti = np.array([
            [ct,        -st*ca,   st*sa,  a*ct],
            [st,         ct*ca,  -ct*sa,  a*st],
            [0,          sa,      ca,     d   ],
            [0,          0,       0,      1   ]
        ])

        T_current = T_current @ Ti
        T_joints.append(T_current.copy())

    return T_current, T_joints


def panda_jacobian(q):
    """
    Geometric Jacobian for Franka Panda (6x7)
    Input:  q - 7 joint angles in radians
    Output: J - 6x7 Jacobian matrix
    """
    T_ee, T_joints = panda_fk(q)
    p_ee = T_ee[:3, 3]
    J = np.zeros((6, 7))

    for i in range(7):
        z = T_joints[i][:3, 2]
        p = T_joints[i][:3, 3]
        J[:3, i] = np.cross(z, p_ee - p)
        J[3:, i] = z

    return J


def panda_elbow_jacobian_y(q):
    """
    1x7 Jacobian of elbow (joint 4 origin) y-coordinate
    """
    T_ee, T_joints = panda_fk(q)
    p_elbow = T_joints[4][:3, 3]  # frame after joint 4

    J_elbow_y = np.zeros(7)
    for i in range(4):  # only joints 1-4 affect elbow
        z = T_joints[i][:3, 2]
        p = T_joints[i][:3, 3]
        col = np.cross(z, p_elbow - p)
        J_elbow_y[i] = col[1]  # y-component

    return J_elbow_y.reshape(1, 7)


def get_elbow_y(q):
    """Returns elbow y-coordinate"""
    _, T_joints = panda_fk(q)
    return T_joints[4][1, 3]


# Quick verification
if __name__ == "__main__":
    q_ready = [0, -np.pi/4, 0, -3*np.pi/4, 0, np.pi/2, np.pi/4]
    T_ee, _ = panda_fk(q_ready)
    print("Ready pose EE position:", T_ee[:3, 3])
    print("Expected approx: [0.307, 0.000, 0.590]")

    J = panda_jacobian(q_ready)
    print("Jacobian shape:", J.shape)
    print("Jacobian rank:", np.linalg.matrix_rank(J))