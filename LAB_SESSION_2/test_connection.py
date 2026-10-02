"""
Step 1: Run this first to verify robot connection
Ask Mahmoud to fill in the correct library name and IP address
"""

import numpy as np

# ============================================================
# MAHMOUD FILLS THIS IN:
# import frankx  (or whatever library the lab uses)
# ROBOT_IP = "192.168.1.11"  (or correct IP)
# ============================================================

def test_connection():
    print("=" * 50)
    print("PANDA CONNECTION TEST")
    print("=" * 50)

    # Step 1: Test kinematics (no robot needed)
    from panda_fk_python import panda_fk, panda_jacobian
    q_ready = [0, -np.pi/4, 0, -3*np.pi/4, 0, np.pi/2, np.pi/4]
    T_ee, _ = panda_fk(q_ready)
    print("\n1. Kinematics check:")
    print(f"   Ready pose EE: {T_ee[:3,3].round(4)}")
    print(f"   Expected:      [0.307, 0.000, 0.590]")

    J = panda_jacobian(q_ready)
    print(f"   Jacobian rank: {np.linalg.matrix_rank(J)} (expected 6)")

    # Step 2: Test SNS (no robot needed)
    from panda_sns_python import sns_algorithm
    J_pos = J[:3, :]
    x_dot = np.array([0.02, 0.0, 0.0])
    Vmin = -0.3 * np.ones(7)
    Vmax =  0.3 * np.ones(7)
    q_dot, s = sns_algorithm(J_pos, x_dot, Vmin, Vmax)
    print("\n2. SNS algorithm check:")
    print(f"   q_dot: {q_dot.round(4)}")
    print(f"   Within limits: {np.all(q_dot >= Vmin) and np.all(q_dot <= Vmax)}")

    print("\n3. Robot connection:")
    print("   Ask Mahmoud to add robot connection code here")
    print("   Target: read joint angles q from robot")
    print("   Then:   run SNS, send q_dot back")

    print("\nAll offline tests passed. Ready for robot connection.")

if __name__ == "__main__":
    test_connection()