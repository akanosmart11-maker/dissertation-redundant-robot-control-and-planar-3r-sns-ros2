#!/usr/bin/env python3
"""
Panda SNS Control Script
Smart Oluwatobi - K2457624
MSc Mechatronics - Kingston University

Connects to Franka Panda and runs SNS velocity control
"""

import pyfranka_interface as franka
import numpy as np
import time

# Import your algorithms
from panda_fk_python import panda_fk, panda_jacobian, get_elbow_y
from panda_sns_python import sns_algorithm

# ===== PARAMETERS =====
ROBOT_IP  = '192.168.1.13'
Vmax_val  = 0.1   # conservative for first test - 10% of max
dt        = 0.01  # 10ms steps (100Hz - safer than 1kHz for first test)
Kp        = 20    # proportional gain

def test_read_joints():
    """Step 1: Just read and print joint angles"""
    print("Connecting to robot...")
    r = franka.Robot_(ROBOT_IP, False, False)
    state = r.getState()
    q = np.array(state.q)
    print("Current joint angles:", q)
    print("EE transform:")
    print(state.T)

    # Verify with your kinematics
    T_ee, _ = panda_fk(q)
    print("Your FK result:", T_ee[:3,3])
    print("Robot FK result:", np.array(state.T)[:3,3])
    print("Difference:", np.linalg.norm(T_ee[:3,3] - np.array(state.T)[:3,3]))
    return r, q

def move_to_joint_position(r, q_target, speed=0.15):
    """Step 2: Move to a joint configuration"""
    print(f"Moving to target joints at speed {speed}...")
    r.move_joints(np.array(q_target), speed)
    state = r.getState()
    print("Reached:", np.array(state.q))

def sns_straight_line(r, delta_position, duration=5.0):
    """
    Step 3: SNS velocity control for straight-line EE motion
    delta_position: [dx, dy, dz] in metres - how far to move
    """
    Vmin = -Vmax_val * np.ones(7)
    Vmax = Vmax_val * np.ones(7)

    # Get start state
    state = r.getState()
    q = np.array(state.q)
    T_ee, _ = panda_fk(q)
    p_start = T_ee[:3, 3]
    p_target = p_start + np.array(delta_position)

    print(f"SNS straight line: {p_start} -> {p_target}")
    print(f"Duration: {duration}s")

    start_time = time.time()
    step = 0

    while time.time() - start_time < duration:
        # Read current state
        state = r.getState()
        q = np.array(state.q)

        # Current EE position from robot FK
        T_robot = np.array(state.T).reshape(4,4)
        p_curr = T_robot[:3, 3]

        # Desired position - quintic timing
        t = time.time() - start_time
        tau = min(t / duration, 1.0)
        gam = 6*tau**5 - 15*tau**4 + 10*tau**3
        gam_dot = (30*tau**4 - 60*tau**3 + 30*tau**2) / duration
        p_des = p_start + gam * np.array(delta_position)
        v_des = gam_dot * np.array(delta_position)

        # Task velocity with feedback
        x_dot = v_des + Kp * (p_des - p_curr)

        # Limit speed
        speed = np.linalg.norm(x_dot)
        if speed > 0.03:  # max 3cm/s
            x_dot = x_dot * 0.03 / speed

        # Jacobian - position rows only
        J = panda_jacobian(q)
        J_pos = J[:3, :]

        # SNS algorithm
        q_dot, s = sns_algorithm(J_pos, x_dot, Vmin, Vmax)

        # Send velocity command
        # NOTE: move_joints takes position not velocity
        # We integrate to get next joint position
        q_next = q + q_dot * dt
        r.move_joints(q_next, 0.1)

        # Print progress every 50 steps
        if step % 50 == 0:
            err = np.linalg.norm(p_des - p_curr) * 1000
            elbow_y = get_elbow_y(q)
            print(f"t={t:.1f}s err={err:.1f}mm s={s:.3f} elbow_y={elbow_y:.4f}m")

        step += 1
        time.sleep(dt)

    # Final state
    state = r.getState()
    q_final = np.array(state.q)
    T_final = np.array(state.T).reshape(4,4)
    print(f"Final EE position: {T_final[:3,3]}")
    print(f"Final error: {np.linalg.norm(T_final[:3,3] - p_target)*1000:.2f} mm")

def gsns_constrained_motion(r, delta_position, y_ceiling=0.30, duration=5.0):
    """
    Step 4: GSNS with elbow y ceiling constraint
    """
    from panda_fk_python import panda_elbow_jacobian_y

    Vmin = -Vmax_val * np.ones(7)
    Vmax_arr = Vmax_val * np.ones(7)

    # GSNS augmented system
    def build_augmented(J_pos, q, y_ceiling, y_floor=-0.5):
        J_elbow = panda_elbow_jacobian_y(q)
        elbow_y = get_elbow_y(q)

        # Shape Cartesian constraint
        cp_Vmax = 0.1
        cp_Amax = 0.5
        v_pos_max = (y_ceiling - elbow_y) / dt
        v_pos_min = (y_floor - elbow_y) / dt
        margin = y_ceiling - elbow_y
        v_stop = np.sqrt(2 * cp_Amax * margin) if margin > 0 else 0
        bmax = min(v_pos_max, cp_Vmax, v_stop)
        bmin = max(v_pos_min, -cp_Vmax)
        if bmax < 0: bmax = 0
        if bmin > 0: bmin = 0

        # Augmented matrix
        A = np.vstack([np.eye(7), J_elbow])
        b_min = np.append(Vmin, bmin)
        b_max = np.append(Vmax_arr, bmax)
        return A, b_min, b_max, elbow_y

    state = r.getState()
    q = np.array(state.q)
    T_start = np.array(state.T).reshape(4,4)
    p_start = T_start[:3,3]
    p_target = p_start + np.array(delta_position)

    print(f"GSNS motion with elbow ceiling at y={y_ceiling}m")
    print(f"Initial elbow y: {get_elbow_y(q):.4f}m")

    start_time = time.time()
    max_elbow_y = get_elbow_y(q)

    while time.time() - start_time < duration:
        state = r.getState()
        q = np.array(state.q)
        T_robot = np.array(state.T).reshape(4,4)
        p_curr = T_robot[:3,3]

        t = time.time() - start_time
        tau = min(t/duration, 1.0)
        gam = 6*tau**5 - 15*tau**4 + 10*tau**3
        gam_dot = (30*tau**4 - 60*tau**3 + 30*tau**2)/duration
        p_des = p_start + gam*np.array(delta_position)
        v_des = gam_dot*np.array(delta_position)
        x_dot = v_des + Kp*(p_des - p_curr)

        speed = np.linalg.norm(x_dot)
        if speed > 0.03:
            x_dot = x_dot * 0.03/speed

        J = panda_jacobian(q)
        J_pos = J[:3,:]

        # GSNS with augmented constraints
        A, b_min, b_max, elbow_y = build_augmented(J_pos, q, y_ceiling)
        max_elbow_y = max(max_elbow_y, elbow_y)

        # Use SNS on augmented system
        q_dot, s = sns_algorithm(A[:3,:], x_dot, b_min[:7], b_max[:7])

        q_next = q + q_dot*dt
        r.move_joints(q_next, 0.1)

        time.sleep(dt)

    print(f"Max elbow y reached: {max_elbow_y:.4f}m")
    print(f"Ceiling was: {y_ceiling}m")
    if max_elbow_y > y_ceiling + 0.005:
        print("CONSTRAINT VIOLATED")
    else:
        print("Constraint satisfied")


# ===== MAIN - RUN IN ORDER =====
if __name__ == "__main__":
    print("="*50)
    print("STEP 1: Read joint angles")
    print("="*50)
    r, q_init = test_read_joints()

    input("\nPress Enter to move to ready pose...")
    ready_pose = [0, -0.785, 0, -2.356, 0, 1.571, 0.785]
    move_to_joint_position(r, ready_pose, speed=0.15)

    input("\nPress Enter to run SNS straight line (2cm in x)...")
    sns_straight_line(r, delta_position=[0.02, 0.0, 0.0], duration=3.0)

    input("\nPress Enter to run GSNS constrained motion...")
    gsns_constrained_motion(r, delta_position=[0.0, 0.05, 0.05],
                           y_ceiling=0.30, duration=5.0)

    print("\nSession complete.")