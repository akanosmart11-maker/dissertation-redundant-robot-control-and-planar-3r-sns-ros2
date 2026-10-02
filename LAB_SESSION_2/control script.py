import time
import numpy as np
from panda_fk_python import panda_fk, panda_jacobian
from panda_sns_python import sns_algorithm

# Parameters
Vmin = -0.1 * np.ones(7)  # conservative limits for first test
Vmax =  0.1 * np.ones(7)
Kp = 30
dt = 0.001  # 1kHz

# Target: move 3cm in x from current position
# (Mahmoud will help fill in the robot API calls)

def control_loop(robot, target_position, duration=3.0):
    start_time = time.time()
    
    while time.time() - start_time < duration:
        # Read current state
        state = robot.get_state()  # Mahmoud fills this in
        q = np.array(state.q)
        
        # Compute kinematics
        T_curr, _ = panda_fk(q)
        p_curr = T_curr[:3, 3]
        J = panda_jacobian(q)
        J_pos = J[:3, :]
        
        # Compute task velocity
        error = target_position - p_curr
        x_dot = Kp * error
        
        # Limit task velocity magnitude
        speed = np.linalg.norm(x_dot)
        if speed > 0.05:  # max 5cm/s
            x_dot = x_dot * 0.05 / speed
        
        # SNS algorithm
        q_dot, s = sns_algorithm(J_pos, x_dot, Vmin, Vmax)
        
        # Send to robot (Mahmoud fills this in)
        robot.send_velocity(q_dot)  # placeholder
        
        time.sleep(dt)
    
    # Stop robot
    robot.send_velocity(np.zeros(7))
    print("Motion complete")
    print(f"Final position: {p_curr}")