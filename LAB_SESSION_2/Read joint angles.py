import frankx  # or whatever library they use
robot = frankx.Robot("192.168.1.11")
state = robot.get_state()
q = state.q  # 7 joint angles in radians
print(q)

# Test script - read joint angles from Panda
# Run this first to verify communication works

print("Connecting to robot...")
# Mahmoud will fill in the correct library and IP
# robot = frankx.Robot("192.168.1.11")
# state = robot.get_state()
# print("Joint angles:", state.q)
print("Replace above with Mahmoud's library syntax")