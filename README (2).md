# Franka Panda Simulation — MSc Project
## Null-Space Saturation with Dynamic Task Priorities
**Student:** Akano Smart Oluwatobi (K2457624)  
**Supervisor:** Dr Claudio Gaz

## Overview
This folder contains MATLAB simulation code for the Franka Emika Panda 7-DOF robot, implementing and validating the SNS, GSNS, and dynamic task priority algorithms from the MSc project.

DH parameters follow the standard DH convention, verified against Gaz et al. (2019) IEEE RA-L kinematics.

## File descriptions

### Kinematics
- **panda_fk.m** Forward kinematics. Input: `q` (7x1 joint angles in radians).  
  Output: `T_ee` (4x4 end-effector transform), `T_joints` (cell array of 8 cumulative transforms, `T_joints{1}`=base, `T_joints{8}`=EE).

- **panda_jacobian.m** Geometric Jacobian (6x7). Top 3 rows: linear velocity.  
  Bottom 3 rows: angular velocity. Verified against finite differences: max error 1.89e-8.

- **panda_elbow_jacobian_y.m** 1x7 Jacobian of elbow (joint 4 origin) y-coordinate.  
  Used as the Cartesian control point for GSNS and dynamic priorities.

### Algorithms (dimension-agnostic, ported from planar_3R)
- **SNS_algorithm.m** Basic SNS (Flacco et al. 2015, Algorithm 1+2).  
  Inputs: `J`, `x_dot`, `Vmin`, `Vmax`, `n`, `m`.

- **SNS_velocity_p.m** Compact SNS implementation for use inside priority stack.

- **GSNS_algorithm.m** Generalised SNS (Kazemipour et al. 2022, Algorithm 1).  
  Augmented constraint system handles joint and Cartesian bounds.

- **dyn_priority_stack.m** Dynamic task priority blending. `alpha=0`: Task 1 only.  
  `alpha=1`: Task 2 only. Intermediate: cosine blend.

### Verification and tests (run in this order)
1. **verify_panda_kinematics.m** — Milestone 1: kinematics verification
2. **panda_SNS_test.m** — Milestone 2: SNS straight-line tracking
3. **panda_GSNS_test.m** — Milestone 3: GSNS elbow ceiling constraint
4. **panda_dynamic_priorities_test.m** — Milestone 4: dynamic priority scheme

## Key results
- **Milestone 1:** Jacobian error 1.89e-8, `||JP|| = 1.25e-14` (machine precision)
- **Milestone 2:** Max EE error 0.001 mm, final error 0.0001 mm
- **Milestone 3:** Elbow held at `y=0.350m`; unconstrained reaches `y=0.384m`
- **Milestone 4:** Ceiling protected (max `y=0.265m`); both other methods violate it (`y=0.333m`); max `alpha=0.315`

## Notes
- Run files from the `franka_panda/` directory
- All algorithms are general and require no changes for different robot geometry
- The SNS and GSNS functions are identical to the planar_3R versions
- Standard DH convention used (not modified DH)
- `T_joints` indexing: `T_joints{i+1}` = frame after joint i
