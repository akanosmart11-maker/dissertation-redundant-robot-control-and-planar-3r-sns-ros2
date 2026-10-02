# dissertation-redundant-robot-control
# MSc Project: Null-Space Saturation with Dynamic Task Priorities
## Complete Code Package — Milestones 1–4
<img width="4284" height="5712" alt="IMG_7330" src="https://github.com/user-attachments/assets/c5b26b79-c4b3-41eb-94c2-975aa33081af" />

Student: Akano Smart
Supervisor: Dr. Claudio Gaz
Date: september 2026

---

## Overview

MATLAB implementations for kinematic control of a planar 3R redundant robot
under hard joint and Cartesian constraints, with a novel dynamic task priority scheme.

## Robot Model

- Type: Planar 3R (3 revolute joints in the XY plane)
- Link lengths: L1 = 1.0 m, L2 = 0.8 m, L3 = 0.6 m
- Task dimension: m = 2 (end-effector x,y position)
- Redundancy: n - m = 1

---

## Milestone 1: Redundancy Resolution

**Files:** planar_3R_forward_kinematics.m, planar_3R_line_tracking.m, planar_3R_animation.m

- Forward kinematics, Jacobian (verified numerically, error ~1e-10)
- Pseudoinverse properties: J*J#*J = J, P*P = P, J*P = 0
- Three null-space strategies produce identical EE paths, different joint trajectories

## Milestone 2: SNS Algorithm

**Files:** planar_3R_saturation_problem_v2.m, planar_3R_SNS.m

Implements Algorithm 1 & 2 from Flacco, De Luca, Khatib (IEEE TRO, 2015).

| Method         | Max error (mm) | Path deviation (mm) | Final error (mm) |
|----------------|---------------|--------------------|--------------------|
| No limits      | 0.01          | 0.002              | 0.000              |
| Naive clamping | 344.3         | 121.3              | 0.012              |
| Task scaling   | 336.8         | 0.033              | 203.8              |
| **SNS**        | **212.2**     | **0.034**          | **0.000**          |

## Milestone 3: GSNS Algorithm

**File:** planar_3R_GSNS.m

Implements Algorithm 1 from Kazemipour, Khatib, Al Khudir, Gaz, De Luca (IEEE RA-L, 2022).

| Method              | Max elbow y (m) | Ceiling (1.55 m) | EE task |
|---------------------|----------------|-------------------|---------|
| No limits           | 1.7935         | VIOLATED          | Perfect |
| SNS (joint only)    | 1.7935         | VIOLATED          | Perfect |
| **GSNS (joint+Cart)**| **1.5500**    | **Satisfied**     | Scaled  |

## Milestone 4: Dynamic Task Priorities (Novel Contribution)

**File:** planar_3R_dynamic_priorities.m

Novel smooth priority transition scheme within the SNS/GSNS framework.

Scenario: EE tracks a circle (Task 1) while elbow y-coordinate is controlled (Task 2).
Priority shifts smoothly when elbow approaches a danger zone using a cosine activation function.

| Method              | Max EE error (mm) | Max elbow y (m) | Ceiling (1.60 m) |
|---------------------|------------------|-----------------|-------------------|
| Task 1 only         | 158.1            | 1.6438          | VIOLATED          |
| Fixed priority      | 256.3            | 1.4221          | Respected         |
| **Dynamic priority**| **233.2**        | **1.4212**      | **Respected**     |

Key finding: Dynamic priority gives better EE tracking than fixed priority (233 vs 256 mm)
while equally respecting the ceiling constraint. Priority shifts only when needed.

---

## File List

| File | Milestone | Description |
|------|-----------|-------------|
| planar_3R_forward_kinematics.m | 1 | FK, Jacobian, pseudoinverse verification |
| planar_3R_line_tracking.m | 1 | Straight-line tracking, null-space strategies |
| planar_3R_animation.m | 1 | Animated visualisation |
| planar_3R_saturation_problem_v2.m | 2 | Saturation problem demonstration |
| planar_3R_SNS.m | 2 | Basic SNS algorithm |
| planar_3R_GSNS.m | 3 | GSNS with Cartesian constraints |
| planar_3R_dynamic_priorities.m | 4 | Dynamic task priorities (novel) |
| dissertation_outline.md | — | Full dissertation structure |
| README.md | — | This file |

---

## Key References

1. Flacco, De Luca, Khatib, "Control of Redundant Robots Under Hard Joint
   Constraints: Saturation in the Null Space", IEEE TRO, vol. 31, no. 3, 2015.
2. Kazemipour, Khatib, Al Khudir, Gaz, De Luca, "Kinematic Control of Redundant
   Robots With Online Handling of Variable Generalized Hard Constraints",
   IEEE RA-L, vol. 7, no. 4, 2022.
