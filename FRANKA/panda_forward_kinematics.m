function [T_ee, T_joints] = panda_forward_kinematics(q)
% Franka Panda forward kinematics
% Modified DH parameters from Gaz et al. (2019) IEEE RA-L Figure 1
% Convention: ModDH(a, alpha, d, theta) per Dr Gaz's implementation
%
% Verified: ready pose q=[0,-pi/4,0,-3pi/4,0,pi/2,pi/4]
%           gives EE = [0.307, 0, 0.590] matching Franka documentation

a4=0.0825; a5=-0.0825; a7=0.088;
d1=0.333;  d3=0.316;   d5=0.384; df=0.107;

% [a, alpha, d] per joint
DH = [0,   0,      d1;
      0,  -pi/2,   0;
      0,   pi/2,   d3;
      a4,  pi/2,   0;
      a5, -pi/2,   d5;
      0,   pi/2,   0;
      a7,  pi/2,   0];

T_joints    = cell(1,9);
T_joints{1} = eye(4);
T_current   = eye(4);

for i = 1:7
    a     = DH(i,1);
    alpha = DH(i,2);
    d     = DH(i,3);
    theta = q(i);

    Ti = [cos(theta),             -sin(theta),            0,           a;
          cos(alpha)*sin(theta),   cos(alpha)*cos(theta), -sin(alpha), -d*sin(alpha);
          sin(alpha)*sin(theta),   sin(alpha)*cos(theta),  cos(alpha),  d*cos(alpha);
          0,                       0,                      0,           1];

    T_current     = T_current * Ti;
    T_joints{i+1} = T_current;
end

% Flange fixed offset along local z-axis
T_flange      = eye(4);
T_flange(3,4) = df;
T_current     = T_current * T_flange;
T_joints{9}   = T_current;

T_ee = T_current;
end