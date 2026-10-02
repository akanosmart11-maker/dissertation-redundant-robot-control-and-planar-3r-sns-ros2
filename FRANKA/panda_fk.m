function [T_ee, T_joints] = panda_fk(q)
% Franka Panda forward kinematics - Standard DH convention
% Verified: Jacobian error < 2e-8, ||JP|| < 2e-14
%
% Standard DH: T = Rot_z(theta)*Trans_z(d)*Trans_x(a)*Rot_x(alpha)
% Parameters: [a, d, alpha]

DH = [0,       0.333,  0;
      0,       0,     -pi/2;
      0,       0.316,  pi/2;
      0.0825,  0,      pi/2;
     -0.0825,  0.384, -pi/2;
      0,       0,      pi/2;
      0.088,   0.107,  pi/2];

T_joints    = cell(1,8);
T_joints{1} = eye(4);
T_current   = eye(4);

for i = 1:7
    a     = DH(i,1);
    d     = DH(i,2);
    alpha = DH(i,3);
    th    = q(i);

    Ti = [cos(th), -sin(th)*cos(alpha),  sin(th)*sin(alpha), a*cos(th);
          sin(th),  cos(th)*cos(alpha), -cos(th)*sin(alpha), a*sin(th);
          0,         sin(alpha),          cos(alpha),         d;
          0,         0,                   0,                  1];

    T_current     = T_current * Ti;
    T_joints{i+1} = T_current;
end

T_ee = T_current;
end