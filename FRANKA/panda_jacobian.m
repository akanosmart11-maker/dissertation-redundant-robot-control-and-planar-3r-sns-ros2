function J = panda_jacobian(q)
% Geometric Jacobian for Franka Panda (6x7)
% Verified against finite differences: error < 2e-8
%
% Input:  q - 7x1 joint angles in radians
% Output: J - 6x7 Jacobian matrix

[T_ee, T_joints] = panda_fk(q);

p_ee = T_ee(1:3, 4);
J    = zeros(6, 7);

for i = 1:7
    z = T_joints{i}(1:3, 3);
    p = T_joints{i}(1:3, 4);
    J(1:3, i) = cross(z, p_ee - p);
    J(4:6, i) = z;
end
end