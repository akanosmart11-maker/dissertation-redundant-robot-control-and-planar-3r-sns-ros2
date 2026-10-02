function J_elbow_y = panda_elbow_jacobian_y(q)
% Returns the 1x7 Jacobian of the elbow y-coordinate
% Elbow = origin of frame after joint 4 = T_joints{5}

[~, T_joints] = panda_fk(q);
p_elbow = T_joints{5}(1:3,4);

J_elbow_y = zeros(1,7);
for i = 1:4
    z = T_joints{i}(1:3,3);
    p = T_joints{i}(1:3,4);
    col = cross(z, p_elbow - p);
    J_elbow_y(i) = col(2);  % y-component only
end
end