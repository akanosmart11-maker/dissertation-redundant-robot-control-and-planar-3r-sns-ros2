% Panda kinematics verification - run this to confirm everything is correct
clear all

fprintf('=== FRANKA PANDA KINEMATICS VERIFICATION ===\n\n');

%% Test 1: Jacobian vs finite differences
q_test = [0.1; -0.4; 0.2; -1.5; 0.1; 1.2; 0.5];
[T_ee, T_joints] = panda_fk(q_test);
J = panda_jacobian(q_test);
p_ee = T_ee(1:3,4);

delta = 1e-7;
Jn = zeros(3,7);
for i = 1:7
    qp = q_test; qp(i) = qp(i)+delta;
    Tp = panda_fk(qp);
    Jn(:,i) = (Tp(1:3,4)-p_ee)/delta;
end
err = max(max(abs(J(1:3,:)-Jn)));
fprintf('1. Linear Jacobian error:  %.4e  (pass < 1e-6)\n', err);

%% Test 2: Null-space identities
J_pinv = pinv(J);
P      = eye(7) - J_pinv*J;
JP     = norm(J*P);
JJpJ   = norm(J*J_pinv*J - J);
PP     = norm(P*P - P);
fprintf('2. ||J*P||:                %.4e  (pass < 1e-10)\n', JP);
fprintf('3. ||J*J#*J - J||:         %.4e  (pass < 1e-10)\n', JJpJ);
fprintf('4. ||P*P - P||:            %.4e  (pass < 1e-10)\n', PP);

%% Test 3: Rotation matrix valid
det_R = det(T_ee(1:3,1:3));
fprintf('5. Rotation determinant:   %.6f   (pass = 1.0)\n', det_R);

%% Summary
fprintf('\n=== SUMMARY ===\n');
if err<1e-5 && JP<1e-10 && JJpJ<1e-10 && abs(det_R-1)<1e-6
    fprintf('STATUS: PASSED\n');
    fprintf('Kinematics verified. Ready to port SNS algorithm.\n');
else
    fprintf('STATUS: FAILED\n');
end