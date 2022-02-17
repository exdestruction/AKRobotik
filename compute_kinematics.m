function [omega, omega_dot, v_dot] = compute_kinematics(q, q_dot, q_2dot, R_W0, R_01, R_12, R_23, R_34, R_45, R_56)
z = [0 0 1]';

R = cat(3, R_01, R_12, R_23, R_34, R_45, R_56);
% initial_v_dot = ;
omega = zeros(3, 6);
omega(:,1) = q_dot(1)*z;

omega_dot = zeros(3, 6);
omega_dot(:,1) = q_2dot(1)*z;

v = zeros(3,6);
v_dot = zeros(3,6);
v_dot(:,1) = inv(R_W0) * [0 -9.81 0]';
% v_dot = 0;

for i = 1:length(q)-1
    invR = inv(R(:,:,i));
    omega(:,i+1) = invR * (omega(:,i) + q_dot(i+1)*z);
    omega_dot(:,i+1) = invR * (omega_dot(:,i) + (z * q_2dot(i+1) + cross(omega(:,i), q_dot(i+1) * z)));
end