dt = 0.001;

if size(q,1) > size(q,2)
	q = q';
	q_dot = q_dot';
	q_2dot = q_2dot';
	tau = tau';
end

sim_steps = length(q(1,:));
n_joints = length(q(:,1));

omega = zeros(3, n_joints, sim_steps);
omega_dot = zeros(3, n_joints, sim_steps);
v_dot = zeros(3, n_joints, sim_steps);
vs_dot = zeros(3, n_joints, sim_steps);

R = cat(3, R_01, R_12, R_23, R_34, R_45, R_56);
p = cat(2, p_01', p_12', p_23', p_34', p_45', p_56');
s = cat(2, s1', s2', s3', s4', s5', s6');



for i = 1:sim_steps
	[omega(:,:,i), omega_dot(:,:,i), v_dot(:,:,i), vs_dot(:,:,i)] = ...
		compute_kinematics(q(:,i), q_dot(:,i), q_2dot(:,i), R, R_W0, p, s);
end


% plot(tout(:,1), squeeze(omega(3,1,:)));
% plot(tout, squeeze(omega_dot(3,1,:)));


function [omega, omega_dot, v_dot, vs_dot] = compute_kinematics(q, q_dot, q_2dot, R, R_W0, p, s)
	z = [0 0 1]';
	n_joints = length(q);
	
	% initial_v_dot = ;
	omega = zeros(3, 6);
	omega(:,1) = q_dot(1)*z;
	
	omega_dot = zeros(3, 6);
	omega_dot(:,1) = q_2dot(1)*z;
	
	v = zeros(3,6);
	v_dot = zeros(3,6);
	v_dot(:,1) = inv(R_W0) * [0 -9.81 0]';
	vs_dot(:,1) = v_dot(:,1);
	% v_dot = 0;
	
	
	for i = 1:n_joints-1
    	invR = inv(R(:,:,i));

    	omega(:,i+1) = invR * (omega(:,i) + q_dot(i+1)*z);
    	omega_dot(:,i+1) = invR * (omega_dot(:,i) + (z * q_2dot(i+1) + cross(omega(:,i), q_dot(i+1) * z)));

		v_dot(:, i+1) = invR * v_dot(:,i) + cross(omega_dot(:, i+1), p(:,i+1)) + cross(omega(:,i+1), cross(omega(:,i+1), p(:,i+1)));
		vs_dot(:,i+1) = v_dot(:,i+1) + cross(omega_dot(:, i+1), s(:, i+1)) + cross(omega(:, i+1), cross(omega(:,i+1), s(:,i+1)));
	end
end