dt = 0.001;

if size(q,1) > size(q,2)
	q = q';
	q_dot = q_dot';
	q_2dot = q_2dot';
	tau_sim = tau_sim';
end

sim_steps = length(q(1,:));
n_joints = length(q(:,1));

omega = zeros(3, n_joints, sim_steps);
omega_dot = zeros(3, n_joints, sim_steps);
v_dot = zeros(3, n_joints, sim_steps);
vs_dot = zeros(3, n_joints, sim_steps);
tau = zeros(1,n_joints, sim_steps);

R = cat(3, R_01, R_12, R_23, R_34, R_45, R_56);
p = cat(2, p_01', p_12', p_23', p_34', p_45', p_56');
s = cat(2, s1', s2', s3', s4', s5', s6');
Ic = cat(3, Ic_1, Ic_2, Ic_3, Ic_4, Ic_5, Ic_6);
m = cat(1, m1, m2, m3, m4, m5, m6);
m = m';



for i = 1:sim_steps
	[omega(:,:,i), omega_dot(:,:,i), v_dot(:,:,i), vs_dot(:,:,i)] = ...
		compute_kinematics(q(:,i), q_dot(:,i), q_2dot(:,i), R, R_W0(:,:,i), p, s);
	tau(:,:,i) = compute_forces_and_torques(omega(:,:,i), omega_dot(:,:,i), vs_dot(:,:,i), m, Ic, p, s, R);
end


% plot(tout(:,1), squeeze(omega(3,1,:)));
% plot(tout, squeeze(omega_dot(3,1,:)));
plot(tout, tau_sim);


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
	v_dot(:,1) = inv(R_W0) * [0 -9.81 0]'; % acceleration of Coordinate System
	vs_dot(:,1) = v_dot(:,1); % acceleration of center of mass (COM)	
	
	for i = 1:n_joints-1
    	invR = inv(R(:,:,i));

    	omega(:,i+1) = invR * (omega(:,i) + q_dot(i+1)*z);
    	omega_dot(:,i+1) = invR * (omega_dot(:,i) + (z * q_2dot(i+1) + cross(omega(:,i), q_dot(i+1) * z)));

		v_dot(:, i+1) = invR * v_dot(:,i) + cross(omega_dot(:, i+1), p(:,i+1)) + cross(omega(:,i+1), cross(omega(:,i+1), p(:,i+1)));
		vs_dot(:,i+1) = v_dot(:,i+1) + cross(omega_dot(:, i+1), s(:, i+1)) + cross(omega(:, i+1), cross(omega(:,i+1), s(:,i+1)));
	end
end

function [tau] = compute_forces_and_torques(omega, omega_dot, vs_dot, m, Ic, p, s, R)
	z = [0 0 1]';
	n_joints = length(omega(1,:));
	f = zeros(3, n_joints+1);
	n = zeros(3, n_joints+1);
	tau = zeros(1, n_joints);

	for i = n_joints:-1:1
		F(:,i) = m(i) * vs_dot(:, i);
		f(:,i) = f(:,i+1) + F(:,i);
		N(:,i) = Ic(:,:,i) * omega_dot(:,i) + cross(omega(:,i), Ic(:,:,i) * omega(:,i));
		n(:,i) = n(:,i+1) + cross(p(:,i)+s(:,i), F(:,i)) + cross(p(:,i), f(:,i+1)) + N(:,i);
		tau(i) = n(:,i)' * inv(R(:,:,i)) * z;
	end

end