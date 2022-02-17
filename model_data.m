clear;

% 1 row - 1 joint
waypoints =[0 35 35 30 30 33 33 30 30 0;
			0 -25 -25 8 8 10 10 8 8 0;
			0 125 125 44 44 42 42 44 44 0;
			0 -65 -65 -58 -58 -58 -58 -58 -58 0;
			0 -85 -85 -116 -116 -120 -120 -116 -116 0;
			0 -70 -70 -185 -185 -185 -185 -185 -185 0];

waypoints = waypoints ./ 180 .* pi; % deg to rad

timepoints = [0 1 2 5 5.5 6 6.5 7 7.5 10];

%calculating q points for PTP
vm_ref = [1, 1, 1, 1, 1, 1]';
bm_ref = [1, 1, 1, 1, 1, 1]';

% se = abs(waypoints(2) - waypoints(1));

for i = 1:length(waypoints(1, :))-1
	[vm(:,i), bm(:,i), tb(:,i), tv(:,i), te(:,i)] = calc_t_ramp_te_ta(waypoints(:,i), waypoints(:,i+1), vm_ref, bm_ref);
end


function [s] = compute_points(tb, tv, te, vm, am, dt)
    s = [];
    s_tmp = [];
    for i = 1:length(vm)
        for t = 0:dt:te
            if t < tb(i)
                point = 0.5 .* am(i) .* t.^2;
            elseif t <= tv(i)
                point = vm(i) .* t - 0.5 .* vm(i).^2 ./ am(i);
            else
                point = vm(i) .* tv(i) - 0.5 .* am(i) .* (te(i) - t).^2;
            end
            s_tmp = [s_tmp point];
        end
        s = [s; s_tmp];
        s_tmp = [];
    end
end