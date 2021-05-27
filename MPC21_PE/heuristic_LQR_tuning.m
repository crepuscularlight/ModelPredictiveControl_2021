% BRRIEF:
%   Template for tuning of Q and R matrices for LQR controller as described 
%   in task 6.
% INPUT:
%   n_samples:  Number of samples considered for the tuning.
%   T0:         Initial condition
%   T_sp:       Set point
%   scen:       Disturbance scenario
% OUTPUT:
%   Q, R: Describes stage cost of LQR controller (x^T Q x + u^T R u)

function [Q, R] = heuristic_LQR_tuning(n_samples, T0, T_sp, scen)

%initialize the figure
figure(2); set(gcf, 'WindowStyle' ,'docked'); grid on; hold on
xlabel('Energy consumption [kWh]'); 
ylabel('Relative norm of steady state deviation');

%set the R as an identical matrix
R=eye(3);
T0_1=T0;
%initialize the Q
Q=zeros(3,3);

dt_best=Inf;
for index = 1:n_samples
    %generate Q randomly
    Q_idx=diag([randi([1 10e6]) randi([1 10e6]) randi([1 10e6])]);
    [T, p, ~,~ , T_v, p_v]=simulate_building(T0_1, @controller_lqr, Q_idx, R, scen, 0);
    T_15=T(:,16);
    
    %dT_relative and total power consumed as the standard.
    dT_relative=norm(T_15-T_sp,2)/norm(T0_1-T_sp,2);
    power_sum = sum(abs(p), 'all')/1000/60;
    RGB='g';
    if T_v==1
        RGB='r';
    end
    if p_v==1
        RGB='b';
    end
    if T_v==1 &p_v==1
        RGB='r';
    end
    if T_v==0 & p_v==0
        RGB='g';
        if(power_sum<16 & dT_relative<dt_best)
            Q=Q_idx;
            Q_pos=[power_sum,dT_relative];
            dt_best=dT_relative;
        end
    end
    scatter(power_sum, dT_relative, [], RGB);
    
end
scatter(Q_pos(1),Q_pos(2),200,'r','filled');
% simulate_building(T0_1, @controller_lqr, Q, R, scen, 1)
end
