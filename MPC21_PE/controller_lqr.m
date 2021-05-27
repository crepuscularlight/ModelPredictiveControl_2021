% BRIEF:
%   Controller function template. Input and output dimension MUST NOT be
%   modified.
% INPUT:
%   Q: State weighting matrix, dimension (3,3)
%   R: Input weighting matrix, dimension (3,3)
%   T: Measured system temperatures, dimension (3,1)
% OUTPUT:
%   p: Heating and cooling power, dimension (3,1)

function p = controller_lqr(Q, R, T, ~, ~)
% controller variables
persistent param;

% initialize controller, if not done already
if isempty(param)
    param = init(Q, R);
end

A=param.A;
B=param.B;
T_sp=param.T_sp;
p_sp=param.p_sp;

% compute control action
%discrete lqr solving function to get the feedback matrix F,u=-Fx
[F,P_inf,~]=dlqr(A,B,Q,R);
delta_T=T-T_sp;
delta_u=-F*delta_T;
p=delta_u+p_sp;
end

function param = init(Q, R)
% get basic controller parameters
param = compute_controller_base_parameters;
% add additional parameters if necessary, e.g.


%calculate P_inf and F


end