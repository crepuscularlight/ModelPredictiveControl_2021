% BRIEF:
%   Controller function template. Input and output dimension MUST NOT be
%   modified.
% INPUT:
%   Q: State weighting matrix, dimension (3,3)
%   R: Input weighting matrix, dimension (3,3)
%   T: Measured system temperatures, dimension (3,1)
%   N: MPC horizon length, dimension (1,1)
% OUTPUT:
%   p: Heating and cooling power, dimension (3,1)

function p = controller_mpc_1_forces(Q, R, T, N, ~)
% controller variables
persistent param forces_optimizer

% initialize controller, if not done already
if isempty(param)
    [param, forces_optimizer] = init(Q, R, N);
end

% evaluate control action by solving MPC problem
[u_mpc,errorcode] = forces_optimizer(T-param.T_sp);
if (errorcode ~= 1)
    warning('MPC1 infeasible');
end
p =u_mpc+param.p_sp;
end

function [param, forces_optimizer] = init(Q, R, N)
% get basic controller parameters
param=compute_controller_base_parameters;

A=param.A;
B=param.B;
Ucons=param.Ucons;
Xcons=param.Xcons;
% get terminal cost
[~,P,~]=dlqr(A,B,Q,R,0);
% implement your MPC using Yalmip2Forces interface here
nx = size(param.A,1);
nu = size(param.B,2);
U = sdpvar(repmat(nu,1,N-1),ones(1,N-1),'full');
X = sdpvar(repmat(nx,1,N),ones(1,N),'full');
T0 = sdpvar(nx,1,'full');
A_x=[eye(3);-eye(3)];
b_x=[Xcons(:,2);-Xcons(:,1)];
A_u=[eye(3);-eye(3)];
b_u=[Ucons(:,2);-Ucons(:,1)];
objective = 0;
constraints = [];
for k = 1:N-1
   constraints=[constraints,X{k+1}==A*X{k}+B*U{k},...
       A_x*X{k+1}<=b_x,...
       A_u*U{k}<=b_u];
   objective = objective + X{k}'*Q*X{k}+U{k}'*R*U{k};
end
objective = objective + X{N}'*P*X{N};
constraints=[constraints,X{1}==T0];
codeoptions = getOptions('solver_name');
forces_optimizer = optimizerFORCES(constraints,objective,codeoptions,T0,U{1});
end
