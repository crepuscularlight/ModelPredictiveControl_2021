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

function p = controller_mpc_1(Q, R, T, N, ~)
% controller variables
persistent param yalmip_optimizer

% initialize controller, if not done already
if isempty(param)
    [param, yalmip_optimizer] = init(Q, R, N);
end

% evaluate control action by solving MPC problem
[u_mpc,errorcode] = yalmip_optimizer(T-param.T_sp);
if (errorcode ~= 0)
    warning('MPC1 infeasible');
end
p = u_mpc+param.p_sp;
end

function [param, yalmip_optimizer] = init(Q, R, N)
% get basic controller parameters
% ...
% get terminal cost
param=compute_controller_base_parameters;

A=param.A;
B=param.B;
Ucons=param.Ucons;
Xcons=param.Xcons;

% implement your MPC using Yalmip here
nx = size(param.A,1);
nu = size(param.B,2);
%define symbolic decision variable
U = sdpvar(repmat(nu,1,N-1),ones(1,N-1),'full');%3xN-1
X = sdpvar(repmat(nx,1,N),ones(1,N),'full');%3xN
T0 = sdpvar(nx,1,'full');
%inialize objective and constraints
objective = 0;
constraints = [];

%construct constraints and cost function
for k = 1:N-1
    if k==1
        constraints=[constraints,X{k+1}==A*X{k}+B*U{k},...
            Ucons(:,1)<=U{k}<=Ucons(:,2)];
    
    else
        constraints = [constraints,...
            X{k+1}==A*X{k}+B*U{k},...
            Xcons(:,1)<=X{k}<=Xcons(:,2),...
            Ucons(:,1)<=U{k}<=Ucons(:,2)];
        
    end
    objective = objective + X{k}'*Q*X{k}+U{k}'*R*U{k};
end

%obtain the P-inf to get the terminal cost
[~,P,~]=dlqr(A,B,Q,R,0);
objective = objective + X{N}'*P*X{N};

%set the solver
ops = sdpsettings('verbose',0,'solver','quadprog');

%ensure starting fron the initial value.
constraints=[constraints,X{1}==T0];
yalmip_optimizer = optimizer(constraints,objective,ops,T0,U{1});
end
