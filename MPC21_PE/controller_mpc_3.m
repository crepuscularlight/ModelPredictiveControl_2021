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

function p = controller_mpc_3(Q, R, T, N, ~)
% controller variables
persistent param yalmip_optimizer

% initialize controller, if not done already
if isempty(param)
    [param, yalmip_optimizer] = init(Q, R, N);
end

% evaluate control action by solving MPC problem
[u_mpc,errorcode] = yalmip_optimizer(T-param.T_sp);
if (errorcode ~= 0)
    warning('MPC3 infeasible');
end
p = u_mpc+param.p_sp;
end

function [param, yalmip_optimizer] = init(Q, R, N)
% get basic controller parameters
param=compute_controller_base_parameters;

A=param.A;
B=param.B;
Ucons=param.Ucons;
Xcons=param.Xcons;

% get terminal cost
[~,P,~]=dlqr(A,B,Q,R,0);

% get terminal set
[A_x,b_x]=compute_X_LQR(Q,R);

% implement your MPC using Yalmip here
nx = size(param.A,1);
nu = size(param.B,2);
U = sdpvar(repmat(nu,1,N-1),ones(1,N-1),'full');
X = sdpvar(repmat(nx,1,N),ones(1,N),'full');
T0 = sdpvar(nx,1,'full');
objective = 0;
constraints = [];
for k = 1:N-1
    if k==1
        constraints=[constraints,
            X{k+1}==A*X{k}+B*U{k},...
            Ucons(:,1)<=U{k}<=Ucons(:,2)];
    
    else
        constraints = [constraints,...
            X{k+1}==A*X{k}+B*U{k},...
            Xcons(:,1)<=X{k}<=Xcons(:,2),...
            Ucons(:,1)<=U{k}<=Ucons(:,2)];
        
    end
    objective = objective + X{k}'*Q*X{k}+U{k}'*R*U{k};
end
constraints = [constraints, A_x*X{N}<=b_x];
objective = objective + X{N}'*P*X{N};

constraints=[constraints,X{1}==T0]
ops = sdpsettings('verbose',0,'solver','quadprog');
yalmip_optimizer = optimizer(constraints,objective,ops,T0,U{1});
end

