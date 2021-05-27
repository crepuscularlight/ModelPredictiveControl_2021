% BRIEF:
%   Controller function template. Input and output dimension MUST NOT be
%   modified.
% INPUT:
%   Q: State weighting matrix, dimension (3,3)
%   R: Input weighting matrix, dimension (3,3)
%   T: Measured system temperatures, dimension (3,1)
%   N: MPC horizon length, dimension (1,1)
%   d: Disturbance matrix, dimension (3,N)
% OUTPUT:
%   p: Heating and cooling power, dimension (3,1)

function p = controller_mpc_5(Q,R,T,N,d)
% controller variables
persistent param yalmip_optimizer

% initialize controller, if not done already
if isempty(param)
    [param, yalmip_optimizer] = init(Q,R,N,d);
end

% evaluate control action by solving MPC problem
[u_mpc,errorcode] = yalmip_optimizer(T-param.T_sp);
if (errorcode ~= 0)
    warning('MPC5 infeasible');
end
p = u_mpc+param.p_sp;
end

function [param, yalmip_optimizer] = init(Q,R,N,d)
load('system/parameters_building');

m_VC=building.m_VC;
m_F1=building.m_F1;
m_F2=building.m_F2;
Bdc = diag([1/m_VC 1/m_F1 1/m_F2]);
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
v = sdpvar(1,1,'full');
T0 = sdpvar(nx,1,'full');
E=sdpvar(repmat(nx,1,N),ones(1,N),'full');

S=eye(nx);
v=100;


objective = 0;
constraints = [];

for k = 1:N-1
    constraints = [constraints,...
        X{k+1}==A*X{k}+B*U{k}+Bdc*d(:,k)...
        Xcons(:,1)-E{k}<=X{k+1}<=Xcons(:,2)+E{k},...
        Ucons(:,1)<=U{k}<=Ucons(:,2),...
        E{k}>=zeros(3,1)];   
    objective = objective + X{k}'*Q*X{k}+U{k}'*R*U{k}+E{k}'*S*E{k}+v*norm(E{k},1);
end
constraints = [constraints, A_x*X{N}<=b_x];
objective = objective + X{N}'*P*X{N}+E{N}'*S*E{N};
constraints=[constraints,X{1}==T0];
ops = sdpsettings('verbose',0,'solver','quadprog');
yalmip_optimizer = optimizer(constraints,objective,ops,T0,U{1});
end