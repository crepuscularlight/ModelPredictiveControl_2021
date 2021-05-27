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

function p = controller_mpc_6(Q,R,T,N,~)
% controller variables
persistent param yalmip_optimizer
global x_hat

% initialize controller, if not done already
if isempty(param)
    [param, yalmip_optimizer] = init(Q,R,T,N);
    param.x_hat(1:3)=T;
    x_hat=[param.x_hat];
end

%calculate the new reference
M=[eye(3)-param.A,-param.B;
   param.C_ref,zeros(3,3)];
reference=pinv(M)*[param.x_hat(4:6);param.b_ref];

% evaluate control action by solving MPC problem
[u_mpc,errorcode] = yalmip_optimizer(T-reference(1:3),reference(1:3),reference(4:6));
if (errorcode ~= 0)
    warning('MPC6 infeasible');
end
p = u_mpc{1}+reference(4:6);

% observer update 
%lecture 7 p29
param.x_hat=param.A_aug*param.x_hat...
    +param.B_aug*p...
    +param.L*param.C_aug*[(T-param.x_hat(1:3));zeros(3,1)];

x_hat=[x_hat,param.x_hat];

% set point update
% ...
end

function [param, yalmip_optimizer] = init(Q,R,T,N)
% get basic controller parameters
param=compute_controller_base_parameters;

A=param.A;
B=param.B;
Pcons=param.Pcons;
Xcons=param.Xcons;
% get terminal cost
[~,P,~]=dlqr(A,B,Q,R,0);

% get terminal set
[A_x,b_x]=compute_X_LQR(Q,R);

% design disturbance observer
% ...
% init state and disturbance estimate variables
% ...
% implement your MPC using Yalmip here
nx = size(param.A,1);
nu = size(param.B,2);
U = sdpvar(repmat(nu,1,N-1),ones(1,N-1),'full');
X = sdpvar(repmat(nx,1,N),ones(1,N),'full');

X_offset=sdpvar(nx,1);
U_offset=sdpvar(nu,1);

Ucons=Pcons-[U_offset,U_offset];
Xcons=param.Tcons-[X_offset,X_offset];

T0 = sdpvar(nx,1,'full');
T_sp = param.T_sp;
p_sp = param.p_sp;
constraints = [];
objective =0;
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

p_in={T0,X_offset,U_offset};
p_out={U{1},objective};

ops = sdpsettings('verbose',0,'solver','quadprog');
yalmip_optimizer = optimizer(constraints,objective,ops,p_in,p_out);

%lecture 7 p29
A_aug=[param.A, eye(3,3);
    zeros(3,3), eye(3,3)];

B_aug=[param.B ;zeros(3,3)];

C_aug=[eye(3,3),zeros(3,3)];

param.A_aug=A_aug; 
param.B_aug=B_aug;
param.C_aug=C_aug;

param.L=place(A_aug',C_aug',[0 0 0 0.3 0.5 0.5]')';
param.p=[0;0;0];

%construct the augment variable.
param.x_hat=[0;0;0;param.Bd*param.d];
end