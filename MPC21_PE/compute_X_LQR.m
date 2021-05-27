% BRRIEF:
%   Template for explicit invariant set computation. You MUST NOT change
%   the output.
% INPUT:
%   Q, R: State and input weighting matrix, dimension (3,3)
% OUTPUT:
%   A_x, b_x: Describes polytopic X_LQR = {x| A_x * x <= b_x}

function [A_x, b_x] = compute_X_LQR(Q, R)
    % get basic controller parameters
    param = compute_controller_base_parameters;
    
    %get the input constraints and the state constraints
    Ucons=param.Ucons;
    Xcons=param.Xcons;
    
    % implement the X_LQR computation and assign the result
    [K,~,~]=dlqr(param.A,param.B,Q,R);
    F=-K;
    
    sys_LQR=LTISystem('A',param.A+param.B*F);
    sys_LQR.x.min=Xcons(:,1);
    sys_LQR.x.max=Xcons(:,2);
    
    Xp=Polyhedron('A',[eye(3);-eye(3);F;-F],...
        'b',[Xcons(:,2);...
           -Xcons(:,1);...
           Ucons(:,2);
           -Ucons(:,1)]);
    %linear constrain matrix 
%     A_x =[eye(3);-eye(3);F;-F]; 
%     b_x = [Xcons(:,2);...
%            -Xcons(:,1);...
%            Ucons(:,2);
%            -Ucons(:,1)];
%   
%      
%     Xp=Polyhedron('A',A_x,'b',b_x);
%     plot(Xp)
    sys_LQR.x.with('setConstraint');
    sys_LQR.x.setConstraint=Xp;
    InvSetLQR=sys_LQR.invariantSet();
    A_x=InvSetLQR.A;
    b_x=InvSetLQR.b;
    InvSetLQR.plot(),alpha(0.5)
end