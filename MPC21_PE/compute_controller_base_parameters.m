function param = compute_controller_base_parameters
    % load truck parameters
    load('system/parameters_building');
    
    m_VC=building.m_VC;
    m_F1=building.m_F1;
    m_F2=building.m_F2;
    
    a_Env_VC=building.a_Env_VC;
    a_F1_VC=building.a_F1_VC;
    a_F2_VC=building.a_F2_VC;
    a_F2_F1=building.a_F2_F1;
    
    b_11=building.b_11;
    b_12=building.b_12;
    b_13=building.b_13;
    b_21=building.b_21;
    b_22=building.b_22;
    b_23=building.b_23;
    b_31=building.b_31;
    b_32=building.b_32;
    b_33=building.b_33;
    
    T_Env=building.T_Env;
    d_VC=building.d_VC;
    d_F1=building.d_F1;
    d_F2=building.d_F2;
    
    
    % Task 1: continuous time dynamics in state space form
    Ac = [-(a_F1_VC+a_F2_VC+a_Env_VC)/m_VC    a_F1_VC/m_VC               a_F2_VC/m_VC; ...
             a_F1_VC/m_F1                    -(a_F1_VC+a_F2_F1)/m_F1     a_F2_F1/m_F1; ...
             a_F2_VC/m_F2                     a_F2_F1/m_F2               -(a_F2_VC+a_F2_F1)/m_F2];
         
    Bc = [b_11/m_VC b_12/m_VC b_13/m_VC;...
          b_21/m_F1 b_22/m_F1 b_23/m_F1;...
          b_31/m_F2 b_32/m_F2 b_33/m_F2];
      
    Bdc = diag([1/m_VC 1/m_F1 1/m_F2]);
    
    %add constant terms to d.
    d = [d_VC d_F1 d_F2]'+[a_Env_VC*T_Env 0 0]';
    
    % Task 2: discretization
    Ts = 60;
%     Euler
%     A = diag([1 1 1])+Ts*Ac
%     B = Bc*Ts
%     Bd = Bdc*Ts

%   Exact
    sysc1=ss(Ac,Bc,eye(3),[]);
    sysd1=c2d(sysc1,Ts);
    A = sysd1.A;
    B = sysd1.B;
    
    sysc2=ss(Ac,Bdc,eye(3),[]);
    sysd2=c2d(sysc2,Ts);
    Bd=sysd2.B;
    
    
    % Task 3: set point computation
    b_ref = [25;-42;-18.5];
    C_ref = eye(3);
    joint_coef=[eye(3)-A -B;...
        C_ref zeros(3,3)];
    
    %pseudo inverse to get the most approximate solution for T and P at
    %stable point.
    result= pinv(joint_coef)*[Bd*d;b_ref];
    T_sp =result(1:3)
    p_sp = result(4:6)
    
    % Task 4: constraints for delta formulation
    Pcons = building.InputConstraints;
    Tcons = building.StateConstraints;
    Ucons = Pcons-[p_sp p_sp];
    Xcons = Tcons-[T_sp T_sp];
    
    % put everything together
    param.A = A;
    param.B = B;
    param.Bd = Bd;
    param.d = d;
    param.b_ref = b_ref;
    param.C_ref = C_ref;
    param.T_sp = T_sp;
    param.p_sp = p_sp;
    param.Ucons = Ucons;
    param.Xcons = Xcons;
    param.Tcons = Tcons;
    param.Pcons = Pcons;
end
