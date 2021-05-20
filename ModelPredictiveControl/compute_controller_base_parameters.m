function param = compute_controller_base_parameters
    % load truck parameters
    load('system/parameters_building');
    a_Env_VC = building.a_Env_VC;
    a_F1_VC = building.a_F1_VC;
    a_F2_VC = building.a_F2_VC;
    a_F2_F1 = building.a_F2_F1;
    m_VC = building.m_VC;
    m_F1 = building.m_F1;
    m_F2 = building.m_F2;
    b_11 = building.b_11;
    b_12 = building.b_12;
    b_13 = building.b_13;
    b_21 = building.b_21;
    b_22 = building.b_22;
    b_23 = building.b_23;
    b_31 = building.b_31;
    b_32 = building.b_32;
    b_33 = building.b_33;
    T_Env = building.T_Env;
    d_VC=building.d_VC;
    d_F1 = building.d_F1;
    d_F2 = building.d_F2;
    m = diag([m_VC, m_F1, m_F2]);
    Ac = [-(a_Env_VC+a_F1_VC+a_F2_VC),a_F1_VC,a_F2_VC;
        a_F1_VC,-a_F1_VC-a_F2_F1,a_F2_F1;
        a_F2_VC,a_F2_F1,-(a_F2_VC+a_F2_F1)];
    Ac = m\Ac;
    % Task 1: continuous time dynamics in state space form
%     Ac = [-(a_Env_VC+a_F1_VC+a_F2_VC)/m_VC,a_F1_VC/m_VC,a_F2_VC/m_VC;
%         -a_F1_VC/m_F1,(a_F1_VC-a_F2_F1)/m_F1,a_F2_F1/m_F1;
%         -a_F2_VC/m_F2,-a_F2_F1/m_F2,(a_F2_VC+a_F2_F1)/m_F2];
    Bc = [b_11/m_VC,b_12/m_VC,b_13/m_VC;
        b_21/m_F1,b_22/m_F1,b_23/m_F1;
        b_31/m_F2,b_32/m_F2,b_33/m_F2];
    Bdc = inv(m);
    d = [d_VC+a_Env_VC*T_Env;d_F1;d_F2];
    
    % Task 2: discretization
    Ts = 60;
    A = Ts*Ac+eye(3);
    B = Ts*Bc;
    Bd = Ts*Bdc;
    
    % Task 3: set point computation
    b_ref = [25;-42;-18.5];
    C_ref = eye(3);
    T_sp = b_ref;
    p_sp = B\(T_sp-A*T_sp-Bd*d);
    
    % Task 4: constraints for delta formulation
    Pcons = building.InputConstraints;
    Tcons = building.StateConstraints;
    Ucons = Pcons-p_sp;
    Xcons = Tcons-T_sp;
    
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
