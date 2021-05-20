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
    % implement the X_LQR computation and assign the result
    A = param.A;
    B = param.B;
    [P,~,~] = dare(A,B,Q,R);
    % Discrete-time LQR controller
    K_E = -1*inv(R + B.'*P*B)*B.'*P*A;
    A_c = A+B*K_E;
    A_x = [eye(3);-eye(3);K_E;-K_E];
    b_x = [param.Xcons(:,2);-param.Xcons(:,1);param.Ucons(:,2);-param.Ucons(:,1)];
    % State constraint set
    X = Polyhedron(A_x,b_x);
    % Discrete-time Algebraic Riccati Equation
    N=30;
    E = X;
    E_cell = cell(N,1);
    E_cell{1} = E;
    % approximate mRPI set
    index=0;
    figure(1)
    for i=1:1:N
        E_last = E_cell{i};
        E_cell{i+1} = intersect(E_last,Polyhedron(E_last.A*A_c,E_last.b));
        if (mldivide(E_cell{i},E_cell{i+1}).volume<0.00001)
            index = i;
            break;
        end
    end
    A_x = E_cell{index}.A; % Important
    b_x = E_cell{index}.b;
    for i=1:1:index
       E_cell{i}.plot('Color',[(1-i/index) (1-i/index) 1])
       alpha 0.5
       hold on;
    end
    T1 = [-2.25;1.75;0.75]
    T2 = [1.5;2.75;-0.25]
    scatter3(T1(1),T1(2),T1(3),'red','filled');
    scatter3(T2(1),T2(2),T2(3),'green','filled');

end