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
persistent param yalmip_optimizer Aaug Baug Caug L M P
% initialize controller, if not done already
if isempty(param)
    [param, yalmip_optimizer] = init(Q,R,T,N);
    Aaug = [param.A,param.Bd; zeros(3), eye(3)]; 
    Caug = [param.C_ref, zeros(3)];
    Baug = [param.B;zeros(3)];
    P = [0.1,0.05,0.1,0.3,0.2,0.1];
    L = (place(Aaug',-Caug',P))';
    M = [param.A-eye(3),param.B;eye(3),zeros(3)];
end

% evaluate control action by solving MPC problem
[u_mpc,errorcode] = yalmip_optimizer(T,param.us,param.dh);
if (errorcode ~= 0)
    warning('MPC6 infeasible');
end
p = u_mpc{1}+param.p_sp;
%param.p_sp = param.B\(param.T_sp-param.A*param.T_sp-param.Bd*(param.d+param.dh));
% observer update
y = T-param.T_sp;
% set point update
aux = Aaug*[param.xh;param.dh] + L*(-y+Caug*[param.xh;param.dh])+Baug*u_mpc{1};
param.xh = aux(1:3);
param.dh =  aux(4:6);
param.us = -param.B\(param.Bd*param.dh);

end

function [param, yalmip_optimizer] = init(Q,R,T,N)
% get basic controller parameters
param = compute_controller_base_parameters;
% get terminal cost
[A_x,b_x]=compute_X_LQR(Q,R);
% get terminal set
[P_f,K,G]=idare(param.A,param.B,Q,R,zeros(3),eye(3));
% design disturbance observer

%param.L = -[eye(3);eye(3)];
% init state and disturbance estimate variables
param.xh = T-param.T_sp;
param.dh = [0;0;0];
param.xs=[0;0;0];
param.us=[0;0;0];
% implement your MPC using Yalmip here
nx = size(param.A,1);
nu = size(param.B,2);
U = sdpvar(repmat(nu,1,N-1),ones(1,N-1),'full');
X = sdpvar(repmat(nx,1,N),ones(1,N),'full');
T0 = sdpvar(nx,1,'full');
us =  sdpvar(nx,1,'full') ;
d =  sdpvar(nx,1,'full') ;
constraints = [X{1}==T0-param.T_sp];
objective = 0;

for k = 1:N-1
    constraints = [constraints,(param.Pcons(:,1)-param.p_sp<=U{k})&(U{k}<=param.Pcons(:,2)-param.p_sp),(param.Tcons(:,1)-param.T_sp<=X{k+1})&(X{k+1}<=param.Tcons(:,2)-param.T_sp),X{k+1}==param.A*X{k}+param.B*U{k}+param.Bd*d];
    objective = objective +X{k}'*Q*X{k}+(U{k}-us)'*R*(U{k}-us);
end
constraints = [constraints, A_x*X{N}<=b_x];
objective = objective +  X{N}'*P_f*X{N};
ops = sdpsettings('verbose',0,'solver','quadprog');
yalmip_optimizer = optimizer(constraints,objective,ops,{T0,us,d},U);
end