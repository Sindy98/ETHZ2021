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
    [param, yalmip_optimizer] = init(Q,R,N);
end

% evaluate control action by solving MPC problem
[u_mpc,errorcode] = yalmip_optimizer(T-param.T_sp,num2cell(d,1));
if (errorcode ~= 0)
     warning('MPC5 infeasible');
end
p = u_mpc{1}+param.p_sp;
end

function [param, yalmip_optimizer] = init(Q,R,N)
% get basic controller parameters
% ...
% get terminal cost
param = compute_controller_base_parameters;
[A_x,b_x]=compute_X_LQR(Q,R);
[P,K,G]=idare(param.A,param.B,Q,R,zeros(3),eye(3));
% get terminal set
% ...
% implement your MPC using Yalmip here
nx = size(param.A,1);
nu = size(param.B,2);
U = sdpvar(repmat(nu,1,N-1),ones(1,N-1),'full');
X = sdpvar(repmat(nx,1,N),ones(1,N),'full');
v = sdpvar(1,1,'full');
T0 = sdpvar(nx,1,'full');
d = sdpvar(repmat(nx,1,N),ones(1,N),'full');
objective = 0;
constraints = [X{1}==T0,v>=0];
for k = 1:N-1
    constraints = [constraints,(param.Ucons(:,1)<=U{k})&(U{k}<=param.Ucons(:,2)),(param.Xcons(:,1)-v<=X{k+1})&(X{k+1}<=param.Xcons(:,2)+v),X{k+1}==param.A*X{k}+param.B*U{k}+d{k}];
    objective = objective +X{k}'*Q*X{k}+U{k}'*R*U{k};
end
constraints = [constraints, A_x*X{N}<=b_x];
objective = objective +  X{N}'*P*X{N};
ops = sdpsettings('verbose',0,'solver','quadprog');
yalmip_optimizer = optimizer(constraints,objective,ops,{T0,d{1:N}},{U{1:N-1},v});
end