function [posEst,linVelEst,oriEst,windEst,driftEst,...
          posVar,linVelVar,oriVar,windVar,driftVar,estState] = ...
    Estimator(estState,actuate,sense,tm,estConst)
% [posEst,linVelEst,oriEst,windEst,driftEst,...
%    posVar,linVelVar,oriVar,windVar,driftVar,estState] = 
% Estimator(estState,actuate,sense,tm,estConst)
%
% The estimator.
%
% The function is initialized for tm == 0, otherwise the estimator does an
% iteration step (compute estimates for the time step k).
%
% Inputs:
%   estState        previous estimator state (time step k-1)
%                   May be defined by the user (for example as a struct).
%   actuate         control input u(k-1), [1x2]-vector
%                   actuate(1): u_t, thrust command
%                   actuate(2): u_r, rudder command
%   sense           sensor measurements z(k), [1x5]-vector, INF entry if no
%                   measurement
%                   sense(1): z_a, distance measurement a
%                   sense(2): z_b, distance measurement b
%                   sense(3): z_c, distance measurement c
%                   sense(4): z_g, gyro measurement
%                   sense(5): z_n, compass measurement
%   tm              time t_k, scalar
%                   If tm==0 initialization, otherwise estimator
%                   iteration step.
%   estConst        estimator constants (as in EstimatorConst.m)
%
% Outputs:
%   posEst          position estimate (time step k), [1x2]-vector
%                   posEst(1): p_x position estimate
%                   posEst(2): p_y position estimate
%   linVelEst       velocity estimate (time step k), [1x2]-vector
%                   linVelEst(1): s_x velocity estimate
%                   linVelEst(2): s_y velocity estimate
%   oriEst          orientation estimate (time step k), scalar
%   windEst         wind direction estimate (time step k), scalar
%   driftEst        estimate of the gyro drift b (time step k), scalar
%   posVar          variance of position estimate (time step k), [1x2]-vector
%                   posVar(1): x position variance
%                   posVar(2): y position variance
%   linVelVar       variance of velocity estimate (time step k), [1x2]-vector
%                   linVelVar(1): x velocity variance
%                   linVelVar(2): y velocity variance
%   oriVar          variance of orientation estimate (time step k), scalar
%   windVar         variance of wind direction estimate(time step k), scalar
%   driftVar        variance of gyro drift estimate (time step k), scalar
%   estState        current estimator state (time step k)
%                   Will be input to this function at the next call.
%
%
% Class:
% Recursive Estimation
% Spring 2021
% Programming Exercise 1
%
% --
% ETH Zurich
% Institute for Dynamic Systems and Control
% Raffaello D'Andrea, Matthias Hofer, Carlo Sferrazza
% hofermat@ethz.ch
% csferrazza@ethz.ch


%% Initialization
if (tm == 0)
    % Do the initialization of your estimator here!
    
    % initial state mean
%     r = rand*estConst.StartRadiusBound;
%     rio = rand*pi;
%     posEst = [r*cos(rio),r*sin(rio)]; % 1x2 matrix
%     linVelEst = [0,0] % 1x2 matrix
%     oriEst = 2*(rand-0.5)*estConst.RotationStartBound; % 1x1 matrix
%     windEst =2*(rand-0.5)*estConst.WindAngleStartBound;% 1x1 matrix
%     driftEst = estConst.GyroDriftStartBound; % 1x1 matrix
    posEst = [0,0]; % 1x2 matrix
    linVelEst = [0,0]; % 1x2 matrix
    oriEst = 0; % 1x1 matrix
    windEst =0.01;% 1x1 matrix
    driftEst =0; % 1x1 matrix
    
    % initial state variance
    posVar = [estConst.StartRadiusBound^2/4,estConst.StartRadiusBound^2/4]; % 1x2 matrix
    linVelVar = [0,0]; % 1x2 matrix
    oriVar = estConst.RotationStartBound^2/3; % 1x1 matrix
    windVar = estConst.WindAngleNoise; % 1x1 matrix
    driftVar = estConst.GyroDriftNoise; % 1x1 matrix
    
    % estimator variance init (initial posterior variance)
    estState.Pm = diag([posVar,linVelVar,oriVar,windVar,driftVar]);
    % estimator state
    estState.xm = [posEst,linVelEst,oriEst,windEst, driftEst];
    % time of last update
    estState.tm = tm;
    estState.sense = zeros(1,5);
else

%% Estimator iteration.
% get time since last estimator update
    dt = tm-estState.tm;
    estState.tm = tm; % update measurement update time
    
    % prior update
    Cdh = estConst.dragCoefficientHydr;
    Cda = estConst.dragCoefficientAir;
    Cw = estConst.windVel;
    Cr = estConst.rudderCoefficient;
    tmp = num2cell(estState.xm);
    [px,py,sx,sy,ori,wind,drif] = deal(tmp{:});
    
%     fun1 = @(p1,p2,sx,sy,ori,wind,drif)[sx,sy,...
%         cos(ori)*(tanh(actuate(1))-Cdh*(sx^2+sy^2))-Cda*(sx-Cw*cos(wind)*sqrt((sx-Cw*cos(wind))^2+(sy-Cw*sin(wind))^2)),...
%         sin(ori)*(tanh(actuate(1))-Cdh*(sx^2+sy^2))-Cda*(sy-Cw*sin(wind)*sqrt((sx-Cw*cos(wind))^2+(sy-Cw*sin(wind))^2)),...
%         Cr*actuate(2),...
%         0,...
%         0];
    xm = ode15s(@odefun,[0,dt],estState.xm);
%      xm = [sx,sy,cos(ori)*(tanh(actuate(1))-Cdh*(sx^2+sy^2))-Cda*(sx-Cw*cos(wind)*sqrt((sx-Cw*cos(wind))^2+(sy-Cw*sin(wind))^2)),...
%         sin(ori)*(tanh(actuate(1))-Cdh*(sx^2+sy^2))-Cda*(sy-Cw*sin(wind)*sqrt((sx-Cw*cos(wind))^2+(sy-Cw*sin(wind))^2)),...
%         Cr*actuate(2),...
%         0,...
%         0]*dt+estState.xm;
    func2 = @(p1,p2,sx,sy,ori,wind,drif)[0,0,1,0,0,0,0;...
        0,0,0,1,0,0,0;...
        0,0,...
        -cos(ori)*Cdh*2*sx-Cda*sqrt((sx-Cw*cos(wind))^2+(sy-Cw*sin(wind))^2)-Cda*(sx-Cw*cos(wind))^2/sqrt((sx-Cw*cos(wind))^2+(sy-Cw*sin(wind))^2),...
        -cos(ori)*Cdh*2*sy-Cda*(sx-Cw*cos(wind))*(sy-Cw*sin(wind))/sqrt((sx-Cw*cos(wind))^2+(sy-Cw*sin(wind))^2),...
        -sin(ori)*(tanh(actuate(1))-Cdh*(sx^2+sy^2)),...
        -Cda*Cw*sin(wind)*sqrt((sx-Cw*cos(wind))^2+(sy-Cw*sin(wind))^2)-Cda*(sx-Cw*cos(wind))^2*Cw*sin(wind)/sqrt((sx-Cw*cos(wind))^2+(sy-Cw*sin(wind))^2)+Cda*(sx-Cw*cos(wind))*(sy-Cw*sin(wind))*Cw*cos(wind)/sqrt((sx-Cw*cos(wind))^2+(sy-Cw*sin(wind))^2),...
        0;0,0,...
        -sin(ori)*Cdh*2*sx-Cda*(sy-Cw*sin(wind))*(sx-Cw*cos(wind))/sqrt((sx-Cw*cos(wind))^2+(sy-Cw*sin(wind))^2),...
        -sin(ori)*Cdh*2*sy-Cda*sqrt((sx-Cw*cos(wind))^2+(sy-Cw*sin(wind))^2)-Cda*(sy-Cw*sin(wind))^2/sqrt((sx-Cw*cos(wind))^2+(sy-Cw*sin(wind))^2),...
        cos(ori)*(tanh(actuate(1))-Cdh*(sx^2+sy^2)),...
        Cda*Cw*cos(wind)*sqrt((sx-Cw*cos(wind))^2+(sy-Cw*sin(wind))^2)+Cda*(sy-Cw*sin(wind))^2*Cw*cos(wind)/sqrt((sx-Cw*cos(wind))^2+(sy-Cw*sin(wind))^2)-Cda*(sx-Cw*cos(wind))*(sy-Cw*sin(wind))*Cw*sin(wind)/sqrt((sx-Cw*cos(wind))^2+(sy-Cw*sin(wind))^2),...
        0;zeros(1,7);zeros(1,7);zeros(1,7)];
    A = func2(px,py,sx,sy,ori,wind,drif)*dt+eye(7);
    L = [zeros(1,4);zeros(1,4);...
        -cos(ori)*Cdh*(sx^2+sy^2),0,0,0;...
        -sin(ori)*Cdh*(sx^2+sy^2),0,0,0;...
        0,Cr*actuate(2),0,0;...
        0,0,1,0;...
        0,0,0,1];
    Pm= A*estState.Pm*A.'+L*diag([estConst.DragNoise,estConst.RudderNoise,estConst.WindAngleNoise,estConst.GyroDriftNoise])*L.';
    % measurement update
    za = sqrt((xm(1)-estConst.pos_radioA(1))^2+(xm(2)-estConst.pos_radioA(2))^2);
    zb = sqrt((xm(1)-estConst.pos_radioB(1))^2+(xm(2)-estConst.pos_radioB(2))^2);
    zc = sqrt((xm(1)-estConst.pos_radioC(1))^2+(xm(2)-estConst.pos_radioC(2))^2);
    zg = xm(5)+xm(7);
    zn = xm(5);
    epsi = diag([estConst.DistNoiseA,estConst.DistNoiseB,estConst.DistNoiseC,estConst.GyroNoise,estConst.CompassNoise]);
    Hk = [(xm(1)-estConst.pos_radioA(1))/za,(xm(2)-estConst.pos_radioA(2))/za,0,0,0,0,0;
        (xm(1)-estConst.pos_radioB(1))/zb,(xm(2)-estConst.pos_radioB(2))/zb,0,0,0,0,0;
        (xm(1)-estConst.pos_radioC(1))/zc,(xm(2)-estConst.pos_radioC(2))/zc,0,0,0,0,0;
        0,0,0,0,1,0,1;
        0,0,0,0,1,0,0;
        ];
    K = Pm*Hk.'/(Hk*Pm*Hk.'+epsi);
    sense =  keepMeasuremnt(sense,[za,zb,zc,zg,zn]);
    estState.xm = xm+(K*(sense - [za,zb,zc,zg,zn]).')';
    estState.pm = (eye(7)-K*Hk)*Pm;
    % Get resulting estimates and variances

    % Output quantities
    posEst = estState.xm(1:2);
    linVelEst = estState.xm(3:4);
    oriEst = estState.xm(5);
    windEst = estState.xm(6);
    driftEst = estState.xm(7);
    var = eig(estState.Pm);
    posVar = var(1:2);
    linVelVar = var(3:4);
    oriVar = var(5);
    windVar = var(6);
    driftVar = var(7);
end
end
function sense = keepMeasuremnt(sense,z)
    for i = 1:5
        if(isnan(sense(i))|isinf(sense(i)))
            sense(i) = z(i);
        end
    end
end
function dydt = odefun(t, estState.xm,estConst)
    tmp = num2cell(estState.xm);
    [px,py,sx,sy,ori,wind,drif] = deal(tmp{:});
    dydt = [sx,sy,...
        cos(ori)*(tanh(actuate(1))-Cdh*(sx^2+sy^2))-Cda*(sx-Cw*cos(wind)*sqrt((sx-Cw*cos(wind))^2+(sy-Cw*sin(wind))^2)),...
        sin(ori)*(tanh(actuate(1))-Cdh*(sx^2+sy^2))-Cda*(sy-Cw*sin(wind)*sqrt((sx-Cw*cos(wind))^2+(sy-Cw*sin(wind))^2)),...
        Cr*actuate(2),...
        0,...
        0];
end