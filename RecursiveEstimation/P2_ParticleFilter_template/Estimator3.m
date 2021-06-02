function [postParticles] = Estimator3(prevPostParticles, sens, act, estConst, km)
% The estimator function. The function will be called in two different
% modes: If km==0, the estimator is initialized. If km > 0, the
% estimator does an iteration for a single sample time interval using the 
% previous posterior particles passed to the estimator in
% prevPostParticles and the sensor measurement and control inputs.
%
% Inputs:
%   prevPostParticles   previous posterior particles at time step k-1
%                       The fields of the struct are [1xN_particles]-vector
%                       (N_particles is the number of particles) 
%                       corresponding to: 
%                       .x_r: x-locations of the robot [m]
%                       .y_r: y-locations of the robot [m]
%                       .phi: headings of the robot [rad]
%                       .kappa: wall offset [m]
%                           
%   sens                Sensor measurement z(k), scalar
%
%   act                 Control inputs u(k-1), [1x2]-vector
%                       act(1): u_f, forward control input
%                       act(2): u_phi, angular control input
%
%   estConst            estimator constants (as in EstimatorConst.m)
%
%   km                  time index k, scalar
%                       corresponds to continous time t = k*Ts
%                       If km==0 initialization, otherwise estimator
%                       iteration step.
%
% Outputs:
%   postParticles       Posterior particles at time step k
%                       The fields of the struct are [1xN_particles]-vector
%                       (N_particles is the number of particles) 
%                       corresponding to: 
%                       .x_r: x-locations of the robot [m]
%                       .y_r: y-locations of the robot [m]
%                       .phi: headings of the robot [rad]
%                       .kappa: wall offset [m]
%
%
% Class:
% Recursive Estimation
% Spring 2021
% Programming Exercise 2
%
% --
% ETH Zurich
% Institute for Dynamic Systems and Control
% Raffaello D'Andrea, Matthias Hofer, Carlo Sferrazza
% hofermat@ethz.ch
% csferrazza@ethz.ch

% Set number of particles:
N_particles = 5000; % obviously, you will need more particles than 10.

%% Mode 1: Initialization
if (km == 0)
    % Do the initialization of your estimator here!
    randAB = randi([1,2],1,N_particles);
    randr = rand(1,N_particles)*estConst.d;
    randrio = rand(1,N_particles)*pi*2;
    x = [estConst.pA(1),estConst.pB(1)];
    y = [estConst.pA(2),estConst.pB(2)];
    postParticles.x_r = x(randAB)+randr.*cos(randrio); % 1xN_particles matrix
    postParticles.y_r = y(randAB)+randr.*sin(randrio); % 1xN_particles matrix
    postParticles.phi = 2*rand(1,N_particles).*estConst.phi_0-estConst.phi_0; % 1xN_particles matrix
    postParticles.kappa = rand(1,N_particles).*estConst.l; % 1xN_particles matrix
    % and leave the function
    return;
end % end init

%% Mode 2: Estimator iteration.
% If km > 0, we perform a regular update of the estimator.
% Implement your estimator here!
xr = prevPostParticles.x_r;
yr = prevPostParticles.y_r;
phi = prevPostParticles.phi;
kappam = prevPostParticles.kappa;

% Prior Update:
vphi = rand(1,N_particles)*estConst.sigma_phi-estConst.sigma_phi/2;
vf = rand(1,N_particles)*estConst.sigma_f-estConst.sigma_f/2;
xm = xr+(act(1)+vf).*cos(phi);
ym = yr+(act(1)+vf).*sin(phi);
phim = phi+act(2)+vphi;

% Posterior Update:
[~,~,zm] = intersection(xm,ym,phim,kappam,estConst);
beta = zeros(1,N_particles);

for i =1:N_particles
    beta(i) = measureNoise(sens-zm(i),estConst.epsilon);
end

sum1 = sum(beta,1);
sum3 = sum(beta,'all');

for k=1:3
    if abs(sum3) < 1e-10
        disp('degenerate');
        % degenerate
        xm = xm + 0.1 * k * rand(1, N_particles) - 0.05 * k;
        ym = ym + 0.1 * k * rand(1, N_particles) - 0.05 * k;
        phim = phim + 0.1 * k * pi * rand(1, N_particles) - 0.05 * k * pi;
        kappam = kappam + 0.2 * k * estConst.l * rand(1, N_particles) - 0.1 * k * estConst.l;
        [xm, ym, kappam] = bound(xm, ym, kappam, N_particles, estConst);
        [~, ~, zm] = intersection(xm,ym,phim,kappam,estConst);
        for i = 1:N_particles
            beta(i) = measureNoise(sens-zm(i), estConst.epsilon);
        end
        sum1 = sum(beta,1);
        sum3 = sum(beta,'all');
    end
end

if abs(sum3) < 1e-10
    % still degenerate
    disp('still degenerate');
    beta = ones(N_particles,1);
    sum1 = sum(beta,1);
    sum3 = sum(beta);
end

betap = sum1/sum3;
tmp =zeros(4,N_particles);
s2 = cumsum(betap);
r2 = rand(1,N_particles);

for i =1:N_particles
    k2 = find(s2>r2(i),1);
    tmp(:,i) = [xm(k2);ym(k2);phim(k2);kappam(k2)];
end

% roughening
k = 0.05;
d = 4;
tmp(1,:) = tmp(1,:) + k*N_particles^(-1/d)*(max(tmp(1,:))-min(tmp(1,:)))*randn(N_particles,1)';
tmp(2,:) = tmp(2,:) + k*N_particles^(-1/d)*(max(tmp(2,:))-min(tmp(2,:)))*randn(N_particles,1)';
tmp(3,:) = tmp(3,:) + k*N_particles^(-1/d)*(max(tmp(3,:))-min(tmp(3,:)))*randn(N_particles,1)';
tmp(4,:) = tmp(4,:) + k*N_particles^(-1/d)*(max(tmp(4,:))-min(tmp(4,:)))*randn(N_particles,1)';

if tmp(4,:) < -estConst.l
    tmp(4,:) = -estConst.l;
elseif tmp(4,:) > estConst.l
    tmp(4,:) = estConst.l;
end

postParticles.x_r = tmp(1,:);
postParticles.y_r = tmp(2,:);
postParticles.phi = tmp(3,:);
postParticles.kappa = tmp(4,:);

end % end estimator

function [x,y,z]= intersection(xm,ym,phim,kappam,p)
    xp = p.contour(:,1);
    yp = p.contour(:,2);
   
    xmin = zeros(size(xp,1),1);
    xmax = zeros(size(xp,1),1);
    ymin = zeros(size(xp,1),1);
    ymax = zeros(size(xp,1),1);
    
    x = Inf(1,size(xm,2));
    y = Inf(1,size(xm,2));
    z = Inf(1,size(xm,2));
    
    for j=1:size(xm,2)
        xp(8:9,1)=kappam(j);
        for i=1:size(xp,1)
            xmin(i) = min(xp(i),xp(mod(i,size(xp,1))+1));
            xmax(i) = max(xp(i),xp(mod(i,size(xp,1))+1));
            ymin(i) = min(yp(i),yp(mod(i,size(xp,1))+1));
            ymax(i) = max(yp(i),yp(mod(i,size(xp,1))+1));
            if(xmin(i)==xmax(i))
                 if abs(tan(phim(j))) > 1e10
                     continue;
                 else
                    tmp = [1,0;tan(phim(j)),-1]\[xmin(i);xm(j)*tan(phim(j))-ym(j)]; 
                    tmp_x = tmp(1);
                    tmp_y = tmp(2);
                    if((tmp_x-xm(j))*cos(phim(j))<0)
                        continue;
                    end
                 end
            else
                if(ymax(i)==ymin(i))
                    if abs(tan(phim(j))) < 1e-10
                        continue;
                    else
                        tmp = [0,1;sin(phim(j)),-cos(phim(j))]\[ymin(i);xm(j)*sin(phim(j))-ym(j)*cos(phim(j))];
                        tmp_x = tmp(1);
                        tmp_y = tmp(2);
                        if((tmp_y-ym(j))*sin(phim(j))<0)
                            continue;
                        end
                     end
                else
                    k = (yp(i)-yp(mod(i+1,size(xp,1))))/(xp(i)-xp(mod(i+1,size(xp,1))));
                    b = yp(i)-k*xp(i);
                    if abs(abs(k) - abs(tan(phim(j)))) < 0.01
                        continue;
                    end
                    tmp = [k,-1;sin(phim(j)),-cos(phim(j))]\[-b;xm(j)*sin(phim(j))-ym(j)*cos(phim(j))];
                    tmp_x = tmp(1);
                    tmp_y = tmp(2);
                    if((tmp_y-ym(j))*sin(phim(j))<0)
                        continue;
                    end
                end
            end
            if (sqrt((tmp_x-xm(j))^2+(tmp_y-ym(j))^2)<z(j)&&(xmin(i)<=tmp_x)&&(tmp_x<=xmax(i))&&(ymax(i)>=tmp_y)&&(ymin(i)<=tmp_y))
                z(j) = sqrt((tmp_x-xm(j))^2+(tmp_y-ym(j))^2);
                x(j) = tmp_x;
                y(j) = tmp_y;
            end

        end
    end
end

function pw = measureNoise(w,epsilon)
    if(isinf(w))
        pw=0;
        return
    end
    w = abs(w);
    if(w<2*epsilon)
        pw = 0.4/epsilon-0.2*w/(epsilon^2);
    else
        if (w<2.5*epsilon)
            pw = 0.4/(epsilon^2)*w-0.8/epsilon;
        else
            if (w<3*epsilon)
                pw = -0.4/(epsilon^2)*w+1.2/epsilon;
            else
                pw=0;
            end
        end
    end
end

function [xm, ym, kappam] = bound(xm, ym, kappam, N_particles, estConst)
% bound and also add noise to anamoly points
    for i = 1:N_particles
        if xm(i)>3
            xm(i) = 3-0.025*rand(1);
        end
        if xm(i)<-estConst.l
            xm(i) = -estConst.l+0.025*rand(1);
        end
        if ym(i)>3
            ym(i) = 3-0.025*rand(1);
        end
        if ym(i)<0
            ym(i) = 0.025*rand(1);
        end
        if kappam(i) < -estConst.l
            kappam(i) = -estConst.l+0.025*rand(1);
        end
        if kappam(i) > estConst.l
            kappam(i) = estConst.l-0.025*rand(1);
        end
    end
end