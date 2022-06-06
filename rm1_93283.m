%% TODO
% make functions out of everything
% pull custom functions into this file
% properly test omni
% start report
% make code presentable and organized
% save values to file

%% start
close all
clear all
clc
addpath ./Beacon_Robot_Localization ./Beacon_Robot_Localization/lib/
%function rm1_93283(N,Dt,r,L,Vn,Wn,v)
% N - número de faróis (4)
% Dt - intervalo de tempo de amostragem dos sensores (1 s)
% r - raio das rodas dos robôs (0.15 m)
% L - separação/afastamento das rodas conforme o modelo cinemático (1 m)
% Vn - incerteza (sigma) na velocidade linear a impor ao robô (input) (0.1 m/s)
% Wn - incerteza (sigma) na velocidade angular a impor ao robô (input) (0.1 m/s)
% v - velocidade linear usada para calcular a segmentação da trajetória (5 m/s)

%if nargin == 0
N = 4;
Dt = 1;
r = 0.15;
L = 1;
Vn = 0.1;
Wn = 0.1;
v= 5;
%end

%% calculate trajectory points
% initial position
P = [0 0 0];

% obtain trajectory
B = BeaconDetection(N, P);

X = P(1):(v*Dt):B(4).X;
Y = pchip([P(1) B.X], [P(2) B.Y], X);

% orientation has to be 0 at the start
thetas = 0;
for i=2:length(X)-1
    thetas = [thetas atan2(Y(i+1)-Y(i), X(i+1)-X(i))];
end

% if the trajectory from pchip ends before the last beacon, add it
if X(end) ~= B(4).X 
    X = [X, B(4).X];
    Y = [Y, B(4).Y];
    thetas = [thetas atan2(Y(end)-Y(end-1), X(end)-X(end-1))];
end

thetas = [thetas thetas(end)];

plot(X, Y);

%% Calculate velocities
%v_l = sqrt((Y(2)-Y(1))^2 + (X(2)-X(1))^2);
%v_a = atan2(Y(2), X(2));

v_l = [];
v_a = [];

for i=2:length(X)
    v_l = [v_l sqrt((Y(i)-Y(i-1))^2 + (X(i)-X(i-1))^2)];
    v_a = [v_a thetas(i)-thetas(i-1)];
end

v_l = v_l/Dt;
v_a = v_a/Dt;

%% EKF
%Define the noise covariances for EKF
% for process
Q=[ Vn^2 0
    0       Wn^2];

% for one landmark observation
R_i=[ B(1).dn^2  0
      0        B(1).an^2];

% for storing the results
xstate_EKF = [0, zeros(1,3)]; % pose at time 0
P_EKF      = 0.01*eye(3);  % initial covariance matrix

% numbers and absolute positions of the landmarks
landmarkxy = [1 B(1).X B(1).Y
              2 B(2).X B(2).Y
              3 B(3).X B(3).Y
              4 B(4).X B(4).Y];

num_steps = size(v_l,2);
for step = 1:num_steps
    % get the data needed for one-step EKF
    % EKF estimate at time t
    xstate_t = xstate_EKF(end,2:4)';
    P_t = P_EKF(end-2:end,:);
    
    % control input at time t
    control_t= [v_l(step) v_a(step)];
    % observation data at time t+1
    B = BeaconDetection(N, [X(step+1) Y(step+1) thetas(step+1)]);
    
    obs_t1 = [];
    needed_landmarks = [];
    for i=1:N
        if ~isnan(B(i).d)
            needed_landmarks = [needed_landmarks;
                                landmarkxy(i,2:3)];
            obs_t1 = [obs_t1 B(i).d B(i).a];
        end
    end
    %obs_t1 = [1 obs_t1(1).d obs_t1(1).a 2 obs_t1(2).d obs_t1(2).a 3 obs_t1(3).d obs_t1(3).a 4 obs_t1(4).d obs_t1(4).a];
    
    %discretization time interval
    R = zeros(size(obs_t1,1), size(obs_t1,1));
    for i=1:size(obs_t1,2)/2
        R(i*2-1:i*2, i*2-1:i*2) = R_i; % adapt to the number of landmarks considered
    end
    %using ekf function
    [xstateT1_T1,PT1_T1] = ekf(xstate_t,P_t,control_t,obs_t1,needed_landmarks,Dt,Q,R);
    
    %update
    xstate_EKF = [xstate_EKF; Dt, xstateT1_T1];
    P_EKF = [P_EKF; PT1_T1];
end

%% draw the estimated robot poses and uncertainty ellipses

figure(1)
arrow_length=2;
axis([-1 5 -3 6])
%if test100==1
%    axis([-25 25 -10 40])
%end

hold on

plot(landmarkxy(:,2),landmarkxy(:,3),'k*','MarkerSize',14); %,'linewidth',6)
text(landmarkxy(:,2)+0.2,landmarkxy(:,3),num2str(landmarkxy(:,1)),'fontweight','bold','fontsize',14)
grid on

for i=0:num_steps
    uncer_p = P_EKF(i*3+1:i*3+2, 1:2);        % get the xy covariance
    
    uncer_x = xstate_EKF(i+1,2);
    uncer_y = xstate_EKF(i+1,3);
    CV=GetCov(uncer_p,uncer_x,uncer_y);  % by wangzhan, make it large on purpose, not now
    plot(CV(1,:),CV(2,:),'-b');
    
    plot(xstate_EKF(i+1,2),xstate_EKF(i+1,3),'bo','linewidth',2);

    % draw the robot heading
    dx = arrow_length*cos(xstate_EKF(i+1,4));
    dy = arrow_length*sin(xstate_EKF(i+1,4));
    quiver(xstate_EKF(i+1,2),xstate_EKF(i+1,3),...
           dx, dy, 0, 'Color', 'b','linewidth',1.2)
    
    %draw the true robot poses for comparison
    
    plot(X(i+1),Y(i+1),'ro','linewidth',2);
    
    dx = arrow_length*cos(thetas(i+1));
    dy = arrow_length*sin(thetas(i+1));
    quiver(X(i+1),Y(i+1),...
           dx, dy, 0, 'Color', 'r','linewidth',1.2)
    
    %pause
end
%end

%% Save estimated positions
fileID = fopen('loc_93283.txt','w');
formatSpec = '%7.4f,%7.4f,%7.4f\n';
for i=1:size(xstate_EKF,1)
    fprintf(fileID,formatSpec,xstate_EKF(i,2),xstate_EKF(i,3),xstate_EKF(i,4))
end

%% DD - wheels angular velocities

% w_l = [];
% w_r = [];
% 
% for i=2:length(X)
%     [VR,VL] = invkinDD(X(i)-X(i-1),NaN,thetas(i)-thetas(i-1),L,Dt);
%     w_l = [w_l VL];
%     w_r = [w_r VR];
% end
% 
% [Vx, Vy, w]=localvels(1,r,L,w_l(1)/r,w_r(1)/r,0);
% 
% w_t = [];
% alphas = [];
% 
% for i=2:length(X)
%     [WT,alpha] = invkinTRI((X(i)-X(i-1))/Dt,(Y(i)-Y(i-1))/Dt,(thetas(i)-thetas(i-1))/Dt,thetas(i-1),L,r,Dt);
%     w_t = [w_t WT];
%     alphas = [alphas alpha];
% end

ws_1 = [];
ws_2 = [];
ws_3 = [];

for i=2:length(X)
    [w1,w2,w3] = invkinOMNI((X(i)-X(i-1))/Dt,(Y(i)-Y(i-1))/Dt,v_a(i-1),thetas(i-1),L,r);
    ws_1 = [ws_1 w1];
    ws_2 = [ws_2 w2];
    ws_3 = [ws_3 w3];
end

% to check if correct
%[vx, vy,w ] = localvels(3,r,L,ws_1(1),ws_2(1),ws_3(1))

fileID = fopen('OMNI_93283.txt','w');
formatSpec = '%7.4f,%7.4f,%7.4f\n';
for i=1:size(ws_1,1)
    fprintf(fileID,formatSpec,ws_1(i),ws_2(i),ws_3(i))
end
