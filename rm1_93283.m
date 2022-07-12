function rm1_93283(N,Dt,r,L,Vn,Wn,v)
% N - número de faróis (4)
% Dt - intervalo de tempo de amostragem dos sensores (1 s)
% r - raio das rodas dos robôs (0.15 m)
% L - separação/afastamento das rodas conforme o modelo cinemático (1 m)
% Vn - incerteza (sigma) na velocidade linear a impor ao robô (input) (0.1 m/s)
% Wn - incerteza (sigma) na velocidade angular a impor ao robô (input) (0.1 m/s)
% v - velocidade linear usada para calcular a segmentação da trajetória (5 m/s)
    
    if nargin == 0
        N = 4;
        Dt = 1;
        r = 0.15;
        L = 1;
        Vn = 0.1;
        Wn = 0.1;
        v= 5;
    end
    
    %% calculate trajectory points
    % initial position
    P = [0 0 0];
    
    % obtain trajectory
    B = BeaconDetection(N, P);
    
    X = P(1):(v*Dt):B(end).X;
    Y = pchip([P(1) B.X], [P(2) B.Y], X);
    
    % orientation has to be 0 at the start
    thetas = 0;
    for i=2:length(X)
        thetas = [thetas atan2(Y(i)-Y(i-1), X(i)-X(i-1))];
    end
    
    % if the trajectory from pchip ends before the last beacon, add it
    if X(end) ~= B(end).X 
        X = [X, B(end).X];
        Y = [Y, B(end).Y];
        thetas = [thetas atan2(Y(end)-Y(end-1), X(end)-X(end-1))];
    end
    
    %% Calculate velocities    
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
    P_EKF = 0.01*eye(3);  % initial covariance matrix
    
    % numbers and absolute positions of the landmarks
    landmarkxy = [];
    for i=1:N
        landmarkxy = [landmarkxy; i B(i).X B(i).Y];
    end
    
    for step = 1:size(v_l,2)
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
    
    drawGraph(landmarkxy, xstate_EKF, P_EKF, X, Y, thetas);
    
    %% Save estimated positions
    fileID = fopen('loc_93283.txt','w');
    formatSpec = '%7.4f,%7.4f,%7.4f\n';
    for i=1:size(xstate_EKF,1)
        fprintf(fileID,formatSpec,xstate_EKF(i,2),xstate_EKF(i,3),xstate_EKF(i,4));
    end
    
    %% DD - wheels angular velocities
    w_l = [];
    w_r = [];
    
    for i=1:size(v_l,2)
        [WR,WL] = invkinDD(v_l(i),v_a(i),L,r);
        w_l = [w_l WL];
        w_r = [w_r WR];
    end
    
    fileID = fopen('DD_93283.txt','w');
    formatSpec = '%7.4f,%7.4f\n';
    for i=1:size(v_l,2)
        fprintf(fileID,formatSpec,w_r(i),w_l(i));
    end
    
    %% TRI - traction wheel angular velocity and wheel angle
    w_t = [];
    alphas = [];
    
    for i=1:size(v_l,2)
        [WT,alpha] = invkinTRI(v_l(i),v_a(i),L,r);
        w_t = [w_t WT];
        alphas = [alphas alpha];
    end
    
    fileID = fopen('TRI_93283.txt','w');
    formatSpec = '%7.4f,%7.4f\n';
    for i=1:size(v_l,2)
        fprintf(fileID,formatSpec,w_t(i),alphas(i));
    end
    
    %% OMNI - wheels angular velocities
    w_1 = [];
    w_2 = [];
    w_3 = [];
    
    for i=1:size(v_l,2)
        [w1,w2,w3] = invkinOMNI(v_l(i),0,v_a(i),L,r);
        w_1 = [w_1 w1];
        w_2 = [w_2 w2];
        w_3 = [w_3 w3];
    end
    
    fileID = fopen('OMNI_93283.txt','w');
    formatSpec = '%7.4f,%7.4f,%7.4f\n';
    for i=1:size(v_l,2)
        fprintf(fileID,formatSpec,w_1(i),w_2(i),w_3(i));
    end

    %% Test inverse kynematics
    disp("Test inverse kynematics")
    for i=1:size(v_l,2)
        disp("Iteration " + i)
        disp("Original velocities: [" + round(v_l(i),4) + " " + 0 + " " + round(v_a(i),4) + "]")

        [Vx, Vy, w]=localvels(1,r,L,w_r(i),w_l(i),0);
        disp("After reverting invkinDD: [" + round(Vx,4) + " " + round(Vy,4) + " " + round(w,4) + "]")

        [Vx, Vy, w]=localvels(2,r,L,w_t(i),alphas(i),0);
        disp("After reverting invkinTRI: [" + round(Vx,4) + " " + round(Vy,4) + " " + round(w,4) + "]")

        [Vx, Vy,w ] = localvels(3,r,L,w_1(i),w_2(i),w_3(i));
        disp("After reverting invkinOMNI: [" + round(Vx,4) + " " + round(Vy,4) + " " + round(w,4) + "]")

        disp("--------------")
    end
end

function [xstate_t1,P_t1] = ekf(xstate_t,P_t,control_t,obs_t1,landmarkxym,Delta_T,Q,R)
% based on the ekf function from classes
    %% prediction step
    %motion model
    pzero =[0 0]; %set noise equal to 0
    
    %motion model estimate
    [xstatet1_t] = motionmodel(xstate_t,control_t,pzero,Delta_T);

    %predicted covariance matrix (uncertainty)
    %jacobian matrix for F
    
    temp = -Delta_T*control_t(1)*sin(xstate_t(3));
    temp2= Delta_T*control_t(1)*cos(xstate_t(3));
    Jfx=[1 0 temp
         0 1 temp2
         0 0 1     
         ];
    
    temp3 =  Delta_T*control_t(1)*cos(xstate_t(3));
    temp4 =  Delta_T*control_t(1)*sin(xstate_t(3));
    
    Jfw=[temp3 0 
         temp4 0
         0     Delta_T
         ];
    
    Pt1_t= Jfx*P_t*Jfx'+Jfw*Q*Jfw';  %uncertainty
    
    %% update step
    
    % get the observations as a column vector
    z_all = obs_t1';
    
    % set noise equal to 0
    nzero = [0 0];   
    
    % predicted observation
    z_pred = [];
    for i=1:size(landmarkxym,1)
        z_pred = [z_pred; sensormodel(landmarkxym(i,:),xstatet1_t,nzero)'];
    end
    
    % innovation
    innov = z_all-z_pred;

    % wrap the angles to [-pi, pi]
    for i=1:size(innov,1)/2
        innov(i*2)=wrap(innov(i*2)); 
    end
    
    %jacobian matrix
    Jh = [];
    for i=1:size(landmarkxym,1)
        Jh = [Jh; jacobi(landmarkxym(i,:),xstatet1_t(1),xstatet1_t(2))];
    end
    
    S = Jh*Pt1_t*Jh'+R;
    K = Pt1_t*Jh'*pinv(S);
    
    %result
    xstatet1_t1 = xstatet1_t'+K*innov;
    Pt1_t1      = Pt1_t - K*Jh*Pt1_t;
    
    xstate_t1 = xstatet1_t1';
    P_t1      = Pt1_t1;
end

function alpha = wrap(alpha)
% based on the wrap function from classes
% make sure angle is in [-pi,pi] radians
	while (alpha > pi)
		alpha = alpha - 2 * pi;
    end

	while (alpha < -pi)
		alpha = alpha + 2 * pi;
    end
end

function drawGraph(landmarkxy, xstate_EKF, P_EKF, X, Y, thetas)
% based on code from classes
% draw the estimated robot poses and uncertainty ellipses
    figure(1)
    
    hold on
    grid on
    
    % mark the landmarks
    plot(landmarkxy(:,2),landmarkxy(:,3),'k*','MarkerSize',14);
    text(landmarkxy(:,2)+0.2,landmarkxy(:,3),num2str(landmarkxy(:,1)),'fontweight','bold','fontsize',14)
    
    arrow_length=2;
    for i=0:size(X, 2)-1
        uncer_p = P_EKF(i*3+1:i*3+2, 1:2); % get the xy covariance
        
        uncer_x = xstate_EKF(i+1,2);
        uncer_y = xstate_EKF(i+1,3);
        CV=GetCov(uncer_p,uncer_x,uncer_y);
        plot(CV(1,:),CV(2,:),'-b');
        
        % draw the estimated robot poses
        plot(xstate_EKF(i+1,2),xstate_EKF(i+1,3),'bo','linewidth',2);

        dx = arrow_length*cos(xstate_EKF(i+1,4));
        dy = arrow_length*sin(xstate_EKF(i+1,4));
        quiver(xstate_EKF(i+1,2),xstate_EKF(i+1,3),dx, dy, 0, 'Color', 'b','linewidth',1.2);
        
        % draw the true robot poses for comparison
        plot(X(i+1),Y(i+1),'ro','linewidth',2);
        
        dx = arrow_length*cos(thetas(i+1));
        dy = arrow_length*sin(thetas(i+1));
        quiver(X(i+1),Y(i+1),dx, dy, 0, 'Color', 'r','linewidth',1.2);
    end
end

function [w_r,w_l] = invkinDD(v_x,w,L,r)
% v_x - speed in x (meters/second)
% w - angular speed (radians/second)
% L - wheel separation (meters)
% r - wheel radius (meters)
    vels = [
        1/r -L/(2*r)
        1/r L/(2*r)
        ] * [v_x w]';
    w_l = vels(1);
    w_r = vels(2);
end

function [w_t,alpha] = invkinTRI(v_x,w,L,r)
% v_x - speed in x (meters/second)
% w - angular speed (radians/second)
% L - wheel separation (meters)
% r - wheel radius (meters)
    alpha = atan2(w*L, v_x);
    w_t = v_x/(r*cos(alpha));
end

function [w1,w2,w3] = invkinOMNI(v_x,v_y,w,L,r)
% v_x - speed in x (meters/second)
% v_y - speed in y (meters/second)
% w - angular speed (radians/second)
% L - wheel separation (meters)
% r - wheel radius (meters)
    vels = [
        0 -1/r L/r
        sqrt(3)/(2*r) 1/(2*r) L/r
        -sqrt(3)/(2*r) 1/(2*r) L/r
    ]*[v_x v_y w]';

    w1 = vels(1);
    w2 = vels(2);
    w3 = vels(3);
end