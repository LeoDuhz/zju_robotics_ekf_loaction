function [] = ekf_localization()
 
% Homework for ekf localization
% Modified by YH on 09/09/2019, thanks to the original open source
% Any questions please contact: zjuyinhuan@gmail.com

    close all;
    clear all;

    disp('EKF Start!')

    time = 0;
    global endTime; % [sec]
    endTime = 60;
    global dt;
    dt = 0.1; % [sec]

    removeStep = 5;

    nSteps = ceil((endTime - time)/dt);

    estimation.time=[];
    estimation.u=[];
    estimation.GPS=[];
    estimation.xOdom=[];
    estimation.xEkf=[];
    estimation.xTruth=[];

    % State Vector [x y yaw]'
    xEkf=[0 0 0]';
    PxEkf = eye(3);

    % Ground True State
    xTruth=xEkf;

    % Odometry Only
    xOdom=xTruth;

    % Observation vector [x y yaw]'
    z=[0 0 0]';

    % Simulation parameter
    global noiseQ
    noiseQ = diag([0.1 0 degreeToRadian(10)]).^2; %[Vx Vy yawrate]

    global noiseR
    noiseR = diag([0.5 0.5 degreeToRadian(5)]).^2;%[x y yaw]
    
    % Covariance Matrix for motion
    convQ=noiseQ.^2;
    
    % Covariance Matrix for observation
    convR=noiseR.^2;
    
    % Other Intial
    sigma_t=eye(3);
    miu_t=[0,0,0]';
    
    % Main loop
    for i=1 : nSteps
        time = time + dt;
        % Input
        u=robotControl(time);
        % Observation
        [z,xTruth,xOdom,u]=prepare(xTruth, xOdom, u);

        % ------ Kalman Filter --------
        % Predict
        %Xn=doMotion(xEkf,u);
        
        %计算两中间变量
        sigma_t_bar=jacobF(miu_t,u)*sigma_t*(jacobF(miu_t,u))'+convQ;%注意这里与课件中的Q、R互换了
        miu_t_bar=doMotion(miu_t,u);
        
        % Update
        Kt=sigma_t_bar*(jacobH(miu_t_bar))'*inv(jacobH(miu_t_bar)*sigma_t_bar*(jacobH(miu_t_bar))'+convR);%计算Kt
        miu_t=miu_t_bar+Kt*(z-doObservation(miu_t_bar));
        sigma_t=(eye(3)-Kt*jacobH(miu_t_bar))*sigma_t_bar;
        
        %得到EKF定位结果
        xEkf=miu_t;
        % -----------------------------

        % Simulation estimation
        estimation.time=[estimation.time; time];
        estimation.xTruth=[estimation.xTruth; xTruth'];
        estimation.xOdom=[estimation.xOdom; xOdom'];
        estimation.xEkf=[estimation.xEkf;xEkf'];
        estimation.GPS=[estimation.GPS; z'];
        estimation.u=[estimation.u; u'];

        % Plot in real time
        % Animation (remove some flames)
        if rem(i,removeStep)==0
            %hold off;
            plot(estimation.GPS(:,1),estimation.GPS(:,2),'*m', 'MarkerSize', 5);hold on;
            plot(estimation.xOdom(:,1),estimation.xOdom(:,2),'.k', 'MarkerSize', 10);hold on;
            plot(estimation.xEkf(:,1),estimation.xEkf(:,2),'.r','MarkerSize', 10);hold on;
            plot(estimation.xTruth(:,1),estimation.xTruth(:,2),'.b', 'MarkerSize', 10);hold on;
            axis equal;
            grid on;
            drawnow;
            %movcount=movcount+1;
            %mov(movcount) = getframe(gcf);
        end 
    end
    close
    
    finalPlot(estimation);
 
end

% control
function u = robotControl(time)
    global endTime;

    T = 10; % sec
    Vx = 1.0; % m/s
    Vy = 0.2; % m/s
    yawrate = 5; % deg/s
    
    % half
    if time > (endTime/2)
        yawrate = -5;
    end
    
    u =[ Vx*(1-exp(-time/T)) Vy*(1-exp(-time/T)) degreeToRadian(yawrate)*(1-exp(-time/T))]';
    
end

% all observations for 
function [z, xTruth, xOdom, u] = prepare(xTruth, xOdom, u)
    global noiseQ;
    global noiseR;

    % Ground Truth
    xTruth=doMotion(xTruth, u);
    % add Motion Noises
    u=u+noiseQ*randn(3,1);
    % Odometry Only
    xOdom=doMotion(xOdom, u);
    % add Observation Noises
    z=xTruth+noiseR*randn(3,1);
end


% Motion Model
function x = doMotion(x, u)
    global dt;
%     x=x+u*dt;
    x=x+[u(1)*cos(x(3))*dt-u(2)*sin(x(3))*dt,u(1)*sin(x(3))*dt+u(2)*cos(x(3))*dt,u(3)*dt]';
end

% Jacobian of Motion Model
function jF = jacobF(x, u)
    global dt;
    jF=[1,0,u(1)*(-sin(x(3)))*dt-u(2)*cos(x(3))*dt;0,1,u(1)*cos(x(3))*dt+u(2)*(-sin(x(3)))*dt;0,0,1];
end

%Observation Model
function z = doObservation(x) %xPred)
    z=x;
 end

%Jacobian of Observation Model
function jH = jacobH(x)
    jH=[1,0,0;0,1,0;0,0,1];
end

% finally plot the results
function []=finalPlot(estimation)
    figure;
    
    plot(estimation.GPS(:,1),estimation.GPS(:,2),'*m', 'MarkerSize', 5);hold on;
    plot(estimation.xOdom(:,1), estimation.xOdom(:,2),'.k','MarkerSize', 10); hold on;
    plot(estimation.xEkf(:,1), estimation.xEkf(:,2),'.r','MarkerSize', 10); hold on;
    plot(estimation.xTruth(:,1), estimation.xTruth(:,2),'.b','MarkerSize', 10); hold on;
    legend('GPS Observations','Odometry Only','EKF Localization', 'Ground Truth');

    xlabel('X (meter)', 'fontsize', 12);
    ylabel('Y (meter)', 'fontsize', 12);
    grid on;
    axis equal;

    % calculate error
    % distance error
    disp('Odemetry distance error:');
    odemetry_distance_error=sqrt(sum((estimation.xTruth(:,1)-estimation.xOdom(:,1)).^2)+sum((estimation.xTruth(:,2)-estimation.xOdom(:,2)).^2));
    disp(odemetry_distance_error);
    disp('My EKF distance error：');
    EKF_distance_error=sqrt(sum((estimation.xTruth(:,1)-estimation.xEkf(:,1)).^2)+sum((estimation.xTruth(:,2)-estimation.xEkf(:,2)).^2));
    disp(EKF_distance_error);
    
    %angle error
    disp('Odemetry angle error:');
    odemetry_angle_error=sum(abs(estimation.xTruth(:,3)-estimation.xOdom(:,3)));
    disp(odemetry_angle_error);
    disp('My EKF angle error:');
    EKF_angle_error=sum(abs(estimation.xTruth(:,3)-estimation.xEkf(:,3)));
    disp(EKF_angle_error);
    
    

end

function radian = degreeToRadian(degree)
    radian = degree/180*pi;
end