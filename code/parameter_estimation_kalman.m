%% initialize system parameters
k = 10;
m = 5;
cf = 0.9;
A = [0 1;-k/m -cf/m];
A2 = [0 1;-8/m -cf/m];
B = [0;0];
C = [1 0];
D = 0;
% process covariance noise
Q_mu = [0 0];
Q_cov = [1 1.5; 1.5 3];
% measurement covariance noise
R_mu = 0;
R_cov = 0.3;

%% system

s1 = ss(A,B,C,D); % LTI system object
[y,t,xs] = step(s1); % plot step response

%% Discrete-time model of the process
dt = 0.1; % sampling time of system
Phi = expm(A*dt);
Phi2 = expm(A2*dt);
H = C;
%% Simulation loop
T = 5000;
n = T/dt;
x = zeros(2,n); % zero the state vector
z = zeros(1,n); % zero the measurement vector
wk = 1; % process noise source
vk = 20; % measurement noise
x(2,1) = 10;
disp(x(1,2));
Gamma = 1;
t = 0:dt:(n-1)*dt; % time vector
for i = 1:n
    if i>1250
        x(:,i+1) = x(:,i)+A2*x(:,i)*dt ;%+ Gamma*wk*mvnrnd(Q_mu,Q_cov,1).';%+ Gamma*wk*mvnrnd(Q_mu,Q_cov,1).';
        
    end
    if i<=1250
        %disp(x(:,i));
        %disp(A*x(:,i));
        x(:,i+1) = x(:,i)+A*x(:,i)*dt ;%+ Gamma*wk*mvnrnd(Q_mu,Q_cov,1).'; %Phi*x(:,i)+ Gamma*wk*mvnrnd(Q_mu,Q_cov,1).'; % generated with input-reffered noise
        %disp('s');
        %disp(Phi2*x(:,i));
        %disp( Gamma*wk*mvnrnd(Q_mu,Q_cov,1).');
    end
    % x(:,i+1) = Phi*x(:,i) + wk*(randn(2,1)); % process noise(direct)
    z(i) = H * x(:,i) + vk*mvnrnd(R_mu,R_cov,1);
end
%%plots

subplot(3,1,1),stairs(t,z)
title('Mass-Spring Damper Kalman Filter Demo')
hold on
stairs(t,x(1,1:end-1),'r')
legend('Measurement','Actual State')
hold off


%% kalman filter
x0 = x(:,1);
k0 = 8;
x_post = [0 10 k0].';
P_post = [0 0 0; 0 0 0;0 0 0];
Q_k = [0.3 0.2 0; 0.4 0.4 0; 0.1 0.1 0 ];
R_k1 = 0.3;
xhat = zeros(3,n);
xhat(:,1) = [x0(1) x0(2) k0 ];
%initial estimate for spring constant

%define initial coefficient matrix
A = [ 0 1 0 ; -k0/m -cf/m  0; 0 0 0 ];

for i=1:n
    x_pri = x_post+ A*x_post*dt ;
    %disp('pri');
    %disp(x_pri);
    %disp(x_post);
    H_k1 = [1 0 0];
    phi_k = A;
    P_pri = phi_k * P_post * (phi_k.') + Q_k;
    K_k1 = P_pri*(H_k1.') /(H_k1*P_pri*(H_k1.') + R_k1);
    Y_k1 = z(i);
    %K_k1 = [0 0 0 ].';
    x_post = x_pri + K_k1*(Y_k1 - H_k1*x_pri);
    disp(Y_k1 - H_k1*x_pri);
    disp(K_k1);
    P_post = (eye(3) - K_k1*H_k1)*P_pri;
    A = [ 0 1 0 ; -x_post(3)/m -cf/m 0; 0 0 0 ]; 
    xhat(:,i+1)=x_post;
end
subplot(3,1,2),stairs(t,z)
hold on
stairs(t,xhat(1,1:end-1),'r')
legend('Measurement','State Estimate')
hold off
subplot(3,1,3),stairs(t,xhat(1,1:end-1),'r')
hold on
stairs(t,x(1,1:end-1),'g')
legend('State Estimate','Actual State')
hold off