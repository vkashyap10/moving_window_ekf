%% initialize system parameters
k = 0.1;
m = 5;
cf = 0.1;
A = [0 1;-k/m -cf/m];
B = [0;1/m];
C = [1 0];
D = 0;

%% system

s1 = ss(A,B,C,D); % LTI system object
[y,t,xs] = step(s1); % plot step response

%% Discrete-time model of the process
dt = 0.1; % sampling time of system
Phi = expm(A*dt);
Gamma = inv(A)*(Phi-eye(2))*B;
H = C;
%% Simulation loop
T = 500;
n = T/dt;
x = zeros(2,n); % zero the state vector
wk = 1; % process noise source
vk = 5; % measurement noise
t = 0:dt:(n-1)*dt; % time vector
for i = 1:n
    x(:,i+1) = Phi*x(:,i) + Gamma*wk*(randn(1,1)); % generated with input-reffered noise
    % x(:,i+1) = Phi*x(:,i) + wk*(randn(2,1)); % process noise(direct)
    z(i) = H * x(:,i) + vk*randn(1);
end
%%plots
subplot(3,1,1),stairs(t,z)
title('Mass-Spring Damper Kalman Filter Demo')
hold on
stairs(t,x(1,1:end-1),'r')
legend('Measurement','Actual State')
hold off
%% Now filter the data with a Kalman filter
Q = [0.031 0;0 0.01]*0.02;
R = 25;
P = [0.5 0; 0 0.5];
xhat = zeros(2,n);
u = 0;
%% Filter loop
for i = 1:n
    [xhat(:,i+1),Pp,K] = kalman(xhat(:,i),z(i),P,Phi,B,u,Q,R,H);
    P = Pp;
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