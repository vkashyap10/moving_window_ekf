%% initialize system parameters
k = 3;
m = 100;
k2 = 0.1;
cf = 0.01;
A = [0 1;-k/m -cf/m];
A2 = [0 1;-k2/m -cf/m];
B = [0;0];
C = [1 0];
D = 0;
% process covariance noise
Q_mu = [0 0];
Q_cov = [0.4 0.4; 0.4 0.4];
% measurement covariance noise
R_mu = 0;
R_cov = 0.4;

%% system

s1 = ss(A,B,C,D); % LTI system object
[y,t,xs] = step(s1); % plot step response

%% Discrete-time model of the process
dt = 0.01; % sampling time of system
Phi = expm(A*dt);
Phi2 = expm(A2*dt);
H = C;
%% Simulation loop
T = 500;
n = T/dt;
x = zeros(2,n); % zero the state vector
z = zeros(1,n); % zero the measurement vector
wk = 1; % process noise source
vk = 5; % measurement noise
x(2,1) = 0.1;
t = 0:dt:(n-1)*dt; % time vector
k_val = zeros(1,n);

for i = 1:n
    if abs(x(1,i))<=0.5
        Gamma = inv(A2)*(Phi2-eye(2))*B;
        x(:,i+1) = x(:,i)+A2*x(:,i)*dt + wk*mvnrnd(Q_mu,Q_cov)*Gamma;%+ Gamma*wk*mvnrnd(Q_mu,Q_cov,1).';
        A2 = [0 1;-k2/m -cf/m];
        k_val(1,i) = k2;
    else
        Gamma = inv(A)*(Phi-eye(2))*B;
        x(:,i+1) = x(:,i)+A*x(:,i)*dt + wk*mvnrnd(Q_mu,Q_cov)*Gamma; %Phi*x(:,i)+ Gamma*wk*mvnrnd(Q_mu,Q_cov,1).'; % generated with input-reffered noise
        k_val(1,i) = k;
        A = [0 1;-k/m -cf/m];

    end
    % x(:,i+1) = Phi*x(:,i) + wk*(randn(2,1)); % process noise(direct)
    z(i) = H * x(:,i) + vk*H*x(:,i)*5.0/100*mvnrnd(R_mu,R_cov,1);
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
k0 = 0;
x_post = [0 0 k0].';
P_post = [0.1 0.1 0.1; 0.1 0.1 0.1;0.1 0.1 0.1];
Q_k = [0.01 0.01 0.01; 0.01 0.01 0.01 ; 0.01 0 0 ];
R_k1 = 21;
xhat = zeros(3,n);
xhat(:,1) = [x0(1) x0(2) k0 ];
%initial estimate for spring constant

%define initial coefficient matrix
A = [ 0 1 0 ; -k0/m -cf/m  0; 0 0 0 ];
lambda = 0.5;
% window for moving window KF
window = 10;
% variables for measurement covariance estimation
e_window = zeros(1,window);
e_exp = zeros(1,window);
% variables for process covariance estimation
e_k = zeros(1,window);
Rk = 0;
for i=1:n
    x_pri = x_post+ A*x_post*dt ;
    H_k1 = [1 0 0];
    phi_k = A;
    P_pri = phi_k * P_post * (phi_k.') + Q_k;
    % measure covariance noise calculation
    if i>1
        Yksmooth = (1-lambda)*z(i-1) + lambda*z(i);
    else
        Yksmooth = z(i);
    end

    if i<=window
        e_window(1,i) = z(i) - Yksmooth;
        e_exp(1,i) = sum(e_window)/i;
        Rk = ((e_window-e_exp)*(e_window-e_exp).')/i;
    else
        e_window = circshift(e_window,-1,1);
        e_window(1,window) = z(i) - Yksmooth;
        e_exp = circshift(e_exp,-1,1);
        e_exp(1,window) = sum(e_window)/window;
        Rk = ((e_window-e_exp)*(e_window-e_exp).')/(window-1);
    end
    
    K_k1 = P_pri*(H_k1.') /(H_k1*P_pri*(H_k1.') + Rk);
    Y_k1 = z(i);
    
    % process covariance estimation
    if i<=window
        e_k(1,i) = z(i) - H_k1*x_pri;
        C0 = sum(e_k)/i;
    else
        e_k = circshift(e_k,-1,1);
        e_k(1,window) = z(i) - H_k1*x_pri;
        C0 = sum(e_k)/window;
    end
    a = [1e-1 1e-1 1e-1].';
    Q_k = Q_k*exp(a.'*H_k1.'*(C0 - Rk - H_k1*P_pri*H_k1.')*H_k1*a);
    
    x_post = x_pri + K_k1*(Y_k1 - H_k1*x_pri);
    P_post = (eye(3) - K_k1*H_k1)*P_pri;
    A = [ 0 1 0 ; -x_post(3)/m -cf/m 0; 0 0 0 ]; 
    xhat(:,i+1)=x_post;
    if i>5
        xhat(3,i+1) = (x_post(3) + xhat(3,i) + xhat(3,i-1) + xhat(3,i-2) + xhat(3,i-3))/5.0;
    end
    
   
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


f1= figure;
figure(f1);
plot(t,xhat(3,2:end))
hold on
plot(t,k_val)
hold off

