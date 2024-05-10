clear all; clc
%% define systems
m1 = 2;
m2 = 1;
g = 10;
L = 0.75;

A = [0 1 0 0;g*(m1+m2)/(m1*L) 0 0 0;0 0 0 1;g*m2/m1 0 0 0];
B = [0;1/(L*m1);0;1/m1];
C = [1 0 0 0; 0 0 1 0];
eig(A)
cm = ctrb(A,B);
rank(cm)
om = obsv(A,C);
rank(om)

%% initialization
x(:,1) = [0.1; 0; 0; 0]; 
y(:,1) = C*x(:,1);
T1 = 20;
T2 = 20;
t = 0:0.01:T1; %0.01 time span of interest
nt1 = length(t); % number of time steps
t = 0:0.01:T2;
nt2 = length(t);
dt = t(2) - t(1);
t = 0:0.01:T1+T2;
nt = length(t);
% reference input
for i = 1:nt
    r(:,i) = 1;
    u(:,i) = 1;
end

%% linear quadratic optimal control 1
Q = C'*C; R = eye (1);
K = lqr(A,B,Q,R)
A3 = A - B*K;
eig(A3)
kg = 0; %-inv(C*inv(A3)*B);
u(:,1) = -K*x(:,1) + kg*r(:,1);

for i = 1:nt1-1
x_dot(:,i) = A*x(:,i) + B*u(:,i);
x(:,i+1) = x(:,i) + x_dot(:,i)*dt;
y(:,i+1) = C*x(:,i+1);
u(:,i+1) = -K*x(:,i+1) + kg*r(:,i+1);
end

% perturbation
m1 = 2;
m2 = m2 + 1;
g = 10;
L = 0.75;

% x(1:,nt1)

A = [0 1 0 0;g*(m1+m2)/(m1*L) 0 0 0;0 0 0 1;g*m2/m1 0 0 0];
B = [0;1/(L*m1);0;1/m1];
C = [1 0 0 0; 0 0 1 0];

for i = nt1:nt1+nt2-2
x_dot(:,i) = A*x(:,i) + B*u(:,i);
x(:,i+1) = x(:,i) + x_dot(:,i)*dt;
y(:,i+1) = C*x(:,i+1);
u(:,i+1) = -K*x(:,i+1) + kg*r(:,i+1);
end

figure
plot(t,x(1,:),'b',t,x(2,:),'r',t,x(3,:),'g',t,x(4,:),'m','linewidth',2)
xline(T1,'k--',{'m2 = m2 + 1'})
set(gca,'fontsize',18)
title('LQR control for $y = [x_1, x_3]^T$','Interpreter', 'latex')
legend({'$x_1$','$x_2$','$x_3$','$x_4$'},'Interpreter', 'latex')
legend boxoff
xlabel('Time (s)')

% %% linear quadratic optimal control 2
% C = [100 0 0 0; 0 0 10 0];
% Q = C'*C; R = eye (1);
% K = lqr(A,B,Q,R)
% A3 = A - B*K;
% eig(A3)
% kg = 0; %-inv(C*inv(A3)*B);
% u(:,1) = -K*x(:,1) + kg*r(:,1);
% for i = 1:nt-1
% x_dot(:,i) = A*x(:,i) + B*u(:,i);
% x(:,i+1) = x(:,i) + x_dot(:,i)*dt;
% y(:,i+1) = C*x(:,i+1);
% u(:,i+1) = -K*x(:,i+1) + kg*r(:,i+1);
% end
% 
% figure
% plot(t,x(1,:),'b',t,x(2,:),'r',t,x(3,:),'g',t,x(4,:),'m','linewidth',2)
% set(gca,'fontsize',18)
% title('LQR control for $y = [100 x_1, 10 x_3]^T$','Interpreter', 'latex')
% legend({'$x_1$','$x_2$','$x_3$','$x_4$'},'Interpreter', 'latex')
% legend boxoff
% xlabel('Time (s)')
% 
% figure
% plot(t,u,'k--','LineWidth',2)
% set(gca,'fontsize',18)
% title('LQR optimal input')
% legend({'$u$'},'Interpreter', 'latex')
% legend boxoff
