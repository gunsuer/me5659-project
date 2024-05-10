clear all; clc
% define systems
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
rank (om)
% initialization
x(:,1) = [0.1;0;0;0]; 
y(:,1) = C*x(:,1);
x_hat(:,1) = [0;0;0;0]; 
y_hat(:,1) = C*x_hat(:,1);
t = 0:0.01:20; %0.01 time span of interest
nt = length(t); % number of time steps
dt = t(2) - t(1);
% reference input
for i = 1:nt
    r(:,i) = 1;
    u(:,i) = 1;
end

%% case 1: state estimation
% eigenvalues of(A - LC) 
p1 = -1; p2 = -2; p3 = -3; p4 = -4;
L = place(A',C',[p1 p2 p3 p4]);L = L'; 
L
A1 = A - L*C;
eig(A1) 
for i = 1:nt-1
x_dot(:,i) = A*x(:,i)+B*u(:,i);
x(:,i+1) = x(:,i) + x_dot(:,i)*dt;
y(:,i+1) = C*x(:,i+1);
x_hat_dot(:,i) = A*x_hat(:,i) + B*u(:,i) + L*(y(:,i)- C*x_hat(:,i));
x_hat(:,i+1) = x_hat(:,i) + x_hat_dot(:,i)*dt;
y_hat(:,i+1) = C*x_hat(:,i+1);
end

figure
plot(t,x(1,:),'r',t,x_hat(1,:),'b-.','linewidth',2)
set(gca,'fontsize',18)
legend({'$x_1$','${\hat x}_1$'},'Interpreter', 'latex')
legend boxoff
title('open loop')
xlabel('Time (s)')
print(gcf,'fig1.png','-dpng','-r300');
figure
plot(t,x(2,:),'r',t,x_hat(2,:),'b-.','linewidth',2)
set(gca,'fontsize',18)
legend({'$x_2$','${\hat x}_2$'},'Interpreter', 'latex')
legend boxoff
title('open loop')
xlabel('Time (s)')
print(gcf,'fig2.png','-dpng','-r300');
figure
plot(t,x(3,:),'r',t,x_hat(3,:),'b-.','linewidth',2)
set(gca,'fontsize',18)
legend({'$x_3$','${\hat x}_3$'},'Interpreter', 'latex')
legend boxoff
title('open loop')
xlabel('Time (s)')
print(gcf,'fig3.png','-dpng','-r300');
figure
plot(t,x(4,:),'r',t,x_hat(4,:),'b-.','linewidth',2)
set(gca,'fontsize',18)
legend({'$x_4$','${\hat x}_4$'},'Interpreter', 'latex')
legend boxoff
title('open loop')
xlabel('Time (s)')
print(gcf,'fig4.png','-dpng','-r300');


%% case 2: observer regulation control
% eigenvalues of(A - BK) 
Mp = 2; % percent overshoot
Ts = 2; % transient time
Mlog = log(Mp/100);
MlogSquared = Mlog^2;
zeta = sqrt(MlogSquared/(pi^2+MlogSquared));
w0 = 4/(Ts*zeta);
P = roots([1 2*zeta*w0 w0^2])
p1 = P(1);
p2 = P(2);
P = [p1 p2 -20 -21]
K = place(A,B,P)
A1 = A - B*K;
eig(A1)
% eigenvalues of(A - LC) 
P2 = eig(A1)-[10;10;10;10];
L = place(A',C',P2);L = L'; 
L
A2 = A - L*C;
eig(A2) 

u(:,1) = -K*x_hat(:,1);
for i = 1:nt-1
x_dot(:,i) = A*x(:,i)+ B*u(:,i);
x(:,i+1) = x(:,i) + x_dot(:,i)*dt;
y(:,i+1) = C*x(:,i+1);
x_hat_dot(:,i) = A*x_hat(:,i)+ B*u(:,i)+ L*(y(:,i)- C*x_hat(:,i));
x_hat(:,i+1) = x_hat(:,i) + x_hat_dot(:,i)*dt;
y_hat(:,i+1) = C*x_hat(:,i+1);
u(:,i+1) = -K*x_hat(:,i+1);
end

figure
plot(t,x(1,:),'r',t,x_hat(1,:),'b-.','linewidth',2)
set(gca,'fontsize',18)
legend({'$x_1$','${\hat x}_1$'},'Interpreter', 'latex')
legend boxoff
title('observer regulation control')
xlabel('Time (s)')
print(gcf,'fig5.png','-dpng','-r300');
figure
plot(t,x(2,:),'r',t,x_hat(2,:),'b-.','linewidth',2)
set(gca,'fontsize',18)
legend({'$x_2$','${\hat x}_2$'},'Interpreter', 'latex')
legend boxoff
title('observer regulation control')
xlabel('Time (s)')
print(gcf,'fig6.png','-dpng','-r300');
figure
plot(t,x(3,:),'r',t,x_hat(3,:),'b-.','linewidth',2)
set(gca,'fontsize',18)
legend({'$x_3$','${\hat x}_3$'},'Interpreter', 'latex')
legend boxoff
title('observer regulation control')
xlabel('Time (s)')
print(gcf,'fig7.png','-dpng','-r300');
figure
plot(t,x(4,:),'r',t,x_hat(4,:),'b-.','linewidth',2)
set(gca,'fontsize',18)
legend({'$x_4$','${\hat x}_4$'},'Interpreter', 'latex')
legend boxoff
title('observer regulation control')
xlabel('Time (s)')
print(gcf,'fig8.png','-dpng','-r300');


%% case 3: observer tracking control
% eigenvalues of(A - BK) 
Mp = 2; % percent overshoot
Ts = 2; % transient time
Mlog = log(Mp/100);
MlogSquared = Mlog^2;
zeta = sqrt(MlogSquared/(pi^2+MlogSquared));
w0 = 4/(Ts*zeta);
P = roots([1 2*zeta*w0 w0^2])
p1 = P(1);
p2 = P(2);
P = [p1 p2 -20 -21]
K = place(A,B,P)
A1 = A - B*K;
eig(A1)
% eigenvalues of(A - LC) 
P2 = eig(A1)-[10;10;10;10];
L = place(A',C',P2);L = L'; 
L
A2 = A - L*C;
eig(A2)

yss = [0; 0.5]; % steady state response
kg = -C*inv(A-B*K)*B
kg = yss(2)/kg(2)

u(:,1) = -K*x_hat(:,1) + kg*r(:,1);
for i = 1:nt-1
x_dot(:,i) = A*x(:,i)+ B*u(:,i);
x(:,i+1) = x(:,i) + x_dot(:,i)*dt;
y(:,i+1) = C*x(:,i+1);
x_hat_dot(:,i) = A*x_hat(:,i)+ B*u(:,i)+ L*(y(:,i)- C*x_hat(:,i));
x_hat(:,i+1) = x_hat(:,i) + x_hat_dot(:,i)*dt;
y_hat(:,i+1) = C*x_hat(:,i+1);
u(:,i+1) = -K*x_hat(:,i+1) + kg*r(:,i+1);
end

figure
plot(t,x(1,:),'r',t,x_hat(1,:),'b--','linewidth',2)
set(gca,'fontsize',18)
legend({'$x_1$','${\hat x}_1$'},'Interpreter', 'latex')
legend boxoff
title('observer tracking control')
xlabel('Time (s)')
print(gcf,'fig9.png','-dpng','-r300');
figure
plot(t,x(2,:),'r',t,x_hat(2,:),'b--','linewidth',2)
set(gca,'fontsize',18)
legend({'$x_2$','${\hat x}_2$'},'Interpreter', 'latex')
legend boxoff
title('observer tracking control')
xlabel('Time (s)')
print(gcf,'fig10.png','-dpng','-r300');
figure
plot(t,x(3,:),'r',t,x_hat(3,:),'b--','linewidth',2)
set(gca,'fontsize',18)
legend({'$x_3$','${\hat x}_3$'},'Interpreter', 'latex')
legend boxoff
title('observer tracking control')
xlabel('Time (s)')
print(gcf,'fig11.png','-dpng','-r300');
figure
plot(t,x(4,:),'r',t,x_hat(4,:),'b--','linewidth',2)
set(gca,'fontsize',18)
legend({'$x_4$','${\hat x}_4$'},'Interpreter', 'latex')
legend boxoff
title('observer tracking control')
xlabel('Time (s)')
print(gcf,'fig12.png','-dpng','-r300');

