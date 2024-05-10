clear all;
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
x(:,1) = [0.1; 0.5; 10; -20]; 
y(:,1) = C*x(:,1);
t = 0:0.01:10; %0.01 time span of interest
nt = length(t); % number of time steps
dt = t(2) - t(1);
% reference input
for i = 1:nt
    r(:,i) = 1;
end

% case 1: open loop
for i = 1:nt-1
x_dot(:,i) = A*x(:,i);
x(:,i+1) = x(:,i) + x_dot(:,i)*dt;
y(:,i+1) = C*x(:,i+1);
end

% figure
% plot(t,y(1,:),'b',t,y(2,:),'r','linewidth',2)
% set(gca,'fontsize',18)
% legend({'$\theta$','$w$'},'Interpreter', 'latex')
% title('open loop')
% legend boxoff
% xlabel('Time (s)')

% case 2: state-feedback
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

yss = [0; 0.5]; % steady state response
kg = -C*inv(A-B*K)*B
kg = yss(2)/kg(2)

u(:,1) = -K*x(:,1) + kg*r(:,1);
for i = 1:nt-1
x_dot(:,i) = A*x(:,i) + B*u(:,i);
x(:,i+1) = x(:,i) + x_dot(:,i)*dt;
y(:,i+1) = C*x(:,i+1);
u(:,i+1) = -K*x(:,i+1) + kg*r(:,i+1);
end

figure
plot(t,r(1,1:nt),'k--',t,y(1,:),'b',t,y(2,:),'r','linewidth',2)
set(gca,'fontsize',18)
legend({'$r$','$\theta$','$w$'},'Interpreter', 'latex')
title('state-feedback')
legend boxoff
xlabel('Time (s)')
% print(gcf,'state_feedback.png','-dpng','-r300');


