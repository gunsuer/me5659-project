clear all;
%% define systems
m1 = 2;
m2 = 1;
g = 10;
L = 0.75;

A = [0 1 0 0;g*(m1+m2)/(m1*L) 0 0 0;0 0 0 1;g*m2/m1 0 0 0];
B = [0;-1/(L*m1);0;1/m1];
C = [1 0 0 0; 0 0 1 0];
eig(A)
cm = ctrb(A,B);
rank(cm)
om = obsv(A,C);
rank(om)

%% initialization
x(:,1) = [0.1; 0; 0; 0]; 
y(:,1) = C*x(:,1);
X(:,1) = x(:,1);
Y(:,1) = C*X(:,1);
t = 0:0.01:10; %0.01 time span of interest
nt = length(t); % number of time steps
dt = t(2) - t(1);

% case 1: linear
for i = 1:nt-1
u(i) = 1;
x_dot(:,i) = A*x(:,i) + B*u(i);
x(:,i+1) = x(:,i) + x_dot(:,i)*dt;
y(:,i+1) = C*x(:,i+1);
end

% case 2: non-linear zero input
for i =1:nt-1
u(i) = 0;
X_dot(1,i) = X(2,i); 
X_dot(2,i) = (g*(m1+m2)/L)*sin(X(1,i))/(m1+m2*sin(X(1,i))^2) - 0.5*m2*sin(2*X(1,i))*X(2,i)^2/(m1+m2*sin(X(1,i))^2) + u(i)*(1/L)*cos(X(1,i))/(m1+m2*sin(X(1,i))^2);
X_dot(3,i) = X(4,i);   
X_dot(4,i) = -L*m2*sin(X(1,i))*X(2,i)^2/(m1+m2*sin(X(1,i))^2) + 0.5*g*m2*sin(2*X(1,i))/(m1+m2*sin(X(1,i))^2) + u(i)/(m1+m2*sin(X(1,i))^2);
X(:,i+1) = X(:,i) + X_dot(:,i)*dt;
Y(:,i+1) = C*X(:,i+1);
end

Y_zero = Y;

X(:,1) = x(:,1);
Y(:,1) = C*X(:,1);

% case 3: non-linear step input
for i =1:nt-1
u(i) = 1;
X_dot(1,i) = X(2,i); 
X_dot(2,i) = (g*(m1+m2)/L)*sin(X(1,i))/(m1+m2*sin(X(1,i))^2) - 0.5*m2*sin(2*X(1,i))*X(2,i)^2/(m1+m2*sin(X(1,i))^2) + u(i)*(1/L)*cos(X(1,i))/(m1+m2*sin(X(1,i))^2);
X_dot(3,i) = X(4,i);   
X_dot(4,i) = -L*m2*sin(X(1,i))*X(2,i)^2/(m1+m2*sin(X(1,i))^2) + 0.5*g*m2*sin(2*X(1,i))/(m1+m2*sin(X(1,i))^2) + u(i)/(m1+m2*sin(X(1,i))^2);
X(:,i+1) = X(:,i) + X_dot(:,i)*dt;
Y(:,i+1) = C*X(:,i+1);
end

Y_step = Y;

%% simulation results
figure
plot(t,Y(1,:),'b',t,y(1,:),'r','linewidth',2)
set(gca,'fontsize',18)
legend({'nonlinear','linear'},'Interpreter', 'latex')
title('pendulum angle')
legend boxoff
xlabel('Time (s)')
ylabel('Angle (rad)')
% print(gcf,'theta_open_loop.png','-dpng','-r300');

figure
plot(t,Y(2,:),'b',t,y(2,:),'r','linewidth',2)
set(gca,'fontsize',18)
legend({'nonlinear','linear'},'Interpreter', 'latex')
title('cart position')
legend boxoff
xlabel('Time (s)')
ylabel('Position (m)')
% print(gcf,'w_open_loop.png','-dpng','-r300');

figure
plot(t,Y_zero(2,:),'b',t,Y_step(2,:),'r','linewidth',2)
set(gca,'fontsize',18)
legend({'zero response','step response'},'Interpreter', 'latex')
title('cart position nonlinear')
legend boxoff
xlabel('Time (s)')
ylabel('Position (m)')
print(gcf,'nonlinear_zero_vs_step_w.png','-dpng','-r300');


