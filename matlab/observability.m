clear;clc

syms g L m1 m2;

% ss model
A = [0 1 0 0;g/L*(m1+m2)/m1 0 0 0;0 0 0 1;g*m2/m1 0 0 0];
B = [0;1/(L*m1);0;1/m1];

C_theta = [1 0 0 0];
C_w = [0 0 1 0];
C_theta_and_w = [1 0 0 0;0 0 1 0];

D = 0;

%% theta
obsv_theta = [C_theta;C_theta*A;C_theta*A^2;C_theta*A^3]
rank(obsv_theta)

inv_obsv_theta = inv(obsv_theta)
q = obsv_theta(4,:)

T = [q;q*A;q*A^2;q*A^3]
% T = cntr
Ac = simplify(inv(T)*A*T)
Bc = simplify(T*B)
Cc = simplify(C_theta*inv(T))

%% w
obsv_w = [C_w;C_w*A;C_w*A^2;C_w*A^3]
rank(obsv_w)

%% theta and w
obsv_theta_and_w = [C_theta_and_w;C_theta_and_w*A;C_theta_and_w*A^2;C_theta_and_w*A^3]
rank(obsv_theta_and_w)

