clear;clc

% syms g L m1 m2;

% params
m1 = 2;
m2 = 1;
g = 10;
L = 0.75;

% ss model
A = [0 1 0 0;g/L*(m1+m2)/m1 0 0 0;0 0 0 1;g*m2/m1 0 0 0]
B = [0;1/(L*m1);0;1/m1]
C = [1 0 0 0]
D = 0

[b,a] = ss2tf(A,B,C,D)

%% ccf
cntr = [B A*B A^2*B A^3*B]
icntr = inv(cntr)
q = icntr(4,:)

T = [q;q*A;q*A^2;q*A^3]
% T = cntr
% Ac = simplify(T*A*inv(T))
% Bc = simplify(T*B)
% Cc = simplify(C*inv(T))

Ac = T*A*inv(T)
Bc = T*B
Cc = C*inv(T)

[b,a] = ss2tf(Ac,Bc,Cc,D)