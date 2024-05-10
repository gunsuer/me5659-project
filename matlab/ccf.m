clear;clc

m1 = 2;
m2 = 1;
g = 9.8;
L = 0.75;

A = [0 1 0 0;g/L*(m1+m2)/m1 0 0 0;0 0 0 1;g*m2/m1 0 0 0];
B = [0;1/(L*m1);0;1/m1];
C = [1 0 0 0];
D = 0;

[E,V] = eig(A)

cntr = [B A*B A^2*B A^3*B]
icntr = inv(cntr)
q = icntr(4,:)
T = [q;q*A;q*A^2;q*A^3]
Ac = T*A*inv(T)
Bc = T*B
Cc = C*inv(T)