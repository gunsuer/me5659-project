clear;clc

syms p11 p12 p13 p14 p22 p23 p24 p33 p34 p44;
m1 = 2;
m2 = 1;
g = 10;
L = 0.75;

A = [0 1 0 0;g*(m1+m2)/(m1*L) 0 0 0;0 0 0 1;g*m2/m1 0 0 0];
P = [p11 p12 p13 p14;p12 p22 p23 p24;p13 p23 p33 p34;p14 p24 p34 p44]

Q = simplify(transpose(A)*P+P*A)

solve(Q == eye(4),[p11 p12 p13 p14 p22 p23 p24 p33 p34 p44])

eig(A)
