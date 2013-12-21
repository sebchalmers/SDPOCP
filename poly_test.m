clear all
close all
clc

syms x y z real

maxDeg = 10;


v1 = x;
v2 = y;
for k = 2:maxDeg
    v1 = [v1;x^k];
    v2 = [v2;y^k];
end
v = [1;v1;v2];
poly = v.'*randn(1+2*maxDeg,1+2*maxDeg)*v;

% poly = 1;
% for k = 1:maxDeg
%     poly = poly + k*x^k;
% end

tic
[newPoly,newVar] = lift(poly,[x,y]);
toc
length(newVar)
newPoly

jacobian(jacobian(newPoly,newVar(2:end)),newVar(2:end))