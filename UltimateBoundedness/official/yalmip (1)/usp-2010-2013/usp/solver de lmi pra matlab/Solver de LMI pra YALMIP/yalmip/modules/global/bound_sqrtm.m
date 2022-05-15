function [Ax, Ay, b] = bound_sqrtm(xL,xU)

% Two lower bounds from tangents
% y < sqrtm(xL) + (x-xL)*1/sqrt(xL)
% y < sqrtm(xU) + (x-xU)*1/sqrt(xU)

% Upper bound from conneting extreme points
% y > sqrt(xU)(x-xL)/(xU-xL) +  sqrt(xL)(xU-x)/(xU-xL)

% can be wrtitten as
% Ax*x + Ay*y < b

xL = max([1e-6 xL]);
Ay = [1;
      1;
      %1;
      -1];
xM = (xL + xU)/2;
b = [sqrt(xL) - xL/(2*sqrt(xL));
     sqrt(xU) - xU/(2*sqrt(xU));
     %sqrt(xM) - xM/(2*sqrt(xM)); 
     sqrt(xU)*(xL)/(xU-xL) -  sqrt(xL)*(xU)/(xU-xL)];
 
Ax = [-1/sqrt(4*xL);
      -1/sqrt(4*xU);
      %-1/sqrt(4*xM);
       sqrt(xU)/(xU-xL) - sqrt(xL)/(xU-xL)];
