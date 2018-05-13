function [Y ] = demoFunctionScale(X,scale)
%this file is a test function general a output value with a input vector
%% Syntax
%  Y=TestFunction(X)  y= = x e^{-x^2 -y^2}
%% Arguments
% * _X_ matrix _(1, 2)_ where _n_ is the number of data points and _d_ is the dimension
%% Output
% * _Y_ vector _(1, 1)_ of output value

%Y=X(:,1).*exp(-X(:,1).^2-.5*X(:,2).^2);

x=X(:,1);
y=X(:,2);
%Y=-20*exp(-0.2*sqrt(0.5*(x.^2+y.^2)))-exp(0.5*(cos(2*pi*x)+cos(2*pi*y)))+exp(1)+20;
%Y=-(y+47)*sin(sqrt(abs(x/2+y+47)))-x*sin(sqrt(abs(x-y-47)))+10*rand(1);
Y=20+x.^2-10*cos(2*pi*x)+y.^2-10*cos(2*pi*y)+10*rand(1);
Y=scale*Y;
%Y=X.^2-6*X+10.5;
end