function [Y ] = TestFunction(X)
%this file is a test function general a output value with a input vector
%% Syntax
%  Y=TestFunction(X)  y= = x e^{-x^2 -y^2}
%% Arguments
% * _X_ matrix _(1, 2)_ where _n_ is the number of data points and _d_ is the dimension
%% Output
% * _Y_ vector _(1, 1)_ of output value

%Y=X(:,1).*exp(-X(:,1).^2-.5*X(:,2).^2);

Y=X.^2-6*X+10.5;
end

