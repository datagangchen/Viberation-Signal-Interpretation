function [y]=RegressModel(Xt,Yt,gam,sig2,alpha,b,x)
%get the regression model
type = 'function estimation';
% [Yp,alpha,b,gam,sig2] = lssvm(Xt,Yt,type);
y= simlssvm({Xt,Yt,type,gam,sig2,'RBF_kernel','preprocess'},{alpha,b},x);
end