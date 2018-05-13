function [ alpha,b] = ParameterSynthesis(X,Y)
disp(' >> gam = 10;');
gam = 50;
disp(' >> sig2 = 0.2;');
sig2 = 0.2;
type = 'classification';
disp(' ');
disp(' >> [alpha,b] = trainlssvm({X,Y,type,gam,sig2,''RBF_kernel'',''preprocess''});');
[alpha,b] = trainlssvm({X,Y,type,gam,sig2,'RBF_kernel','preprocess'});

Xt = 10.*rand(50,2)-5;
Ytest = simlssvm({X,Y,type,gam,sig2,'RBF_kernel','preprocess'},{alpha,b},Xt);
disp(' ');
disp(' >> plotlssvm({X,Y,type,gam,sig2,''RBF_kernel'',''preprocess''},{alpha,b});');
figure; plotlssvm({X,Y,type,gam,sig2,'RBF_kernel','preprocess'},{alpha,b});
end

