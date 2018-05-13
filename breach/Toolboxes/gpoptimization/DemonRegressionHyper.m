function DemonRegressionHyper()
%This function demonstrate the optimization of Hyper parameters
% Programmed by Gang Chen
%Date 06/08/2017 (MM/DD/YY)
fhandle=@TestFunction;
lb=[1];
ub=[5];
[ Xs, Xt, Yt] = gp_sample_Hyper(fhandle,lb,ub,'plot', false);

[ sigma ] =OptHyper(Xt, Yt);
sigma
any2vec(sigma)
end

