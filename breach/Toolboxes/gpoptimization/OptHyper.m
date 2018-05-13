function [ sigma ] =OptHyper(X , Y )
%This function optimizae the hyperparameters
%By Gang Chen
%date 06/08/2017
% (X,Y) the training data pairs
%Xte the text point

%% SAY WHICH CODE WE WISH TO EXERCISE
id = [1,1]; % use Gauss/Exact
id = [1,2; 3,2; 4,2]; % compare Laplace
id = [1,3; 2,3; 3,3]; % study EP
id = [1,5; 2,5]; % look into KL (takes quite a while)
id = [1,4; 2,4; 3,4; 4,4]; % deal with VB
mu = 1.0; s2 = 0.01^2; nu = 3;
pg = {@priorGauss,mu,s2};                          % Gaussian prior
pc = {@priorClamped};                         % equivalent to above
p1 = {@priorSmoothBox1,0,3,15};  % smooth box constraints lin decay
sn = 1;
cov = {@covSEiso}; sf = 1; ell = 0.4;                             % setup the GP
hyp0.cov  = log([ell;sf]);
mean = {@meanSum,{@meanLinear,@meanConst}}; a = 1/5; b = 1;       % m(x) = a*x+b
hyp0.mean = [a;b];

lik_list = {'likGauss','likLaplace','likSech2','likT'};   % possible likelihoods
inf_list = {'infGaussLik','infLaplace','infEP','infVB','infKL'};% inference algs

Ncg = 20;                                   % number of conjugate gradient steps

idSelect=3;
  lik = lik_list{id(idSelect,1)};                                % setup the likelihood
  if strcmp(lik,'likT')
    nu = 4;
    hyp0.lik  = log([nu-1;sqrt((nu-2)/nu)*sn]);
  else
    hyp0.lik  = log(sn);
  end
  inf = inf_list{id(idSelect,2)};
  fprintf('OPT: %s/%s\n',lik_list{id(idSelect,1)},inf_list{id(idSelect,2)})
 % prior.mean = {pg;pc};  % Gaussian prior for first, clamp second par
%prior.cov  = {p1;[]}; % box prior for first, nothing for second par
%im = {@infPrior,@infExact,prior};                % inference method
hyp = minimize(hyp0,'gp', -Ncg, inf, mean, cov, lik, X, Y); % opt hypers


  sigma=hyp;
  
  
end

