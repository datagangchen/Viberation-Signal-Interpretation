% demonstrate usage of regression
%
% Copyright (c) by Gang Chen 2-17-6-9.

function [mu, sigma] = gpRegression(Xt,Yt,Xs,iselect, Ncg,plotresult, opthyper)

% Sequential regression using GP algorithms
%% Syntax
%   [queries, Yt] = gpRegression(Xt, Yt, Xs, T)
%   [queries, Yt] = gpRegression(..., 'Name',Value)
%% Arguments

% * _Xt_ matrix _(nt, d)_ of initial data
% * _Yt_ vector _(nt, 1)_ of initial observations
% * _Xs_ matrix _(ns, d)_ of search data
% * _iselect_ integer for the index of optimization method
%*_Ncg_  interger for the number of conjugate gradient steps
%*_plot_ bool value for the whether plot the result
%% SAY WHICH CODE WE WISH TO EXERCISE
id = [1,1]; % use Gauss/Exact
id = [1,2; 3,2; 4,2]; % compare Laplace
id = [1,3; 2,3; 3,3]; % study EP
id = [1,5; 2,5]; % look into KL (takes quite a while)
id = [1,4; 2,4; 3,4; 4,4]; % deal with VB

seed = 197; randn('seed',seed), rand('seed',seed)

%ntr = size(Xt,1); nte = size(Xs,1);                        % number of training and test points

sn = 0.2;

cov = {@covSEiso}; sf = 1; ell = 0.4;                             % setup the GP
hyp0.cov  = log([ell;sf]);
mean = {@meanSum,{@meanLinear,@meanConst}};

a = zeros(size(Xt,2),1);
a(end)=1;
b = 0;       % m(x) = a*x+b
hyp0.mean = [a;b];
lik_list = {'likGauss','likLaplace','likSech2','likT'};   % possible likelihoods
inf_list = {'infGaussLik','infLaplace','infEP','infVB','infKL'};% inference algs



  lik = lik_list{id(iselect,1)};                                % setup the likelihood
  if strcmp(lik,'likT')
    nu = 4;
    hyp0.lik  = log([nu-1;sqrt((nu-2)/nu)*sn]);
  else
    hyp0.lik  = log(sn);
  end
  inf = inf_list{id(iselect,2)};
 % fprintf('OPT: %s/%s\n',lik_list{id(iselect,1)},inf_list{id(iselect,2)})

  if   opthyper % strcmp(opthyper,'true')

        hyp = minimize(hyp0,'gp', -Ncg, inf, mean, cov, lik, Xt, Yt); % opt hypers

  else
       hyp = hyp0;

   % hyp = minimize(hyp0,'gp', -Ncg, inf, mean, cov, lik, Xt, Yt); % opt hypers
  end
  [ymu, ys2] = gp(hyp, inf, mean, cov, lik, Xt, Yt, Xs);  % predict
  [nlZ(1)] = gp(hyp, inf, mean, cov, lik, Xt, Yt);
  mu=ymu;
  sigma=ys2;

  %% output the result
if plotresult
        figure;
  col = {'k',[.8,0,0],[0,.5,0],'b',[0,.75,.75],[.7,0,.5]};                % colors  
sdscale = 0.5;                  % how many sd wide should the error bars become?
      [Z, IZ] = sort(Xs(:,1));
     plot(Z,mu(IZ),'Color',col{iselect},'LineWidth',2)
    leg{1} = sprintf('%s/%s -lZ=%1.2f',...
                                lik_list{id(iselect,1)},inf_list{id(iselect,2)},nlZ);
    hold on
      ysd = sdscale*sqrt(ys2);
  fill([Z;flipud(Z)],[ymu(IZ)+ysd(IZ);flipud(ymu(IZ)-ysd(IZ))],...
       col{iselect},'EdgeColor',col{iselect},'FaceAlpha',0.1,'EdgeAlpha',0.3);
   hold on
    plot(Z,ymu(IZ) ,'Color',col{iselect},'LineWidth',2)
    hold on
    plot(Xt,Yt,'k+'),
    plot(Xt,Yt,'ko'), legend(leg)
end

%plot(xtr,ytr,'k+'), plot(xtr,ytr,'ko'), legend(leg)