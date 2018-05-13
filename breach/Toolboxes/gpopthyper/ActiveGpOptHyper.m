% This file is part of GpOptimization.

% Copyright (c) by Gang Chen, 2017

function [Xtnew,Ytnew, Sigma] = ActiveGpOptHyper( fhandle,Xt, Yt, Xs, T,varargin)
%%
% Sequential active sample using GP algorithms
%% Syntax
%   [queries, Yt] = ActiveGpOptHyper(fhandle, Xt, Yt, Xs, T)
%   [queries, Yt] = ActiveGpOptHyper(..., 'Name',Value)
%% Arguments
% * _fhandle_ function's handle which function returns observations with input
% variables
% * _Xt_ matrix _(nt, d)_ of initial data
% * _Yt_ vector _(nt, 1)_ of initial observations
% * _Xs_ matrix _(ns, d)_ of search data
% * _T_ integer for the number of iterations
%% Name-Value Pair Arguments
% * _algo_ string code for the query algorithm
%    possible values are 'gpucb' (default), 'chaining', 'ei' for purely sequential optimization
%    and 'gpucbpe', 'greedyucb' for batch sequential optimization
% * _noise_ scalar for the noise standard deviation (default: 1e-2)
% * _u_ scalar for the negative logarithm of the upper bound probability (default: 3)
% * _B_ integer for the size of the batch in batch sequential optimization (default: 1)
% * _plot_ boolean to visualize the optimization (default: false)
% * _verbose_ boolean to monitor the optimization (default: true)
%% Outputs
% * _queries_ vector _(1,T)_ of queries indices
% * _Yt_ vector _(T,1)_ of observations
%% See also
% <gpucb.html gpucb> | <chaining_ucb.html chaining_ucb>

ip = inputParser;
ip.addOptional('algo', 'gpucb');
ip.addOptional('u', 3);
ip.addOptional('B', 1);
ip.addOptional('plot', false);
ip.addOptional('verbose', false);
ip.addOptional('errorhold',0.001);
ip.addOptional('iselect',1);
ip.addOptional('Scale',1);
ip.addOptional('opthyper',false);
ip.parse(varargin{:});
opt = ip.Results;


if opt.B>1 & ~(strcmp(opt.algo,'gpucbpe') | strcmp(opt.algo,'gpbucb') | strcmp(opt.algo,'greedyucb'))
    error('The algorithm is not adapted for batch optimization');
end
fprintf('\n Sampling the system\n');
for iter=1:T
[mu, sigma] = gpRegression(Xt,Yt,Xs,1, 10,opt.plot,opt.opthyper);
[Xt,Yt] = ActiveGpHyper( fhandle,Xt, Yt, Xs, mu,sigma, iter,'algo', opt.algo,'Scale', opt.Scale);
end
tempzeros=ones(length(Xs(:,1)),1)*Xt(:,1)'-Xs(:,1)*ones(1,length(Xt(:,1)));
tempzeros=ceil(abs(tempzeros)./max(max(abs(tempzeros))))';
tempzeros=floor(sum(tempzeros)/length(Xt(:,1)));


Xtnew=Xt;
Ytnew=Yt;

Sigma=zeros(size(Xt,1),1);
for ind=1:size(Xt,1)
   index=find(Xs(:,2)==Xt(ind,2));
   Sigma(ind)=sigma(index(end));
    
end



 if opt.plot
        figure;
        iselect=opt.iselect;
        id = [1,4; 2,4; 3,4; 4,4]; % deal with VB
        lik_list = {'likGauss','likLaplace','likSech2','likT'};   % possible likelihoods
inf_list = {'infGaussLik','infLaplace','infEP','infVB','infKL'};% inference algs
  col = {'k',[.8,0,0],[0,.5,0],'b',[0,.75,.75],[.7,0,.5]};                % colors  
sdscale = 3;                  % how many sd wide should the error bars become?

      [Z, IZ] = sort(Xs(:,1));
     plot(Z,mu(IZ),'Color',col{iselect},'LineWidth',2)
    leg{1} = sprintf('%s/%s -lZ=%1.2f',...
                                lik_list{id(iselect,1)},inf_list{id(iselect,2)},1);
    hold on
      ysd = sdscale*sqrt(sigma);
  fill([Z;flipud(Z)],[mu(IZ)+ysd(IZ);flipud(mu(IZ)-ysd(IZ))],...
       col{iselect},'EdgeColor',col{iselect},'FaceAlpha',0.1,'EdgeAlpha',0.3);
   hold on
    plot(Z,mu(IZ) ,'Color',col{iselect},'LineWidth',2)
    hold on
    plot(Xt,Yt,'k+'),
    plot(Xt,Yt,'ko'), legend(leg)
 end

end

