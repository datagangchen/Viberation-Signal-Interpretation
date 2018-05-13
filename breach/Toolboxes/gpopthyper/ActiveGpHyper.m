% This file is part of GpOptimization.

% Copyright (c) by Gang Chen, 2017

function [Xt,Yt] = ActiveGpHyper( fhandle,Xt, Yt, Xs, mu,sigma, iter,varargin)
%%
% Sequential active sample using GP algorithms
%% Syntax
%   [queries, Yt] = ActiveGpHyper(fhandle, Xt, Yt, Xs, T)
%   [queries, Yt] = ActiveGpHyper(..., 'Name',Value)
%% Arguments
% * _fhandle_ function's handle which function returns observations with input
% variables
% * _Xt_ matrix _(nt, d)_ of initial data
% * _Yt_ vector _(nt, 1)_ of initial observations
% * _Xs_ matrix _(ns, d)_ of search data
%*_mu_estimated mean value 
%*_sigma_estimated convariance
% * _iter_ integer for the number of iterations
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
ip.addOptional('verbose', true);
ip.addOptional('errorhold',0.001);
ip.addOptional('Scale',1);
ip.parse(varargin{:});
opt = ip.Results;


if opt.B>1 & ~(strcmp(opt.algo,'gpucbpe') | strcmp(opt.algo,'gpbucb') | strcmp(opt.algo,'greedyucb'))
    error('The algorithm is not adapted for batch optimization');
end


[ns,d] = size(Xs);


    % compute ucb
    u = opt.u + log(iter^2*pi^2/6);

    switch opt.algo
      case 'gpucb'
        % GP-UCB - Srinivas et al. (2012)
        ucb = gpucb(sigma, u, ns);
        target = mu + ucb;
       [target_val,xt] = max(target);

      case 'gpacb'
        % GP-ACB -
        muadp=(mu-min(mu))/(max(mu)-min(mu));
        acb=muadp.*sigma;
        uacb = gpucb(acb, u, ns);
        target = mu + uacb;
        [target_val,xt] = max(target);  

    end
    
    % query point
    yt =-feval(fhandle,Xs(xt(end),:));
   if ~isempty(yt)
    Xt= [Xt; Xs(xt(end),:)]; % if batch, Xt is already partially updated
    Yt= [Yt; yt];
 
   else
       Xs(xt(end),:)=[];
   end

   
end