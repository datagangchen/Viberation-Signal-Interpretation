% This file is part of GpOptimization.

% Copyright (c) by Gang Chen, 2017

function [Xs, Xt, Yt] = gp_sample_hyper(fhandle,lb,ub,varargin)
%%
% Sample a Gaussian process in a hypercube and perform Bayesian inferance
%% Syntax
%   [f, Xs, Fs, Xt, Yt] = gp_sample(..., 'Name',Value)
%% Input
% -handle : function handle 
% -boundery:_(2,d)_maxtrix_ boundery of the variables 
%% Name-Value Pair Arguments
% * _d_ integer for the dimension (default: 1)
% * _size_ scalar for the length of the hypercube (default: 40)
% * _ns_ integer for the number of sampled points (default: 1000)
% * _nt_ integer for the number of training points (default: 10)
% * _kernel_ kernel function (default: <kernel_se_normiso.html @kernel_se_normiso>)
% * _basis_ basis function (default: <basis_none.html @basis_none>)
% * _noise_ scalar for the noise standard deviation (default: 0.01)
% * _plot_ boolean to visualize the optimization (default: true)
% * _posterior_ boolean to compute the Bayesian inferance (default: true)
% * _verbose_ boolean to monitor the process (default: true)
%% Outputs
% * _f_ function which returns noisy observations
% * _Xs_ matrix _(ns, d)_ of sampled data
% * _Fs_ matrix _(ns, 1)_ of sampled values
% * _Xt_ matrix _(nt, d)_ of training data
% * _Yt_ vector _(nt, 1)_ of training observations
Xt=[];
Xs=[];
Yt=[];

ip = inputParser;
ip.addOptional('d', size(lb,1));
ip.addOptional('ns', 1000);
ip.addOptional('nt', 5);
ip.addOptional('kernel', @kernel_se_normiso);
ip.addOptional('basis', @basis_none);
ip.addOptional('noise', 1e-2);
ip.addOptional('plot', true);
ip.addOptional('posterior', true);
ip.addOptional('verbose', true);
ip.addOptional('Scale',1);
ip.parse(varargin{:});
opt = ip.Results;

if opt.verbose; fprintf('generating fn...\n'); end
Xs = rand(opt.ns, opt.d);
Xs = rand(opt.ns, opt.d);
% linsample=linspace(0,1,opt.ns);
% Xs=linsample'*ones(1,opt.d);
range=ub-lb;
Trans_Xs=diag(range)*Xs'+lb*ones(1,opt.ns);
Xs=Trans_Xs';
% Xt=rand(opt.nt, opt.d);
% % randsample=opt.ns*rand(1,opt.nt);
% Trans_Xt=diag(range)*Xt'+lb*ones(1,opt.nt);
Xt=Xs(1:opt.nt,:);

Yt=zeros(size(Xt,1),1);
for index=1:size(Xt,1)
Yt(index,1) =-feval(fhandle,Xt(index,:));
end

end