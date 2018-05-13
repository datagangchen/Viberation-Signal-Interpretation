% This file is part of GpOptimization.

% Copyright (c) by Gang Chen, 2015

function [Xs, Xt, Yt] = gp_sample(fhandle,lb,ub,varargin)
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


ip = inputParser;
ip.addOptional('d', size(lb,1));
ip.addOptional('ns', 1500);
ip.addOptional('nt', 40);
ip.addOptional('kernel', @kernel_se_normiso);
ip.addOptional('basis', @basis_none);
ip.addOptional('noise', 1e-2);
ip.addOptional('plot', true);
ip.addOptional('posterior', true);
ip.addOptional('verbose', true);
ip.parse(varargin{:});
opt = ip.Results;

if opt.verbose; fprintf('generating fn...\n'); end
Xs = rand(opt.ns, opt.d);
% linsample=linspace(0,1,opt.ns);
% Xs=linsample'*ones(1,opt.d);
range=ub-lb;
Trans_Xs=diag(range)*Xs'+lb*ones(1,opt.ns);
Xs=Trans_Xs';
Xt=rand(opt.nt, opt.d);
% randsample=opt.ns*rand(1,opt.nt);
Trans_Xt=diag(range)*Xt'+lb*ones(1,opt.nt);
Xt=Trans_Xt';

% Xt=Xs(round(randsample),:);
Kss = opt.kernel(Xs,Xs);

% Fs = cholpsd(Kss)' * randn(opt.ns, 1);
% f = @(X) Fs(X) + opt.noise * randn(size(X,1), 1);
% f=TestFunction(Xs);
Yt=zeros(length(Xt),1);
for index=1:length(Xt)
Yt(index,1) =-feval(fhandle,Xt(index,:));
end



if opt.posterior
%     plot(Xt, Yt, 'r+');
%     hold on

    if opt.verbose; fprintf('Bayesian inferance...\n'); end
    Ktt = Kss(1:opt.nt,1:opt.nt);
    Ht = opt.basis(Xt);
    BayesInv = gp_inf(Ktt, Yt, opt.noise, Ht);
    
    if opt.verbose; fprintf('Prediction...\n'); end
    Kts = Kss(1:opt.nt,:);
    Hs = opt.basis(Xs);
    dKss = diag(Kss);
    [mu, s2] = gp_pred(Kts, dKss, BayesInv, Ht, Hs);

    if opt.plot
        [Z,IZ] = sort(Xs(:,1));
        plot(Z, mu(IZ), 'k');
        fill([Z; flipdim(Z,1)], [mu(IZ)+2*s2(IZ); flipdim(mu(IZ)-2*s2(IZ),1)], ...
             [7 8 7]/8, 'EdgeColor','None', 'FaceAlpha',.8);
    end
else
    if opt.plot
        if opt.d == 1
            [Z,IZ] = sort(Xs(:,1));
            plot(Z, Fs(IZ), 'k');
        elseif opt.d == 2
            scatter(Xs(:,1),Xs(:,2),30,Fs,'fill');
        end
    end
end

end