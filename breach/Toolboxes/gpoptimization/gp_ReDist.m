% This file is part of GpOptimization.
%
% GpOptimization is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, version 3 of the License.
%
% GpOptimization is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with GpOptimization. If not, see <http://www.gnu.org/licenses/>.
%
% Copyright (c) by Emile Contal, 2015

function [rub,rlb] = gp_ReDist(fhandle,lb,ub,varargin)
%%
% Posterior mean and variance of GP given the kernel matrices and the Bayesian inferance
%% Syntax
%   [mu, sigma2] = gp_pred(Kts, dKss, BayesInv)
%   [mu, sigma2] = gp_pred(Kts, dKss, BayesInv, Ht, Hs)
%% Arguments
% * _Kts_ matrix _(nt, ns)_ of kernel between the points of _Xt_ and _Xs_
% * _dKss_ matrix _(ns, 1)_ of diagonal kernel between the points of _Xs_
% * _BayesInv_ structure array returned by _<gp_inf.html gp_inf>(Ht, Ktt, Yt, noise)_
% * _Ht_ matrix _(nt, b)_ of basis data for _Xt_
% * _Hs_ matrix _(ns, b)_ of basis data for _Xs_
%% Outputs
% * _mu_ matrix _(ns, 1)_ of posterior mean $E[f(X_s) \mid X_t, Y_t]$
% * _sigma2_ matrix _(ns, 1)_ of posterior variance $V[f(X_s) \mid X_t, Y_t]$
%% See also
% <gp_inf.html gp_inf> | <gp_dist.html gp_dist>



ip = inputParser;
ip.addOptional('d', size(lb,1));
ip.addOptional('ns', 1000);
ip.addOptional('nt', 20);
ip.addOptional('kernel', @kernel_se_normiso);
ip.addOptional('basis', @basis_none);
ip.addOptional('noise', 1e-2);
ip.addOptional('plot', true);
ip.addOptional('posterior', true);
ip.addOptional('verbose', true);
ip.addOptional('scale', 0.064);
ip.parse(varargin{:});
opt = ip.Results;

if opt.verbose; fprintf('generating fn...\n'); end
range=ub-lb;
temp=zeros(1,size(lb,1));
while 1
    
Xt=rand(opt.nt, opt.d);
% randsample=opt.ns*rand(1,opt.nt);
Trans_Xt=diag(range)*Xt'+lb*ones(1,opt.nt);
Xt=Trans_Xt';
% f = @(X) Fs(X) + opt.noise * randn(size(X,1), 1);
Yt=zeros(length(Xt),1);
for index=1:length(Xt)
Yt(index,1) =-feval(fhandle,Xt(index,:));
end
[maxF,IF]=max(Yt);
max_x=Xt(IF,:);
tempx=max_x./range';
  distemp=(tempx-temp)*(tempx-temp)';
  if distemp<0.2*opt.scale
      break;
  end
 temp=tempx;
end

nXt=Xt*diag((1./max(Xt)));
n_Xt=size(Xt,1);
X_temp=[];
for index=1:n_Xt
    temp=nXt(IF,:)-nXt(index,:);
    Dis=temp*temp';
    if Dis<=opt.scale 
    X_temp=[X_temp;Xt(index,:)];
    end
end


%calculate the bonds
   rub=max(X_temp)';
   rlb=min(X_temp)';
   

end

