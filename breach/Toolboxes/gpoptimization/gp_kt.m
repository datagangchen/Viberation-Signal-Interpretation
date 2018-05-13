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

function [mu, sigma2,kt] = gp_kt(Kts, dKss, BayesInv, Ht, Hs,Xs)
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

if nargin<5; Hs = []; end
if nargin<4; Ht = []; end

% mu
if Ht
    mu = Hs*BayesInv.bet + Kts'*BayesInv.invCY; % (ns x 1)
else
    mu = Kts'*BayesInv.invCY;
end
% [Z, IZ] = sort(Xs(:,1));
% X=linspace(min(Xs(:,1)),max(Xs(:,1)),length(Z));
% Yinter=interp1(Z,mu(IZ),X,'spline');

% mul=[flip(Yinter(1:5)) Yinter flip(Yinter(end-5:end))];

% Ykt=zeros(size(Yinter));
% for iter=1:length(Yinter)
%     Ykt(iter)=kurtosis(mul(iter:iter+10));
% end
% Ykt=diff(Yinter);
% InterKt=interp1(X(1:end-1),Ykt,Z,'spline');
% for iter=1:length(mu)
%    mukt(IZ(iter))=abs(InterKt(iter));
% end
% figure(3)
% plot(Z,mukt(IZ));
% sigma2
if nargout > 1
    Vf = BayesInv.RC'\Kts; % (nt x ns)
    covf = dKss - sum(Vf.*Vf, 1)';

    if Ht
        Rs = Hs' - Ht'*solve_chol(BayesInv.RC, Kts); % (b x ns)
        LHCH = cholpsd(Ht'*(solve_chol(BayesInv.RC, Ht))); % (b x b)
        Vb = LHCH'\Rs; % (b x ns)
        covb = sum(Vb.*Vb, 1)';
        sigma2 = covb + covf;
    else
        sigma2 = covf;
    end
    nmukt=(mu-min(mu))/(max(mu)-min(mu));

%     kt=sigma2;
%     kt=(sigma2+1-nmukt')/2;
    
    kt=nmukt.*sigma2;
%     kt=mukt/max(mukt)+sigma2;
end

end


