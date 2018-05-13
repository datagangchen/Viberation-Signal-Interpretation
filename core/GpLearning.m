function [Xt,Yt] = GpLearning(obj,Xs,Xt,Yt, varargin)
% Gp learning for the parameters

ip = inputParser;
ip.addOptional('algo', 'gpucb');
ip.addOptional('u', 3);
ip.addOptional('B', 1);
ip.addOptional('kernel', @kernel_se_normiso);
ip.addOptional('basis', @basis_none);
ip.addOptional('noise', 1e-2);
ip.addOptional('plot', true);
ip.addOptional('verbose', true);
ip.addOptional('Kss', []);
ip.addOptional('errorhold',0.001);
ip.addOptional('class','robust')
ip.parse(varargin{:});
opt = ip.Results;
% tempmax=0;
% check input


opt.Kss=opt.kernel(Xs,Xs);
if strcmp(opt.algo, 'chaining') & isempty(opt.Kss)
    error('You must provide Kss in the optional arguments for using the chaining algorithm');
end
if opt.B>1 & ~(strcmp(opt.algo,'gpucbpe') | strcmp(opt.algo,'gpbucb') | strcmp(opt.algo,'greedyucb'))
    error('The algorithm is not adapted for batch optimization');
end
target_plot = [];

[ns,d] = size(Xs);

queries = [];
Ht = opt.basis(Xt);
Hs = opt.basis(Xs);
Ktt = opt.kernel(Xt,Xt);
Kts = opt.kernel(Xt,Xs);
dKss = opt.kernel(Xs,'diag');

BayesInv = gp_inf(Ktt, Yt, opt.noise, Ht);

for iter=1:obj.Episode
    % GP pred
    [mu, s2] = gp_pred(Kts, dKss, BayesInv, Ht, Hs);

    % compute ucb
    u = opt.u + log(iter^2*pi^2/6);

    switch opt.algo
      case 'gpucb'
        % GP-UCB - Srinivas et al. (2012)
        ucb = gpucb(s2, u, ns);
        target = mu + ucb;
        [target_val,xt] = max(target);
      case 'gpkt'
        % GP-UCB - Srinivas et al. (2012)
        [mu, s2,kt] = gp_kt(Kts, dKss, BayesInv, Ht, Hs,Xs);
        ktb = gpucb(kt, u, ns);
        target = mu + ktb;
        [target_val,xt] = max(target);  
      case 'chaining'
        % Chainig-UCB - Contal et al. (2015)
        D2 = gp_dist(opt.Kss, Kts, Kts, dKss, dKss, BayesInv);
        ucb = chaining_ucb(D2, s2, opt.u + log(iter^2*pi^2/6));
        target = mu + ucb;
        [target_val,xt] = max(target);
      case 'ei'
        % Expected Improvement
        fmax = max(Yt);
        ni = (mu-fmax) ./ sqrt(s2);
        target = (mu-fmax).*normcdf(ni) + sqrt(s2) .* normpdf(ni);
        [target_val,xt] = max(target);
      case 'greedyucb'
        % Batch-greedy ucb policy, similar to GP-BUCB without the initialization phase
        xt = [];
        for b=1:opt.B
            ucb = gpucb(s2, u, ns);
            target = mu + ucb;
            [target_val,xtb] = max(target);
            xt = [xt; xtb];
            if b<opt.B
                Ht = [Ht; opt.basis(Xs(xtb,:))];
                Ktt12 = opt.kernel(Xt,Xs(xtb,:));
                Ktt22 = opt.kernel(Xs(xtb,:),Xs(xtb,:));
                BayesInv = gp_inf_update(Ktt12, Ktt22, [Yt; mu(xt)], opt.noise, BayesInv, Ht);
                Kts = [Kts; opt.kernel(Xs(xtb,:),Xs)];
                [~, s2] = gp_pred(Kts, dKss, BayesInv, Ht, Hs);
                Xt = [Xt; Xs(xtb,:)];
            end
        end
      case 'gpucbpe'
        % GP-UCB-PE - Contal et al. (2013)
        ucb = gpucb(s2, u, ns);
        target = mu + ucb;
        [target_val,xtb] = max(target);
        xt = [xtb];
        lcb = mu - ucb;
        ydot = max(lcb);
        is_in_Rt = target >= ydot;
        for b=2:opt.B
            % update
            Ht = [Ht; opt.basis(Xs(xtb,:))];
            Ktt12 = opt.kernel(Xt,Xs(xtb,:));
            Ktt22 = opt.kernel(Xs(xtb,:),Xs(xtb,:));
            BayesInv = gp_inf_update(Ktt12, Ktt22, [Yt; mu(xt)], opt.noise, BayesInv, Ht);
            Kts = [Kts; opt.kernel(Xs(xtb,:),Xs)];
            [~, s2] = gp_pred(Kts, dKss, BayesInv, Ht, Hs);
            Xt = [Xt; Xs(xtb,:)];
            % next query
            s2_in_Rt = s2 .* is_in_Rt;
            [~, xtb] = max(s2_in_Rt);
            xt = [xt; xtb];
        end
    end
    
    % query point
    queries = [queries xt];
    
    %% assign parameter
        Xtt=Xs(xt(end),:);
        scale= zeros(size(obj.scale_param));
        time = zeros(size(obj.time_param));
        [ii,jj]=size(scale);

        parfor index_i =1:ii
        scale(index_i,:)=Xtt(jj*(index_i-1)+1:jj*index_i);
        time(index_i,:) =Xtt(jj*ii+jj*(index_i-1)+1: jj*ii +jj*index_i);  
        end
        obj.Nimf = Xtt(2*ii*jj+1:2*ii*jj+length(obj.Nimf));
        obj.Nfreq= Xtt(2*ii*jj+length(obj.Nimf)+1:end);
        obj.time_param=time;
        obj.scale_param = scale; 
    
    switch opt.class
        case 'robust'
            yt=obj.get_class();
            
        case 'class'
            yt=obj.get_robust();
    end

 

    Xt = [Xt; Xs(xt(end),:)]; % if batch, Xt is already partially updated
    Yt = [Yt; yt];

    

    % GP update
    Ht = [Ht; opt.basis(Xt(end,:))];
    Ktt12 = opt.kernel(Xt(1:end-1,:),Xt(end,:));
    Ktt22 = opt.kernel(Xt(end,:),Xt(end,:));
    BayesInv = gp_inf_update(Ktt12, Ktt22, Yt, opt.noise, BayesInv, Ht);
    Kts = [Kts; opt.kernel(Xt(end,:),Xs)];
end

end

