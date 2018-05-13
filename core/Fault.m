classdef Fault
     properties(SetAccess = private, Hidden)
         imfbin      % number of imfs used
         Fs          % sampling frequency
         data        %vibration data
         label       % label of the system
         index_imf        % the index of chosen imf
         energy_param % param for predicate
         time_param  % param for time
         freq_param
         alph        %regulator
         robustness 
         misclass  
         history    % the history of learning
         ins_freq_all
         ins_energy_all
         Episode
         bound     % bound for frequency calcaculation
         ub
         lb
         parl

     end
  %% the formula is with template 
  %% Alw_Ev_[0,t](freq =pi_1) and (energy >pi_2)
    
    methods   
         function obj = initialize(obj,arg1,arg2,arg3,arg4,varargin)
             obj.data = arg1;
             obj.label= arg2;
             obj.Fs = arg3;
             obj.Episode=arg4;
             if length(varargin)==1
                 obj.imfbin =varargin{1};
                  obj.alph = 0.1;
                 obj.bound =0.05;
                 obj.parl = false;
             elseif length(varargin)==2
                 obj.imfbin =varargin{1};
                 obj.alph = varargin{2};                
                 obj.bound =0.05;
                 obj.parl = false;     
             elseif length(varargin)==3
                 obj.imfbin =varargin{1};
                 obj.alph = varargin{2};
                 obj.bound =varargin{3}; 
                 obj.parl=false;
             elseif length(varargin)==4                    
                 obj.imfbin =varargin{1};
                 obj.alph = varargin{2};
                 obj.bound =varargin{3};
                 obj.parl =varargin{4};
             else 
                 obj.imfbin = 5;
                 obj.alph = 0.1;
                 obj.bound =0.05;
                 obj.parl= false;
             end
             
             LA= size(obj.data,1);
             signal=obj.data;
             Fss=obj.Fs;
             Imfbin=obj.imfbin;
             Ins_freq_all=[];
             Ins_energy_all=[];
             parfor index =1:LA
                 [imf,residual,info] = emd(signal(index,:),'Interpolation','pchip','Display',false);
                 [hs,f,t,imfinsf,imfinse] = hht(imf,Fss);
                 Ins_freq_all(:,:,index)=imfinsf(:,1:Imfbin)';
                 Ins_energy_all(:,:,index)=imfinse(:,1:Imfbin)';
             end
%              Ins_freq_all=Ins_freq_all./max(max(max(Ins_freq_all)));
%              Ins_energy_all=Ins_energy_all./max(max(max(Ins_energy_all)));
             obj.ins_freq_all=Ins_freq_all;
             obj.ins_energy_all=Ins_energy_all;
         
%              save freq_all.mat Ins_freq_all
%              save energy_all.mat Ins_energy_all
             obj.index_imf = ones(obj.imfbin,1);  
 
             obj.time_param = abs(randn(obj.imfbin,1));
             obj.freq_param = abs(randn(obj.imfbin,1));
             obj.energy_param = abs(randn(obj.imfbin,1));
             obj.lb=zeros(3*obj.imfbin,1);
             obj.ub=zeros(3*obj.imfbin,1);
             obj.ub(1:obj.imfbin)=length(obj.ins_freq_all(1,:,1))/obj.Fs;
             obj.ub(obj.imfbin+1:2*obj.imfbin)=max(max(max(obj.ins_freq_all)));
             obj.ub(2*obj.imfbin+1:end)=max(max(max(obj.ins_energy_all)));
 
         end
         
         function yt =init_param(obj, Xtt,optobj)    
             
                Np=obj.imfbin;
                obj.index_imf = Xtt(3*Np+1:end)';
                obj.time_param=Xtt(1:Np)';
                obj.freq_param = Xtt(Np+1:2*Np)';
                obj.energy_param=Xtt(2*Np+1:3*Np)'; 
                
                LD=size(obj.data,1);
                switch optobj
                    case 'robust'
                        yt=obj.get_robust(LD);

                    case 'class'
                        yt=obj.get_class(LD);
                end
         end
         function [Xs,Xt,Yt] = init_robust(obj)
          %assing the training set (Xs, Xt, Yt)  
                ns =2000;
                nt=10;
                dim=size(obj.ub,1);
                Xs = rand(ns, dim);
                range=obj.ub-obj.lb;
                Trans_Xs=diag(range)*Xs'+obj.lb*ones(1,ns);
                Xs=Trans_Xs';         
                Nindex = randsrc(ns,length(obj.index_imf),[0 1]);
                for index =1:ns
                    if sum(Nindex(index,:))==0 
                        Nindex(index,1)=1;
                    end     
                end
                Xs=[Xs, Nindex];
                
                
                
                out = randsrc(nt,1,0:1:ns);
                Xt=Xs(out,:);
                Yt=zeros(size(Xt,1),1);
                LT=length(out);
                 ii=obj.imfbin;
                for index = 1: LT    
                    Xtt= Xt(index,:);            
                    obj.index_imf= Xtt(3*ii+1:end)';
                    obj.time_param=Xtt(1:ii)';
                    obj.freq_param = Xtt(ii+1:2*ii)';  
                    obj.energy_param = Xtt(3*ii+1:end)';
                    LD=size(obj.data,1);
                    Yt(index) = obj.get_robust(LD);
                end
               
         end
         function [Xs,Xt,Yt] = init_class(obj)
          %assing the training set (Xs, Xt, Yt)  
    %assing the training set (Xs, Xt, Yt)  
                ns =2000;
                nt=10;
                dim=size(obj.ub,1);
                Xs = rand(ns, dim);
                range=obj.ub-obj.lb;
                Trans_Xs=diag(range)*Xs'+obj.lb*ones(1,ns);
                Xs=Trans_Xs';         
                Nindex = randsrc(ns,length(obj.index_imf),[0 1]);
                for index =1:ns
                    if sum(Nindex(index,:))==0 
                        Nindex(index,1)=1;
                    end     
                end
                Xs=[Xs, Nindex];
                
                
                
                out = randsrc(nt,1,0:1:ns);
                Xt=Xs(out,:);
                Yt=zeros(size(Xt,1),1);
                LT=length(out);
                 ii=obj.imfbin;
                for index = 1: LT    
                    Xtt= Xt(index,:);            
                    obj.index_imf= Xtt(3*ii+1:end)';
                    obj.time_param=Xtt(1:ii)';
                    obj.freq_param = Xtt(ii+1:2*ii)';  
                    obj.energy_param = Xtt(3*ii+1:end)';
                    LD=size(obj.data,1);
                    Yt(index) = obj.get_class(LD);
                end
                  
               
         end         

         
         function [y] = get_robust(obj,LD)
             % calculate the robustness for the set of trajectories
             robust=zeros(LD,1);
            if obj.parl
                 parfor index = 1:LD
                    ins_freqs =obj.ins_freq_all(:,:,index);
                    ins_energys=obj.ins_energy_all(:,:,index);
                    robust(index) = Robustness(obj,ins_freqs,ins_energys); 
                 end
            else
                for index = 1:LD
                    ins_freqs =obj.ins_freq_all(:,:,index);
                    ins_energys=obj.ins_energy_all(:,:,index);
                    robust(index) = Robustness(obj,ins_freqs,ins_energys); 
                end     
            end
                robust = robust.*obj.label;
            
              y= min(robust) -obj.alph*sum(obj.index_imf);
         end
         
         function y = get_class(obj,LD)
             % calculate the miss class rate for the set of trajectories.
             robust=zeros(LD,1);
             if obj.parl
                 parfor index = 1:LD
                    ins_freqs =obj.ins_freq_all(:,:,index);
                    ins_energys=obj.ins_energy_all(:,:,index);
                    robust(index,1) = Robustness(obj,ins_freqs,ins_energys); 
                 end   
             else
                 for index = 1:LD
                    ins_freqs =obj.ins_freq_all(:,:,index);
                    ins_energys=obj.ins_energy_all(:,:,index);
                    robust(index,1) = Robustness(obj,ins_freqs,ins_energys); 
                 end                  
                 
             end
                   
             class1 =robust>=0; 
             class2 = robust<0;
             class = class1-class2;
             miss_cla = sum(abs(class - obj.label))/(2*length(robust));         
             y = 1-miss_cla -obj.alph*sum(obj.index_imf);          
         end
         
         function [Xt,Yt] = get_act_robust(obj,Xs,Xt,Yt,varargin)
             %active learning for the chosen Xt+1
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
                ip.addOptional('optobj','robust')
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
                     %   fprintf('sample: %f\n',Xtt);
                        Np=obj.imfbin;

                        obj.index_imf = Xtt(3*Np+1:end)';
                        obj.time_param=Xtt(1:Np)';
                        obj.freq_param = Xtt(Np+1:2*Np)';
                        obj.energy_param=Xtt(2*Np+1:3*Np)';
                        LD=size(obj.data,1);

                   
                    switch opt.optobj
                        case 'robust'
                            yt=obj.get_robust(LD);

                        case 'class'
                            yt=obj.get_class(LD);
                    end
                    fprintf('the value is %f\n',yt);

                    if ~isnan(yt)

                    Xt = [Xt; Xs(xt(end),:)]; % if batch, Xt is already partially updated
                    Yt = [Yt; yt];
                    else
                     Xs(xt(end),:)
                    end



                    % GP update
                    Ht = [Ht; opt.basis(Xt(end,:))];
                    Ktt12 = opt.kernel(Xt(1:end-1,:),Xt(end,:));
                    Ktt22 = opt.kernel(Xt(end,:),Xt(end,:));
                    BayesInv = gp_inf_update(Ktt12, Ktt22, Yt, opt.noise, BayesInv, Ht);
                    Kts = [Kts; opt.kernel(Xt(end,:),Xs)];
                    
                    fprintf('iteration number is : %d\n',iter);
                end

                          
         end
         
         function [Xt,Yt] = get_act_class(obj,Xs,Xt,Yt)
             % the overall learning process
             fprintf('hello\n');
         end
    end
end