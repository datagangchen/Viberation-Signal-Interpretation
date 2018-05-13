classdef FaultDetection
     properties(SetAccess = private, Hidden)
         Nbin        % number of frequency bins 
         imfbin      % number of imfs used
         Fs          % sampling frequency
         data        %vibration data
         label       % label of the system
         ins_freq    % instantenous frequency
         ins_energy  % instantaneous energy
         FreqIndex   % frequecy index of the fault signal
         Nfreq       % the index of chosen frequency bins
         Nimf        % the index of chosen imf
         scale_param % param for predicate
         time_param  % param for time
         alph        %regulator
         belta       %regulator
         robustness 
         misclass  
         history    % the history of learning
         ins_freq_all
         ins_energy_all
         Episode
         ub
         lb

     end
    
    
    methods   
         function obj = initialize(obj,arg1,arg2,arg3,arg4,varargin)
             obj.data = arg1;
             obj.label= arg2;
             obj.Fs = arg3;
             obj.Episode=arg4;
             if length(varargin)>1
                 obj.Nbin= varargin{1};
                 obj.imfbin =varargin{2};
                 obj.alph = varargin{3};
                 obj.belta = varargin{4};
                 
             else 
                 obj.imfbin = 5;
                 obj.Nbin =10;
                 obj.alph = 0.1;
                 obj.belta = 0.1;
             end
             
             LA= size(obj.data,1);
             signal=obj.data;
             Fss=obj.Fs;
             Imfbin=obj.imfbin;
             Ins_freq_all=obj.ins_freq_all;
             Ins_energy_all=obj.ins_energy_all;
             parfor index =1:LA
                 [imf,residual,info] = emd(signal(index,:),'Interpolation','pchip','Display',false);
                 [hs,f,t,imfinsf,imfinse] = hht(imf,Fss);
                 Ins_freq_all(:,:,index)=imfinsf(:,1:Imfbin)';
                 Ins_energy_all(:,:,index)=imfinse(:,1:Imfbin)';
             end
             obj.ins_freq_all=Ins_freq_all;
             obj.ins_energy_all=Ins_energy_all;
             obj.ins_freq =Ins_freq_all(:,:,1);
             obj.ins_energy=Ins_energy_all(:,:,1);
             obj.FreqIndex = obj.Fs/8:obj.Fs/(4*obj.Nbin):obj.Fs*3/8;
             
             obj.Nimf = ones(obj.imfbin,1);
             obj.Nfreq = ones(obj.Nbin,1);   
             obj.scale_param = randn(obj.imfbin,obj.Nbin);
             obj.time_param = abs(randn(obj.imfbin,obj.Nbin));
             obj.lb=zeros(2*obj.imfbin*obj.Nbin,1);
             obj.ub=zeros(2*obj.imfbin*obj.Nbin,1);
             obj.ub(1:obj.imfbin*obj.Nbin)=length(obj.ins_freq(1,:))/obj.Fs;
             obj.ub(obj.imfbin*obj.Nbin+1:end)=max(max(obj.ins_energy));
             
         end
         
         function [Xs,Xt,Yt] = init_robust(obj)
          %assing the training set (Xs, Xt, Yt)  
                ns =10000;
                nt=10;
                dim=size(obj.ub,1);
                Xs = rand(ns, dim);
                range=obj.ub-obj.lb;
                Trans_Xs=diag(range)*Xs'+obj.lb*ones(1,ns);
                Xs=Trans_Xs';         
                Nindex = randsrc(ns,length(obj.Nimf)+length(obj.Nfreq),[0 1]);
                for index =1:ns
                    if sum(Nindex(index,1:length(obj.Nimf)))==0 || sum(Nindex(index,length(obj.Nimf)+1:end))==0
                        Nindex(index,1)=1;
                        Nindex(index,end)=1;
                    end
                    
                end
                Xs=[Xs, Nindex];
                
                
                
                out = randsrc(nt,1,0:1:ns);
                Xt=Xs(out,:);
                Yt=zeros(size(Xt,1),1);
                scale= zeros(size(obj.scale_param));
                time = zeros(size(obj.time_param));
                
                LT=length(out);
                for index = 1: LT
                    [ii,jj]=size(scale);
                    Xtt= Xt(index,:);
                    parfor index_i =1:ii
                        scale(index_i,:)=Xtt(jj*(index_i-1)+1:jj*index_i);
                        time(index_i,:) =Xtt(jj*ii+jj*(index_i-1)+1: jj*ii +jj*index_i);  
                    end
                    obj.Nimf = Xtt(2*ii*jj+1:2*ii*jj+length(obj.Nimf));
                    obj.Nfreq= Xtt(2*ii*jj+length(obj.Nimf)+1:end);
                    obj.time_param=time;
                    obj.scale_param = scale;
                    LD=size(obj.data,1);
                    Yt(index) = obj.get_robust(LD);
                end
               
         end
         function [Xs,Xt,Yt] = init_class(obj, ub,lb)
          %assing the training set (Xs, Xt, Yt)  
                ns =1000;
                nt=10;
                dim=size(ub,1);
                Xs = rand(ns, dim);
                range=ub-lb;
                Trans_Xs=diag(range)*Xs'+lb*ones(1,ns);
                Xs=Trans_Xs';
                Nfreqs = randsrc(ns,length(obj.Nfreq),[0 1]);
                Nimfs = randsrc(ns,length(obj.Nimf),[0 1]);
                Xs=[Xs, Nfreqs, Nimfs];
                out = randsrc(nt,1,0:1:ns);
                Xt=Xs(out,:);
                Yt=zeros(size(Xt,1),1);
                scale= zeros(size(obj.scale_param));
                time = zeros(size(obj.time_param));
                LT=length(out);
                for index = 1: LT
                    [ii,jj]=size(scale);
                    Xtt= Xt(index,:);
                    parfor index_i =1:ii
                        scale(index_i,:)=Xtt(jj*(index_i-1)+1:jj*index_i);
                        time(index_i,:) =Xtt(jj*ii+jj*(index_i-1)+1: jj*ii +jj*index_i);  
                    end
                    obj.Nimf = Xtt(2*ii*jj+1:2*ii*jj+length(obj.Nimf));
                    obj.Nfreq= Xtt(2*ii*jj+length(obj.Nimf)+1:end);
                    obj.time_param=time;
                    obj.scale_param = scale;
                    Yt(index) = obj.get_class();
                end
               
         end         

         
         function [y] = get_robust(obj,LD)
             % calculate the robustness for the set of trajectories
         
             parfor index = 1:LD
                ins_freqs =obj.ins_freq_all(:,:,index);
                ins_energys=obj.ins_energy_all(:,:,index);
                robust(index) = Cal_Robustness(obj,ins_freqs,ins_energys); 
             end       
              y= min(robust) -obj.alph*sum(obj.Nfreq)-obj.belta*sum(obj.Nimf);
         end
         
         function y = get_class(obj)
             % calculate the miss class rate for the set of trajectories.
             robust= zeros(size(obj.data,1),1);
             LD=length(robust);
             parfor index = 1:LD
                ins_freqs =obj.ins_freq_all(:,:,index);
                ins_energys=obj.ins_energy_all(:,:,index);
                robust(index) = Cal_Robustness(obj,ins_freqs,ins_energys); 
             end               
             class =robust>0;             
             miss_cla = abs(class - obj.label)/(2*length(robust));            
             y = 1-miss_cla -obj.alph*sum(obj.Nfreq)-obj.belta*sum(obj.Nimf);             
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
                        LD=size(obj.data,1);

                   
                    switch opt.optobj
                        case 'robust'
                            yt=obj.get_robust(LD);

                        case 'class'
                            yt=obj.get_class(LD);
                    end

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