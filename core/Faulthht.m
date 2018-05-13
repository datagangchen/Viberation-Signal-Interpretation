classdef Faulthht
     properties(SetAccess = private, Hidden)
         Fbin      % number of freqency components used
         Fs          % sampling frequency
         data        %vibration data
         label       % label of the system
         index_Fbin       % the index of chosen imf
         energy_param % param for predicate
         time_param  % param for time
         freq_param
         alph        %regulator
         robustness 
         misclass  
         history    % the history of learning
         freq_vector
         time_freq_all
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
                 obj.Fbin =varargin{1};
                  obj.alph = 0.1;
                 obj.bound =0.05;
                 obj.parl = false;
             elseif length(varargin)==2
                 obj.Fbin =varargin{1};
                 obj.alph = varargin{2};                
                 obj.bound =0.05;
                 obj.parl = false;     
             elseif length(varargin)==3
                 obj.Fbin =varargin{1};
                 obj.alph = varargin{2};
                 obj.bound =varargin{3}; 
                 obj.parl=false;
             elseif length(varargin)==4                    
                 obj.Fbin =varargin{1};
                 obj.alph = varargin{2};
                 obj.bound =varargin{3};
                 obj.parl =varargin{4};
             else 
                 obj.Fbin = 5;
                 obj.alph = 0.1;
                 obj.bound =0.05;
                 obj.parl= false;
             end
             
             LA= size(obj.data,1);
             signal=obj.data;
             Fss=obj.Fs;
  
             Time_freq_all=[];
             Freq_vector =[];
             parfor index =1:LA
                 [tf,f] = HHTspectrum(signal(index,:),Fss);
                 Time_freq_all(:,:,index)=abs(tf);  
                 Freq_vector(index,:) =f';
             end
             obj.time_freq_all=Time_freq_all;
             obj.freq_vector = Freq_vector(1,:);

             obj.index_Fbin = ones(obj.Fbin,1);  
 
             obj.time_param = abs(randn(obj.Fbin,1));
             obj.freq_param = abs(randn(obj.Fbin,1));
             obj.energy_param = abs(randn(obj.Fbin,1));
             obj.lb=zeros(5*obj.Fbin,1);
             obj.ub=zeros(5*obj.Fbin,1);
             obj.ub(1:3*obj.Fbin)=length(obj.time_freq_all(1,:,1))/obj.Fs;
             obj.ub(3*obj.Fbin+1:4*obj.Fbin)=max(size(Time_freq_all,1));
             obj.ub(4*obj.Fbin+1:end)=max(max(max(Time_freq_all)));
 
         end
         
         function yt =init_param(obj, Xtt,optobj)    
             
                Np=obj.Fbin;
                
                obj.time_param=Xtt(1:3*Np)';
                obj.freq_param = Xtt(3*Np+1:4*Np)';
                obj.energy_param=Xtt(4*Np+1:5*Np)'; 
                obj.index_Fbin = Xtt(5*Np+1:end)';
                
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
                nt=30;
                dim=size(obj.ub,1);
                Xs = rand(ns, dim);
                range=obj.ub-obj.lb;
                Trans_Xs=diag(range)*Xs'+obj.lb*ones(1,ns);
                Xs=Trans_Xs';         
                Nindex = randsrc(ns,obj.Fbin,[0 1]);
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
                for index = 1: LT    
                    Xtt= Xt(index,:) ;       
                    Yt(index) = obj.init_param(Xtt, 'robust');
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
                Nindex = randsrc(ns,obj.Fbin,[0 1]);
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
                for index = 1: LT    
                    Xtt= Xt(index,:);               
                    Yt(index) = obj.init_param(Xtt, 'class');
                end
         end         

         
         function [y] = get_robust(obj,LD)
             % calculate the robustness for the set of trajectories
             robust=zeros(LD,1);
             tf =obj.time_freq_all;
            if obj.parl
                 parfor index = 1:LD
                    robust(index) = ImageRobustness(obj,tf(:,:,index)); 
                 end
            else
                for index = 1:LD
                    robust(index) = ImageRobustness(obj,tf(:,:,index)); 
                end     
            end
                robust = robust.*obj.label;
            
              y= min(robust) -obj.alph*sum(obj.index_Fbin);
         end
         
         function y = get_class(obj,LD)
             % calculate the miss class rate for the set of trajectories.
             robust=zeros(LD,1);
             tf =obj.time_freq_all;
             if obj.parl
                 parfor index = 1:LD
                    robust(index,1) = ImageRobustness(obj,tf(:,:,index)); 
                 end   
             else
                 for index = 1:LD
                    robust(index,1) = ImageRobustness(obj,tf(:,:,index)); 
                 end                  
                 
             end
                   
             class1 =robust>=0; 
             class2 = robust<0;
             class = class1-class2;
             miss_cla = sum(abs(class - obj.label))/(2*length(robust));         
             y = 1-miss_cla -obj.alph*sum(obj.index_Fbin);          
         end
         
         function [Xt,Yt] = cal_act_robust(obj,Xs,Xt,Yt)
             %active learning for the chosen Xt+1
              % Gp learning for the parameters
              [Xt,Yt] = get_act_robust(obj,Xs,Xt,Yt,'optobj','robust');
                  
         end
         
         function [Xt,Yt] = cal_act_class(obj,Xs,Xt,Yt)
             % the overall learning process
             [Xt,Yt] = get_act_robust(obj,Xs,Xt,Yt,'optobj','class');
         end
    end
end