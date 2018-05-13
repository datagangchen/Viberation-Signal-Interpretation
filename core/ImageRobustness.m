function [robustness] = ImageRobustness(obj, TimeFreqency)
% This function calculate the robustness with given formula 
% obj is the object for formula 
%         obj.ins_imf : the instanouse imf 
%         obj.inf_energy: the instanouse energy 
%         obj.label   : the label for signal 
%         obj.param   : the parameter for predicate
%         obj.Nimf    : the index for chosen imf
%         obj.Nfreq   : the index for chosen frequency bins
%         obj.FreqIndex: the frequecy index 

%% assemble the system
Num_Fbin =obj.Fbin;


  % create signal names list
  name =[]; 
  
%   for index_j = 1: Num_Fbin
%       str= [ 'imffr_', num2str(index_j)];
%       name = [name, {str}];
%   end
  
  for index_i =1:Num_Fbin
      str = ['imfen_', num2str(index_i)];
      name =[name, {str}];
  end
  
%% create formula tree
            comma=',';           
            al='alw_[';
            ev='ev_[';
            rb=']';
            leq='<';
            geq='>=';
            lbr='(';
            rbr=')';   
            ts='[t]';
            imp = '=>';
Tree.ftree=tree;

%% create formula
 
   
        
        for index_i = 1:Num_Fbin
                if obj.index_Fbin(index_i)==1
                    %% inital the trees struct
                    tree1=tree;
                    Atree1.ftree=tree1;
                    Atree1.rule=[];
                    Atree1.time=[];
                    Atree1.dir=[];
                    Atree1.param_name=[];
                    Atree1.param_value=[];
                    Atree1.op=[];
                    Atree2=Atree1;


                    t1=num2str(0);
 
                    t2=num2str((length(TimeFreqency(1,:))-1)/obj.Fs-obj.time_param(index_i));
                    f1=[al t1 comma t2 rb ]; % alw_[t1,t2]
                    %% generate the first Atree
                    Atree1.ftree=tree(f1);
                    Atree1.rule=[1];
                    Atree1.time=[0,(length(TimeFreqency(1,:,1))-1)/obj.Fs];
                    Atree1.dir=[];
                    Atree1.param_name=[];
                    Atree1.param_value=[];
                    Atree1.op=[];
                    %% generate the second atree
                    param_name2= ['t_1_', num2str(index_i)];

                    f21=[ev t1 comma param_name2 rb]; % ev_[t1, ti]
                    Atree2.ftree=tree(f21);
                    Atree2.rule=[2];
                    Atree2.time=[0, obj.time_param(index_i)];
                    Atree2.dir=[ ];
                    Atree2.param_name={param_name2};
                    Atree2.param_value=[obj.time_param(index_i)];
                    Atree2.op=[];    
                    %% generate the second atree
                    param_name22= ['t_2_', num2str(index_i+Num_Fbin)];

                    f22=[al t1 comma param_name22 rb]; % ev_[t1, ti]
                    Atree22.ftree=tree(f22);
                    Atree22.rule=[2];
                    Atree22.time=[0, obj.time_param(index_i)];
                    Atree22.dir=[ ];
                    Atree22.param_name={param_name22};
                    Atree22.param_value=[obj.time_param(index_i)];
                    Atree22.op=[];  
                     %% generate the second atree
                    param_name23= ['t_3_', num2str(index_i+2*Num_Fbin)];

                    f23=[ev t1 comma param_name23 rb]; % ev_[t1, ti]
                    Atree23.ftree=tree(f23);
                    Atree23.rule=[2];
                    Atree23.time=[0, obj.time_param(index_i)];
                    Atree23.dir=[ ];
                    Atree23.param_name={param_name23};
                    Atree23.param_value=[obj.time_param(index_i)];
                    Atree23.op=[];                     
 
                  %% generate the thrid atree
                    name2 = name{index_i};
                    value = num2str(0.3*obj.energy_param(index_i));
                    f4=[lbr name2 ts leq value rbr];% (e[t]> e1)
                    Atree4=Atree2;
                    Atree4.ftree=tree(f4);
                    Atree4.rule=[];
                    Atree4.time=[];
                    Atree4.dir=[-1];
                    Atree4.param_name={};
                    Atree4.param_value=[];
                    Atree4.op=[];  
                    %% generate the thrid atree
                    value = num2str(obj.energy_param(index_i));
                    f5=[lbr name2 ts geq value rbr];% (e[t]> e1)
                    Atree5=Atree2;
                    Atree5.ftree=tree(f5);
                    Atree5.rule=[];
                    Atree5.time=[];
                    Atree5.dir=[-1];
                    Atree5.param_name={};
                    Atree5.param_value=[];
                    Atree5.op=[];  
             %% (w>p1) ^(w<p2) ^ (fe> p_i_j)
                    temp1=CombAtree(Atree22,Atree4,1);
                    temp2=CombAtree(Atree23,Atree5, 2);
                    temp = CombAtree(temp1, temp2,5);
                    temp=CombAtree(Atree2, temp,2);
                    temp=CombAtree(Atree1,temp,1);

                    
                    if isempty(Tree.ftree.get(1))   
                             Tree = temp;
                    else
                             Tree = CombAtree(Tree, temp,3);
                    end                    
                  
                end
        %% assign inner formula
        end



 
    if isempty(Tree.ftree.get(1))
        robustness=NaN;
        return;
        fprintf('error: The parameter leads to no formula!\n')
    end

%% calculate robustness
  vars = [name];  % variables for signals values
  params = [Tree.param_name];  % parameters related to signal or to be  % used in temporal logic formula
  init=zeros(1,length(name));
  p0 = [ init ...      % default for initial values of signals 
         Tree.param_value' ];    % default for parameters p0,p1 and p2
  Sys = CreateSystem(vars,params, p0); % creates the Sys structure
  P = CreateParamSet(Sys);
  P1 = Refine(P,2);
%% The formula 
  formula=Tree.ftree;
  phi=formula.get(1);
%   Tree.param_value
 [phi1, phistruct] = QMITL_Formula('phi_tmp__', phi);

  traj.time=0:1/obj.Fs:(size(TimeFreqency,2)-1)/obj.Fs;
  Traj = zeros(Num_Fbin, size(TimeFreqency,2));
  for index = 1 : Num_Fbin
      freq = obj.freq_param(index); 
      lbf =max(1, round((1- obj.bound)*freq));
      ubf = min(size(TimeFreqency,1), round((1+obj.bound)*freq));
      ubf = max(ubf,lbf);
      Traj(index,:)=  max(TimeFreqency(lbf:ubf,:));
  end

  traj.X = Traj;
  traj.param=p0;  
  P1.traj=traj;
  P1.Xf=traj.X(:,end);
 
robustness = QMITL_Eval(Sys, phi1, P1, P1.traj,0);
end

