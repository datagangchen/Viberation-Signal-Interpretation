function [robustness] = Robustness(obj,ins_freq,ins_energy)
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
Num_imf =obj.Fbin;


  % create signal names list
  name =[]; 
  
  for index_j = 1: Num_imf
      str= [ 'imffr_', num2str(index_j)];
      name = [name, {str}];
  end
  
  for index_i =1:Num_imf
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
Tree.ftree=tree;

%% create formula
 
   
        
        for index_i = 1:Num_imf
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
                    t2=num2str((length(ins_freq(1,:)-1)/2)/obj.Fs-obj.time_param(index_i));
                    f1=[al t1 comma t2 rb ];
                    %% generate the first Atree
                    Atree1.ftree=tree(f1);
                    Atree1.rule=[1];
                    Atree1.time=[0,(length(ins_freq(1,:))-1)/obj.Fs];
                    Atree1.dir=[];
                    Atree1.param_name=[];
                    Atree1.param_value=[];
                    Atree1.op=[];
                    %% generate the second atree
                    param_name2= ['t_', num2str(index_i)];

                    f2=[ev t1 comma param_name2 rb];
                    Atree2.ftree=tree(f2);
                    Atree2.rule=[2];
                    Atree2.time=[0, obj.time_param(index_i)];
                    Atree2.dir=[ ];
                    Atree2.param_name={param_name2};
                    Atree2.param_value=[obj.time_param(index_i)];
                    Atree2.op=[];           
                    %% generate the thrid atree
                    name1 = name{index_i};
                    name2 = name{Num_imf+index_i};
                    freq = obj.freq_param(index_i);
                    name1_v1 = num2str((1+obj.bound)*freq);
                    name1_v2 = num2str((1-obj.bound)*freq);
                    f3=[lbr name1 ts leq name1_v1 rbr];
                    Atree3=Atree1;
                    Atree3.ftree=tree(f3);
                    Atree3.rule=[];
                    Atree3.time=[];
                    Atree3.dir=[-1];
                    Atree3.param_name=[];
                    Atree3.param_value=[];
                    Atree3.op=[];

                    %% generate the forth atree


                    f4=[lbr name1 ts geq name1_v2 rbr];
                    Atree4=Atree3;
                    Atree4.ftree=tree(f4);
                    Atree4.rule=[];
                    Atree4.time=[];
                    Atree4.dir=[1];
                    Atree4.param_name=[];
                    Atree4.param_value=[];
                    Atree4.op=[];            
                    %% generate the fifth atree

                    name2_v2  = ['p_', num2str(index_i) ];
                    f5=[lbr name2 ts geq name2_v2 rbr];
                    Atree5=Atree3;
                    Atree5.ftree=tree(f5);
                    Atree5.rule=[];
                    Atree5.time=[];
                    Atree5.dir=[-1];
                    Atree5.param_name={name2_v2};
                    Atree5.param_value=[obj.energy_param(index_i)];
                    Atree5.op=[];  
             %% (w>p1) ^(w<p2) ^ (fe> p_i_j)
                    temp=CombAtree(Atree3,Atree4,4);
                    temp=CombAtree(temp,Atree5,4);
                    temp=CombAtree(Atree2,temp, 2);
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
 
  traj.time=0:1/obj.Fs:length(ins_freq(1,1:end-1))/obj.Fs;
  traj.X = [ins_freq; ins_energy];
  traj.param=p0;  
  P1.traj=traj;
  P1.Xf=traj.X(:,end);
 
robustness = QMITL_Eval(Sys, phi1, P1, P1.traj,0);
end

