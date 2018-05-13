function [Atree] =CombAtree(Atree1,Atree2,rule)
%This function combine two formula into one and save it in a struct data
%type
%atree1, atree2: the two formula that need combining
%rule          : the production rule used to combine the two formula. which
%                 use index to denotes the rule as follows:
%                 1:A-->  \square A
%                 2:A-->   \diamond A
%                 3:A-->A v B
%                 4:A--> A ^ B
%                 5:A--> \mu
%Atree         : the combined formula which is a data struct has fields as
               %follows:
%               -tree: the structure of the formula
%               -rule: the sequence of rule used to generate the formula
%                  which uses index to denotes the rule as follows:
%                      1:A-->  \square A
%                      2:A-->   \diamond A
%                      3:A-->A v B
%                      4:A--> A ^ B
%                      5:A--> \mu
%               -time: the time bound sequence in the formula
%               -dir:  indicate < />= in predicate 
%                            -1---  <
%                            1---  >=
%               -mag:  the magnitude sequence of the predicate
%               -operator: the temporal operator sequence in the formula which use             
%                          index to denote  the operator
%                         1---  \square
%                         2---- \diamond
                         
%By:  Gang Chen
%Date: 2/11/2018



%% example 




%% start
      lb='(';
      rb=')';
      and=' and ';
      or=' or ';
      imp = '=>';
switch rule
    
    case 1  %\square
        tr1=Atree1.ftree;
        tr2=Atree2.ftree;
        f1=tr1.get(1);
        f2=tr2.get(1);
  
        formula=[ f1 lb f2 rb ];
   %% add the two tree      
        newtree=tree(formula);
        newtree=newtree.graft(1,tr1);
        newtree=newtree.graft(1,tr2);
        
        rule=[Atree1.rule; Atree2.rule];
        time=[Atree1.time; Atree2.time];
        dir=[Atree1.dir; Atree2.dir];
        param_name=[Atree1.param_name; Atree2.param_name];
        param_value=[Atree1.param_value; Atree2.param_value];     
        if isempty(Atree2.op)
         op=[ Atree1.op;3;Atree2.op];    
        else
         op=[ Atree1.op;0;Atree2.op]; 
        end
 %% Create the combined tree       
        Atree.ftree=newtree;
        Atree.rule=rule;
        Atree.time=time;
        Atree.dir=dir;
        Atree.param_name=param_name;
        Atree.param_value=param_value;
        Atree.op=op;
        
        
    case 2
        tr1=Atree1.ftree;
        tr2=Atree2.ftree;
        f1=tr1.get(1);
        f2=tr2.get(1);
        formula=[ f1 lb f2 rb ];
   %% add the two tree      
        newtree=tree(formula);
        newtree=newtree.graft(1,tr1);
        newtree=newtree.graft(1,tr2);
        
        rule=[Atree1.rule; Atree2.rule];
        time=[Atree1.time; Atree2.time];
        dir=[Atree1.dir; Atree2.dir];
        param_name=[Atree1.param_name; Atree2.param_name];
        param_value=[Atree1.param_value; Atree2.param_value];  
        if isempty(Atree2.op)
         op=[ Atree1.op;3;Atree2.op];    
        else
         op=[ Atree1.op;0;Atree2.op]; 
        end   
 %% Create the combined tree       
        Atree.ftree=newtree;
        Atree.rule=rule;
        Atree.time=time;
        Atree.dir=dir;
        Atree.param_name=param_name;
        Atree.param_value=param_value;
        Atree.op=op;      
        
        
    case 3 %or
        tr1=Atree1.ftree;
        tr2=Atree2.ftree;
        f1=tr1.get(1);
        f2=tr2.get(1);
        formula=[lb f1 rb or lb f2 rb];
   %% add the two tree      
        newtree=tree(formula);
        newtree=newtree.graft(1,tr1);
        newtree=newtree.graft(1,tr2);
        
        rule=[Atree1.rule; 3; Atree2.rule];
        time=[Atree1.time; Atree2.time];
        dir=[Atree1.dir; Atree2.dir];
        param_name=[Atree1.param_name; Atree2.param_name];
        param_value=[Atree1.param_value; Atree2.param_value];  
        op=[Atree1.op;0; Atree2.op];    
 %% Create the combined tree       
        Atree.ftree=newtree;
        Atree.rule=rule;
        Atree.time=time;
        Atree.dir=dir;
        Atree.param_name=param_name;
        Atree.param_value=param_value;
        Atree.op=op;        
        
        
        
    case 4%and
        tr1=Atree1.ftree;
        tr2=Atree2.ftree;
        f1=tr1.get(1);
        f2=tr2.get(1);
        formula=[lb  f1 rb and lb f2 rb];
   %% add the two tree      
        newtree=tree(formula);
        newtree=newtree.graft(1,tr1);
        newtree=newtree.graft(1,tr2);
        
        rule=[Atree1.rule; 4; Atree2.rule];
        time=[Atree1.time; Atree2.time];
        dir=[Atree1.dir; Atree2.dir];
        param_name=[Atree1.param_name; Atree2.param_name];
        param_value=[Atree1.param_value; Atree2.param_value];  
        op=[Atree1.op;0; Atree2.op];    
 %% Create the combined tree       
        Atree.ftree=newtree;
        Atree.rule=rule;
        Atree.time=time;
        Atree.dir=dir;
        Atree.param_name=param_name;
        Atree.param_value=param_value;
        Atree.op=op;       
     case 5% =>
        tr1=Atree1.ftree;
        tr2=Atree2.ftree;
        f1=tr1.get(1);
        f2=tr2.get(1);
        formula=[ lb f1 rb imp lb f2 rb];
   %% add the two tree      
        newtree=tree(formula);
        newtree=newtree.graft(1,tr1);
        newtree=newtree.graft(1,tr2);
        
        rule=[Atree1.rule; 5; Atree2.rule];
        time=[Atree1.time; Atree2.time];
        dir=[Atree1.dir; Atree2.dir];
        name = Atree1.param_name;
%         for index =1:length(name)
%             name{index} =  [ name{index},'_do'];
%         end
        param_name=[name; Atree2.param_name];
        param_value=[Atree1.param_value; Atree2.param_value];  
        op=[Atree1.op;0; Atree2.op];    
 %% Create the combined tree       
        Atree.ftree=newtree;
        Atree.rule=rule;
        Atree.time=time;
        Atree.dir=dir;
        Atree.param_name=param_name;
        Atree.param_value=param_value;
        Atree.op=op;          
    otherwise
        
        disp('non existing rule!')
        
end
        
        
        
        





end

