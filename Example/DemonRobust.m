InitBreach
t=0:0.01:10;
x1= linspace(0,10,1001);
formulas = QMITL_ReadFile('spec.stl');
formula = phiss2;
name ={'x1'};
param_name={'p1'};
param_value=2.5;
  vars = [name];  % variables for signals values
  
  params = [param_name];  % parameters related to signal or to be  

  % used in temporal logic formula
  init=zeros(1,length(name));
  p0 = [ init ...      % default for initial values of signals 
         param_value' ];    % default for parameters p0,p1 and p2

  Sys = CreateSystem(vars,params, p0); % creates the Sys structure
  P = CreateParamSet(Sys);
  P1 = Refine(P,3);

 %[phi1, phistruct] = QMITL_Formula('phi_tmp__', formu);%
  traj.X=x1;
  traj.time=t;
  traj.param=p0;  
  P1.traj=traj;
  P1.Xf=traj.X(:,end);

   
  robust = QMITL_Eval(Sys, formula, P1, P1.traj,0)
   