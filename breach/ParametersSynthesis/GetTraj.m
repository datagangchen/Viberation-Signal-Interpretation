function [ Pf ] =GetTraj( Sys,  falsif_opt, prop_opt)
%COMPUTETRAJ computes trajectories for a system given initial conditions
% and parameters
% 
% Synopsis:   Pf = ComputeTraj(Sys,P0,tspan [,u])
% 
% Inputs:
% -  Sys   : System (needs to be compiled)
% -  P0    : Initial conditions and params given in a parameter set or in
%            an array of size Sys.DimP x size(P0,2). size(P0,2) is the
%            number of trajectories that will be computed.
% -  tspan : Interval of the form [t0, tf]. Fixed time instants can also
%            be specified tspan = [t0 t1 ... tN];
% -  falsif_opt : the parameter for input and range of the input
%            
%           
% 
% Output:
%  -  Pf   Parameter set augmented with the field traj containing
%          computed trajectories if the input is a param set. The field
%          traj_ref is filled. If the P0 is an  array of parameter values,
%          then Pf is an array of trajectories.
  if isfield(prop_opt, 'order')
    order = prop_opt.order;
  else
    order = 1:numel(prop_opt.params);
  end
  params_prop.names = prop_opt.params(order);
  Pprop = CreateParamSet(Sys, params_prop.names);
  
maxTraj=falsif_opt.nb_max_call;
Bound=falsif_opt.ranges;
lb=Bound(:,1);
ub=Bound(:,2);
Xs =rand( maxTraj,size(Bound,1));
range=ub-lb;
Trans_Xs=diag(range)*Xs'+lb*ones(1,maxTraj);
sample=Trans_Xs';
input=sample(1,:);
Pp=SetParam(Pprop,[23,24],input);
Pf=ComputeTraj(Sys, Pp, Sys.tspan);
for index=2:maxTraj;
    input=sample(index,:);
    Pp=SetParam(Pprop,[23,24],input);
    Pr = ComputeTraj(Sys, Pp, Sys.tspan);% Pr stores the falsifying trajectories;
    Pf=SConcat(Pf,Pr);
end

end

