function [p, rob, Pr] = ReqMining(Sys, phi, falsif_opt, prop_opt, iter_max)
%REQMINING requirement mining algorithm
%
% Synopsis: [p, rob, Pfals] = Req_Mining(Sys, template_phi, falsif_opt, prop_opt, iter_max)
%
%  Example:
%  
%  mdl = 'Autotrans_shift';
%  Sys = CreateSimulinkSystem(mdl, {}, {}, [], 'UniStep1');
%  Sys.tspan = 0:.01:50;
%  
%  prop_opt.params = {'vmax', 'rpm_max'};
%  
%  prop_opt.monotony   = [1 1];
%  prop_opt.p_interval = [0 200   ;...  % for vmax
%                         0 6000 ];     % for rmp_max
%  prop_opt.p_tol      = [1 1];
%  
%  falsif_opt.params = {'throttle_u0'};
%  falsif_opt.range = [0 100];
%
%  falsif_opt.iter_max = 100;
%  falsif_opt.nb_init = 10;
%  falsif_opt.nb_iter = 100;
%  
%  [p, rob,Pr] = ReqMining(Sys, phi1, falsif_opt, prop_opt);
%  
%  See also GetPropParamBin 
  

%% Create system and input strategy

%  params_u = falsif_opt.params;
%  ranges = falsif_opt.ranges;
%  nb_init = falsif_opt.nb_init;
  
%  Pu = CreateParamSet(Sys,params_u, ranges);
  
  if isfield(prop_opt, 'order')
    order = prop_opt.order;
  else
    order = 1:numel(prop_opt.params);
  end
  
  % options for parameter synthesis 
  
  params_prop.names = prop_opt.params(order);
  monotony = prop_opt.monotony(order);
  p_interval = prop_opt.p_interval(order,:);
  p_tol = prop_opt.p_tol(order);
 
  Pprop = CreateParamSet(Sys, params_prop.names);
 
  if (~exist('iter_max','var'))
      iter_max = falsif_opt.iter_max;
  end 
        
  %% Arguments for falsification
  
  %% First trajectory
  
  i = 0;
  Pr = ComputeTraj(Sys, Pprop, Sys.tspan); % Pr stores the falsifying trajectories

  % Mine parameters for first simulation
  [p, rob] = GetPropParamBin(Sys, phi, Pr,  params_prop.names, monotony, p_interval, ...
                            p_tol, Pr.traj);
  
  params_prop.values = p;
  i = 1;
  Pr = SetParam(Pr, params_prop.names, params_prop.values');
  
  iter = 1;
  fprintf('\n');
  
  %% Main loop
  while iter<=iter_max
    
    % update display
    clc;    
    fprintf('Iter %d/%d\n\n', iter, iter_max);
    fprintf('Attempts falsifying PSTL formula:\n\n');
    fprintf('%s\n\n', disp(phi)); 
    fprintf('with:\n\n');
    for i = 1:numel(params_prop.names)
      fprintf('%s= %g\n', params_prop.names{i}, p(i))
    end
    fprintf('\n');

    % Falsification step 
    [Propt, rob, nb_call] = Falsify(Sys, phi, falsif_opt, params_prop);
    
    ifalse = find(rob<0);
    
    if ~isempty(ifalse)
        Pfals = Sselect(Propt, ifalse);
        fprintf('\nGetting new parameters for formula:\n');
        Pr = SConcat(Pr,Pfals);
        % Parameter inferrence step
        [p, rob] = GetPropParamBin(Sys, phi, Pr,  params_prop.names, monotony, p_interval, ...
                                   p_tol, Pr.traj);
        params_prop.values = p;
    else
        fprintf(['\np:' num2str(p) '\n']) ;  
        return;
    end
    iter=iter+1;
  end
  fprintf('\nReq mining stopped after max number of iterations.\n');
  fprintf(['\np: ' num2str(p) '\n' ]);
  
 