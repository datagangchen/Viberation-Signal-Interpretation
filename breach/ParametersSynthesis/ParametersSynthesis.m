function [hyperplane] =ParametersSynthesis(Sys, phi, falsif_opt, prop_opt, iter_max)
%
% Synopsis: [p, rob, Pfals] = ParametersSynthesis(Sys, template_phi, falsif_opt, prop_opt, iter_max)
% 
%hyperplane : the hyperplane for the parameters;
%
%  Example:
%  
%  mdl = 'epasCar';
%  Sys = CreateSimulinkSystem(mdl, {}, {}, [], 'UniStep1');
%  Sys.tspan = 0:.01:10;
%  
% prop_opt.params = {'wy', 'alat'};
% prop_opt.monotony   = [1 1];
% prop_opt.p_interval = [0 0.8   ;...  % for wy
%                        0 6 ];     % for alat
% prop_opt.p_tol      = [1 1];
%   
% System and falsification parameters
% falsif_opt.params = {'u0_u0' ... ,
%                      'bL_u0'... 
%                     };
% falsif_opt.ranges = [0 50; ...
%                      0 100; ...   
% ];
% falsif_opt.nb_init = 6;
% falsif_opt.nb_iter = 1000;
% falsif_opt.nb_max_call = 1000;
%  
%  [p, rob,Pr] = ParametersSynthesis(Sys, phi1, falsif_opt, prop_opt);
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
  
  %% This part is for EPASCar model
%   input=[0.002,140];
%   Pprop=SetParam(Pprop,[21,22],input);
%   Pf = ComputeTraj(Sys, Pprop, Sys.tspan);% Pr stores the falsifying trajectories;
%   
%   Traj_value=Pr.traj.X;
%   tspan=Pr.traj.time;
%   wy=Traj_value(1,:);
%   alat=Traj_value(2,:);
%   [wy_OS,wy_RT,wy_PT,yawGain,alat_OS,alat_RT,alat_PT]=FeatureExtract(tspan,wy,alat);
%   plot(tspan,wy);
%   figure
%   plot(tspan,alat)
%% This part is for automation model
  [ Pf ] =GetTraj( Sys,falsif_opt, prop_opt);
  a=Pf.traj
  
  ubound=p_interval(:,2);
  lbound=p_interval(:,1);
  %% function estimation
  %bais=0;
 
  Sdata=[];
  bai=-1:0.5:1;
 for index=1:5
    bais=bai(index);
    fun = @(theta) QMITL_Robustness(Sys,Pprop, phi, params_prop.names,Pf.traj,theta,bais); 
   [Xs,Xt, Yt] = gp_sample(fun,lbound,ubound,'plot',false);
   [queries, label,SYt0,SXt0] = GpOpt(fun,Xt,Yt,Xs, iter_max,'algo','gpucb','plot',false,'scale',4,'threshold',0.05);
   samples=[SYt0, SXt0];
   save sample_gr_Jm.mat samples
%    Frontier
%    hold on
   delete sample_gr_Jm.mat
   Sdata=[Sdata;samples];
 end
 save gr_Jm.mat Sdata
%   X=Samples(:,1);
%   Y=Samples(:,2);
% type = 'function estimation';
% [Yp,alpha,b,gam,sig2] = lssvm(X,Y,type);
% xt0=0:0.1:ubound(1);
% yt0 = simlssvm({X,Y,type,gam,sig2,'RBF_kernel','preprocess'},{alpha,b},xt');
% temp=[xt0',yt0];
% plot(xt0',yt0);
% figure; plotlssvm({X,Y,type,gam,sig2,'RBF_kernel','preprocess'},{alpha,b});
%   %bais=40;
%   bais=40;
%   fun = @(theta) QMITL_Robustness(Sys,Pprop, phi, params_prop.names,Pf.traj,theta,bais);
%   [Xs,Xt, Yt] = gp_sample(fun,lbound,ubound,'plot',false);
%   [queries, label,Yt1] = GpOpt(fun,Xt,Yt,Xs, iter_max,'algo','gpkt','plot',true,'scale',0.05,'threshold',10);
%   Samples=Xs(queries,:);
%   X=Samples(:,1);
%   Y=Samples(:,2);
% type = 'function estimation';
% [Yp,alpha,b,gam,sig2] = lssvm(X,Y,type);
% xt40=0:0.1:ubound(1);
% yt40 = simlssvm({X,Y,type,gam,sig2,'RBF_kernel','preprocess'},{alpha,b},xt');
% 
%   %bais=80;
%   bais=80;
%   fun = @(theta) QMITL_Robustness(Sys,Pprop, phi, params_prop.names,Pf.traj,theta,bais);
%   [Xs,Xt, Yt] = gp_sample(fun,lbound,ubound,'plot',false);
%   [queries, label,Yt2] = GpOpt(fun,Xt,Yt,Xs, iter_max,'algo','gpkt','plot',true,'scale',0.05,'threshold',10);
%   Samples=Xs(queries,:);
%   X=Samples(:,1);
%   Y=Samples(:,2);
% type = 'function estimation';
% [Yp,alpha,b,gam,sig2] = lssvm(X,Y,type);
% xt80=0:0.1:ubound(1);
% yt80 = simlssvm({X,Y,type,gam,sig2,'RBF_kernel','preprocess'},{alpha,b},xt');
% 
% 
%   %bais=-40;
%   bais=-40;
%   fun = @(theta) QMITL_Robustness(Sys,Pprop, phi, params_prop.names,Pf.traj,theta,bais);
%   [Xs,Xt, Yt] = gp_sample(fun,lbound,ubound,'plot',false);
%   [queries, label,Yt3] = GpOpt(fun,Xt,Yt,Xs, iter_max,'algo','gpkt','plot',true,'scale',0.005,'threshold',10);
%   Samples=Xs(queries,:);
%   X=Samples(:,1);
%   Y=Samples(:,2);
% type = 'function estimation';
% [Yp,alpha,b,gam,sig2] = lssvm(X,Y,type);
% xt_40=0:0.1:ubound(1);
% yt_40 = simlssvm({X,Y,type,gam,sig2,'RBF_kernel','preprocess'},{alpha,b},xt');
% 
%   %bais=-80;
%   bais=-80;
%   fun = @(theta) QMITL_Robustness(Sys,Pprop, phi, params_prop.names,Pf.traj,theta,bais);
%   [Xs,Xt, Yt] = gp_sample(fun,lbound,ubound,'plot',false);
%   [queries, label,Yt4] = GpOpt(fun,Xt,Yt,Xs, iter_max,'algo','gpkt','plot',true,'scale',0.005,'threshold',10);
%   Samples=Xs(queries,:);
%   X=Samples(:,1);
%   Y=Samples(:,2);
% type = 'function estimation';
% [Yp,alpha,b,gam,sig2] = lssvm(X,Y,type);
% xt_80=0:0.1:ubound(1);
% yt_80 = simlssvm({X,Y,type,gam,sig2,'RBF_kernel','preprocess'},{alpha,b},xt');
% 
% Bounds=zeros(length(xt),6);
% Bounds(:,1)=xt';
% Bounds(:,2)=yt40;
% Bounds(:,3)=yt80;
% Bounds(:,4)=yt0;
% Bounds(:,5)=yt_40;
% Bounds(:,6)=yt_80;
%   save Bounds.mat Bounds
%   SampleY=zeros(length(Yt1),5);
%   SampleY(:,1)=Yt0;
%   SampleY(:,2)=Yt1;
%   SampleY(:,3)=Yt2;
%   SampleY(:,4)=Yt3;
%   SampleY(:,5)=Yt4;
%     save SampleY.mat SampleY
% plot(Bounds(:,1),Bounds(:,2),Bounds(:,1),Bounds(:,3),Bounds(:,1),Bounds(:,4),Bounds(:,1),Bounds(:,5),Bounds(:,1),Bounds(:,6));
%% classification type
%   Samples_0=Xs(queries,:);
%   label_0=label;
%   bais=30;
%    %bais=50
%   fun = @(theta) QMITL_Robustness(Sys,Pprop, phi, params_prop.names,Pf.traj,theta,bais);
%   [Xs,Xt, Yt] = gp_sample(fun,lbound,ubound,'plot',false);
%   [queries, label] = GpOpt(fun,Xt,Yt,Xs, iter_max,'algo','gpkt','plot',false,'scale',1,'threshold',100);
%   Samples_1=Xs(queries,:);
%   label_1=label;
%   %bias=100;
%    bais=60; 
%   fun = @(theta) QMITL_Robustness(Sys,Pprop, phi, params_prop.names,Pf.traj,theta,bais);
%   [Xs,Xt, Yt] = gp_sample(fun,lbound,ubound,'plot',false);
%   [queries, label] = GpOpt(fun,Xt,Yt,Xs, iter_max,'algo','gpkt','plot',false,'scale',1,'threshold',100);
%   Samples_2=Xs(queries,:);
%   label_2=label;
%      bais=-60; 
%   fun = @(theta) QMITL_Robustness(Sys,Pprop, phi, params_prop.names,Pf.traj,theta,bais);
%   [Xs,Xt, Yt] = gp_sample(fun,lbound,ubound,'plot',false);
%   [queries, label] = GpOpt(fun,Xt,Yt,Xs, iter_max,'algo','gpkt','plot',false,'scale',1,'threshold',50);
%   Samples_4=Xs(queries,:);
%   label_4=label;
%      bais=-30; 
%   fun = @(theta) QMITL_Robustness(Sys,Pprop, phi, params_prop.names,Pf.traj,theta,bais);
%   [Xs,Xt, Yt] = gp_sample(fun,lbound,ubound,'plot',false);
%   [queries, label] = GpOpt(fun,Xt,Yt,Xs, iter_max,'algo','gpkt','plot',false,'scale',1,'threshold',100);
%   Samples_3=Xs(queries,:);
%   label_3=label;
%   gam = 10;
%   sig2 = 0.2;
%   type ='classification';
%   colormap cool;
%   map = colormap;
%   [alpha,b] = trainlssvm({Samples_0,label_0',type,gam,sig2,'RBF_kernel','preprocess'}); 
%   [XX,YY,ZZd1,model]= GetContour({Samples_0,label_0',type,gam,sig2,'RBF_kernel','preprocess'},{alpha,b});
%   [alpha,b] = trainlssvm({Samples_1,label_1',type,gam,sig2,'RBF_kernel','preprocess'}); 
%   [XX,YY,ZZd2,model]= GetContour({Samples_1,label_1',type,gam,sig2,'RBF_kernel','preprocess'},{alpha,b});
%   [alpha,b] = trainlssvm({Samples_2,label_2',type,gam,sig2,'RBF_kernel','preprocess'}); 
%   [XX,YY,ZZd3,model]= GetContour({Samples_2,label_2',type,gam,sig2,'RBF_kernel','preprocess'},{alpha,b});
%    [alpha,b] = trainlssvm({Samples_3,label_3',type,gam,sig2,'RBF_kernel','preprocess'}); 
%   [XX,YY,ZZd4,model]= GetContour({Samples_3,label_3',type,gam,sig2,'RBF_kernel','preprocess'},{alpha,b});
%    [alpha,b] = trainlssvm({Samples_4,label_4',type,gam,sig2,'RBF_kernel','preprocess'}); 
%   [XX,YY,ZZd5,model]= GetContour({Samples_4,label_4',type,gam,sig2,'RBF_kernel','preprocess'},{alpha,b});
%   ZZd=ZZd1+ZZd2+ZZd3+ZZd4+ZZd5;
%    eval('[C,h]=contourf(XX,YY,ZZd);','warning(''no surface plot feasable'');'); 
 hyperplane=0;

end

