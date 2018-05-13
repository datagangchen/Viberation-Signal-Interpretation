function [ r ] =QMITL_RobustnessSyn(Sys,P, phi, params,theta,bias)
%Given a set of trajector and formula phi, calculte the robustness
%  
Ptmp = SetParam(P, params,[2,100]);
Ptmp=SetParam(Ptmp,[23,24],theta);

Pf=ComputeTraj(Sys, Ptmp, Sys.tspan);
traj=Pf.traj;
rtmp=[];
for index=1:size(traj,1)
[val,tau]= QMITL_Eval(Sys, phi, Ptmp, traj(index,:),0);
rtmp(index)=min(val);
end
r=min(rtmp)-bias;
end

