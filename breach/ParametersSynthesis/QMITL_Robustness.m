function [ r ] =QMITL_Robustness(Sys,P, phi, params,traj,theta,bias)
%Given a set of trajector and formula phi, calculte the robustness
%  
Ptmp = SetParam(P, params, theta);
rtmp=[];
for index=1:size(traj,1)
[val,tau]= QMITL_Eval(Sys, phi, Ptmp, traj(index,:),0);
rtmp(index)=min(val);
end
r=min(rtmp)-bias;
end

