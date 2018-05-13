function TestDemon()
%This file is used to test the gpotimization method
%   Detailed explanation goes here
%Sample data and 10 training points from a GP with isotropic squared-exponential kernel:
fhandle=@TestFunction;
lb=[1];
ub=[5];
[ Xs, Xt, Yt] = gp_sample_Hyper(fhandle,lb,ub);

%Run the optimization algorithm for 20 iterations:


[queries, observations] = GpOpt_Hyper(fhandle,Xt,Yt,Xs,10,'plot',true,'algo', 'gpucb');

end

