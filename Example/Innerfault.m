%% artifical bearing data
function [results] = Innerfault(arg1)


if nargin>0
    method = arg1;
else
    method = 'robust';
end

dataInner=  load('209.mat');
xInner = dataInner.X209_DE_time;
fsInner = dataInner.X209RPM/60;

dataNorm  = load('97.mat');
xNorm =dataNorm.X097_DE_time;
fsNorm = dataNorm.X097RPM/60;

dataOuter=  load('234.mat');
xOuter = dataOuter.X234_DE_time;
fsOuter = dataOuter.X234RPM/60;

dataBall=  load('222.mat');
xBall = dataBall.X222_DE_time;
fsBall = dataBall.X222RPM/60;


LS=1000;
LN=60;
%% Generate data set
data=[];
parfor index =1:round(LN/3)
    data=[data; xNorm(LS*(index-1)+1:LS*index,1)']
end
parfor index =1:round(LN/3)
    data=[data; xOuter(LS*(index-1)+1:LS*index,1)']
end
parfor index =1:round(LN/3)
    data=[data; xBall(LS*(index-1)+1:LS*index,1)']
end
parfor index =1:LN
    data=[data; xInner(LS*(index-1)+1:LS*index,1)']
end

label=-ones(LN+3*round(LN/3),1);
label(3*round(LN/3)+1:end)=-label(LN+1:end);


Fs=12000;
Fbin=5;
Episode=500;
alpha = 0.001;
bound = 0.1;
Bear= Faulthht;
Bear=Bear.initialize(data,label,Fs,Episode, Fbin,alpha,bound,true);


switch method
    case 'class'
    [Xs,Xt,Yt] = Bear.init_class();
    [Xt,Yt] = Bear.cal_act_class(Xs,Xt,Yt);
    case 'robust'
     [Xs,Xt,Yt] = Bear.init_robust();
    [Xt,Yt] = Bear.cal_act_robust(Xs,Xt,Yt);       
        
end

[a,b] = max(Yt);
params = Xt(b,:);
index =find(params(5*Fbin+1:end) ==1);
time1 = params(index);
time2 = params(index+Fbin);
time3 = params(index + 2*Fbin);
freq = params(index+3*Fbin);
energy = params(index+4*Fbin);
results.performance = a;
results.time1 = time1;
results.time2 = time2;
results.time3 = time3;
results.freq = freq;
results.energy = energy;
results.methods = method;
results.type = 'Inner';
results.Xt = Xt;
results.Yt = Yt;
end
% results = [{Xt}, {Yt},{param}];

