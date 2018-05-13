% clear all
close all

dataInner=  load('209.mat');

xInner = dataInner.X209_DE_time;
fsInner = dataInner.X209RPM/60;
tInner = (0:length(xInner)-1)/fsInner;

dataNorm  = load('97.mat');
xNorm =dataNorm.X097_DE_time;
fsNorm = dataNorm.X097RPM/60;

LS=1000;
LN=50;
%% Generate data set
data=[];
parfor index =1:LN
    data=[data; xNorm(LS*(index-1)+1:LS*index,1)']
end

parfor index =1:LN
    data=[data; xInner(LS*(index-1)+1:LS*index,1)']
end

label=-ones(2*LN,1);
label(LN+1:end)=-label(LN+1:end);


Fs=12000;
imfbin=1;
Episode=500;
Bear= Fault;
Bear=Bear.initialize(data,label,Fs,Episode, imfbin,0.001,0.4,true);
[Xs,Xt,Yt] = Bear.init_class();
Xtt=[0.04,5000,1,1];
Yt=Bear.init_param(Xtt,'robust');

