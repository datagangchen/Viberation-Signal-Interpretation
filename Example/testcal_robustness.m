t=1:1:1000;
t=t/20;
signal=sin(2*pi*t)+sin(5*pi*t)+randn(1);

obj.Fbin=1;
obj.Fs=20;
obj.data=signal;
% [imf,residual,info] = emd(obj.data(1,:),'Interpolation','pchip','Display',false);
% [hs,f,t,imfinsf,imfinse] = hht(imf,obj.Fs);

% obj.IMF=imf(:,1:obj.imfbin)';
% ins_freq =imfinsf(:,1:obj.Fbin)';
% ins_energy=imfinse(:,1:obj.imfbin)';

obj.index_Fbin = ones(obj.Fbin,1);
obj.bound=0.05;

obj.time_param = abs(randn(3*obj.Fbin,1));
obj.freq_param = 100*rand(obj.Fbin,1);
obj.energy_param = obj.freq_param;
[robustness] = ImageRobustness(obj,abs(round(10*randn(100,1000))))



