function [spectrum, frequency] = HHTspectrum(x,Fs)
% Get the spectrum for input imfs

 [imf,residual,info] = emd(x,'Interpolation','pchip','Display',false);
 [hs,f,t,imfinsf,imfinse] = hht(imf,Fs);
 
 imfnum =10;
 freqbin = 100;
 if imfnum > size(imf,2)
     imfnum = size(imf,2);
 end
 
 tf= [];
 frequency = linspace(0,Fs/2,freqbin);
 delta = frequency(2)-frequency(1);
 for index =1:imfnum
    f = round(abs(imfinsf(:,index))/delta);
   [m,n] = size(imfinsf(:,index));
    One = ones(freqbin,m);
    diags = diag( freqbin:-1:1);
    Increase = diags*One;
    Zeros = Increase - ones(freqbin,1)*f';
    Logic = Zeros == 0;
    temptf = Logic*diag(imfinse(:,index));
    if isempty(tf)
        tf = temptf;
    else
        tf = tf+temptf;
    end
    
    
 end
     
     
     
     
 spectrum = abs(tf);
end

