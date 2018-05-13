function [tf, freq] = hhtcwt(x,fs, isemd, display,AX)
%UNTITLED3 use imf to get time frequency represention
if nargin==4 
    hf = gcf;
    hf.NextPlot = 'replace';
    AX = axes('parent',hf);
end

imfnum = 5;
 [imf] = emd(x);

 if imfnum> size(imf,2)
     imfnum = size(imf,2);
 end
 tf=[];
if isemd
    for index = 1:imfnum
        [wt,freq] = cwt(imf{index},fs);
        if isempty(tf)
         tf = wt;
        else
          tf = tf+wt;
        end
    end
else
    
   [tf,freq] = cwt(x,fs);  
end
tf = abs(tf);

if display

   t = 0:1/fs:(length(x)-1)/fs;

%    [t,coiweightt,ut] = engunits(t,'unicode','time');
%    xlbl = ['Time (',ut,')'];
%     [freq,coiweightf,uf] = engunits(freq,'unicode');
%     ylbl = ['Frequency (',uf,'Hz)']; 

%     hf = gcf;
%     hf.NextPlot = 'replace';
%     AX = axes('parent',hf);
    imagesc(AX,t,log2(freq),tf);

    Yticks = 2.^(round(log2(min(freq))):round(log2(max(freq))));
    AX.YLim = log2([min(freq), max(freq)]);
    AX.YTick = log2(Yticks);
    AX.YDir = 'normal';
    set(AX,'YLim',log2([min(freq),max(freq)]), ...
        'layer','top', ...
        'YTick',log2(Yticks(:)), ...
        'YTickLabel',num2str(sprintf('%g\n',Yticks)), ...
        'layer','top')

%     title('Magnitude Scalogram');
%     xlabel(xlbl);
%     ylabel(ylbl);
end

end

