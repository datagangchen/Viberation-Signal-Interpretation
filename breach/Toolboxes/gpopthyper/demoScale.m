function demoScale()
%This example demo the regression process
clear all;
close all;
fhandle=@demoFunctionScale;
lb=[-5.12;-5.12];
ub=[5.12;5.12];

%Max loop number
scale=.1:1:10;
           Loop=1;
           inital=20;
           iter=100;
mindata=zeros(2, length(scale),4);
alldata=zeros(2,length(scale),4,iter+inital);
method={'false', 'true'};
active={'gpucb','gpacb'};
for i=1:2
    for j=1:2
        
        for iid=1:length(scale)

            OptimialValue=zeros(Loop, iter);
            OptimialSigma=zeros(Loop, iter);
            
            for index=1:Loop;
                [ Xs, Xt, Yt] = gp_sample_hyper(fhandle,lb,ub,'nt',inital, 'ns',1500,'Scale',scale(iid) );
                 [Xtnew,Ytnew,Sigma] = ActiveGpOptHyper( fhandle,Xt, Yt, Xs,iter,'plot',false,'opthyper',method{i},'Scale',scale(iid),'algo', active{j});
                for ind=1:iter
                    [tmax,tsigma]=max(Ytnew(1:inital+ind));
                    OptimialValue(index,ind)=-tmax;
                    if tsigma>length(Sigma)
                        tsigma=length(Sigma);
                    end
                    
                    OptimialSigma(index,ind)=Sigma(tsigma(1));
                   
               end
            end
             Optimal=OptimialValue;
             Osigma=OptimialSigma;
             
             [tmin,tindex]=min(Optimal);
              mindata(1,iid,j+2*(i-1))=tmin(1)/scale(iid)*1.25;
              mindata(2,iid,j+2*(i-1))=Osigma(tindex(1));
              alldata(1,iid,j+2*(i-1),1:length(Optimal))=Optimal;
              alldata(2,iid,j+2*(i-1),1:length(Osigma))=Osigma;
        end
    end
  
    
end

save mindata_11_7.mat mindata

save alldata_11_7.mat alldata

 
end

