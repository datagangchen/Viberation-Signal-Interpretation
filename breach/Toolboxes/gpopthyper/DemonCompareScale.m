function DemonCompareScale()
%Investigate the effect of scale facter to the optimal result for the
%process.

% 2017-11-7

clear all;
close all;
fhandle=@demoFunctionScale;
lb=[-5.12;-5.12];
ub=[5.12;5.12];

%Max loop number
scale=.1:0.2:10;
method={'true', 'false'};
active={'gpucb','gpacb'};

           inital=20;
           iter=300;
           loop=5;
 scale_prediction=zeros(4,length(scale),loop);          

    for i=1:2
        for j=1:2
            for iid=1:length(scale)
              for ind=1:loop   
                    [ Xs, Xt, Yt] = gp_sample_hyper(fhandle,lb,ub,'nt',inital, 'ns',1500,'Scale',scale(iid) );
                    [Xtnew,Ytnew,Sigma] = ActiveGpOptHyper( fhandle,Xt, Yt, Xs,iter,'plot',false,'opthyper',method{i},'Scale',scale(iid),'algo', active{j});
                    [tmax,tsigma]=max(Ytnew);
                       scale_prediction(j+2*(i-1),iid,ind)=-tmax/scale(iid);
             end          
            end
        end
    end
  save scale_prediction.mat scale_prediction  
 end



