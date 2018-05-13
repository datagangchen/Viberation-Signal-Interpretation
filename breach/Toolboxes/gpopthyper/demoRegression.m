function demoRegression()
%This example demo the regression process
clear all;
close all;
 fhandle=@demoFunction;
lb=[-5.12];
ub=[5.12];

%Max loop number
Loop=2;
inital=10;
iter=30;
OptimialValue=zeros(Loop, iter);
for index=1:Loop;
[ Xs, Xt, Yt] = gp_sample_hyper(fhandle,lb,ub,'nt',inital, 'ns',1500);

 [Xtnew,Ytnew] = ActiveGpOptHyper( fhandle,Xt, Yt, Xs,iter,'plot',false,'opthyper',false);
for ind=1:iter
    OptimialValue(index,ind)=-max(Ytnew(1:inital+ind));
end
end
Optimal=sum(OptimialValue)/Loop;

plot(Optimal);


% figure
% x=Xtnew(:,1);
% y=Xtnew(:,2);
% z=Ytnew;
% scatter3(x,y,z,'*');
% view(-30,10)
 
end

