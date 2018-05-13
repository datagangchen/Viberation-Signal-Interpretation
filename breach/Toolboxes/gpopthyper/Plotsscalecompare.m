function Plotsscalecompare()
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
        % steps=0.1;
        % Inter=3;
        % degree=6;
        % %%GP-ACB (Hyper Opt)
        % optimal=load('Optimal200-2true_0.5_gpacb.mat');
        % optimal=optimal.OptimialValue;
        % optimal=sum(optimal)/2;
        % Diff=diff(optimal);
        % index=find(Diff<0);
        % x=1:1:length(optimal);
        % y = optimal;
        % p = polyfit(x,y,degree);
        % x1=1:Inter:300;
        % y1 = polyval(p,x1);
        % yy1 = smooth(x1,y1,steps,'loess')*1.25;
        % 
        % optimal=load('Optimal200-2truegpacb.mat');
        % optimal=optimal.OptimialValue;
        % optimal=sum(optimal)/2;
        % Diff=diff(optimal);
        % index=find(Diff<0);
        % x=1:1:length(optimal);
        % y = optimal;
        % p = polyfit(x,y,degree);
        % x1=1:Inter:300;
        % y1 = polyval(p,x1);
        % yy2= smooth(x1,y1,steps,'loess')*1.25;
        % 
        % yy2(50:end)=yy2(50,1)*ones(length(yy2(50:end)),1);
        % 
        % optimal=load('Optimal200-2true_2_gpacb.mat');
        % optimal=optimal.OptimialValue;
        % optimal=sum(optimal)/2;
        % Diff=diff(optimal);
        % index=find(Diff<0);
        % x=1:1:length(optimal);
        % y = optimal;
        % p = polyfit(x,y,degree);
        % x1=1:Inter:300;
        % y1 = polyval(p,x1);
        % yy3 = smooth(x1,y1,steps,'loess')*1.25;
        % 
        % optimal=load('Optimal200-2true_4_gpacb.mat');
        % optimal=optimal.OptimialValue;
        % optimal=sum(optimal)/2;
        % Diff=diff(optimal);
        % index=find(Diff<0);
        % x=1:1:length(optimal);
        % y = optimal;
        % p = polyfit(x,y,degree);
        % x1=1:Inter:300;
        % y1 = polyval(p,x1);
        % yy4 = smooth(x1,y1,steps,'loess')*1.25;
        % 
        % acbopt=[min(yy1)*2,min(yy2),min(yy3)/2,min(yy4)/4];
        % %% GP-UCB  (Hyper Opt)
        % optimal=load('Optimal200-2true_0.5_gpucb.mat');
        % optimal=optimal.OptimialValue;
        % optimal=sum(optimal)/2;
        % Diff=diff(optimal);
        % index=find(Diff<0);
        % x=1:1:length(optimal);
        % y = optimal;
        % p = polyfit(x,y,degree);
        % x1=1:Inter:300;
        % y1 = polyval(p,x1);
        % yy1 = smooth(x1,y1,steps,'loess')*1.25;
        % 
        % optimal=load('Optimal200-2true.mat');
        % optimal=optimal.OptimialValue;
        % optimal=sum(optimal)/2;
        % Diff=diff(optimal);
        % index=find(Diff<0);
        % x=1:1:length(optimal);
        % y = optimal;
        % p = polyfit(x,y,degree);
        % x1=1:Inter:300;
        % y1 = polyval(p,x1);
        % yy2= smooth(x1,y1,steps,'loess')*1.25;
        % 
        % yy2(50:end)=yy2(50,1)*ones(length(yy2(50:end)),1);
        % 
        % optimal=load('Optimal200-2true_2_gpucb.mat');
        % optimal=optimal.OptimialValue;
        % optimal=sum(optimal)/2;
        % Diff=diff(optimal);
        % index=find(Diff<0);
        % x=1:1:length(optimal);
        % y = optimal;
        % p = polyfit(x,y,degree);
        % x1=1:Inter:300;
        % y1 = polyval(p,x1);
        % yy3 = smooth(x1,y1,steps,'loess')*1.25;
        % 
        % optimal=load('Optimal200-2true_4_gpucb.mat');
        % optimal=optimal.OptimialValue;
        % optimal=sum(optimal)/2;
        % Diff=diff(optimal);
        % index=find(Diff<0);
        % x=1:1:length(optimal);
        % y = optimal;
        % p = polyfit(x,y,degree);
        % x1=1:Inter:300;
        % y1 = polyval(p,x1);
        % yy4 = smooth(x1,y1,steps,'loess')*1.25;
        % ucbopt=[min(yy1)*2,min(yy2),min(yy3)/2,min(yy4)/4];
        % %% GP-ACB
        % optimal=load('Optimal200-2false_0.5_gpacb.mat');
        % optimal=optimal.OptimialValue;
        % optimal=sum(optimal)/2;
        % Diff=diff(optimal);
        % index=find(Diff<0);
        % x=1:1:length(optimal);
        % y = optimal;
        % p = polyfit(x,y,degree);
        % x1=1:Inter:300;
        % y1 = polyval(p,x1);
        % yy1 = smooth(x1,y1,steps,'loess')*1.25;
        % 
        % optimal=load('Optimal200-2falsegpacb.mat');
        % optimal=optimal.OptimialValue;
        % optimal=sum(optimal)/2;
        % Diff=diff(optimal);
        % index=find(Diff<0);
        % x=1:1:length(optimal);
        % y = optimal;
        % p = polyfit(x,y,degree);
        % x1=1:Inter:300;
        % y1 = polyval(p,x1);
        % yy2= smooth(x1,y1,steps,'loess')*1.25;
        % 
        % yy2(50:end)=yy2(50,1)*ones(length(yy2(50:end)),1);
        % 
        % optimal=load('Optimal200-2false_2_gpacb.mat');
        % optimal=optimal.OptimialValue;
        % optimal=sum(optimal)/2;
        % Diff=diff(optimal);
        % index=find(Diff<0);
        % x=1:1:length(optimal);
        % y = optimal;
        % p = polyfit(x,y,degree);
        % x1=1:Inter:300;
        % y1 = polyval(p,x1);
        % yy3 = smooth(x1,y1,steps,'loess')*1.25;
        % 
        % optimal=load('Optimal200-2false_4_gpacb.mat');
        % optimal=optimal.OptimialValue;
        % optimal=sum(optimal)/2;
        % Diff=diff(optimal);
        % index=find(Diff<0);
        % x=1:1:length(optimal);
        % y = optimal;
        % p = polyfit(x,y,degree);
        % x1=1:Inter:300;
        % y1 = polyval(p,x1);
        % yy4 = smooth(x1,y1,steps,'loess')*1.25;
        % acb=[min(yy1)*2,min(yy2),min(yy3)/2,min(yy4)/4];
        % %% GP-UCB
        % optimal=load('Optimal200-2false_0.5_gpucb.mat');
        % optimal=optimal.OptimialValue;
        % optimal=sum(optimal)/2;
        % Diff=diff(optimal);
        % index=find(Diff<0);
        % x=1:1:length(optimal);
        % y = optimal;
        % p = polyfit(x,y,degree);
        % x1=1:Inter:300;
        % y1 = polyval(p,x1);
        % yy1 = smooth(x1,y1,steps,'loess')*1.25;
        % 
        % optimal=load('Optimal200-2false.mat');
        % optimal=optimal.OptimialValue;
        % optimal=sum(optimal)/2;
        % Diff=diff(optimal);
        % index=find(Diff<0);
        % x=1:1:length(optimal);
        % y = optimal;
        % p = polyfit(x,y,degree);
        % x1=1:Inter:300;
        % y1 = polyval(p,x1);
        % yy2= smooth(x1,y1,steps,'loess')*1.25;
        % 
        % yy2(50:end)=yy2(50,1)*ones(length(yy2(50:end)),1);
        % 
        % optimal=load('Optimal200-2false_2_gpucb.mat');
        % optimal=optimal.OptimialValue;
        % optimal=sum(optimal)/2;
        % Diff=diff(optimal);
        % index=find(Diff<0);
        % x=1:1:length(optimal);
        % y = optimal;
        % p = polyfit(x,y,degree);
        % x1=1:Inter:300;
        % y1 = polyval(p,x1);
        % yy3 = smooth(x1,y1,steps,'loess')*1.25;
        % 
        % optimal=load('Optimal200-2false_4_gpucb.mat');
        % optimal=optimal.OptimialValue;
        % optimal=sum(optimal)/2;
        % Diff=diff(optimal);
        % index=find(Diff<0);
        % x=1:1:length(optimal);
        % y = optimal;
        % p = polyfit(x,y,degree);
        % x1=1:Inter:300;
        % y1 = polyval(p,x1);
        % yy4 = smooth(x1,y1,steps,'loess')*1.25;
        % ucb=[min(yy1)*2,min(yy2),min(yy3)/2,min(yy4)/4];
%% plot the data with Demoncomparescale()
optimal=load('scale_prediction.mat');
optimal=optimal.scale_prediction;
y_ucb=sum(optimal(1,:,:)')/length(optimal(1,1,:));
y_acb=sum(optimal(2,:,:)')/length(optimal(1,1,:));
y_ucbopt=sum(optimal(3,:,:)')/length(optimal(1,1,:));
y_acbopt=sum(optimal(4,:,:)')/length(optimal(1,1,:));
x1=.1:0.2:10;
steps=0.1;
ucb= smooth(x1,y_ucb,steps,'loess');
acb= smooth(x1,y_acb,steps,'loess');
ucbopt= smooth(x1,y_ucbopt,steps,'loess');
acbopt= smooth(x1,y_acbopt,steps,'loess');
%% Plot results
x1=[0.5, 1,2,4];
plot(x1,acbopt,'-.', x1,ucbopt,'-', x1,acb,'g-+', x1,ucb,'--', 'LineWidth',2)

legend('GP-ACB (Hyper Opt)','GP-UCB (Hyper Opt)','GP-ACB','GP-UCB')
axis([0 5 0 12])
xlabel('\xi');
ylabel('Prediction Error (%)');

set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [2 1 5 5]);
set(gcf, 'Alphamap',0.1);
set(gcf, 'Colormap', hsv);
set(gcf,'Units', 'inches');
set(gcf,'Position',[3, 3, 6, 6]);
set(gcf,'OuterPosition',[2.5,2.5,6,6])
set(gcf,'Color','white')
set(gca,'FontUnits','points')
set(gca,'FontSize',20);
set(gca,'Color','none');
set(gca,'Box','on');

end

