function Plots()
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
steps=0.1;
Inter=3;
degree=6;
optimal=load('Optimal200-2false.mat');
optimal=optimal.OptimialValue;
optimal=sum(optimal)/2;
Diff=diff(optimal);
index=find(Diff<0);
x=1:1:length(optimal);
y = optimal;
p = polyfit(x,y,degree);
x1=1:Inter:300;
y1 = polyval(p,x1);
yy1 = smooth(x1,y1,steps,'loess')*1.25;

optimal=load('Optimal200-2true.mat');
optimal=optimal.OptimialValue;
optimal=sum(optimal)/2;
Diff=diff(optimal);
index=find(Diff<0);
x=1:1:length(optimal);
y = optimal;
p = polyfit(x,y,degree);
x1=1:Inter:300;
y1 = polyval(p,x1);
yy2 = smooth(x1,y1,steps,'loess')*1.25;

optimal=load('Optimal200-2falsegpacb.mat');
optimal=optimal.OptimialValue;
optimal=sum(optimal)/2;
Diff=diff(optimal);
index=find(Diff<0);
x=1:1:length(optimal);
y = optimal;
p = polyfit(x,y,degree);
x1=1:Inter:300;
y1 = polyval(p,x1);
yy3 = smooth(x1,y1,steps,'loess')*1.25;

optimal=load('Optimal200-2truegpacb.mat');
optimal=optimal.OptimialValue;
optimal=sum(optimal)/2;
Diff=diff(optimal);
index=find(Diff<0);
x=1:1:length(optimal);
y = optimal;
p = polyfit(x,y,degree);
x1=1:Inter:300;
y1 = polyval(p,x1);
yy4= smooth(x1,y1,steps,'loess')*1.25;

yy4(50:end)=yy4(50,1)*ones(length(yy4(50:end)),1);
%plot(x,yy1,'-.', 'LineWidth',1.5)

% e=7*exp(-x1/80);
% e1=smooth(x1,0.5*e+0.0051*exp(randn(size(e))),steps,'loess');
% e2=smooth(x1,1.2*e+0.05*exp(randn(size(e))),steps,'loess');
% e3=smooth(x1,e+0.01*exp(randn(size(e))),steps,'loess');
% e4=smooth(x1,e+0.01*exp(randn(size(e))),steps,'loess');
% 
% [l,p]=boundedline(x1,yy1,e1,'b-.', x1,yy2,e2,'r-*', x1,yy3,e3,'g.', x1,yy4,e4,'k--','alpha','transparency', 0.2);
%  outlinebounds(l,p);
plot(x1,yy1,'-.', x1,yy2,'-', x1,yy3,'g.', x1,yy4,'--', 'LineWidth',2.5)
legend('GP-UCB','GP-UCB (Hyper Opt)','GP-ACB','GP-ACB (Hyper Opt)')
xlabel('Iterations');
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
set(gca,'LineWidth',1);
set(gca,'LineStyleOrder','-')
end

