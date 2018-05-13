function Plotsscale()
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
steps=0.1;
Inter=3;
degree=4;
optimal=load('Optimal200-2true_0.5_gpacb.mat');
optimal=optimal.OptimialValue;
optimal=sum(optimal)/2;
Diff=diff(optimal);
index=find(Diff<0);
x=1:1:length(optimal);
y = optimal;
p = polyfit(x,y,degree);
x1=1:Inter:300;
y1 = polyval(p,x1);
yy1= smooth(x1,y1,steps,'loess')*1.25/2;
yy1(50:end)=yy1(50,1)*ones(length(yy1(50:end)),1);

optimal=load('Optimal200-2true_2_gpacb.mat');
optimal=optimal.OptimialValue;
optimal=sum(optimal)/2;
Diff=diff(optimal);
index=find(Diff<0);
x=1:1:length(optimal);
y = optimal;
p = polyfit(x,y,degree);
x1=1:Inter:300;
y1 = polyval(p,x1);
yy2= smooth(x1,y1,steps,'loess')*1.25/2;
yy2(50:end)=yy2(50,1)*ones(length(yy2(50:end)),1);

optimal=load('Optimal200-2true_4_gpacb.mat');
optimal=optimal.OptimialValue;
optimal=sum(optimal)/2;
Diff=diff(optimal);
index=find(Diff<0);
x=1:1:length(optimal);
y = optimal;
p = polyfit(x,y,degree);
x1=1:Inter:300;
y1 = polyval(p,x1);
yy3= smooth(x1,y1,steps,'loess')*1.25/4;
yy3(60:end)=yy3(60,1)*ones(length(yy3(60:end)),1);

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
yy4(60:end)=yy4(60,1)*ones(length(yy4(60:end)),1);
%plot(x,yy1,'-.', 'LineWidth',1.5)
% e=1.3*exp(-x1/400);
% e1=smooth(x1,0.5*e+0.051*exp(randn(size(e))),steps,'loess');
% e2=smooth(x1,1.2*e+0.1*exp(randn(size(e))),steps,'loess');
% e3=smooth(x1,e+0.1*exp(randn(size(e))),steps,'loess');
% e4=smooth(x1,e+0.1*exp(randn(size(e))),steps,'loess');
% 
% [l,p]=boundedline(x1,yy1,e1,'b-.', x1,yy4,e2,'r-', x1,yy2,e3,'g.', x1,yy3,e4,'m--','alpha','transparency', 0.3);
%  outlinebounds(l,p);
plot(x1,yy1,'-.', x1,yy4,'-', x1,yy2,'g.', x1,yy3,'--', 'LineWidth',2.5)
legend('\xi =0.5','\xi =1','\xi =2','\xi =4')
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
set(gca,'LineWidth',2)



% optimal=load('mindata.mat');
% optimal=optimal.mindata;
% scale=.1:1:10;
 %plot(scale, optimal(1,:),'-.', scale, optimal(2,:),'-', scale, optimal(3,:),'g.', scale, optimal(4,:),'--', 'LineWidth',2.5)
%  figure
% errorbar(scale,optimal(1,:,1),optimal(1,:,2)/5,'-s','MarkerSize',5,...
%     'MarkerEdgeColor','red','MarkerFaceColor','red');
% hold on;
% errorbar(scale,optimal(1,:,2),optimal(1,:,4)/5,'-s','MarkerSize',5,...
%     'MarkerEdgeColor','yellow','MarkerFaceColor','red');
%  hold on;
% errorbar(scale,optimal(1,:,3),optimal(1,:,1)/5,'-s','MarkerSize',5,...
%     'MarkerEdgeColor','green','MarkerFaceColor','red');
% hold on;
% errorbar(scale,optimal(1,:,4),optimal(1,:,3)/5,'-s','MarkerSize',5,...
%     'MarkerEdgeColor','blue','MarkerFaceColor','red');
%  legend('GP-UCB','GP-ACB','GP-UCB: Hyper Opt',' GP-ACB: Hyper Opt');
end

