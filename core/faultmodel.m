%% Define the parameters for the bearing.

n = 8;         % Number of rolling element bearings
d = 0.02;      % Diameter of rolling elements 
p = 0.12;      % Pitch diameter of bearing
thetaDeg = 15; % Contact angle in degrees
fPin =22.5;
fs = 20E3;          % Sample Rate (Hz)

%% The impacts occur at the ball pass frequency-inner race (BPFI);
bpfi = n*fPin/2*(1 + d/p*cosd(thetaDeg));
t = 0:1/fs:0.05-1/fs;
fImpact = 3000;
noise =  randn(size(t))/20+sin(4*pi*fImpact*t)/20;


tImpact = 0:1/fs:5e-3-1/fs;
xImpact = sin(2*pi*fImpact*tImpact).*kaiser(length(tImpact),40)';

xComb = zeros(size(t));
xComb(1:round(fs/bpfi):end) = 1;
xBper_bpfi = 0.33*conv(xComb,xImpact,'same')+noise;

%%  Outer-race defects

xComb = zeros(size(t));
bpfo =  n*fPin/2*(1 - d/p*cosd(thetaDeg));
xComb(1:round(fs/bpfo):end) = 1;
xBper_bpfo = 0.33*conv(xComb,xImpact,'same')+noise;

%% all-Passing Frequency Roller ~BPFR!
bpfr=  p/d*fPin*(1 - d/p*cosd(thetaDeg)^2);
xComb = zeros(size(t));
xComb(1:round(fs/bpfr):end) = 1;
xBper_bpfr = 0.33*conv(xComb,xImpact,'same')+noise;

figure
subplot(2,2,1)
plot(t,xBper_bpfi)
xlim([0 0.05])
ylim([-0.4,0.4])
xlabel('Time (s)')
ylabel('Acceleration')
subplot(2,2,2)
plot(t,xBper_bpfo)
xlim([0 0.05])
ylim([-0.4,0.4])
xlabel('Time (s)')
ylabel('Acceleration')
subplot(2,2,3)
plot(t,xBper_bpfr)
xlim([0 0.05])
ylim([-0.4,0.4])
xlabel('Time (s)')
ylabel('Acceleration')
subplot(2,2,4)
plot(t,noise)
xlim([0 0.05])
ylim([-0.4,0.4])
xlabel('Time (s)')
ylabel('Acceleration')

set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [3 3 8 6]);
set(gcf, 'Alphamap',0.01);
set(gcf, 'Colormap', cool);
set(gcf,'Units', 'inches');
set(gcf,'Position',[2, 2, 6, 6]);
set(gcf,'OuterPosition',[1,1,7,7])
set(gcf,'Color','white')
% set(gca,'LineWidth',1)
% set(gca,'FontUnits','points')
% set(gca,'FontSize',20);
% set(gca,'Color','none');
% set(gca,'Box','off');
set(findall(gcf,'-property','FontSize'),'FontSize',16)

