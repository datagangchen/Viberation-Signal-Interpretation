function [wy_OS,wy_RT,wy_PT,yawGain,alat_OS,alat_RT,alat_PT]=FeatureExtract(time,wy,alat);
%% wy
% values for yaw-related metrics
% for a certain initial setup,
t_reference = 1.22; %s, 50% value time
t = time;
% yaw overshoot
wy_OS = (max(wy) - wy(end))/wy(end);
t90wy = t(find(wy < 0.9*wy(end),1));
wy_RT = t90wy - t_reference;
[v_wy, i_wy] = max(wy);
tPTwy = t(i_wy);
wy_PT = tPTwy - t_reference;
yawGain = wy(end)/45.8*pi/180;

% latacc overshoot
alat_OS = (max(alat) - alat(end))/alat(end);
t90alat = t(find(alat < 0.9*alat(end),1));
alat_RT = t90alat - t_reference;
[v_alat, i_alat] = max(alat);
tPTalat = t(i_alat);
alat_PT = tPTalat - t_reference;
end