%{
EEG_BEH_Relationship
Author: Tom Bullock
Date: 10.14.22

%}

clear
close all

% set dir 
%sourceDir = '/home/waldrop/Desktop/WTF_EYE/Data_Compiled'; % for SCAM (EXP1)
sourceDir = '/home/waldrop/Desktop/SCAMS/Data_Compiled'; % for SCAMS (EXP2)

% load beh data
load([sourceDir '/' 'Modelling_Data.mat'])

% calculate delta between spatial fix and spatial move precision AND color
% fix and color move precision
spatial_delta = modelSD(:,1) - modelSD(:,3);
color_delta = modelSD(:,2) - modelSD(:,4);


% load slope data (sjs x conds x freqs [just use 1] x times)
load([sourceDir '/' 'IEM_Slopes_Within_Fixed.mat'])

% find mean slope delta 0.5-1.5 s (i.e. during the eye-movement window)
% between the fixation and move conditions for spaital and color

spatial_delta_slope = squeeze(mean(rSl_total(:,1,1,258:513) - rSl_total(:,3,1,258:513),4));
color_delta_slope = squeeze(mean(rSl_total(:,2,1,258:513) - rSl_total(:,4,1,258:513),4));


[r_spatial,p_spatial] = corr(spatial_delta, spatial_delta_slope);
[r_color, p_color] = corr(color_delta, color_delta_slope);


