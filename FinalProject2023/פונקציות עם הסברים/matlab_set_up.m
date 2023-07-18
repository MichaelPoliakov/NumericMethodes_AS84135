%% General Set Up
clc; close all; 
format shortG; format compact;

set(0,'defaultAxesXAxisLocation','origin','defaultAxesYAxisLocation','origin')
set(0,'defaultAxesFontWeight','bold')
set(groot,'defaultAxesLabelFontSizeMultiplier',1.5);

% Tick 
set(groot,'defaultAxesTickDir','both','defaultAxesTickLength',4e-3.*[1 1])
set(groot,'defaultAxesTickDirMode','manual')
set(groot,'defaultAxesXMinorTick','on','defaultAxesYMinorTick','on')

set(groot,'defaultAxesTitleFontSizeMultiplier',1.5);


set(groot, 'defaultLineLineWidth', 3,'defaultLineMarkerSize',8);
