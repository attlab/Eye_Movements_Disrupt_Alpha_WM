%{
EYE_Upsample_250Hz_Files
Author: Tom Bullock
Date: 05.26.22

%}

clear
close all

load('s8_0112.mat')

edfStruct.FSAMPLE.time_test = repelem(edfStruct.FSAMPLE.time,1,4);
edfStruct.FSAMPLE.pa_test = repelem(edfStruct.FSAMPLE.pa,1,4);
edfStruct.FSAMPLE.gx_test = repelem(edfStruct.FSAMPLE.gx,1,4);
edfStruct.FSAMPLE.gy_test = repelem(edfStruct.FSAMPLE.gx,1,4);