%{
Author: Tom Bullock
Date: 11.14.19

List bad channels for specific subjects
%}

function badChannels = badChannelInfo(sjNum)

if sjNum==1
    badChannels = {'FT7','T7','F7','F5','FT8','T8','TP8','P6','P8','F6'};
elseif sjNum==2
    badChannels = {'CP6','P4'}; 
elseif sjNum==3
    badChannels = {'T8','TP8','C1','F1','P5','P2','P1','P3','P4','CP2','F8','C3','O1','O2','P9','T7'}; % participant excluded due to excessive chan rej
elseif sjNum==4
    badChannels = {'T7','FC4','P6','T8','F4','P10','P9','FT7','AF8','CP6','F5','F7','Fp2','F8','Iz','C3'}; %
elseif sjNum==5
    badChannels = {'T8','TP8','P6','P8','F7','FC3'};
elseif sjNum==6
    badChannels = {'P10','TP7','P10','P9','AF8','AF7','T7','F8'};
elseif sjNum==7
    badChannels = {'P10','F6','P9','P8','AF7','F3','AF8','F7','Iz','FP2','F8'};
elseif sjNum==8
    badChannels = {'F8','FC6','FP2'};
elseif sjNum==9
    badChannels = {'FT8','F5','C6','FP2','AF8','P3','P10','P2','IZ','Fp1','PO4'};
elseif sjNum==10
    badChannels = {'FT8','F8','F6','FC6','C6','T8','AF7','FP1','FC4','F5','O1','P3'};
elseif sjNum==11
    badChannels = {'T7','CP3','C5','CP5','FT7','AF7','F7'};
elseif sjNum==12
    badChannels = {'CP5','FC5','P8','F8','F6','F7'};
elseif sjNum==13
    badChannels = {'FP2','AF8','FT7','T7','PO4'};
elseif sjNum==14
    badChannels = {};
elseif sjNum==15
    badChannels = {'F8','F6','TP7','T7','AF7','PO8','P5','F5','P9','FT7','F4'};
elseif sjNum==16
    badChannels = {'F4','AF8'};
elseif sjNum==17
    badChannels = {'P9'};
elseif sjNum==21
    badChannels = {'TP7','CP2'};
elseif sjNum==22
    badChannels = {'T8'};
elseif sjNum==23
    badChannels = {'T8','F5','P7'};
elseif sjNum==24
    badChannels = {'TP7','F8','FT8'};
elseif sjNum==27
    badChannels = {'F5','FT8'};
elseif sjNum==28
    badChannels = {'TP7','C5'};    
else
    badChannels = {};
end


return
