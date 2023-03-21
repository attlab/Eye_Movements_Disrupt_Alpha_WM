%{
EEG_Merge_Broken_Files
Author: Tom Bullock
Date: 11.21.19
Updated: 11.11.21

Notes: saves a merged .mat file
%}

clear

% add EEGLAB to path
cd('/home/waldrop/Desktop/WTF_EYE/Dependancies/eeglab14_1_2b')
eeglab
close all
cd ..

% set dirs
parentDir = '/home/waldrop/Desktop/SCAMS';
sourceDir = [parentDir '/EEG_Raw'];

sjNum=27;
session=1;

% set file parts [file names may vary depending on what they were saved as]
filePart1 = sprintf('sj%02d_se%02d_SCAMS.bdf',sjNum,session); 
filePart2 = sprintf('sj%02d_se%02d_Part2_SCAMS.bdf',sjNum,session); 
%filePart3 = sprintf('sj%02d_se%02d_part3_SCAMS.bdf',sjNum,session); 


% load parts
EEG1=pop_biosig([sourceDir '/' filePart1]);
EEG2=pop_biosig([sourceDir '/' filePart2]);
%EEG3=pop_biosig([sourceDir '/' filePart3]);

% special case - sj11se01part1 - remove events from 1668 onwards because
% crashy crashy (manually found latency of first block of second cond)
%EEG1 = pop_select(EEG1,'nopoint',(1572017-1):1817290);

% % % special case three files merged
% filePart3 = sprintf('sj%02d_se%02d_wtfEye_2.bdf',sjNum,session);
% EEG3=pop_biosig([sourceDir '\\' filePart3]);
% % EEG = pop_mergeset(EEG,EEG3);

% remove some events from sj27se01part1 coz crash
EEG1 = pop_select(EEG1,'point',[1,3066085-1]);

% merge
EEG = pop_mergeset(EEG1,EEG2);
%EEG = pop_mergeset(EEG,EEG3);

% save
save([sourceDir '/' sprintf('sj%02d_se%02d_SCAMS.mat',sjNum,session)],'EEG','-v7.3')