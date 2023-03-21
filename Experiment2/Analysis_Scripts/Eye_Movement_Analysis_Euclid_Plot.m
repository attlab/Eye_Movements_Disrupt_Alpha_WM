%{
Eye_Movement_Analysis_Plot
Author: Tom Bullock
Date: 06.05.20

plot euclidian errors for conds 3 and 4 for ems 1 and 2

SOMETHING IS UP WITH THIS SCRIPTO

%}

clear
close all

sourceDir = '/home/waldrop/Desktop/SCAMS/EYE_Euclidian_Error_Mats';
destDir = '/home/waldrop/Desktop/SCAMS/Data_Compiled';
destDirPlot = '/home/waldrop/Desktop/SCAMS/Plots';

% select subjects (multiple)
subjects = [1:2,4:17,21:24,27:28]; % kicked out subjects 03 (bad EEG) and 19 (bad BEH) AND 20 coz missing eye-tracking file

h=figure('units','normalized','OuterPosition',[0,0,1,1]);
for iSub=1:length(subjects)
    
    sjNum = subjects(iSub);
    
    % load data
    load([sourceDir '/' sprintf('sj%02d_euclid_error.mat',sjNum)])
    
    % get means for each em and each cond
    cond3(iSub,1) = nanmean([euclidErrorStruct(3).em1]);
    %cond3(iSub,2) = nanmean([euclidErrorStruct(3).em2]);
    
    cond4(iSub,1) = nanmean([euclidErrorStruct(4).em1]);
    %cond4(iSub,2) = nanmean([euclidErrorStruct(4).em2]);
    
    % get times
    cond3_idx(iSub,1) = nanmean([euclidErrorStruct(3).em1_idx]);
    %cond3_idx(iSub,2) = nanmean([euclidErrorStruct(3).em2_idx]);
    
    cond4_idx(iSub,1) = nanmean([euclidErrorStruct(4).em1_idx]);
    %cond4_idx(iSub,2) = nanmean([euclidErrorStruct(4).em2_idx]);
    
end

cond3_visAngle = visAngleCalculate_pix(cond3,120);
cond4_visAngle = visAngleCalculate_pix(cond4,120);






% convert lat index to wide format for R
for iData=1:2
    
    % get data 
    if iData==1
        theseData = [cond3_visAngle,cond4_visAngle]; % euclid error
    else
        theseData = [cond3_idx,cond4_idx]; % euclid error time
    end
    
    nSubs = size(theseData,1);

    % convert into R "long" format
%     dumMem = [zeros(nSubs,1), ones(nSubs,1), zeros(nSubs,1), ones(nSubs,1)];
%     dumEyes = [zeros(nSubs,1), zeros(nSubs,1), ones(nSubs,1), ones(nSubs,1)];
    dumMem = [zeros(nSubs,1), ones(nSubs,1)];
    dumEyes = [ones(nSubs,1), ones(nSubs,1)];
    
    % create column vector of subjects
    sjNums = 1:nSubs;
    %sjNums = [sjNums, sjNums, sjNums, sjNums]';
    sjNums = [sjNums, sjNums]';
    
    % use colon operator function to get into long format
    if iData==1
        euclid_error_LF = [theseData(:), sjNums, dumMem(:), dumEyes(:)];
    elseif iData==2
        euclid_time_LF = [theseData(:), sjNums, dumMem(:), dumEyes(:)];
    end
    
end


save([destDir '/' 'EYE_Euclid_Error.mat'],'cond3_visAngle','cond4_visAngle','cond3_idx','cond4_idx','euclid_error_LF','euclid_time_LF','-v7.3')





% generate plots for dist
subplot(1,2,1);
for iPlot=1:2
    
    if iPlot==1; theseData = cond3; thisTitle = 'Spatial Move'; thisColor = [0,204,0];thisX=1;
    elseif iPlot==2; theseData = cond4; thisTitle = 'Color Move'; thisColor = [204,0,204]; thisX=2;
    end
    
    % compute mean and SEM in pixels
    thisMean = mean(theseData);
    thisSEM = std(theseData,0,1)./sqrt(size(theseData,1));
    
    % convert to visual angles
    thisMean = visAngleCalculate_pix(thisMean,120);
    thisSEM = visAngleCalculate_pix(thisSEM,120);
    
    %subplot(1,2,iPlot);
    bar(thisX,thisMean,'FaceColor',thisColor./255); hold on
    

    
    errorbar(thisX+.25,thisMean,thisSEM,...
        'LineStyle','none',...
        'LineWidth',4,...
        'color','k',...
        'CapSize',0)
    
    set(gca,...
        'box','off',...
        'linewidth',1.5,...
        'xlim',[.5,2.5],...
        'fontsize',36,...
        'xTick',[],...
        'XTickLabel',{' ',' ',' ',' '})
    
    pbaspect([1,1,1]);
    
    %title(thisTitle);
        
    
end

%plot individual data points using plotSpread
theseData = [cond3,cond4];
theseData = visAngleCalculate_pix(theseData,120);
plotSpread(theseData,'distributionMarkers',{'.'},'distributionColors',{'k'});
set(findall(1,'type','line','color','k'),'markerSize',24) %Change marker size




% generate plots for time (sample)
subplot(1,2,2);
for iPlot=1:2
    
    if iPlot==1; theseData = cond3_idx; thisTitle = 'Spatial Move'; thisColor = [0,204,0];thisX=1;
    elseif iPlot==2; theseData = cond4_idx; thisTitle = 'Color Move'; thisColor = [204,0,204]; thisX=2;
    end
    
    % compute mean and SEM in pixels
    thisMean = mean(theseData);
    thisSEM = std(theseData,0,1)./sqrt(size(theseData,1));
    
    % convert from samples (500 Hz) to milliseconds
    thisMean = thisMean*2;
    
%     % convert to visual angles
%     thisMean = visAngleCalculate_pix(thisMean,120);
%     thisSEM = visAngleCalculate_pix(thisSEM,120);
    
    %subplot(1,2,iPlot);
    bar(thisX,thisMean,'FaceColor',thisColor./255); hold on
    
    errorbar(thisX+.25,thisMean,thisSEM,...
        'LineStyle','none',...
        'LineWidth',4,...
        'color','k',...
        'CapSize',0)
    
    set(gca,...
        'box','off',...
        'linewidth',1.5,...
        'xlim',[.5,2.5],...
        'ylim',[0,120],...
        'fontsize',36,...
        'xTick',[],...
        'XTickLabel',{' ',' ',' ',' '});
    
    
    
    pbaspect([1,1,1]);
    
    %title(thisTitle);
        
    
end

%plot individual data points using plotSpread
theseData = [cond3_idx,cond4_idx];
theseData = theseData*2;
%theseData = visAngleCalculate_pix(theseData,120);
plotSpread(theseData,'distributionMarkers',{'.'},'distributionColors',{'k'});
set(findall(1,'type','line','color','k'),'markerSize',24) %Change marker size

saveas(h,[destDirPlot '/' 'Eye_Error_Bar_Plot.eps'],'epsc')



