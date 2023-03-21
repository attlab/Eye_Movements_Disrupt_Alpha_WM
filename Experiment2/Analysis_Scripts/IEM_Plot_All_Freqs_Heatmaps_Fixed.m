%{
IEM_Plot_Cross_TT_BF_Thresholded (for SCAM)
Author: Tom Bullock
Date: 04.16.20

Load in the slope data from the crossTT analyses and either apply a t-test
based threshold (p<.05) or a bayes factor threshold (BF>xx) or no threshold

%}

clear 
close all

% select thresholding type (0=none, 1=p<.05,2=BF)
thresholdType=2;

% set dirs
sourceDir = '/home/waldrop/Desktop/SCAMS';
slopeDataDir = [sourceDir '/Data_Compiled'];
%bfTestDir = '/home/waldrop/Desktop/WTF_EYE/R_Results';
bfTestDir = [sourceDir '/R_Results'];
%plotDir = '/home/waldrop/Desktop/WTF_EYE/Plots';
plotDir = [sourceDir '/Plots'];

% load matrices of real and permuted slopes from cross TT analysis (from
% matlab cluster)
load([slopeDataDir '/' 'IEM_Slopes_All_Freqs_Fixed.mat'])

% average across subs
rSl = squeeze(mean(rSl_total,1));


% loop through real slope mat and zero out any p<.05 (frequentist)
if thresholdType==1
    for iCond=1:size(rSl_total,2)
        for iFreq=1:size(rSl_total,3)
            for iTime=1:size(rSl_total,4)
                
                if ttest(rSl_total(:,iCond,iFreq,iTime),pSl_total(:,iCond,iFreq,iTime))
                    rSl_Adj(iCond,iFreq,iTime) = rSl(iCond,iFreq,iTime);
                else
                    rSl_Adj(iCond,iFreq,iTime) = 0;
                end
                
            end
        end
    end
end

% loop through real slope mat and zero out any BF<3 (bayes)
if thresholdType==2
    load([bfTestDir '/' 'BF_Results_IEM_All_Freqs_RvP_Fixed.mat'])
    for iCond=1:size(rSl_total,2)
        for iFreq=1:size(rSl_total,3)
            for iTime=1:size(rSl_total,4)
                if bfGrandOutput(iTime,iCond,iFreq)>3
                    rSl_Adj(iCond,iFreq,iTime) = rSl(iCond,iFreq,iTime);
                else
                    rSl_Adj(iCond,iFreq,iTime) = 0;
                end
            end
        end
    end
end



% loop through and plot heatmaps for all freqs results


for iPlot=1:4
    
    h=figure('Units','normalized','Position',[0,0,.25,.2]);
   
    %subplot(2,2,iPlot);
    
    ih = imagesc(squeeze(rSl_Adj(iPlot,1:30,:)),[0,.004]);
    
    set(gca,'xtick',[0.1,8,16,24,32,40],...
        'xTickLabel',{'-0.5',' 0  ', ' 0.5',' 1  ',' 1.5',' 2  '},...
        'ytick',[1,10,20,30],...
        'YTickLabel',[0,10,20,30],...
        'fontsize',24,...
        'ydir','normal',...
        'LineWidth',1.5,...
        'box','off')

   % lines for events
    vline(128/16,'w--') % targ on
    vline(192/16,'w--') % targ off
    if iPlot>2
        vline(256/16,'w:') % EM1
        %vline(384/16,'w:') % EM2
        vline(512/16,'w:') % EM3
    end
   
%     set(gca,...
%         'ytick',[1,10,20,30],...
%         'YTickLabel',[0,10,20,30],...
%         'fontsize',24,...
%         'ydir','normal',...
%         'LineWidth',1.5,...
%         'box','off')
    
    pbaspect([2,1,1])
    colormap('jet')
    %colorbar

    saveas(h,[plotDir '/' 'IEM_All_Freqs_Cond' num2str(iPlot) '.eps'],'epsc')
    
end

%'xTickLabel',{'-0.5',' 0  ', ' 0.5',' 1  ',' 1.5',' 2  '},...








% 
% 
% 
% 
% 
% for iCond=1:size(all_rSl,2)
%     for iTr=1:size(all_rSl,3)
%         for iTe=1:size(all_rSl,4)
%             if ttest(all_rSl(:,iCond,iTr,iTe),all_pSl(:,iCond,iTr,iTe))==1
%                 rSl_Adj(iCond,iTr,iTe) = rSl(iCond,iTr,iTe);
%             else
%                 rSl_Adj(iCond,iTr,iTe) = 0;
%             end
%         end
%     end
% end
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% % average real and perm slope data over subjects
% rSl = squeeze(mean(all_rSl,1));
% pSl = squeeze(mean(all_pSl,1));
% 
% % load matrix of BFs (outputted from R script)
% load([bfTestDir '/' 'Cross_TT_Slope_BF_TTESTS.mat'])
% 
% % reorganize bfSummary mat to be same dims as slopes
% bfSummary = permute(bfSummary, [3,2,1]);
% 
% % select threshold
% if thresholdType==0
%     % no thresholding
%     rSl_Adj = rSl;
% elseif thresholdType==1
%     % loop through real slope mat and zero out any p<.05 (frequentist)
%     for iCond=1:size(all_rSl,2)
%         for iTr=1:size(all_rSl,3)
%             for iTe=1:size(all_rSl,4)
%                 if ttest(all_rSl(:,iCond,iTr,iTe),all_pSl(:,iCond,iTr,iTe))==1
%                     rSl_Adj(iCond,iTr,iTe) = rSl(iCond,iTr,iTe);
%                 else
%                     rSl_Adj(iCond,iTr,iTe) = 0;
%                 end
%             end
%         end
%     end
% elseif thresholdType==2 
%     % loop through real slope mat and zero out any BF<3 (bayesian)
%     rSl_Adj = [];
%     for i=1:size(bfSummary,1)
%         for j=1:size(bfSummary,2)
%             for k=1:size(bfSummary,3)    
%                 if bfSummary(i,j,k)>3
%                     rSl_Adj(i,j,k) = rSl(i,j,k);
%                 else
%                     rSl_Adj(i,j,k) = 0;
%                 end               
%             end
%         end
%     end
% end
% 
% 
% % plot all train/test combos
% h=figure('units','normalized','outerposition',[0 0 1 1]);
% xNew=-500:2501/40:2000;
% yNew=xNew;
% cLims = [0 .004];
% 
% % plotting loop
% plotIdx = 0;
% % for k=[ 1,6,11,16,...
% %         2,5,9,13,...
% %         3,7,10,14,...
% %         4,8,12,15]
% for k=[1,6,11,16,3,9,8,14]
%     
%     plotIdx = plotIdx+1;
%     
%     % set up plot titles
%     if      k==1; trainCond=1; testCond=1;
%     elseif  k==2; trainCond=1; testCond=2;
%     elseif  k==3; trainCond=1; testCond=3;
%     elseif  k==4; trainCond=1; testCond=4;
%     elseif  k==5; trainCond=2; testCond=1;
%     elseif  k==6; trainCond=2; testCond=2;
%     elseif  k==7; trainCond=2; testCond=3;
%     elseif  k==8; trainCond=2; testCond=4;
%     elseif  k==9; trainCond=3; testCond=1;
%     elseif  k==10; trainCond=3; testCond=2;
%     elseif  k==11; trainCond=3; testCond=3;
%     elseif  k==12; trainCond=3; testCond=4;
%     elseif  k==13; trainCond=4; testCond=1;
%     elseif  k==14; trainCond=4; testCond=2;
%     elseif  k==15; trainCond=4; testCond=3;
%     elseif  k==16; trainCond=4; testCond=4;
%     end
%     
%     if      trainCond==1; trainLabel = 'S/F';
%     elseif  trainCond==2; trainLabel = 'C/F';
%     elseif  trainCond==3; trainLabel = 'S/M';
%     elseif  trainCond==4; trainLabel = 'C/M';
%     end
%     
%     if      testCond==1; testLabel = 'S/F';
%     elseif  testCond==2; testLabel = 'C/F';
%     elseif  testCond==3; testLabel = 'S/M';
%     elseif  testCond==4; testLabel = 'C/M';
%     end
%     
%     tl = ['Train: ' trainLabel ' --> ' 'Test: ' testLabel];
%     
%     % create subplots
%     subplot(2,4,plotIdx)
%     
%     % select condition for plotting
%     dataForPlot = squeeze(rSl_Adj(k,:,:));
%     
%     % generate plot
%     imagesc(xNew,yNew,dataForPlot,cLims);
%     
%     % plot settings
%     axisLims = [-500,2000];
%     axisTicks = -500:500:2000;
%     axisLabels = {'-.5',' 0 ',' .5',' 1 ','1.5',' 2 '};
%     
%     set(gca,...
%         'XLim',axisLims,...
%         'YLim',axisLims,...
%         'XTick',axisTicks,...
%         'YTick',axisTicks,...
%         'XTickLabel',axisLabels,...
%         'YTickLabel',axisLabels,...
%         'FontSize',20,...
%         'LineWidth',1.5)
%     
%     ylabel('Training (s)')
%     xlabel('Testing (s)')
%     pbaspect([1,1,1])
%     title(tl)
%     colormap('jet')
%     %colorbar % comment to remove c-bars for pub
%     
%     % target onset
%     vline(0,'w--')
%     hline(0,'w--')
%     
%     % target offset
%     vline(250,'w--')
%     hline(250,'w--')
%     
%     if ~ismember(k,[1,6])
%         % first eye movement (R or L)
%         vline(500,'w:')
%         hline(500,'w:')
%         
%         % second eye movement (R or L)
%         vline(1000,'w:')
%         hline(1000,'w:')
%         
%         % third eye movement (back to fix)
%         vline(1500,'w:')
%         hline(1500,'w:')
%     end
%     
%     clear dataForPlot
%     
% end
% 
% 
% 
% 





















% %%%%%%%%%%%%%%%%% ORIGINAL %%%%%%%%%%%%%%%%%%
% 
% % set up plot
% xNew=-500:2501/40:2001;
% yNew=xNew;
% cLim=[0 .005];
% 
% 
% 
% %% do all plots
% h=figure;
% for i=1:16
%     
%     subplot(4,4,i)
%     %pause(.5)
%     
% %    % sets up plot, gets correct stuff etc.
% %     if i==1; tl='Train EC, Test EC'; m=13;
% %     elseif i==2; tl='Train EC, Test EO'; m=1;
% %     elseif i==3; tl='Train EC, Test ECM'; m=2;
% %     elseif i==4; tl='Train EC, Test EOM'; m=3;
% %         
% %     elseif i==5; tl='Train EO, Test EO'; m=14;
% %     elseif i==6; tl='Train EC, Test EO'; m=4;
% %     elseif i==7; tl='Train EC, Test ECM'; m=5;
% %     elseif i==8; tl='Train EC, Test EOM'; m=6;
% %         
% %     elseif i==9; tl='Train ECM, Test ECM'; m=15;
% %     elseif i==10; tl='Train ECM, Test EC'; m=7;
% %     elseif i==11; tl='Train ECM, Test EO'; m=8;
% %     elseif i==12; tl='Train ECM, Test EOM'; m=9;   
% %       
% %     elseif i==13; tl='Train EOM, Test EOM'; m=16;
% %     elseif i==14; tl='Train EOM, Test EC'; m=10;
% %     elseif i==15; tl='Train EOM, Test EO'; m=11;
% %     elseif i==16; tl='Train EOM, Test ECM'; m=12;  
% %     end  
%     
%     
%     % sets up plot, gets correct stuff etc.
%     if i==1; tl='Train EC, Test EC'; m=13;
%     elseif i==2; tl='Train EC, Test EO'; m=1;
%     elseif i==3; tl='Train EC, Test ECM'; m=2;
%     elseif i==4; tl='Train EC, Test EOM'; m=3;
%         
%     elseif i==5; tl='Train EO, Test EO'; m=14;
%     elseif i==6; tl='Train EO, Test EC'; m=4; 
%     elseif i==7; tl='Train EO, Test ECM'; m=5;
%     elseif i==8; tl='Train EO, Test EOM'; m=6;
%         
%     elseif i==9; tl='Train ECM, Test ECM'; m=15;
%     elseif i==10; tl='Train ECM, Test EC'; m=7;
%     elseif i==11; tl='Train ECM, Test EO'; m=8;
%     elseif i==12; tl='Train ECM, Test EOM'; m=9;   
%       
%     elseif i==13; tl='Train EOM, Test EOM'; m=16;
%     elseif i==14; tl='Train EOM, Test EC'; m=10;
%     elseif i==15; tl='Train EOM, Test EO'; m=11;
%     elseif i==16; tl='Train EOM, Test ECM'; m=12;  
%     end  
%     
%     
%     
%         
%     imagesc(xNew,yNew,squeeze(rSl_Adj(m,:,:)),cLim)
%     pbaspect([1,1,1])
%     
%     L=vline(0);
%     set(L,'linewidth',2,'linestyle','--','color','w')
%     L=hline(0);
%     set(L,'linewidth',2,'linestyle','--','color','w')
%     L=vline(250);
%     set(L,'linewidth',2,'linestyle','--','color','w')
%     L=hline(250);
%     set(L,'linewidth',2,'linestyle','--','color','w')
%     
%     if i==3||i==4||i==7||i==8
%         L=vline(300);
%         set(L,'linewidth',2,'linestyle','--','color','w')
%         L=hline(300);
%         set(L,'linewidth',2,'linestyle','--','color','w')
%     end
%     
%     %cbar
%     colormap('jet')
%         %pbaspect([1 2 1])
%     set(gca,'fontsize',6, 'linewidth',2)
%     title(tl)
%     
%     %saveas(h,[plotDir '/' 'TTcross_combo_' num2str(i) '.eps',],'epsc')     
%     
% end
