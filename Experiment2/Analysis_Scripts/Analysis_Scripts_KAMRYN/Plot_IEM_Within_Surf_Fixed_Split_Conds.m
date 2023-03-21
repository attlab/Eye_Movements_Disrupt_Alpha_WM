%{
Plot_IEM_Within_Surf_Fixed
Author: Tom Bullock
Date: 10.23.20

Plot fixed model CRFs as surf plots.

%}

clear
close all

% set dirs
sourceDir = '/home/waldrop/Desktop/WTF_EYE/IEM_Results_TT_Within_Fixed'; 
plotDir = '/home/waldrop/Desktop/WTF_EYE/Plots';


% select subjects
subjects = [1:7,9,14,16:20,22:27,31]; 

for iSub=1:length(subjects)
    sjNum = subjects(iSub);
    
   load([sourceDir '/' sprintf('sj%02d_fixed_IEM.mat',sjNum)])
   
   allData(iSub,1,:,:) = em.tfs_cond1.total;
   allData(iSub,2,:,:) = em.tfs_cond2.total;
   allData(iSub,3,:,:) = em.tfs_cond3.total;
   allData(iSub,4,:,:) = em.tfs_cond4.total;
   
   %allData(iSub,:,:,:) = squeeze(mean(em.allTF,3)); % avg across iterations, only 1 freq (alpha) currently (allTF = cond x freq x iter x sample x chan)
   
end


% plot data
% set up plot
thisSR=256;
xNewTick = [1, thisSR*.5, thisSR*1,  thisSR*1.5, thisSR*2, thisSR*2.5];
xNewTickLabel = [];%['-0.5';'  0 ';' 0.5';'  1 ';' 1.5';'  2 '];
xNewLim = [1 640];
yNewTick = [1 2 3 4 5 6 7 8];
zNewTick = [.2 .4 .6 .8 1];
zNewTickLabel = [];%zNewTick;
yNewTickLabel = [];%['-180';'    ';'-90 ';'    ';' 0  ';'    ';' 90 ';'    '];
yNewLim = [1 8];
colorBar = [0 .8];
thisAxisFontsize = 20;

for iCond=1:4
    
    h=figure('units','normalized','outerposition',[0 0 .65 1]);

    if      iCond==1; mc='Spatial'; ec='Fix';
    elseif  iCond==2; mc='Color';   ec='Fix';
    elseif  iCond==3; mc='Spatial'; ec='Move';
    elseif  iCond==4; mc='Color';   ec='Move';
    end
    
    %subplot(2,2,iCond);
    %surf(squeeze(mean(all_real_total(:,iCond,f,:,:),1)),'linestyle','none','FaceAlpha',1); 
    surf(squeeze(mean(allData(:,iCond,:,:),1)),'linestyle','none','FaceAlpha',1); 

    
    %ylabel('Time(s)')
    %zlabel('Chan. Resp.(uV2)')
    %title([mc ' / ' ec])
    view(70,30)
    zlim([0 1])
    colormap('jet')
    caxis([.2,1]);
    pbaspect([1 2 1])
    set(gca,'yTick',xNewTick,'yticklabel',xNewTickLabel,'xTick',yNewTick,...
        'xticklabel',yNewTickLabel,'ylim',xNewLim,'fontsize',thisAxisFontsize,...
        'LineWidth',4,'ztick',zNewTick,'zticklabel',zNewTickLabel);
    grid('off')
    
    saveas(h,[plotDir '/' 'IEM_Surf_Fixed_' mc '_' ec '.eps'],'epsc')
    saveas(h,[plotDir '/' 'IEM_Surf_Fixed_' mc '_' ec '.jpg'],'jpeg')
end












% h=figure('units','normalized','outerposition',[0 0 1 1]);
% 
% for iPlot=1:4
%     
%     subplot(2,2,iPlot)
%     
%     surf(squeeze(mean(allData(:,iPlot,:,:),1)))
%     
%     set(gca,...
%         'ytick',linspace(1,640,6),...
%         'yticklabel',-500:500:2000);
%     
%     title(['Cond' num2str(iPlot)])
%     
% end
