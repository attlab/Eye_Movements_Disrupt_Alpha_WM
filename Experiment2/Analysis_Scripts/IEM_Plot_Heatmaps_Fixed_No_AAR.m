%{
Plot_IEM_Within_HEAT
Author: Tom Bullock
Date: 10.23.20

Plot fixed model CRFs as surf plots.

%}

clear
close all

% set dirs
sourceDir = '/home/waldrop/Desktop/SCAMS/IEM_Results_TT_Within_Fixed_No_AAR'; 
plotDir = '/home/waldrop/Desktop/SCAMS/Plots';


% select subjects
subjects = [1,2,4:17,20:24,27:28]; 

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
xNewTickLabel = ['-0.5';'  0 ';' 0.5';'  1 ';' 1.5';'  2 '];
xNewLim = [1 640];
yNewTick = [1 2 3 4 5 6 7 8];
zNewTick = [.2 .4 .6 .8 1];
zNewTickLabel = [];%zNewTick;
yNewTickLabel = ['-180';'    ';'-90 ';'    ';' 0  ';'    ';' 90 ';'    '];
yNewLim = [1 8];
colorBar = [0 .8];
thisAxisFontsize = 16;

for iCond=1:4
    
    h=figure('units','normalized','outerposition',[0    0.5067    0.4850    0.4467]);

    if      iCond==1; mc='Spatial'; ec='Fix';
    elseif  iCond==2; mc='Color';   ec='Fix';
    elseif  iCond==3; mc='Spatial'; ec='Move';
    elseif  iCond==4; mc='Color';   ec='Move';
    end
    
    %subplot(2,2,iCond);
    %surf(squeeze(mean(all_real_total(:,iCond,f,:,:),1)),'linestyle','none','FaceAlpha',1); 
    %surf(squeeze(mean(allData(:,iCond,:,:),1)),'linestyle','none','FaceAlpha',1); 
    
    imagesc(squeeze(mean(allData(:,iCond,:,:),1))');  
    pcolor(squeeze(mean(allData(:,iCond,:,:),1))'); 
    shading interp

    
    %xlabel('Time(s)','FontSize',14)
    %ylabel('Chan. Resp.(uV2)','FontSize',14)
    %title([mc ' / ' ec],'FontSize',12)
    %view(70,30)
    %zlim([0 1])
    
    colormap('jet')
    caxis([0,.8]);
    pbaspect([2 1 1])
    
    set(gca,'xlim',[1 640],...
        'XTick',xNewTick,...
        'xTickLabel',xNewTickLabel,...
        'YTick',1:8,...
        'YTickLabel',yNewTickLabel,...
        'FontSize',40,...
        'LineWidth',2,...
        'layer','top')
    
    
    % stim onset line
   line([128,128],[.5,8.5],...
       'linewidth',5,...
       'color','w',...
       'linestyle','--'); 
   
   %text(128,8.2,'On','fontsize',17);
   
       % stim offset line
   line([192,192],[.5,8.5],...
       'linewidth',6,...
       'color','w',...
       'linestyle','--'); 
   
   %text(192,8.2,'Off','fontsize',17);

   if ismember(iCond,[3,4])
       
       for iEyeLine=1:4
           
           if iEyeLine==1; xLoc = 256; thisColor = 'w'; thisText = 'Cue'; % fix-dot moves
           %elseif iEyeLine==2; xLoc = 350;  thisColor = 'r'; thisText = 'Eyes lim'; %256+(256/2); % eyes must be on fix dot by this point
           elseif iEyeLine==3; xLoc = 512;  thisColor = 'w'; thisText = 'RF'; % fix-dot moves back
           %elseif iEyeLine==4; xLoc = 627;  thisColor = 'r'; thisText = 'Eyes lim'; % eyes must be back to center by this point
           end
           
           % stim offset line
           line([xLoc,xLoc],[.5,8.5],...
               'linewidth',6,...
               'color',thisColor,...
               'linestyle',':');
           
           % text on top of line
           %text(xLoc,8.2,thisText,'FontSize',17);
           
       end
       
   end
    
    
%     set(gca,'yTick',xNewTick,'yticklabel',xNewTickLabel,'xTick',yNewTick,...
%         'xticklabel',yNewTickLabel,'ylim',xNewLim,'fontsize',thisAxisFontsize,...
%         'LineWidth',4,'ztick',zNewTick,'zticklabel',zNewTickLabel);
%     grid('off')
%     
    saveas(h,[plotDir '/' 'IEM_Heat_Fixed_No_AAR_' mc '_' ec '.eps'],'epsc')
    saveas(h,[plotDir '/' 'IEM_Heat_Fixed_No_AAR_' mc '_' ec '.jpg'],'jpeg')
end