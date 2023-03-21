%{
Plot_IEM_Within_HEAT
Author: Tom Bullock
Date: 10.23.20

Plot fixed model CRFs as surf plots.

%}

clear
close all

% set dirs
sourceDir = '/home/waldrop/Desktop/WTF_EYE/IEM_Results_TT_Within_Fixed_Color'; 
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
xNewTickLabel = ['-0.5';'  0 ';' 0.5';'  1 ';' 1.5';'  2 '];
xNewLim = [1 640];
yNewTick = [1 2 3 4 5 6 7 8];
zNewTick = [.2 .4 .6 .8 1];
zNewTickLabel = [];%zNewTick;
yNewTickLabel = ['-180';'    ';'-90 ';'    ';' 0  ';'    ';' 90 ';'    '];
yNewLim = [1 8];
colorBar = [0 .8];
thisAxisFontsize = 20;

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

    
    %ylabel('Time(s)')
    %zlabel('Chan. Resp.(uV2)')
    %title([mc ' / ' ec])
    %view(70,30)
    %zlim([0 1])
    
    colormap('jet')
    caxis([0,.4]);
    pbaspect([2 1 1])
    
    set(gca,'xlim',[1 640],...
        'XTick',xNewTick,...
        'xTickLabel',xNewTickLabel,...
        'YTick',1:8,...
        'YTickLabel',yNewTickLabel,...
        'FontSize',24,...
        'LineWidth',2,...
        'layer','top')
    
    
    % stim onset line
   line([128,128],[.5,8.5],...
       'linewidth',3,...
       'color','w',...
       'linestyle','--'); 
   
       % stim offset line
   line([192,192],[.5,8.5],...
       'linewidth',3,...
       'color','w',...
       'linestyle','--'); 
   

   if ismember(iCond,[3,4])
       
       for iEyeLine=1:3
           
           if iEyeLine==1; xLoc = 256;
           elseif iEyeLine==2; xLoc = 256+(256/2);
           elseif iEyeLine==3; xLoc = 512;
           end
           
           % stim offset line
           line([xLoc,xLoc],[.5,8.5],...
               'linewidth',3,...
               'color','w',...
               'linestyle',':');
           
       end
       
   end
    
    
    
%     set(gca,'yTick',xNewTick,'yticklabel',xNewTickLabel,'xTick',yNewTick,...
%         'xticklabel',yNewTickLabel,'ylim',xNewLim,'fontsize',thisAxisFontsize,...
%         'LineWidth',4,'ztick',zNewTick,'zticklabel',zNewTickLabel);
%     grid('off')
%     
    saveas(h,[plotDir '/' 'IEM_Heat_Fixed_COLOR_' mc '_' ec '.eps'],'epsc')
    saveas(h,[plotDir '/' 'IEM_Heat_Fixed_COLOR_' mc '_' ec '.jpg'],'jpeg')
end





