%{
Purpose: loops through raw data, gets ITC, plot stuff, do ITC conversion
into long format for R studio processing
Author: Tom
Date: 06.13.17
%}

clear 
close all

%sourceFolder = '/home/waldrop/Desktop/WTF_EYE/EEG_Bandpassed';
sourceFolder = '/home/waldrop/Desktop/WTF_EYE/EEG_Processed';
%destFolder = '/home/bullock/WTF/Data_Compiled';
destFolder = '/home/waldrop/Desktop/WTF_EYE/Data_Compiled';

%% create plots?
suppressPlots = 'on';

subjects = [1:7,9:14,16:20,22:27,31];

elects = [25:27 29 30 62:64];

% preallocate "theseData" mat and fill with NANs
%%%allRawData = NaN(18,4,640,960);

for i=1:length(subjects)
    
    sjNum = subjects(i)
    
    load([sourceFolder '/' sprintf('sj%02d_eeg.mat',sjNum)])
    
    for j=1:4
        
        if j==1; title = ['Sj ' num2str(sjNum) ' Eyes Closed PO and O Elects']; thisVert = 250;
        elseif j==2; title = ['Sj ' num2str(sjNum) ' Eyes Open PO and O Elects']; thisVert = 250;
        elseif j==3; title = ['Sj ' num2str(sjNum) ' Eyes Closed Masked PO and O Elects']; thisVert = [250 300];
        elseif j==4; title = ['Sj ' num2str(sjNum) ' Eyes open Masked PO and O Elects']; thisVert = [250 300];
        end
        
        %load([sourceFolder '/' sprintf('sj%02d_cond%02d_Alpha.mat',sjNum,j)])
        
        % reformat to fit existing script
        clear eeg
        %eeg.data = permute(band.eeg,[3,1,2]);
        eeg.data = permute(allConds(j).eegs,[3,1,2]);
        
        theseData = squeeze(mean(eeg.data(:,elects,1:768),2))';
        times = [-500,768,256];%-500:2501/640:2000;
        avewidth = 50; % 1=no smoothing, increase number to increase vertical smoothing
        decimate = 1; % 1=default
        
        % set axes
        theseLimits = [-500, 2000, -30, 30, -3,13, 0, 1, 10];
        theseCbarLims = [-25 25];
        
        % % phasesort options
        % ms_center = 0;
        % prct = 0; % percent trials to reject for low amp (default = 25)
        % freq = 8; % lower end of freq range (sorts by phase at freq of max power in this range)
        % maxfreq = 12; % higher end of freq range
        
        % % phase sort at target onset (0ms)
        % figure; erpimage(theseData,[],times, title, avewidth, decimate, 'vert',[250,300], 'phasesort', [0, 25, 8, 12]) % plots phase sorted trials
        
        
        
        
        
        % get phase coherence info (unsorted trials)
        if ~strcmp(suppressPlots,'on')
            figure;
        end
        
        [outdata,outvar,outtrials,limits,axhndls, ...
            erp,amps,cohers,cohsig,ampsig,outamps,...
            phsangls,phsamp,sortidx,erpsig] =erpimage(theseData,[],times, title, avewidth, decimate, 'vert',thisVert, 'coher',[10, 10, 0.05],'limits',theseLimits,'cbar','on','caxis',[-25 25], 'NoShow','on');
        
        % create a var with missingn trials filled with NaNs 9theseData =
        % 640*960)
        
        % if trials missing, replace with NaNs before storing in mat.
        if size(theseData,2)~=960
            for t=size(theseData,2):960
                theseData(:,t) = NaN; 
            end
        end
        
        allRawData(i,j,:,:) = theseData;
        %allSmoothData(i,j,:,:) = outData;
        allITC(i,j,:) = cohers;
        allITCsig(i,j,:) = cohsig;
        nRawTrial(i,j) = size(eeg.data,1);
        nSmootTrials(i,j) = size(outdata,2);
    end
end

% save ITC
save([destFolder '/' 'ITC_Results.mat'],'allITC')


% plot all mean coherence
% h=figure;
% plot(-500:2501/640:2000,squeeze(mean(allITC(:,1,:),1)),'k'); hold on
% plot(-500:2501/640:2000,squeeze(mean(allITC(:,2,:),1)),'c'); hold on
% plot(-500:2501/640:2000,squeeze(mean(allITC(:,3,:),1)),'g'); hold on
% plot(-500:2501/640:2000,squeeze(mean(allITC(:,4,:),1)),'m'); hold on

h=figure('OuterPosition',[119    77   713   501],'Units','normalized');

for iCond=1:4
    
    if iCond==1; thisLineStyle='-'; thisColor = 'r'; % S/F
    elseif iCond==2; thisLineStyle='-'; thisColor = 'b'; % C/F
    elseif iCond==3; thisLineStyle='-'; thisColor = 'g'; % S/M
    elseif iCond==4; thisLineStyle='-'; thisColor = 'm'; % C/M
    end
    
    plot(linspace(-250,2000,577),squeeze(mean(allITC(:,iCond,64:640),1)),'color',thisColor,'LineWidth',3); hold on
    
end


set(gca,...
    'xlim',[-250,2000],...
    'XTick',[-250, 0:500:2000],...
    'XTickLabel',[{'-.25','  0 ',' 0.5',' 1  ',' 1.5',' 2  '}],...
    'box','off',...
    'LineWidth',1.5,...
    'FontSize',24);

pbaspect([2,1,1])


% stim offset line
for iLine=1:5
    
    if      iLine==1; thisX=0;
    elseif  iLine==2; thisX=250;
    elseif  iLine==3; thisX=500;
    elseif  iLine==4; thisX=1000;
    elseif  iLine==5; thisX=1500;
    end
    
    line([thisX,thisX],[0,.3],...
        'linewidth',3,...
        'color','k',...
        'linestyle','--');
    
end

ylabel('ITC (Normalized)','fontsize',24)
xlabel('Time (s)','fontsize',24)

legend('Spatial/Fix','Color/Fix','Spatial/Move','Color/Move','fontsize',16)


saveas(h,['/home/waldrop/Desktop/WTF_EYE/Plots/' 'ITC_Plot.jpg'],'jpeg')
saveas(h,['/home/waldrop/Desktop/WTF_EYE/Plots/' 'ITC_Plot.eps'],'epsc')


% vline(thisSR*.5,'k--')
% vline(thisSR*.75,'k--')
% vline(thisSR*1,'k--')
% vline(thisSR*1.5,'k--')
% vline(thisSR*2,'k--')













% % plot all ITCs
% h=figure;
% for i=1:4
%     subplot(1,4,i)
%     imagesc(-500:2501/640:200,1:500,squeeze(nanmean(allRawData(:,i,:,:)))',[-20 20]);
%     pbaspect([1,1,1])
% end


% convert data from matlab (18x4x640) to LF (72x4x640) for R studio processing
%% SLOPE CONVERTSION TO LONG FORMAT FOR R

% convert into R "long" format (72x4)
dumEyes = [zeros(18,1), ones(18,1), zeros(18,1), ones(18,1)];
dumMask = [zeros(18,1), zeros(18,1), ones(18,1), ones(18,1)];

% create column vector of subjects
sjNums = 1:size(allITC,1);
sjNums = [sjNums, sjNums, sjNums, sjNums]';


%%% REAL TOTAL POWER %%%
realMat = [];
for i=1:size(allITC,3) % loop through times

    % get subs at each timepoint
    realMat = allITC(:,:,i);

    %use colon operator and create a long format matrix (72*4*27*128)...4 cols are DV, ID, EYES, MASKING
    allITC_LF(:,:,i) = [realMat(:), sjNums, dumEyes(:), dumMask(:)];

end

parsave([destFolder '/' 'ITC_DATA.mat'],allITC,allITCsig,nRawTrial,nSmootTrials,allRawData,allITC_LF)

parsave([destFolder '/' 'ITC_DATA_LF.mat'],allITC_LF)

%% SAVE ITC MAT ONLY (for easy export to local machine
parsave([destFolder '/' 'ITC_DATA_FOR_SCATTER.mat'],allITC)