% ASSAY 1 
% WHOLE SAMPLE ELECTRICAL ACTIVITY QUANTIFICATION

%% variables
clear;
close all;
% Path- choose a date
dateFolder = uigetdir('/Volumes/Dev_Electrophys/Zengyou_Ephys_data/SCZ_MEA/SCZ_MaxOne/');

%%%%%batch analysis: choose a date, and analyze all the activity
% Get a list of all files and folders in this folder.
subFolders = getsub(dateFolder);


for i= 6:length(subFolders)
    %good for only one activity
    activityScan = fullfile(subFolders(i).folder,subFolders(i).name,'Activity_Scan');
    runfolder = getsub(activityScan);
    if ~isempty(runfolder)
        for j = 1:length(runfolder)
            runID = runfolder(j).name;% if more runID for one chip
            fprintf(">>>>>>>>starting file %d.%d of %d: %s \n ",i,j,length(subFolders),runID);
            folder_name = fullfile(activityScan,runID);
            assay_1(folder_name);
            desFolder = '/Users/zye10/Documents/Maxone/Figures/';
            savepng(gcf,folder_name,desFolder);
            close
        end
    else
        continue
    end
end
fprintf("All activity done for now \n");

function savepng(fig,folder_name,desFolder)
    %%save figure_ Zengyou
    folders = regexp(folder_name,filesep,'split');
    runID = folders{end};
    chipID = folders{end-2};
    dateR = folders{end-3};
    typeR = folders{end-1}(1); % A or N
    
    supTitle = strcat(chipID,'_',dateR,'_',typeR,runID);
    sgtitle(supTitle, 'FontSize', 14', 'FontWeight', 'Bold','Interpreter', 'latex')
    
    outPath = strcat(desFolder,supTitle,'.png');
    saveas(fig,outPath)
    fprintf("\n file saved: %s \n",outPath);
    % clear
end
function assay_1(folder_name)
    % ASSAY_1
    fullScanFolder = mxw.fileManager(folder_name);

    % Thresholds
    thr_spike_rate = 0.05;
    thr_amp = 15;
    
    % Compute spikes features
    spikeRate = mxw.activityMap.computeSpikeRate(fullScanFolder);
    amplitude90perc = abs(mxw.activityMap.computeAmplitude90percentile(fullScanFolder));
    
    %% Plot
    
    fr_factor = 5;
    fr_min = thr_spike_rate; % 0.2;
    fr_max = fr_min + (max(spikeRate)-fr_min)/fr_factor;
    
    amp_factor = 4;
    amp_min = thr_amp; % 20;
    amp_max = amp_min + (max(amplitude90perc)-amp_min)/amp_factor;
    
    
    % firing rate
    figure('color',[1 1 1],'position',[100 100 1300 700]);
    subplot(2,3,1);
    mxw.plot.activityMap(fullScanFolder, spikeRate, 'Ylabel', 'Hz', 'CaxisLim', [fr_min fr_max],'Figure',false,'Title','Whole-Sample Mean Firing Rate','colormap','parula');
    line([300 800],[2000+400 2000+400],'Color','k','LineWidth',5); axis off;
    text(340,2100+500,'0.5 mm','color','k');xlim([200 3750]);ylim([150 2500])
    % amplitude
    subplot(2,3,2);
    mxw.plot.activityMap(fullScanFolder, amplitude90perc,  'CaxisLim', [amp_min amp_max], 'Ylabel', 'Amplitude \muV','Figure',false,'Title','Whole-Sample Spike Amplitude','RevertColorMap', true ,'colormap','parula' );
    line([300 800],[2000+400 2000+400],'Color','k','LineWidth',5); axis off;
    text(340,2100+500,'0.5 mm','color','k');xlim([200 3750]);ylim([150 2500])
    
    % active area
    idx = (spikeRate>thr_spike_rate & amplitude90perc>thr_amp);
    % active electrodes
    subplot(2,3,3);
    mxw.plot.activityMap(fullScanFolder, double(idx), 'Ylabel', '','Figure',false,'Title',['Active Electrodes = ',num2str((sum(idx)/length(idx))*100,'%.2f'),' %'],'colormap','parula');
    line([300 800],[2000+400 2000+400],'Color','k','LineWidth',5); axis off; cbh=colorbar;set(cbh,'YTick',[0 1])
    text(340,2100+500,'0.5 mm','color','k');xlim([200 3750]);ylim([150 2500])
    
    %idx = (amplitude90perc>thr_amp);
    % firing rate distribution
    subplot(2,3,4);h = histogram(spikeRate(idx),0:.05:ceil(max(spikeRate)));
    ylabel('Counts');xlabel('Mean Firing Rate [Hz]');box off;
    h.FaceColor = 'b'; h.EdgeColor = 'b'; h.FaceAlpha = 1;
    legend(['MFR = ',num2str(mean(spikeRate(idx)),'%.2f'),' Hz,  sd = ',num2str(std(spikeRate(idx)),'%.2f')])
    xlim([0 mxw.util.percentile(spikeRate,99)])
    
    idx = (spikeRate>thr_spike_rate);
    % amplitude distribution
    subplot(2,3,5);h = histogram(amplitude90perc(idx),ceil(0:1:max(amplitude90perc)));ylabel('counts');
    ylabel('Counts');xlabel('Mean Spike Amplitude [\muV]');box off;
    h.FaceColor = [0 0.8 0]; h.EdgeColor = [0 0.8 0]; h.FaceAlpha = 1;
    legend(['MSA = ',num2str(mean(amplitude90perc(idx)),'%.2f'),' \muV,  sd = ',num2str(std(amplitude90perc(idx)),'%.2f')])
    
    % Electrode Percentage with MFA and MSA
    subplot(2,3,6); hold on
    x = [0 0.5 1 10 100];
    y = histcounts(spikeRate(spikeRate>thr_spike_rate),x);
    y = y/sum(y)*100;
    c = [0 0 1;0 0 0.7;0 0 0.5;0 0 0.3];
    for i = 1:length(y)
        bar(i,y(i), 'FaceColor',c(i,:));
    end
    
    x = [0 30 100 1000];
    y = histcounts(abs(amplitude90perc),x);
    y = y/sum(y)*100;
    c = [0 1 0;0 0.7 0;0 0.5 0;0 0.3 0];
    for i = 1:length(y)
        bar(i+6,y(i), 'FaceColor',c(i,:));
    end
    
    set(gca,'XTick',[1:4 7:9],'XTicklabel',[{'0-0.5'},{'0.5-1'},{'1-10'},{'>10'},{'0-30'},{'30-100'},{'>100'}],'fontsize',7)
    
    
    xlabel('MFR [Hz]                           MSA [\muV]','fontsize',10)
    ylabel(['Electrodes Percentage'],'fontsize',10)


end

function subFolders = getsub(x)
%get subfolder
    if isfolder(x)
        files = dir(x);
        % Get a logical vector that tells which is a directory.
        dirFlags = [files.isdir];
        % Extract only those that are directories.
        dfolders = files(dirFlags);
        % remove '.' and '..'
        subFolders = dfolders(~ismember({dfolders(:).name},{'.','..'}));
    else
        subFolders = [];
    end
end