% ZY: V3 for 6 wells
% 1. need to save image inside assay one for loop of well (done)
% 2. need to filter the electrode to ensure each electrode will be one cell
% (done 20200213 with el2cell function and modified code)
% ALl work on 20200211
% ASSAY 1 
% WHOLE SAMPLE ELECTRICAL ACTIVITY QUANTIFICATION

% filename1 = 'existingData.csv';
% % Read the CSV as a table
% t = readtable(filename1);
% % Add a new column to the end of the table
% numOfColumn = size(t, 2);
% newCol = num2cell(d1(1:end,7)); % Your new column
% t.(numOfColumn+1) = newCol;
% % Change column name if needed
% t.Properties.VariableNames{numOfColumn+1} = 'newCol';
% % Write to CSV file
% writetable(t, 'new.csv')


%% variables
clear;
close all;
Assay = [0,1,1]; %activity,network,filter
% Path- choose a date
dateFolder = uigetdir('/Users/zye10/Documents/Maxone/MaxData/Britt_ZY/');
% dateFolder = uigetdir('/Volumes/Dev_Electrophys/Zengyou_Ephys_data/SCZ_MEA/SCZ_MaxOne/');
[~,dateR] = fileparts(dateFolder);
fprintf(">>>>>>>>starting analyzing activity and network for date %s<<<<<<<<<<\n",dateR);
%%%%%batch analysis: choose a date to analyze all activity and network
global plotPath
global tablePath
project = 'ZY_MaxOne';
outRoot = '/Users/zye10/Documents/Maxone/Projects/';
outPath = fullfile(outRoot,project,'DayBatch/');

plotPath = fullfile(outPath,'plots');
tablePath = fullfile(outPath,'tables');
if ~isfolder(plotPath)
    mkdir(plotPath);
end
if ~isfolder(tablePath)
    mkdir(tablePath);
end
% Get a list of chip folders in this date.
subFolders = getsub(dateFolder);
fprintf("A total of %d chips found:\n",length(subFolders));
fprintf("%s\n",subFolders.name);
% ask to choose chip name with index list
[chip_indx,ckp] = listdlg('ListString',{subFolders.name});

%analyze activity and network for each chip
for i= chip_indx%1:length(subFolders)
    %good for activity
    chipID = subFolders(i).name;
    % dateFolder is subFolder(i).folder
    if Assay(1)
        activityScan = fullfile(dateFolder,chipID,'Activity_Scan');
        runfolderA = getsub(activityScan);
        if ~isempty(runfolderA)
            for j = 1:length(runfolderA)
                runID = runfolderA(j).name;% if more runID for one chip
                fprintf(">>>>>>>>starting activity file %d.%d of %d: chip%s_%s \n ",i,j,length(subFolders),chipID,runID);
                folder_nameA = fullfile(activityScan,runID);
                % save_exl(folder_nameA)
                try
                    assay_1(folder_nameA);
                catch
                    continue
                end
                % savepng(gcf,folder_nameA);
                close
            end
        else
            %continue
            fprintf(">>>>>>>>No activityScan for %d of %d: %s \n ",i,length(subFolders),chipID);
        end
    end
    if Assay(2) || Assay(3)
        %good for network
        networkScan = fullfile(dateFolder,chipID,'Network');
        runfolderN = getsub(networkScan);
        if ~isempty(runfolderN)
            for j = 1:length(runfolderN)
                runID = runfolderN(j).name;% if more runID for one chip;
                fprintf(">>>>>>>>starting network file %d.%d of %d: chip%s_%s \n ",i,j,length(subFolders),chipID,runID);
                folder_nameN = fullfile(networkScan,runID);
                % save_exl(folder_nameN)
                if Assay(2)
                    try
                        assay_3(folder_nameN);
                    catch
                        continue
                    end
                end
                if Assay(3)
                    try
                        assay_3_filter(folder_nameN)
                    catch
                        continue
                    end                        
                end
                % savepng(gcf,folder_nameN);
                close
            end
        else
            fprintf(">>>>>>>>No networkScan for %d of %d: %s \n ",i,length(subFolders),chipID);
            continue
        end
    end
end
fprintf("All done for now \n");

function save_exl(folder_name)
    %% create exel with chipnane
    % create sheetnames with A/N
    % collect the infortable here with max function
    % transpose struct: structfun(@(fold) fold(:), networkAct, 'UniformOutput', false)
    % struct2table, need to change the name
    % T.Properties.VariableNames([1 3]) = {'Gender' 'Height'}
    % T.Properties.VariableNames{'Var4'} = 'Weight'

    global tablePath
    folders = regexp(folder_name,filesep,'split');
    runID = folders{end};
    chipID = folders{end-2};
    dateR = folders{end-3};
    typeR = folders{end-1}(1);
 
    if typeR == 'A' % activity
        %get infortable
        fullScanFolder = mxw.fileManager(folder_name);
        % Compute spikes features
        spikeRate = mxw.activityMap.computeSpikeRate(fullScanFolder);
        amplitude90perc = abs(mxw.activityMap.computeAmplitude90percentile(fullScanFolder));
        infortable = array2table([spikeRate,amplitude90perc]);
        infortable.Properties.VariableNames = {strcat('spikeRate_',chipID,'_',dateR),strcat('amplitude90perc_',chipID,'_',dateR)};
        
    elseif typeR == 'N' % network
        networkAnalysisFile = mxw.fileManager(folder_name);
        %%all these could be parameters in function
        % Bin size for spike counts
        bin_size = 0.02;
        % Gamma of Gaussian to convolve
        gaussian_gamma = 0.3; % in seconds
        % Compute network activity
        networkAct = mxw.networkActivity.computeNetworkAct(networkAnalysisFile, 'BinSize', bin_size, 'file', 1,'GaussianSigma', gaussian_gamma);
        infortable = array2table([transpose(networkAct.time),networkAct.firingRate]);
        infortable.Properties.VariableNames = {'time',strcat('firingRate_',chipID,'_',dateR)};
    end
    % save file
    filename=fullfile(tablePath,strcat(chipID,'.xlsx'));
    writetable(infortable,filename,'sheet',strcat(dateR,'_',typeR,runID));

end
function savepng(fig,folder_name)
    %% save figure_ Zengyou
    global plotPath
    folders = regexp(folder_name,filesep,'split');
    runID = folders{end};
    chipID = folders{end-2};
    dateR = folders{end-3};
    typeR = folders{end-1}(1); % A or N
    
    supTitle = strcat(chipID,'_',dateR,'_',typeR,runID); % set in assay
    % sgtitle(supTitle, 'FontSize', 14', 'FontWeight', 'Bold','Interpreter', 'latex')
    
    filename = fullfile(plotPath,strcat(supTitle,'.png'));
    saveas(fig,filename);
    fprintf("\n file saved: %s \n",filename);
    % clear
end
function assay_1(folder_name)
    %% ASSAY_1
    global tablePath plotPath
    folders = regexp(folder_name,filesep,'split');
    runID = folders{end};
    chipID = folders{end-2};
    dateR = folders{end-3};
    typeR = folders{end-1}(1); % A or N
    %check single well or multi-well plate
    testFileMan = mxw.fileManager(folder_name);
    try
        wellsInfo = h5info(testFileMan.fileNameList{1},'/wells/');
    end
    
    clear fullScanFolder %because this is a cell, need to clear old data
    if exist('wellsInfo')
        
        well_labels = {wellsInfo.Groups.Name};
        if length(well_labels)>1
            % [indx,tf] = listdlg('ListString',well_labels);
            % ZY: make it defaut for 6 well
            indx = 1:6;
            tf = 1;
            if tf == 1
                for j= 1:length(indx)
                    fullScanFolder{j} = mxw.fileManager(folder_name,indx(j));
                end
            end
        else
            fullScanFolder{1} = mxw.fileManager(folder_name,1);
        end
    else
        fullScanFolder{1} = mxw.fileManager(folder_name);
    end
    infor_A = table();
    infor_sum = table(); %save the average value
%     m_Hz = []; %original code is []
%     m_uV = [];
%     m_pct = [];
    for i= 1:length(fullScanFolder)
    
        % Thresholds
        thr_spike_rate = 0.05;
        thr_amp = 5;% change from 15 to 5
    
        % Compute spikes features
        spikeRate = mxw.activityMap.computeSpikeRate(fullScanFolder{i});
        amplitude90perc = abs(mxw.activityMap.computeAmplitude90percentile(fullScanFolder{i}));
        
        % ZY if one well do not have any signal, go to next
        if ~sum(spikeRate)
            continue
        end
    
        %% Plot
        
        % set color bar scale for firing rate and amplitude activity maps
        fr_factor = 5;
        fr_min = thr_spike_rate;
        fr_max = fr_min + (max(spikeRate)-fr_min)/fr_factor;
        
        amp_factor = 4;
        amp_min = thr_amp;
        amp_max = amp_min + (max(amplitude90perc)-amp_min)/amp_factor;
        
        % create figure
        if isprop(fullScanFolder{i}.fileObj,'wellID')
            figName = [chipID '_' dateR '_Well' num2str(fullScanFolder{i}.wellID) '_' typeR runID]; % ZY: add chipID and date
        else
            % figName=''; % ZY: change from '' to chipID-Date
            figName = strcat(chipID,'_',dateR,'_',typeR,runID);
        end
        
        figure('color', 'w','position',[0 100 1300 700]);
        % sgtitle(figName);
        
        % ZY set size
        sgtitle(figName,'FontSize', 14', 'FontWeight', 'Bold','Interpreter', 'latex');
        % firing rate
        subplot(2,3,1);
        mxw.plot.activityMap(fullScanFolder{i}, spikeRate, 'Ylabel', 'Hz', 'CaxisLim', [fr_min fr_max],'Figure',false,'Title','Whole-Sample Mean Firing Rate','colormap','parula');
        line([300 800],[2000+400 2000+400],'Color','k','LineWidth',5); axis off;
        text(340,2100+500,'0.5 mm','color','k');xlim([200 3750]);ylim([150 2500])
        % amplitude
        subplot(2,3,2);
        mxw.plot.activityMap(fullScanFolder{i}, amplitude90perc,  'CaxisLim', [amp_min amp_max], 'Ylabel', 'Amplitude \muV','Figure',false,'Title','Whole-Sample Spike Amplitude','RevertColorMap', true ,'colormap','parula' );
        line([300 800],[2000+400 2000+400],'Color','k','LineWidth',5); axis off;
        text(340,2100+500,'0.5 mm','color','k');xlim([200 3750]);ylim([150 2500])
        
        % active area
        idx = (spikeRate>thr_spike_rate & amplitude90perc>thr_amp);
        % active electrodes
        subplot(2,3,3);
        mxw.plot.activityMap(fullScanFolder{i}, double(idx), 'Ylabel', '','Figure',false,'Title',['Active Electrodes = ',num2str((sum(idx)/length(idx))*100,'%.2f'),' %'],'colormap','parula');
        line([300 800],[2000+400 2000+400],'Color','k','LineWidth',5); axis off; cbh=colorbar;set(cbh,'YTick',[0 1])
        text(340,2100+500,'0.5 mm','color','k');xlim([200 3750]);ylim([150 2500])
        
        %idx = (amplitude90perc>thr_amp);
        % firing rate distribution
        subplot(2,3,4);h = histogram(spikeRate(idx),0:.05:ceil(max(spikeRate)));
        ylabel('Counts');xlabel('Mean Firing Rate [Hz]');box off;
        h.FaceColor = 'b'; h.EdgeColor = 'b'; h.FaceAlpha = 1;
        legend(['MFR = ',num2str(mean(spikeRate(idx)),'%.2f'),' Hz,  sd = ',num2str(std(spikeRate(idx)),'%.2f')])
        
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
        ylabel(['Electrodes Percentage'],'fontsize',10);
        %% ZY: save information
        m_Hz = mean(spikeRate(idx)); %original code is []
        m_uV = mean(amplitude90perc(idx));
        m_pct = sum(idx)/length(idx)*100;
        % save each well to outPath as png
        filename = fullfile(plotPath,strcat(figName,'.png'));
        saveas(gcf,filename);
        fprintf("\n file saved: %s \n",filename);
        close
        %ZY: save excel
        infortable = array2table([spikeRate,amplitude90perc]);
        infortable.Properties.VariableNames = {strcat('spikeRate_',figName),strcat('amplitude90perc_',figName)};
        infor_summary = array2table([m_Hz;m_uV;m_pct]);
        infor_summary.Properties.VariableNames = {strcat('m_Hz_',figName),strcat('m_uV_',figName),strcat('m_pct_',figName)};
        
        if isempty(infor_A)
            infor_A = infortable;
            infor_sum = infor_summary;
        else
            infor_A = [infor_A infortable];
            infor_sum = [infor_sum infor_summary];
        end

    end
    sheetName = strcat(dateR,'_',typeR,runID);
    sheetName_sum = strcat(dateR,'_',typeR,runID,'_sum');
    filename=fullfile(tablePath,strcat(chipID,'.xlsx'));
    writetable(infor_A,filename,'sheet',sheetName);
    writetable(infor_sum,filename,'sheet',sheetName_sum);
    
end

function assay_3(folder_name) %ZY: bad plot for burst duration
    %ZY assay 3
    global tablePath plotPath
    folders = regexp(folder_name,filesep,'split');
    runID = folders{end};
    chipID = folders{end-2};
    dateR = folders{end-3};
    typeR = folders{end-1}(1); % A or N
    %ZY: excel file
    sheetName = strcat(dateR,'_',typeR,runID);
    exl_name=fullfile(tablePath,strcat(chipID,'.xlsx'));
    % ZY remove 1
    testFileMan = mxw.fileManager(folder_name);
    try
        wellsInfo = h5info(testFileMan.fileNameList{1},'/wells/');
    end
    
    clear networkAnalysisFile
    if exist('wellsInfo')
        
        well_labels = {wellsInfo.Groups.Name};
        if length(well_labels)>1
            % [indx,tf] = listdlg('ListString',well_labels);
            % ZY: make it defaut for 6 well
            indx = 1:6;
            tf = 1;
            if tf == 1
                for j= 1:length(indx)
                    %ZY: in case no recording
                    try 
                        networkAnalysisFile{j} = mxw.fileManager(folder_name,indx(j));
                    catch
                        warning('Wrong with chipID well %d: no signal',j);
                        networkAnalysisFile{j}=[]; %ZY: remove the one without signal
                        continue;
                    end
                end
            end
        else
            networkAnalysisFile{1} = mxw.fileManager(folder_name,1);
        end
    else
        networkAnalysisFile{1} = mxw.fileManager(folder_name);
    end
    
    
    %%
    
    % Bin size for spike counts
    bin_size = 0.02;
    % Threshold to detect bursts
    thr_burst = 1.2; % in rms
    % Gamma of Gaussian to convolve
    gaussian_gamma = 0.3; % in seconds
    % threshold to find the start and stop time of the bursts,
    thr_start_stop = 0.3; % 0.3 means 30% value of the burst peak
    
    %% execute
    
    % Load data information
    % infor_N = table();
    for j=1:length(networkAnalysisFile)
        %ZY need to jump out of loop if networkAnalusisFils is empty or not
        %recorded
        if isempty(networkAnalysisFile{j})
            continue
        end
        % Compute network activity
        networkAct = mxw.networkActivity.computeNetworkAct(networkAnalysisFile{j}, 'BinSize', bin_size, 'file', 1,'GaussianSigma', gaussian_gamma);
        % ZY if one well do not have any signal, go to next
        if ~sum(networkAct.firingRate)
            continue
        end
        networkStats = mxw.networkActivity.computeNetworkStats(networkAct, 'Threshold', thr_burst);

        % Plotting:
         if isprop(networkAnalysisFile{j}.fileObj,'wellID')
            wellID = ['Well' num2str(networkAnalysisFile{j}.wellID) '_'];
            figName = [chipID '_' dateR '_' wellID typeR runID]; % ZY: add chipID and date
            
         else
            % figName=''; % ZY: change from '' to chipID-Date
            figName = strcat(chipID,'_',dateR,'_',typeR,runID);
        end
        
        figure('color', 'w','position',[0 100 1300 700]);
        % sgtitle(figName);
        
        % ZY set size
        sgtitle(figName,'FontSize', 14', 'FontWeight', 'Bold','Interpreter', 'latex');
                
        % Raster Plot
        ax(1)=subplot(2, 3, 1);
        mxw.plot.rasterPlot(networkAnalysisFile{j}, 'file', 1, 'Figure', false);box off;
        xlim([0 floor(max(networkAct.time))-0.5])
        ylim([0 length(networkAnalysisFile{j}.fileObj.map.channel)])
        % Histogram gaussian convolution
        ax(2)=subplot(2, 3, 4);
        mxw.plot.networkActivity(networkAct, 'Threshold', thr_burst, 'Figure', false);box off;
        hold on;plot(networkStats.maxAmplitudesTimes,networkStats.maxAmplitudesValues,'or')
        xlim([0 floor(max(networkAct.time))-0.5])
        linkaxes(ax, 'x')
        
        if length(networkStats.maxAmplitudesTimes)>2
            
            
            % Burst Peak
            subplot(2, 3, 2);
            mxw.plot.networkStats(networkStats, 'Option', 'maxAmplitude',  'Figure', false, ...
                'Ylabel', 'Counts', 'Xlabel', 'Burst Peak [Hz]', 'Title', 'Burst Peak Distribution','Bins',20,'Color','b'); box off;
            legend(['Mean Burst Peak = ',num2str(mean(networkStats.maxAmplitudesValues),'%.2f'), ' sd = ',num2str(std(networkStats.maxAmplitudesValues),'%.2f')])
            
            % IBI
            subplot(2, 3, 3);
            mxw.plot.networkStats(networkStats, 'Option', 'maxAmplitudeTimeDiff',  'Figure', false,...
                'Ylabel', 'Counts', 'Xlabel', 'Interburst Interval [s]', 'Title', 'Interburst Interval Distribution','Bins',20,'Color','b'); box off;
            legend(['Mean Interburst Interval = ',num2str(mean(networkStats.maxAmplitudeTimeDiff),'%.2f'),' sd = ',num2str(std(networkStats.maxAmplitudeTimeDiff),'%.2f')])
            
            % Synchrony, Percentage Spikes within burst
            subplot(2, 3, 5);
            % Burst Amplitudes
            amp = networkStats.maxAmplitudesValues';
            % Burst Times
            peak_times = networkStats.maxAmplitudesTimes;
            
            edges = [];
            for i = 1:length(amp)
                % take a sizeable (±6 s) chunk of the network activity curve
                % around each burst peak point
                % ZY: 3 rather than 6 works better
                idx = networkAct.time>(peak_times(i)-3) & networkAct.time<(peak_times(i)+3);
                t1 = networkAct.time(idx);
                a1 = networkAct.firingRate(idx)';
                
                % get the amplitude at the desired peak width
                peakWidthAmp = (amp(i)-round(amp(i)*thr_start_stop));
                % get the indices of the peak edges
                idx1 = find(a1<peakWidthAmp & t1<peak_times(i));
                idx2 = find(a1<peakWidthAmp & t1>peak_times(i));
                
                %         if ~isempty(idx1) && ~isempty(idx2)
                %         t_before = [];
                %         t_after = [];
                %
                %         else
                
                if ~isempty(idx1) && ~isempty(idx2)
                    t_before = t1(idx1(end));
                    t_after = t1(idx2(1));
                    %         end
                    edges = [edges; t_before t_after];
                end
            end
            
            subplot(2, 3, 1);
            hold on;
            for i = 1:length(edges)
                line([edges(i,1),edges(i,1)],[0 length(networkAnalysisFile{j}.fileObj.map.channel)],'Color','b')
                line([edges(i,2),edges(i,2)],[0 length(networkAnalysisFile{j}.fileObj.map.channel)],'Color','r')
            end
            
            % identify spikes that fall within the bursts
            ts = ((double(networkAnalysisFile{j}.fileObj(1).spikes.frameno) - double(networkAnalysisFile{j}.fileObj(1).firstFrameNum))/networkAnalysisFile{j}.fileObj(1).samplingFreq)';
            ch = networkAnalysisFile{j}.fileObj(1).spikes.channel;
            spikes_per_burst = [];
            ts_within_burst = [];
            ch_within_burst = [];
            
            for i = 1:length(edges)
                
                idx = (ts>edges(i,1) & ts<edges(i,2));
                spikes_per_burst = [spikes_per_burst sum(idx)];
                
                ts_within_burst = [ts_within_burst ts(idx)];
                ch_within_burst = [ch_within_burst ch(idx)'];
                
            end
            
            % Synchrony, Percentage Spikes within burst
            subplot(2, 3, 5);
            h = histogram(spikes_per_burst,20);
            h.FaceColor = 'b'; h.EdgeColor = 'b'; h.FaceAlpha = 1;
            box off;ylabel('Counts');xlabel('Number of Spikes Per Burst')
            title(['Spikes Within Burst = ', num2str(sum(spikes_per_burst/length(ts))*100,'%.1f'),' %'])
            legend(['Mean Spikes Per Burst = ',num2str(mean(spikes_per_burst),'%.2f'), ' sd = ',num2str(std(spikes_per_burst),'%.2f')])
            
            % Burst Duration
            subplot(2, 3, 6);
            h = histogram(abs(edges(:,1) - edges(:,2)),20);
            h.FaceColor = 'b'; h.EdgeColor = 'b'; h.FaceAlpha = 1;
            box off;ylabel('Counts');xlabel('Time [s]')
            title(['Burst Duration'])
            legend(['Mean Burst Duration = ',num2str(mean(abs(edges(:,1) - edges(:,2))),'%.2f'), ' s sd = ',num2str(std(abs(edges(:,1) - edges(:,2))),'%.2f')])
            
        end
        %% ZY: save information
        % save each well to outPath as png
        filename = fullfile(plotPath,strcat(figName,'.png'));
        saveas(gcf,filename);
        fprintf("\n file saved: %s \n",filename);
        close
        %ZY: save excel
        infortable = array2table([transpose(networkAct.time),networkAct.firingRate]);
        infortable.Properties.VariableNames = {'time',strcat('firingRate_',figName)};
        colNum = (j-1)*2+1;
        colName=num2col(colNum);
        % sheetName = strcat(dateR,'_',typeR,runID);
        % filename=fullfile(tablePath,strcat(chipID,'.xlsx'));
        writetable(infortable,exl_name,'sheet',sheetName,'range',colName);
        % pause(3)
    end

end
function assay_3_filter(folder_name)
    %ZY assay 3 with el2cell to pick peaked electrode as cell
    global tablePath plotPath
    folders = regexp(folder_name,filesep,'split');
    runID = folders{end};
    chipID = folders{end-2};
    dateR = folders{end-3};
    typeR = folders{end-1}(1); % A or N
    %ZY: excel file
    sheetName = strcat(dateR,'_',typeR,runID,'_filtered');
    exl_name=fullfile(tablePath,strcat(chipID,'.xlsx'));
    % ZY remove 1
    testFileMan = mxw.fileManager(folder_name);
    try
        wellsInfo = h5info(testFileMan.fileNameList{1},'/wells/');
    end
    
    clear networkAnalysisFile
    if exist('wellsInfo')
        
        well_labels = {wellsInfo.Groups.Name};
        if length(well_labels)>1
            % [indx,tf] = listdlg('ListString',well_labels);
            % ZY: make it defaut for 6 well
            indx = 1:6;
            tf = 1;
            if tf == 1
                for j= 1:length(indx)
                    %ZY: in case no recording
                    try 
                        networkAnalysisFile{j} = mxw.fileManager(folder_name,indx(j));
                    catch
                        warning('Wrong with chipID well %d: no signal',j);
                        networkAnalysisFile{j}=[]; %ZY: remove the one without signal
                        continue;
                    end
                end
            end
        else
            networkAnalysisFile{1} = mxw.fileManager(folder_name,1);
        end
    else
        networkAnalysisFile{1} = mxw.fileManager(folder_name);
    end
    
    
    %%
    
    % Bin size for spike counts
    bin_size = 0.02;
    % Threshold to detect bursts
    thr_burst = 1.2; % in rms
    % Gamma of Gaussian to convolve
    gaussian_gamma = 0.3; % in seconds
    % threshold to find the start and stop time of the bursts,
    thr_start_stop = 0.3; % 0.3 means 30% value of the burst peak
    
    %% execute
    
    % Load data information

    for j=1:length(networkAnalysisFile)
        %ZY need to jump out of loop if networkAnalusisFils is empty or not
        %recorded
        if isempty(networkAnalysisFile{j})
            continue
        end
        % Compute network activity
        % check multiwell or not
        if isprop(networkAnalysisFile{j}.fileObj,'wellID')
            wellID = ['Well' num2str(networkAnalysisFile{j}.wellID) '_'];
            figName = [chipID '_' dateR '_' wellID typeR runID '_filtered'];
        else
            figName=strcat(chipID,'_',dateR,'_',typeR,runID,'_filtered');
        end
        %ZY: new add as suggested
        el_list = el2cell(networkAnalysisFile{j},figName);% the list we want to filter out, save figure with wellID
        times_channels = mxw.util.computeRelativeSpikeTimes(networkAnalysisFile{j});
        idx = ismember(networkAnalysisFile{j}.fileObj.map.electrode,el_list);
        selected_channels = unique(networkAnalysisFile{j}.fileObj.map.channel(idx),'stable');
        spx.time = times_channels.time(ismember(times_channels.channel, selected_channels));
        spx.channel = times_channels.channel(ismember(times_channels.channel, selected_channels));
        networkAct = mxw.networkActivity.computeNetworkAct(spx, 'BinSize', bin_size, 'file', 1,'GaussianSigma', gaussian_gamma);
        if ~sum(networkAct.firingRate)
            continue
        end
        %ZY delete as suggested: networkAct = mxw.networkActivity.computeNetworkAct(networkAnalysisFile{j}, 'BinSize', bin_size, 'file', 1,'GaussianSigma', gaussian_gamma);
        networkStats = mxw.networkActivity.computeNetworkStats(networkAct, 'Threshold', thr_burst);
        
        % Plotting:

        
        figure('color', 'w','position',[0 100 1300 700]);
        sgtitle(figName,'FontSize', 14', 'FontWeight', 'Bold','Interpreter', 'latex');
        
        % Raster Plot
        ax(1)=subplot(2, 3, 1);
        %ZY replace mxw.plot.rasterPlot(networkAnalysisFile{j}, 'file', 1, 'Figure', false);box off;
        mxw.plot.rasterPlot(spx,'Figure', false);box off;
        xlim([0 floor(max(networkAct.time))-0.5])
        ylim([0 length(networkAnalysisFile{j}.fileObj.map.channel)])
        % Histogram gaussian convolution
        ax(2)=subplot(2, 3, 4);
        mxw.plot.networkActivity(networkAct, 'Threshold', thr_burst, 'Figure', false);box off;
        hold on;plot(networkStats.maxAmplitudesTimes,networkStats.maxAmplitudesValues,'or')
        xlim([0 floor(max(networkAct.time))-0.5])
        linkaxes(ax, 'x')
        
        if length(networkStats.maxAmplitudesTimes)>2
            
            
            % Burst Peak
            subplot(2, 3, 2);
            mxw.plot.networkStats(networkStats, 'Option', 'maxAmplitude',  'Figure', false, ...
                'Ylabel', 'Counts', 'Xlabel', 'Burst Peak [Hz]', 'Title', 'Burst Peak Distribution','Bins',20,'Color','b'); box off;
            legend(['Mean Burst Peak = ',num2str(mean(networkStats.maxAmplitudesValues),'%.2f'), ' sd = ',num2str(std(networkStats.maxAmplitudesValues),'%.2f')])
            
            % IBI
            subplot(2, 3, 3);
            mxw.plot.networkStats(networkStats, 'Option', 'maxAmplitudeTimeDiff',  'Figure', false,...
                'Ylabel', 'Counts', 'Xlabel', 'Interburst Interval [s]', 'Title', 'Interburst Interval Distribution','Bins',20,'Color','b'); box off;
            legend(['Mean Interburst Interval = ',num2str(mean(networkStats.maxAmplitudeTimeDiff),'%.2f'),' sd = ',num2str(std(networkStats.maxAmplitudeTimeDiff),'%.2f')])
            
            % Synchrony, Percentage Spikes within burst
            subplot(2, 3, 5);
            % Burst Amplitudes
            amp = networkStats.maxAmplitudesValues';
            % Burst Times
            peak_times = networkStats.maxAmplitudesTimes;
            
            edges = [];
            for i = 1:length(amp)
                % take a sizeable (±6 s) chunk of the network activity curve
                %ZY: why 6? it doesn't work after filter, try 3
                % around each burst peak point
                idx = networkAct.time>(peak_times(i)-3) & networkAct.time<(peak_times(i)+3);
                t1 = networkAct.time(idx);
                a1 = networkAct.firingRate(idx)';
                
                % get the amplitude at the desired peak width
                peakWidthAmp = (amp(i)-round(amp(i)*thr_start_stop));
                % get the indices of the peak edges
                idx1 = find(a1<peakWidthAmp & t1<peak_times(i));
                idx2 = find(a1<peakWidthAmp & t1>peak_times(i));
                
                %         if ~isempty(idx1) && ~isempty(idx2)
                %         t_before = [];
                %         t_after = [];
                %
                %         else
                
                if ~isempty(idx1) && ~isempty(idx2)
                    t_before = t1(idx1(end));
                    t_after = t1(idx2(1)); %ZY wrong
                    %         end
                    edges = [edges; t_before t_after]; % ZY wrong at column 66 and others
                end
            end
            
            subplot(2, 3, 1);
            hold on;
            for i = 1:length(edges)
                line([edges(i,1),edges(i,1)],[0 length(networkAnalysisFile{j}.fileObj.map.channel)],'Color','b')
                line([edges(i,2),edges(i,2)],[0 length(networkAnalysisFile{j}.fileObj.map.channel)],'Color','r')
            end
            
            % identify spikes that fall within the bursts
            ts = ((double(networkAnalysisFile{j}.fileObj(1).spikes.frameno) - double(networkAnalysisFile{j}.fileObj(1).firstFrameNum))/networkAnalysisFile{j}.fileObj(1).samplingFreq)';
            ch = networkAnalysisFile{j}.fileObj(1).spikes.channel;
            spikes_per_burst = [];
            ts_within_burst = [];
            ch_within_burst = [];
            
            for i = 1:length(edges)
                
                idx = (ts>edges(i,1) & ts<edges(i,2));
                spikes_per_burst = [spikes_per_burst sum(idx)];
                
                ts_within_burst = [ts_within_burst ts(idx)];
                ch_within_burst = [ch_within_burst ch(idx)'];
                
            end
            
            % Synchrony, Percentage Spikes within burst
            subplot(2, 3, 5);
            h = histogram(spikes_per_burst,20);
            h.FaceColor = 'b'; h.EdgeColor = 'b'; h.FaceAlpha = 1;
            box off;ylabel('Counts');xlabel('Number of Spikes Per Burst')
            title(['Spikes Within Burst = ', num2str(sum(spikes_per_burst/length(ts))*100,'%.1f'),' %'])
            legend(['Mean Spikes Per Burst = ',num2str(mean(spikes_per_burst),'%.2f'), ' sd = ',num2str(std(spikes_per_burst),'%.2f')])
            
            % Burst Duration
            subplot(2, 3, 6);
            h = histogram(abs(edges(:,1) - edges(:,2)),20);
            h.FaceColor = 'b'; h.EdgeColor = 'b'; h.FaceAlpha = 1;
            box off;ylabel('Counts');xlabel('Time [s]')
            title(['Burst Duration'])
            legend(['Mean Burst Duration = ',num2str(mean(abs(edges(:,1) - edges(:,2))),'%.2f'), ' s sd = ',num2str(std(abs(edges(:,1) - edges(:,2))),'%.2f')])
            
        end
        %% ZY: save information
        % save each well to outPath as png
        filename = fullfile(plotPath,strcat(figName,'.png'));
        saveas(gcf,filename);
        fprintf("\n file saved: %s \n",filename);
        close
        %ZY: save excel
        infortable = array2table([transpose(networkAct.time),networkAct.firingRate]);
        infortable.Properties.VariableNames = {'time',strcat('firingRate_',figName)};
        colNum = j*2-1;
        colName=num2col(colNum);
        % sheetName = strcat(dateR,'_',typeR,runID);
        % filename=fullfile(tablePath,strcat(chipID,'.xlsx'));
        writetable(infortable,exl_name,'sheet',sheetName,'range',colName);
        pause(3)
    end
end
function el_list = el2cell(networkAnalysisFile,figName)
    global plotPath
    %% default filter   
%     switch nargin
%         case 1
%             thr_spike_rate = 0.05;
%             thr_amp = 15;
%         case 2
%             thr_amp = 15;
%         otherwise
%             error('min 1 and max 3 inputs are accepted.')           
%     end        
    thr_spike_rate = 0.02;
	thr_amp = 5;
    % get chipID, D=dateR,runID
%     folder_name = networkAnalysisFile.referencePath;
%     folders = regexp(folder_name,filesep,'split');
%     runID = folders{end};
%     chipID = folders{end-2};
%     dateR = folders{end-3};
%     typeR = folders{end-1}(1); % A or N
    png_name = strcat(figName,'_file'); %ZY need well info for plate
    raw=zeros(220,120);
    map=zeros(220,120); % transpose after loading value
    % the function MEA to get the 90% amplitude, will not have NaN value
    spikeRate = mxw.activityMap.computeSpikeRate(networkAnalysisFile);
    amplitude90perc = abs(mxw.activityMap.computeAmplitude90percentile(networkAnalysisFile));
    electrode = networkAnalysisFile.processedMap.electrode;
    % idx1 = amplitude90perc>thr_amp;
    idx1 = (spikeRate>thr_spike_rate & amplitude90perc>thr_amp);
    amplitude90perc_f = amplitude90perc(idx1);% filter by threshold
    % spikeRate_f = spikeRate(idx1);
    electrode_f=electrode(idx1);
    for i=1:length(electrode)
        idx2 = electrode(i);
        raw(idx2+1) = amplitude90perc(i);
    end
    for i=1:length(electrode_f)
        idx2 = electrode_f(i);
        map(idx2+1) = amplitude90perc_f(i);
    end
    colormap('parula'); % hot or parula
    x=[1 220];y=[1 120];
    clims = [0 max(amplitude90perc)];
    imagesc(x,y,raw',clims);
    
    set(gcf,'Position',[100 100 880 480]) % electrode matrix 220x120
    colorbar;
    %% get the new index with peak amp, save the peak electrode information
    % need to add para for distance nearby
    hold on
    regmax = imregionalmax(map);
    spy(regmax','ms',5)
    % save png
    sgtitle(png_name,'FontSize', 14', 'FontWeight', 'Bold','Interpreter', 'latex')
    filename = fullfile(plotPath,strcat(png_name,'.png'));
    saveas(gcf,filename);
    fprintf("\n filter file saved: %s \n",filename);
    close
    % electrode picked: output
    el_list = find(regmax==1)-1;
    %would be the index for maxwell function
%     final_idx = double.empty(0,length(el_list)); % the position in the raw electrode list
%     for i=1:length(el_list)
%         final_idx(i)=find(fullScanFolder.processedMap.electrode==el_list(i));
%     end
    
end
function colName = num2col(colNum)
    % create excel number to column
    z = 'A':'Z';
    colLetter = arrayfun(@(x)z(rem(floor(x*26.^(1-floor(log(x)/log(26)+1):0)),26)),colNum,'un',0);
    colName = sprintf('%s1', colLetter{1});
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