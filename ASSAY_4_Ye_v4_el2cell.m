%% ASSAY 4
%% Version
% V2: works well at 2/3/2020
% 
% for the case of multi-scan in one day, like before and after dugs
% we will need to count the total file number first for subplot size
% MULTI-DAY BURST METRICS
% V3: date 2/17/2020, make it work for 6-well place, ongoing (finish on
%   2/17/2020
% Next: add el2cell (finished on 2/17/2020)
% 2/20/2020: final improve for excel and add electrode plot
clear;
close all;

global plotPath
% global tablePath
project = 'ZY_MaxOne';
outRoot = '/Users/zye10/Documents/Maxone/Projects/';
outPath = fullfile(outRoot,project,'Network_Batch_filtered/');
% outPath = '/Users/zye10/Documents/Maxone/tst/';
plotPath = fullfile(outPath,'plots');
tablePath = fullfile(outPath,'tables');
if ~isfolder(plotPath)
    mkdir(plotPath);
end
if ~isfolder(tablePath)
    mkdir(tablePath);
end
exl_name=fullfile(tablePath,'MaxOne_batch_network.xlsx');
%% variables
% Bin size for spike counts
bin_size = 0.02;
% set threshold to detect bursts
thr_burst = 1.2; % in rms 1.5
% set bin size
gaussian_bin_size = 0.3; % in seconds 0.1(this is gaussian_gamma in ASSAY_3)
thr_start_stop = 0.3; % 0.3 means 30% value of the burst peak


% Set path to folder containing the recording
% main_path = uigetdir('/Volumes/Dev_Electrophys/Zengyou_Ephys_data/SCZ_MEA/SCZ_MaxOne/');
main_path = '/Volumes/Dev_Electrophys/Zengyou_Ephys_data/SCZ_MEA/SCZ_MaxOne/';
% main_path = '/Users/zye10/Documents/Maxone/MaxData/Britt_ZY/';
% cd(main_path)
% f = dir;
% c = 0;
% dates =[];
% for i = 3:length(f)
%     c = c+1;
%     dates{c} = [f(i).name];
% end
%f = dir([main_path,'/',dates{1}]);

title4 = 'Assay_4_plate ready-filtered';
prompt = {'Enter MEA chip number sep with ,:','Enter the dates'};
dims = [1 35];
%you can type the default here
def_Meas = '4310';%'M02078,M02089,M02035,M02088';
def_dates = '200103,200116,200129';%'200129,200214';
definput = {def_Meas,def_dates};
answer = inputdlg(prompt,title4,dims,definput); %ask for chips and dates

%% Network

% MEAs_id = str2double(answer);
MEAs_id = strsplit(answer{1},',');
dates = strsplit(answer{2},',');

%% matlab could not tell cap or not
if exist(fullfile(main_path,dates{1},MEAs_id{1},'Network'),'dir')
    exp_type = 'Network';  
else
    error('No Network recording found')
end




for k = 1:length(MEAs_id)
    chip_id = MEAs_id{k};    
    % ZY: store infor
    fprintf('>>>>>start the chip %s\n',chip_id);
    miltiFile=0;
    
    files = {}; % all full file path
    for i = 1:length(dates)
        p = fullfile(main_path,dates{i},chip_id,exp_type);%day+type
        
        %ZY: need to add code for multi-scan for the chip in one day
        if exist(p,'dir')==7
            % cd(p)
            f = getsub(p); % get runID folder
            %check if multi
            if length(f)>1 && miltiFile ==0
                miltiFile = 1;
            end
            if ~isempty(f) %the case of more than 3? need to count the total files number for fig size
                for j =1:length(f)
                    files{end+1} = fullfile(p,f(j).name);
                end
            end
        end
    end
    %% zy: check multiwell or not
    try
        testFileMan = mxw.fileManager(files{1});
        wellsInfo = h5info(testFileMan.fileNameList{1},'/wells/');
    end
    
    clear networkAnalysisFiles %because this is a cell, need to clear old data
    % applly el2cell inside the loop

    for m = 1:length(files)
        if exist('wellsInfo')       
            well_labels = {wellsInfo.Groups.Name};
            if length(well_labels)>1
            % [indx,tf] = listdlg('ListString',well_labels);
            % ZY: make it defaut for 6 well
                indx = 1:6;
                tf = 1;
                if tf == 1
                    for n= 1:length(indx)
                        try
                            networkAnalysisFiles{m,n} = mxw.fileManager(files{m},indx(n));%date,well
                        catch
                            warning('Wrong with chipID well %d: no signal',n);
                            networkAnalysisFiles{m,n}=[]; %ZY: remove the one without signal
                            continue;
                        end
                    end
                end
           
            else
                networkAnalysisFiles{m,1} = mxw.fileManager(files{m},1);
            end
        else
            networkAnalysisFiles{m,1} = mxw.fileManager(files{m});
        end
    end
    % go along wells
    for j = 1:size(networkAnalysisFiles,2)
        ncol = size(networkAnalysisFiles,1);
        Bu_amp = double.empty(0,ncol);
        Bu_IBI = double.empty(0,ncol);
        row_names = cell(1,ncol);
        if isempty(networkAnalysisFiles{1,j})
            continue
        elseif isprop(networkAnalysisFiles{1,j}.fileObj,'wellID')
            wellID = networkAnalysisFiles{1,j}.wellID;
            figName = ['Chip_',chip_id,'_Well',num2str(wellID),'_filtered'];
            VariableNames = {strcat('Bu_amp_Hz_W',num2str(wellID)),...
                strcat('Bu_IBI_s_W',num2str(wellID))};
        else
            figName=['Chip_',chip_id ',_filtered'];
            VariableNames = {'Bu_amp_Hz','Bu_IBI_s'};
        end
        max_FR = 0;
        % need to chang size for nrow and n col
        plot_filter = 1;% 0 if not plot
        nrow = 3+plot_filter;
        width = 350* ncol; height = 200*nrow;
        figure('color', 'w','position',[0 100 width height]);
        for i = 1:ncol %each day      
            if isempty(networkAnalysisFiles{i,j})
                continue
            end
            networkAnalysisFile = networkAnalysisFiles{i,j};
            %ZY: new add as suggested
            %% change to new seperate function, 1 means do the plot
%             plot_filter = 1;% 0 if not plot
%             nrow = 3+plot_filter;
            subplot(nrow,ncol,i+ncol*3); %add pick piture to 4th row if plotkey =1
            try
                el_list = el2cell_v3(networkAnalysisFile,plot_filter);% the list we want to filter out, save figure with wellID
            catch
                continue
            end
            if i<ncol
                colorbar.Ticks = [];
            end
            %png_filter = strcat(figName,'_file');
            %title(png_filter,'FontSize', 14', 'FontWeight', 'Bold','Interpreter', 'latex')
            %   filename = fullfile(plotPath,strcat(png_name,'.png'));
            %   saveas(gcf,filename);
            
            %%
            
            times_channels = mxw.util.computeRelativeSpikeTimes(networkAnalysisFile);
            idx = ismember(networkAnalysisFile.fileObj.map.electrode,el_list);
            selected_channels = unique(networkAnalysisFile.fileObj.map.channel(idx),'stable');
            spx.time = times_channels.time(ismember(times_channels.channel, selected_channels));
            spx.channel = times_channels.channel(ismember(times_channels.channel, selected_channels));
            networkAct = mxw.networkActivity.computeNetworkAct(spx, 'BinSize', bin_size, 'file', 1,'GaussianSigma', gaussian_bin_size);
            if ~sum(networkAct.firingRate)
                continue
            end
            %% zy date+chip+runID
            folders = regexp(networkAnalysisFile.referencePath,filesep,'split');
            runID = folders{end};
            % chipID = folders{end-2};
            dateR = folders{end-3};
            % if multi-scan in one day: m, add runID
            if miltiFile
                row_names{i}=strcat(dateR,'-',runID); 
            else                 
                row_names{i} = dateR;
            end

            %% plot
            % networkAnalysisFile = mxw.fileManager(files{i});
        
            % Compute network activity
            % zy:change file value to 1 as ASSAY_3
            % networkAct = mxw.networkActivity.computeNetworkAct(networkAnalysisFile, 'BinSize', 0.02, 'file', 1,'GaussianSigma', gaussian_bin_size);
%             if ~sum(networkAct.firingRate)
%                 continue
%             end
            networkStats = mxw.networkActivity.computeNetworkStats(networkAct, 'Threshold', thr_burst);
            max_FR = max([max_FR max(networkAct.firingRate)]);
            
            % POpulation Firing Rate
            subplot(nrow,ncol,i)
            % zy: change MEAS_id to chip_id
            mxw.plot.networkActivity(networkAct, 'Threshold', thr_burst, 'Figure', false,'Title',[row_names{i} '-',chip_id]);box off;%'\_'
            % Title.Interpreter = 'latex';
            hold on;plot(networkStats.maxAmplitudesTimes,networkStats.maxAmplitudesValues,'or')
            xlim([0 floor(max(networkAct.time))-0.5])
            % ylim([0 ylim_burst]) % needs to be fixed
            
            subplot(nrow,ncol,i+ncol)
            mxw.plot.networkStats(networkStats, 'Option', 'maxAmplitude',  'Figure', false, ...
                'Ylabel', 'Counts', 'Xlabel', 'Burst Peak [Hz]', 'Title', '','Bins',20,'Color','g'); box off;
            legend(['Mean BP = ',num2str(mean(networkStats.maxAmplitudesValues),'%.1f'), ' HZ - sd = ',num2str(std(networkStats.maxAmplitudesValues),'%.1f')])
            % bug here for some without value
            if ~networkStats.maxAmplitudesValues
                xlim([mean(networkStats.maxAmplitudesValues)-std(networkStats.maxAmplitudesValues)*3 mean(networkStats.maxAmplitudesValues)+std(networkStats.maxAmplitudesValues)*3])
            end
            subplot(nrow,ncol,i+ncol*2)
            mxw.plot.networkStats(networkStats, 'Option', 'maxAmplitudeTimeDiff',  'Figure', false,...
                'Ylabel', 'Counts', 'Xlabel', 'Interburst Interval [s]', 'Title', '','Bins',20,'Color','r'); box off;
            legend(['Mean IBI = ',num2str(mean(networkStats.maxAmplitudeTimeDiff),'%.1f'),' s - sd = ',num2str(std(networkStats.maxAmplitudeTimeDiff),'%.1f')])
            % bug here for some without value
            if ~networkStats.maxAmplitudeTimeDiff
                xlim([mean(networkStats.maxAmplitudeTimeDiff)-std(networkStats.maxAmplitudeTimeDiff)*3 mean(networkStats.maxAmplitudeTimeDiff)+std(networkStats.maxAmplitudeTimeDiff)*3])
            end
            Bu_amp(i) = mean(networkStats.maxAmplitudesValues);
            Bu_IBI(i) = mean(networkStats.maxAmplitudeTimeDiff);
            
        end
        
        for i=1:ncol
            subplot(nrow,ncol,i);
            ylim([0 max_FR]);
        end
        %ZY: save first figure
        if ~isempty(Bu_amp)
            savepng2(gcf,[figName,'_N1']);
            close
            
            %2nd plot B
            
            figure('color','w','position',[100 100 400 480]);hold on %change high from 400
            subplot(2,1,1)
            bar(Bu_amp,'g')
            ylim([0 max(Bu_amp)*1.2]) %ZY: original +10
            ylabel('Burst Peak [Hz]');box off
            set(gca,'xtick',[])
            subplot(2,1,2)
            bar(Bu_IBI,'r')
            ylim([0 max(Bu_IBI)*1.2]) %ZY: original +10
            set(gca,'xtick',1:ncol,'xticklabel',row_names)
            ylabel('Interburst Interval [s]');box off
            
            %savepng for figure2
            savepng2(gcf,[figName,'_N2']);
            %save_exl
            infortable = array2table(transpose([Bu_amp;Bu_IBI]),'VariableNames',VariableNames);
            infortable.Properties.RowNames = row_names(~cellfun('isempty',row_names));
            infortable.Properties.DimensionNames{1}='dates';
            % exl_name=fullfile(tablePath,'MaxOne_batch_network.xlsx');
            sheet_name = ['Chip_',chip_id '_filter5'];
            colNum = (j-1)*3+1;
            colName=num2col(colNum);
            writetable(infortable,exl_name,'sheet',sheet_name,'WriteRowNames',true, 'range',colName);
            close
        else
            close
        end
    end           
end
fprintf("\n All done!\n");
% clear
function colName = num2col(colNum)
    % create excel number to column
    z = 'A':'Z';
    colLetter = arrayfun(@(x)z(rem(floor(x*26.^(1-floor(log(x)/log(26)+1):0)),26)),colNum,'un',0);
    colName = sprintf('%s1', colLetter{1});
end
function savepng2(fig,chipID)
    %% save figure_ Zengyou
    global plotPath
    % folders = regexp(folder_name,filesep,'split');
    % runID = folders{end};
    % chipID = folders{end-2};
    % dateR = folders{end-3};
    % typeR = folders{end-1}(1); % A or N
    
    supTitle = chipID;%onlyput chipID here
    sgtitle(supTitle, 'FontSize', 14', 'FontWeight', 'Bold','Interpreter', 'latex')
    
    filename = fullfile(plotPath,strcat(supTitle,'.png'));
    saveas(fig,filename);
    fprintf("\n file saved: %s \n",filename);
    
end

function el_list = el2cell_list(networkAnalysisFile)
%    global plotPath
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
    thr_spike_rate = 0.05;
	thr_amp = 15;
    % get chipID, D=dateR,runID
%     folder_name = networkAnalysisFile.referencePath;
%     folders = regexp(folder_name,filesep,'split');
%     runID = folders{end};
%     chipID = folders{end-2};
%     dateR = folders{end-3};
%     typeR = folders{end-1}(1); % A or N
%    png_name = strcat(figName,'_file'); %ZY need well info for plate
    map=zeros(220,120); % transpose after loading value
    % the function MEA to get the 90% amplitude, will not have NaN value
    spikeRate = mxw.activityMap.computeSpikeRate(networkAnalysisFile);
    amplitude90perc = abs(mxw.activityMap.computeAmplitude90percentile(networkAnalysisFile));
    
    % idx1 = amplitude90perc>thr_amp;
    idx1 = (spikeRate>thr_spike_rate & amplitude90perc>thr_amp);
    amplitude90perc_f = amplitude90perc(idx1);% filter by threshold
    % spikeRate_f = spikeRate(idx1);
    electrode_f=networkAnalysisFile.processedMap.electrode(idx1);
    for i=1:length(electrode_f)
        idx2 = electrode_f(i);
        map(idx2) = amplitude90perc_f(i);
    end
%     colormap('parula'); % hot or parula
%     x=[1 220];y=[1 120];
%     clims = [thr_amp max(amplitude90perc)];
%     imagesc(x,y,map',clims);
%     
%     set(gcf,'Position',[100 100 880 480]) % electrode matrix 220x120
%     colorbar;
%     %% get the new index with peak amp, save the peak electrode information
%     % need to add para for distance nearby
%     hold on
    regmax = imregionalmax(map);
%     spy(regmax','ms',5)
%     % save png
%     sgtitle(png_name,'FontSize', 14', 'FontWeight', 'Bold','Interpreter', 'latex')
%     filename = fullfile(plotPath,strcat(png_name,'.png'));
%     saveas(gcf,filename);
%     fprintf("\n filter file saved: %s \n",filename);
%     close
    % electrode picked: output
    el_list = find(regmax==1);
    %would be the index for maxwell function
%     final_idx = double.empty(0,length(el_list)); % the position in the raw electrode list
%     for i=1:length(el_list)
%         final_idx(i)=find(fullScanFolder.processedMap.electrode==el_list(i));
%     end
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