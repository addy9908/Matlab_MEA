%% ASSAY 2
% ASSAY_2_Ye_V2: single well works well at 2/3/2020
% V3: date 2/14/2020, make it work for 6-well place, ongoing (finished on
%       2/17/2020)
% 2/21/2020: V4 with bugs fixed
% for the case of multi-scan in one day, like before and after dugs
% we will need to count the total file number first for subplot size
% MULTI-DAY ACTIVITY QUANTIFICATION
%%ZY
clear;
close all;
global plotPath
% global tablePath
project = 'Britt_ZY';
outRoot = '/Users/zye10/Documents/Maxone/Projects/';
outPath = fullfile(outRoot,project,'Activity_Batch');
% outPath = '/Users/zye10/Documents/Maxone/tst/';
plotPath = fullfile(outPath,'plots');
tablePath = fullfile(outPath,'tables');
if ~isfolder(plotPath)
    mkdir(plotPath);
end
if ~isfolder(tablePath)
    mkdir(tablePath);
end
exl_name=fullfile(tablePath,'MaxOne_batch_activity.xlsx');

%% variables

% Thresholds
thr_spike_rate = 0.02;% Hz from 0.05 to 0.02
thr_amp = 5; %defaut 25 uV, 1st at 15, now at 5


% Paths
% main_path = uigetdir('/Volumes/Dev_Electrophys/Zengyou_Ephys_data/SCZ_MEA/SCZ_MaxOne/');
main_path = '/Users/zye10/Documents/Maxone/MaxData/Britt_ZY/';
% main_path = '/Volumes/Dev_Electrophys/Zengyou_Ephys_data/SCZ_MEA/SCZ_MaxOne/';
% cd(main_path)


title2 = 'Assay_2_plate ready';
prompt = {'Enter MEA chip number sep with ,:','Enter the dates'};
dims = [1 35];
%you can type the default here
def_Meas = 'M02078,M02089,M02035,M02088';
def_dates = '200129,200214';
definput = {def_Meas,def_dates};
answer = inputdlg(prompt,title2,dims,definput); %ask for chips and dates

%% SCAN


% MEAs_id = str2double(answer); 
% MEAs_id = str2num(answer{1}); %ZY3 for plate:
MEAs_id = strsplit(answer{1},',');
dates = strsplit(answer{2},',');

%ZY3: for plate need to use MEAS_id{1}
if exist(fullfile(main_path,dates{1},MEAs_id{1},'Activity_Scan'),'dir')
    exp_type = 'Activity_Scan';
else
    error('No Activity Scan found')
end


% make files for each chip, and check well or not
for k = 1:length(MEAs_id)
    chip_id = MEAs_id{k};
    fprintf('>>>>>start the chip %s\n',chip_id);
    miltiFile=0;%multi file per day?
    
    files = {};% all full file path

    for i = 1:length(dates)
        p = fullfile(main_path,dates{i},chip_id,exp_type);%day+type      
      
        %% ZY: need to add code for multi-scan for the chip in one day
        if exist(p,'dir')==7
            % cd(p)
            f = getsub(p); % get runID folder
            %ZY: check if multi file per day
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
    
    clear fullScanFolders %because this is a cell, need to clear old data
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
                        fullScanFolders{m,n} = mxw.fileManager(files{m},indx(n));%date,well
                    end
                end
           
            else
                fullScanFolders{m,1} = mxw.fileManager(files{m},1);
            end
        else
            fullScanFolders{m,1} = mxw.fileManager(files{m});
        end
    end

    % go along wells
    for j = 1:size(fullScanFolders,2) %ZY3: need the code to go through wells here, simply repeat by idx1:6 well first, then days       
        %% ZY: well
        ncol = size(fullScanFolders,1);
        m_Hz = double.empty(0,ncol); %original code is []
        m_uV = double.empty(0,ncol);
        m_pct = double.empty(0,ncol);
        row_names = cell(1,ncol);  
   
        %infor_A = table();
        if isprop(fullScanFolders{1,j}.fileObj,'wellID')
            wellID = fullScanFolders{1,j}.wellID;
            figName = ['Chip_',chip_id,'_Well' num2str(wellID)]; %,'_filtered'];
            VariableNames = {strcat('m_Hz_W',num2str(wellID)),...
                            strcat('m_uV_W',num2str(wellID)),...
                            strcat('m_pct_W',num2str(wellID))};
        else
            figName=['Chip_',chip_id]; %',_filtered');
            VariableNames = {'m_Hz','m_uV','m_pct'};
        end
        nrow = 3;
        width = 420* ncol; height = 240*nrow;
        figure('color','w','position',[10 100 width height]);   
        c = 0;  
        
        for i= 1:ncol %ZY file per well
            
            fullScanFolder = fullScanFolders{i,j};
            %% old code
            % fullScanFolder = mxw.fileManager(files{i},k);%% ZY: need runID, use fullfile, add k for wellN
            
            spikeRate = mxw.activityMap.computeSpikeRate(fullScanFolder);
            amplitude90perc = abs(mxw.activityMap.computeAmplitude90percentile(fullScanFolder));
            idx = (spikeRate>thr_spike_rate & amplitude90perc>thr_amp);
            if ~sum(spikeRate)
                continue
            end
            subplot(3,ncol,i) %should change all length to size(fullScanFolders,1)
            mxw.plot.activityMap(fullScanFolder, spikeRate,'Ylabel', 'Hz', 'CaxisLim', [0 2],'Figure',false,'colormap','parula');
            if i==1
                line([300 800],[2000+400 2000+400],'Color','k','LineWidth',5);
                text(340,2100+500,'0.5 mm','color','k');
            end
            axis off; xlim([200 3750]);ylim([150 2600])
            title(['MFR = ', num2str(mean(spikeRate(idx)),'%.2f'),'Hz']);
            m_Hz(i) = mean(spikeRate(idx));
            
            
            if i<ncol
                c = colorbar;c.Ticks = [];
            end
            
            subplot(3,ncol,i+ncol*2)
            % ZY: add CaxisLim
            mxw.plot.activityMap(fullScanFolder, double(idx),'Ylabel', '','CaxisLim', [0 1], 'Figure',false,'Title',['Active Elctrodes = ',num2str((sum(idx)/length(idx))*100,'%.2f'),' %'],'colormap','parula');
            xlim([200 3750]);ylim([150 2500])
            m_pct(i) = (sum(idx)/length(idx))*100;
            if i<ncol
                c = colorbar;c.Ticks = [];
            else
                c = colorbar;c.Ticks = [0 1];
            end
            set(gca,'Xtick',[],'Ytick',[])
            ax = gca; ax.XColor = [0 0 0];ax.YColor = [1 1 1];
            %% zy date+chip+runID
            folders = regexp(fullScanFolder.referencePath,filesep,'split');
            runID = folders{end};
            % chipID = folders{end-2};
            dateR = folders{end-3};
            % if multi-scan in one day: m, add runID
            if miltiFile
                row_names{i}=strcat(dateR,'_',runID);
            else
                row_names{i} = dateR;
            end
            % typeR = folders{end-1}(1);
            xlabel([row_names{i},'_',figName],'Interpreter', 'latex');
            ax.XLabel.Visible = 'on'; % zy: default is off
            idx = (spikeRate>thr_spike_rate); % might need to be changed ZY:this is amp not fre
            subplot(3,ncol,i+ncol)
            % ZY: change the heatmap range
            mxw.plot.activityMap(fullScanFolder, amplitude90perc, 'Ylabel', '\muV', 'CaxisLim', [thr_amp 60], 'Figure',false,'RevertColorMap', true ,'colormap','parula' );
            xlim([200 3750]);ylim([150 2500]);axis off;
            title(['MSA = ', num2str(mean(amplitude90perc(idx)),'%.2f'),'\muV']);
            m_uV(i) = mean(amplitude90perc(idx)); % collect for figure 2
            if i<ncol
                c = colorbar;c.Ticks = [];
            end
            
        end
        %ZY: save first figure
        if ~isempty(m_Hz)% ZY: ignore if no signal
            savepng2(gcf,[figName '_A1']);% j is not real wellid
            close
            figure('color','w','position',[10 100 400 720]);hold on
            %% summary figure B
        
            subplot(3,1,1)
            bar(m_Hz,'b')
            ylim([0 max(m_Hz)+0.1])
            ylabel('MFR [Hz]');box off
            set(gca,'xtick',[])
            subplot(3,1,2)
            bar(m_uV,'g')
            ylim([0 max(m_uV)+3])
            ylabel('MSA [\muV]');box off
            set(gca,'xtick',[])
            subplot(3,1,3)
            bar(m_pct,'r')
            ylim([0 max(m_pct)*1.2]) %ZY: reset scale
            set(gca,'xtick',1:length(row_names),'xticklabel',row_names)
            ylabel('Active Area [%]');box off
            
            %savepng for figure2
            savepng2(gcf,[figName '_A2']);
            %save_exl
            infortable = array2table(transpose([m_Hz;m_uV;m_pct]),'VariableNames',VariableNames);
            infortable.Properties.RowNames = row_names(~cellfun('isempty',row_names));
            infortable.Properties.DimensionNames{1}='dates';
            infortable;
            % exl_name=fullfile(tablePath,'MaxOne_batch_activity.xlsx');
            
            sheet_name = ['Chip_',chip_id,'_Th_A',num2str(thr_amp)];
            colNum = (j-1)*4+1;
            colName=num2col(colNum);
            writetable(infortable,exl_name,'sheet',sheet_name,'WriteRowNames',true,'range',colName);
            close
        else
            close
        end
    end
end
close
fprintf("\n all done!\n");


% clear

function savepng2(fig,figName)
    %% save figure_ Zengyou
    global plotPath
    % folders = regexp(folder_name,filesep,'split');
    % runID = folders{end};
    % chipID = folders{end-2};
    % dateR = folders{end-3};
    % typeR = folders{end-1}(1); % A or N
    
    supTitle = figName;%onlyput chipID here
    sgtitle(supTitle, 'FontSize', 14', 'FontWeight', 'Bold','Interpreter', 'latex')
    
    filename = fullfile(plotPath,strcat(supTitle,'.png'));
    saveas(fig,filename);
    fprintf("\n file saved: %s \n",filename);
    % clear
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