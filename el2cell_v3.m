% electrode2cell sorting
% Zengyou Ye, date: 02/05/2020

%% test
% data_root = '/Users/zye10/Documents/Maxone/MaxData/';
% project = 'Britt_ZY';
% folder_name=fullfile(data_root,project,'200129/M02078/Network/000002');
% networkAnalysisFile = mxw.fileManager(folder_name,2);
% el_list = el2cell_v3(networkAnalysisFile,1);
function el_list = el2cell_v3(networkAnalysisFile,plot_filter)
    if ~isempty(networkAnalysisFile)
        %filter
        thr_spike_rate = 0.02; %from 0.05 to 0.02
        thr_amp = 5;
        % infor for the title (remake title outside)    
%         folders = regexp(networkAnalysisFile.referencePath,filesep,'split');
%         runID = folders{end};
%         chipID = folders{end-2};
%         dateR = folders{end-3};
        % get recording infor
        spikeRate = mxw.activityMap.computeSpikeRate(networkAnalysisFile);
        amplitude90perc = abs(mxw.activityMap.computeAmplitude90percentile(networkAnalysisFile));
        electrode = networkAnalysisFile.processedMap.electrode;

        if sum(spikeRate) && sum(amplitude90perc) 
            % initialize electrode matrix, transpose after loading value
            raw = zeros(220,120);
            raw_filter=zeros(220,120); %
            
            % idx1 = amplitude90perc>thr_amp;
            idx1 = (spikeRate>thr_spike_rate & amplitude90perc>thr_amp);
            amplitude90perc_f = amplitude90perc(idx1);% filter by threshold
            % spikeRate_f = spikeRate(idx1);
            electrode_f=electrode(idx1);
            % attention, electrode starts from 0
            for i=1:length(electrode)
                idx2 = electrode(i);
                raw(idx2+1) = amplitude90perc(i);
            end
            for i=1:length(electrode_f)
                idx3 = electrode_f(i);
                raw_filter(idx3+1) = amplitude90perc_f(i);
            end
            regmax = imregionalmax(raw_filter);
            % electrode picked
            el_list = find(regmax==1)-1; % number elec back
            % surf(Z)
            % hold on
            % imagesc(Z)
            if nargin>1 && plot_filter
                colormap('parula'); % hot or parula
                x=[1 220];y=[1 120];
                clims = [0 max(amplitude90perc)];
                imagesc(x,y,raw',clims);
                
                % set(gcf,'Position',[0 0 1760 960]) % electrode matrix 220x120
                pbaspect([22 12 1]); % the raio of figure
                set(gca,'xticklabel',[],'yticklabel',[]);
                
                colorbar;
                %% get the new index with peak, save the peak electrode information
                hold on
                
                %would be the index for maxwell function
                % final_idx = double.empty(0,length(elec_picked));
                % for i=1:length(elec_picked)
                %     final_idx(i)=find(electrode==elec_picked(i));
                % end
                spy(regmax','ms',2)
                
                msg1 = sprintf('thr_fre=%s, thr_amp=%d',num2str(thr_spike_rate),thr_amp);
                msg2 = sprintf('total=%d, filtered=%d, picked=%d',...
                    length(electrode),length(electrode_f), length(el_list));
                xlabel({msg1;msg2},'Interpreter', 'latex');
                % title_name = [chipID,'_',dateR, '_',runID];
                % title(title_name,'FontSize', 14', 'FontWeight', 'Bold','Interpreter', 'latex')
            end
        end
    end
end


