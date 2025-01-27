% Importing TTL data into MatlabRewardOnsetStim_data
close all; clear variables; clc;
set(0,'defaultfigurecolor',[1 1 1])

%% Define paths and data to analyze
path2data = 'C:\Users\Anacker1\Desktop\Fiber-Data\Fibers\MH\DIGGING TASK\REVERSAL\'; %Location of the data
path2savefolder = 'C:\Users\Anacker1\Desktop\Fiber-Data\Results\Digging Task\REVERSAL\'; %Path to save

% path2data = '/Volumes/TOSHIBA ART/Backup Server New York/Lauren/FP scripts/Examples/Female coh 7/Day 1 context A/'; %Location of the data
% path2savefolder = '/Volumes/TOSHIBA ART/Backup Server New York/Lauren/FP scripts/Results/Female coh 7/Day 1 context A/'; %Path to save

%% Define the list of mice to analyzed
% Determine it by the folders included in the path2data
mice_list = dir(path2data);
for o = length(mice_list):-1:1
    if mice_list(o).isdir == 0
        mice_list(o) = [];        
    else
        if strcmp(mice_list(o).name,'.') == 1 || strcmp(mice_list(o).name,'..') == 1
            mice_list(o) = [];
        end
    end
end
Nmice = length(mice_list);

Animals_cinnamon_correct = [2 3 4]; %% 5 and 6 had it on the left side. Mouse 1 from the first cohort had reward on the right
if  contains(path2data,'Reversal')
    Animals_cinnamon_correct = setdiff(1:Nmice,Animals_cinnamon_correct);
end

FPrig = 1; % 1 for new rig, 2 for old rig (Pedro's data)
Behav_data = 'Digging';%'Digging';'Barnes',Ymaze
IdChannel = {'dCA1GCamp6f','vCA1GCamp6f'}; %{'Psychlight','jRGECO1b','GCamp6f'};{'GCamp6f',[], 'Psychlight'} 'dCA1GCamp6f','vCA1GCamp6f'
color2plotchannels = {'r','g'};%{'r',[],'g'}'r','g'
MiceperFile = 1;
min2remove = 0.1; % Minutes to remove from the beginning of the recording
TRANGE = [-5 30]; % window size [start time relative to epoc onset, window duration]TRANGE = [-5 30]; 
BASELINE_PER.Onset = [-5 0]; % baseline period within our window
BASELINE_PER.Offset = [-5 0]; % baseline period within our window
BASELINE_PER.ShockOnset = [-5 0]; % baseline period within our window
% BASELINE_PER.Reward = [-3 -1]; % baseline period within our window

LimitXaxis = [-3 3]; % Window to visualize in the plot. If empty, full length of TRANGE
Limits_Selected_Epochs.short = [2 3]; % Duration criteria for selected Choice epochs (in seconds)
Limits_Selected_Epochs.long = [4 180];
Limits_IEI = 3;
Epoch2analyze = 3; %Time in seconds to get measurements aligned to behavior:Nmice

for i = 1:length(IdChannel)
    if ~isempty(IdChannel{i})
        limits2plot.(IdChannel{i}) = []; %If empty, automatically adjusted to the data.
    end
end
limits2plot.all = [];
Color_scale = []; % For the heatmaps. If empty, automatically adjusted to the data.

show_plot = 1; % If 0, plots are not displayed##### if you want the plots to pop up : show_plot= 1
save_plot = 1; % If 0, plots are not saved
reanalysis = 0; % If 1, the code runs for sessions already analyze. If 0, session excluded if analysed
overwrite = 1; % If reanalyzing data, 1 if you want to save them and overwrite results and 0 if not

for m = 1:Nmice
    % % Define the path to the mouse and find the folders to%%%%
    % analyze:Nmice///1
    if contains(path2data,'\')
        path2mouse = [path2data,mice_list(m).name,'\'];
    else
        path2mouse = [path2data,mice_list(m).name,'/'];
    end
    sessions = dir(path2mouse);
    for o = length(sessions):-1:1
        if strcmp(sessions(o).name,'.') == 1 || strcmp(sessions(o).name,'..') == 1  ||sessions(o).isdir == 0
            sessions(o)=[];
        end
    end
    
    % Create a folder to save the data for this mouse and define the path
    if exist([path2savefolder,mice_list(m).name],'dir') == 0
        mkdir([path2savefolder,mice_list(m).name])
    end
    if contains(path2data,'\')
        path2save_mouse = [path2savefolder,mice_list(m).name,'\'];
    else
        path2save_mouse = [path2savefolder,mice_list(m).name,'/'];
    end
   
    %% Loop for all the sessions for the mouse
   
    for s = 1:length(sessions)
        if MiceperFile == 2
            splitname = split(sessions(s).name,'_');
            Mouse_name{1} = splitname{1};
            Mouse_name{2} = splitname{2};
            for o = 1:length(Mouse_name)
                if strcmp(mice_list(m).name,Mouse_name{o}) == 1
                    mouse2analyze = o;
                end
            end
        else
            mouse2analyze = 1;
        end

        % Define the path to the session and create folder to save if needed:
        PATH2SESSION = [path2mouse,sessions(s).name];
        if exist([path2save_mouse,sessions(s).name],'dir') == 0
            mkdir([path2save_mouse,sessions(s).name])
        end
        if contains(path2data,'\')
            PATH2SAVE = [path2save_mouse,sessions(s).name,'\'];
        else
            PATH2SAVE = [path2save_mouse,sessions(s).name,'/'];
        end

        tmp = dir(PATH2SESSION);
        for o = length(tmp):-1:1
            if contains(tmp(o).name,'.')
                tmp(o) = [];
            end
        end
        counter = 0;
        for o = 1:length(tmp)
            if tmp(o).isdir == 1
                counter = counter+1;
            end
        end
        if counter >= 1
            if contains(path2data,'\')
                PATH2SESSION = [PATH2SESSION,'\',tmp(1).name];
            else
                PATH2SESSION = [PATH2SESSION,'/',tmp(1).name];
            end
        end

        % Check if results are already saved for this session
        done = exist([PATH2SAVE,'Session_analysis.mat'],'file');
        if done == 0 || reanalysis == 1
            if exist([PATH2SAVE,'figures'],'dir') == 0
                mkdir([PATH2SAVE,'figures'])
            end

            %% Read the data
            % Now read the specified data from our block into a Matlab structure.
            data = TDTbin2mat(PATH2SESSION, 'TYPE', {'epocs', 'scalars', 'streams'});         

            % Extract the names of the STREAMS in data (A405A or X05A)
            STREAM_STORE1 = cell(1,length(IdChannel));
            STREAM_STORE2 = cell(1,length(IdChannel));
            stream_Names = fieldnames(data.streams);
            for o = 1:length(stream_Names)
                if contains(stream_Names{o},'05') && contains(stream_Names{o},'A')
                    STREAM_STORE1{1} = stream_Names{o};
                elseif contains(stream_Names{o},'05') && contains(stream_Names{o},'C')
                    STREAM_STORE1{2} = stream_Names{o};
                elseif (contains(stream_Names{o},'65') || contains(stream_Names{o},'70')) && contains(stream_Names{o},'A')
                    STREAM_STORE2{1} = stream_Names{o};
                elseif (contains(stream_Names{o},'65') || contains(stream_Names{o},'70')) && contains(stream_Names{o},'C')
                    STREAM_STORE2{2} = stream_Names{o};
                end
            end

            %Adjust the lengths of the channels
            minLength = ones(length(STREAM_STORE1),1)*nan;
            for o = 1:length(STREAM_STORE1)
                if ~isempty(STREAM_STORE1{o})
                    minLength(o) = length(data.streams.(STREAM_STORE1{o}).data);
                end
            end
            minLength = min(minLength);

            for o = 1:length(STREAM_STORE1)
                if ~isempty(STREAM_STORE1{o})
                    if length(data.streams.(STREAM_STORE1{o}).data) > minLength
                        data.streams.(STREAM_STORE1{o}).data...
                            (minLength+1:length(data.streams.(STREAM_STORE1{o}).data)) = [];
                    end
                    if length(data.streams.(STREAM_STORE2{o}).data) > minLength
                        data.streams.(STREAM_STORE2{o}).data...
                            (minLength+1:length(data.streams.(STREAM_STORE2{o}).data)) = [];
                    end
                end
            end

            % Gets time vector for our data
            Fs = data.streams.(STREAM_STORE1{1}).fs;
            Max_idx = minLength;
            dt = 1/Fs;
            time = 0:dt:dt*(Max_idx-1);

            dummie = 1:length(time);
            dummie = dummie(time > min2remove*60);
            idx_remove = dummie(1);
            clear dummie

            dwn_rate = 10;
            Fs = Fs/dwn_rate;
            %             time_dwn = arrayfun(@(i) mean(time(i:i+dwn_rate-1)),1:dwn_rate:length(time)-dwn_rate+1);
            time_dwn = downsample(time(idx_remove:end),dwn_rate);
            N = length(time_dwn);
            clear time

            %% Get the dFF for both sensors
            % Pre_processing parameters:

            detrend_method = 1; % 0, no detrending; 1, using global detrending; 2, using local detreding
            filter_flag = 1; % 0, no highpass filter. 1, highpass filter.
            zeroing_flag = 0; % 0, no zeroing; 1, zeroing using moving window and percentile substraction, even while it does not do zeroing it is not going through the very short ones< 1 min therefore you need to adjust following lines
            zeroing_moving_wnd = 30; % Size of the moving window for the zeroing. In seconds.The original timing was 60
            zeroing_dt = 10; % Time steps for the zeroing. In seconds. The original timing is 30
            percentile = 8; % Percetile used for the substraction during zeroing. This is the original timing
            plot_flag = 0; % If 1, plot the different steps of the preprocessing. This is the original timing

            for channel = 1:length(IdChannel)
                if ~isempty(IdChannel{channel})
                    % downsample 10x
                    data405 = data.streams.(STREAM_STORE1{channel}).data(idx_remove:end);
                    data405 = downsample(data405,10);
                    if contains(STREAM_STORE2{channel},'560')
                        data560 = data.streams.(STREAM_STORE2{channel}).data(idx_remove:end);
                        data560 = downsample(data560,10);
                        data465 = [];
                    else
                        data465 = data.streams.(STREAM_STORE2{channel}).data(idx_remove:end);
                        data465 = downsample(data465,10);
                        data560 = [];
                    end

                    figure
                    if show_plot == 0
                        set(gcf,'visible','off')
                    end
                    plot(time_dwn,data405,'k')
                    hold on
                    if ~isempty(data465)
                        plot(time_dwn,data465,'b')
                        legend('405 Channel','465 Channel')
                    elseif ~isempty(data560)
                        plot(time_dwn,data560,'c')
                        legend('405 Channel','560 Channel')
                    end
                    xlim([0 time_dwn(end)])
                    xlabel('Time (s)')
                    ylabel('Fluorescence')
                    title([IdChannel{channel},' Raw channels'])
                    saveas(gcf,[PATH2SAVE,'figures\',IdChannel{channel},' individual raw channels.jpg'])
                    saveas(gcf,[PATH2SAVE,'figures\',IdChannel{channel},' individual raw channels.fig'])

                    %Baseline the 465 using the 405
                    if ~isempty(data465)
                        bls = polyfit(data405, data465, 1);
                        Y_fit = bls(1) .* data405 + bls(2);
                        Y_dF = data465 - Y_fit;
                        %dFF using 405 fit as baseline
                        dFF.(IdChannel{channel}).raw.uncorrected = 100*(Y_dF)./Y_fit;
                    elseif ~isempty(data560)
                        bls = polyfit(data405, data560, 1);
                        Y_fit = bls(1) .* data405 + bls(2);
                        Y_dF = data560 - Y_fit;
                        %dFF using 405 fit as baseline
                        dFF.(IdChannel{channel}).raw.uncorrected = 100*(Y_dF)./Y_fit;
                    end

                    if strcmp(IdChannel{channel},'Psychlight')
                        dFF.(IdChannel{channel}).raw.uncorrected = smooth(dFF.(IdChannel{channel}).raw.uncorrected,102)';
                    end

% %                     % Filter instead of smoothing
%                     ftype = 'low';
%                     n = 2; % 2nd order filter
%                     Wn = 1/((Fs)/2);
%                     %     Wn = 1/((sampling_rate_ds)/2);
%                     [a,b] = butter(n,Wn,ftype);
%                     dFF.(IdChannel{channel}).smoothed.uncorrected = filtfilt(a,b,double(dFF.(IdChannel{channel}).raw.uncorrected));
% %                     %

                    figure
                    if show_plot == 0
                        set(gcf,'visible','off')
                    end
                    plot(time_dwn,Y_dF,'k')
                    hold on
                    plot(time_dwn,dFF.(IdChannel{channel}).raw.uncorrected,'r')
                    xlim([0 time_dwn(end)])
                    xlabel('Time (s)')
                    ylabel('Fluorescence')
                    if ~isempty(data465)
                        legend('Corrected 465','raw dFF')
                    elseif ~isempty(data560)
                        legend('Corrected 560','raw dFF')
                    end
                    title([IdChannel{channel},' Raw channels'])
                    saveas(gcf,[PATH2SAVE,'figures\',IdChannel{channel},' corrected and dFF.jpg'])
                    saveas(gcf,[PATH2SAVE,'figures\',IdChannel{channel},' corrected and dFF.fig'])

                    [dFF.(IdChannel{channel}).unfilt,dFF.(IdChannel{channel}).filt] = preprocessing_Marie(dFF.(IdChannel{channel}).raw.uncorrected,...
                        time_dwn,Fs,detrend_method,zeroing_moving_wnd,zeroing_dt,percentile,filter_flag,zeroing_flag,plot_flag);

                    figure
                    if show_plot == 0
                        set(gcf,'visible','off')
                    end
                    plot(time_dwn,dFF.(IdChannel{channel}).raw.uncorrected+6)
                    hold on
                    plot(time_dwn,dFF.(IdChannel{channel}).filt.uncorrected)
                    plot(time_dwn,dFF.(IdChannel{channel}).filt.corrected-6)
                    xlim([0 time_dwn(end)])
                    xlabel('Time (s)')
                    ylabel('dFF')
                    legend('Raw + 6','High-pass Filtered','Zero corrected (8th percentile) - 6')
                    title([IdChannel{channel},' Processing Steps'])
                    saveas(gcf,[PATH2SAVE,'figures\',IdChannel{channel},' processing steps.jpg'])
                    saveas(gcf,[PATH2SAVE,'figures\',IdChannel{channel},' processing steps.fig'])

                end
            end

            channel_count = 0;
            for channel = 1:length(IdChannel)
                if ~isempty(IdChannel{channel}) 
                    channel_count = channel_count + 1;
                end
            end
            %% Z-score dFF
            for channel = 1:length(IdChannel)
                if ~isempty(IdChannel{channel})
                    fld1 = fieldnames(dFF.(IdChannel{channel}));
                    for flt = 1:length(fld1)
                        fld2 = fieldnames(dFF.(IdChannel{channel}).(fld1{flt}));
                        for correction = 1:length(fld2)
                            zscored_dFF.(IdChannel{channel}).(fld1{flt}).(fld2{correction}) = zscore(dFF.(IdChannel{channel}).(fld1{flt}).(fld2{correction}));
                        end
                    end
                end
            end
            %% Define behavioral data
            filt_mode = 'unfilt'; % unfilt or filt. Filt referes to the highpass.
            zero_mode = 'uncorrected';
            for get_behav_data = 1
                if ~isempty(Behav_data)
                    behav_files = dir(fullfile(path2mouse,'*.xlsx'));
                    if length(behav_files) > 1
                        for i = length(behav_files):-1:1
                            if strcmp(behav_files(i).name(1),'.') || strcmp(behav_files(i).name(1),'~')
                                behav_files(i) = [];
                            end
                        end
                    end
                    if length(behav_files) > 1
                        for i = 1:length(behav_files)
                            if contains(behav_files(i).name,'Reward')
                                reward_file = i;
                            else
                                behav_file = i;
                            end
                        end
                    end
                    if strcmp(Behav_data,'Ymaze')
                        Behavior_code = {'Start','Decision','Left','Right','PostDec'};
                        Color_code = {'c','r','g','b','Y'};
                    elseif strcmp(Behav_data,'Digging')
                        Behavior_code = {'Start','Investcorrect','Investwrong','cinnamon','garlic','PostDec'};%%'IW','DC','DW'
                        Color_code = {'c','r','g','b','y','m'};%
                    elseif strcmp(Behav_data,'Barnes')
                        Behavior_code = {'Start','Decision','Left','Right'};
                        Color_code = {'c','r','g','b'};
                    end

                    %%% For Matlab 2020 or later %%%
                    [Behav_time,Choice_score,Choice_epochs] = ExtractBehavData_Maryam...
                        (behav_files(behav_file),mice_list(m).name,sessions(s).name,time_dwn(1),Behavior_code, Behav_data);
                    [Reward_times] = ExtractRewardData_Maryam...
                        (behav_files(reward_file),mice_list(m).name,sessions(s).name,time_dwn(1),Behavior_code);
                    %%% For Matlab 2019 or earlier %%%
                    %                     [Behav_time,Choice_score,Choice_epochs] = ExtractBehavData_Lauren_v2019...
                    %                         (behav_files,mice_list(m).name,time_dwn(1),concatenation_threshold,duration_threshold);


                    for o = 1:length(Behavior_code)
                        if ~isempty(Choice_epochs.(Behavior_code{o}))
                            Choice_epochs.(Behavior_code{o})(:,2) = Choice_epochs.(Behavior_code{o})(:,2) + 0.99;
                        end
                    end
                    if ~isempty(Reward_times)
                        Reward_times(2) = Reward_times(2) + 0.99;
                        Choice_epochs.Reward = Reward_times;
                    else
                        Choice_epochs.Reward = [];
                    end

                    % Raw dFF
                    figure
                    if show_plot == 0
                        set(gcf,'visible','off')
                    end
                    count = 0;
                    for channel = 1:length(IdChannel)
                        if ~isempty(IdChannel{channel})
                            tmp_max = max(dFF.(IdChannel{channel}).(filt_mode).(zero_mode));
                            tmp_min = min(dFF.(IdChannel{channel}).(filt_mode).(zero_mode));
                            count = count + 1;
                            subplot(channel_count,1,count)
                            hold on
                            for o = 1:length(Behavior_code)
                                for i = 1:size(Choice_epochs.(Behavior_code{o}),1)
                                    h1(i) = area([Choice_epochs.(Behavior_code{o})(i,1) Choice_epochs.(Behavior_code{o})(i,2)],...
                                        [tmp_max tmp_max],tmp_min,'FaceColor',Color_code{o},'EdgeColor','none','FaceAlpha',0.25);
                                end
                            end
                            plot(time_dwn,dFF.(IdChannel{channel}).(filt_mode).(zero_mode),color2plotchannels{channel})
                            if ~isempty(Reward_times)
                                xline(Reward_times(1),'k')
                                xline(Reward_times(2),'k')
                            end
                            xlabel('Time (s)')
                            ylabel('dFF')
                            xlim([time_dwn(1) time_dwn(end)])
                            ylim([tmp_min tmp_max])
                            title(IdChannel{channel})
                        end
                    end
                    if contains(PATH2SAVE,'\')
                        saveas(gcf,[PATH2SAVE,'figures\All_session_Choice_bouts.jpg'])
                        saveas(gcf,[PATH2SAVE,'figures\All_session_Choice_bouts.fig'])

                    else
                        saveas(gcf,[PATH2SAVE,'figures/All_session_Choice_bouts.jpg'])
                        saveas(gcf,[PATH2SAVE,'figures/All_session_Choice_bouts.fig'])

                    end

                    zoom_limits = [0 30]; % Window to plot zoom of the session (In seconds)
                    figure
                    if show_plot == 0
                        set(gcf,'visible','off')
                    end
                    count = 0;
                    for channel = 1:length(IdChannel)
                        if ~isempty(IdChannel{channel})
                            tmp_max = max(dFF.(IdChannel{channel}).(filt_mode).(zero_mode));
                            tmp_min = min(dFF.(IdChannel{channel}).(filt_mode).(zero_mode));
                            count = count + 1;
                            subplot(channel_count,1,count)
                            hold on
                            for o = 1:length(Behavior_code)
                                for i = 1:size(Choice_epochs.(Behavior_code{o}),1)
                                    h1(i) = area([Choice_epochs.(Behavior_code{o})(i,1) Choice_epochs.(Behavior_code{o})(i,2)],...
                                        [tmp_max tmp_max],tmp_min,'FaceColor',Color_code{o},'EdgeColor','none','FaceAlpha',0.25);
                                end
                            end
                            plot(time_dwn,dFF.(IdChannel{channel}).(filt_mode).(zero_mode),color2plotchannels{channel})
                            if ~isempty(Reward_times)
                                xline(Reward_times(1),'k')
                                xline(Reward_times(2),'k')
                            end
                            xlabel('Time (s)')
                            ylabel('dFF')
                            xlim(zoom_limits)
                            ylim([tmp_min tmp_max])
                            title(IdChannel{channel})
                        end
                    end
                    if contains(PATH2SAVE,'\')
                        saveas(gcf,[PATH2SAVE,'figures\All_session_Choice_bouts_zoomIn.jpg'])
                        saveas(gcf,[PATH2SAVE,'figures\All_session_Choice_bouts_zoomIn.fig'])
                    else
                        saveas(gcf,[PATH2SAVE,'figures/All_session_Choice_bouts_zoomIn.jpg'])
                        saveas(gcf,[PATH2SAVE,'figures/All_session_Choice_bouts_zoomIn.fig'])
                    end

                    % Zscored dff 
                    figure
                    if show_plot == 0
                        set(gcf,'visible','off')
                    end
                    count = 0;
                    for channel = 1:length(IdChannel)
                        if ~isempty(IdChannel{channel})
                            tmp_max = max(zscored_dFF.(IdChannel{channel}).(filt_mode).(zero_mode));
                            tmp_min = min(zscored_dFF.(IdChannel{channel}).(filt_mode).(zero_mode));
                            count = count + 1;
                            subplot(channel_count,1,count)
                            hold on
                            for o = 1:length(Behavior_code)
                                for i = 1:size(Choice_epochs.(Behavior_code{o}),1)
                                    h1(i) = area([Choice_epochs.(Behavior_code{o})(i,1) Choice_epochs.(Behavior_code{o})(i,2)],...
                                        [tmp_max tmp_max],tmp_min,'FaceColor',Color_code{o},'EdgeColor','none','FaceAlpha',0.25);
                                end
                            end
                            plot(time_dwn,zscored_dFF.(IdChannel{channel}).(filt_mode).(zero_mode),color2plotchannels{channel})
                            if ~isempty(Reward_times)
                                xline(Reward_times(1),'k')
                                xline(Reward_times(2),'k')
                            end
                            xlabel('Time (s)')
                            ylabel('dFF')
                            xlim([time_dwn(1) time_dwn(end)])
                            ylim([tmp_min tmp_max])
                            title(IdChannel{channel})
                        end
                    end
                    if contains(PATH2SAVE,'\')
                        saveas(gcf,[PATH2SAVE,'figures\All_session_zscored_Choice_bouts.jpg'])
                        saveas(gcf,[PATH2SAVE,'figures\All_session_zscored_Choice_bouts.fig'])
                    else
                        saveas(gcf,[PATH2SAVE,'figures/All_session_zscored_Choice_bouts.jpg'])
                        saveas(gcf,[PATH2SAVE,'figures/All_session_zscored_Choice_bouts.fig'])
                    end

                    zoom_limits = [0 30]; % Window to plot zoom of the session (In seconds)time_dwn(1) 30[time_dwn(1) 10]
                    figure
                    if show_plot == 0
                        set(gcf,'visible','off')
                    end
                    count = 0;
                    for channel = 1:length(IdChannel)
                        if ~isempty(IdChannel{channel})
                            tmp_max = max(zscored_dFF.(IdChannel{channel}).(filt_mode).(zero_mode));
                            tmp_min = min(zscored_dFF.(IdChannel{channel}).(filt_mode).(zero_mode));
                            count = count + 1;
                            subplot(channel_count,1,count)
                            hold on
                            for o = 1:length(Behavior_code)
                                for i = 1:size(Choice_epochs.(Behavior_code{o}),1)
                                    h1(i) = area([Choice_epochs.(Behavior_code{o})(i,1) Choice_epochs.(Behavior_code{o})(i,2)],...
                                        [tmp_max tmp_max],tmp_min,'FaceColor',Color_code{o},'EdgeColor','none','FaceAlpha',0.25);
                                end
                            end
                            plot(time_dwn,zscored_dFF.(IdChannel{channel}).(filt_mode).(zero_mode),color2plotchannels{channel})
                            if ~isempty(Reward_times)
                                xline(Reward_times(1),'k')
                                xline(Reward_times(2),'k')
                            end
                            xlabel('Time (s)')
                            ylabel('dFF')
                            xlim(zoom_limits)
                            ylim([tmp_min tmp_max])
                            title(IdChannel{channel})
                        end
                    end
                    if contains(PATH2SAVE,'\')
                        saveas(gcf,[PATH2SAVE,'figures\All_session_zscored_Choice_bouts_zoomIn.jpg'])
                        saveas(gcf,[PATH2SAVE,'figures\All_session_zscored_Choice_bouts_zoomIn.fig'])
                    else
                        saveas(gcf,[PATH2SAVE,'figures/All_session_zscored_Choice_bouts_zoomIn.jpg'])
                        saveas(gcf,[PATH2SAVE,'figures/All_session_zscored_Choice_bouts_zoomIn.fig'])
                    end
                end
            end
            %% Cross-Correlation
            if channel_count > 1
                Max_lag = 10; % Max time lag in seconds to test for correlation
                corr_type = 'Pearson'; % Pearson or Spearman
                filt_mode = 'raw'; % unfilt or filt. Filt referes to the highpass.
                zero_mode = 'uncorrected'; % uncorrected or corrected. Corrected referes to the 8th percentile correction.
                comb = nchoosek(1:length(IdChannel),2);
                for o = 1:size(comb,1)
                    if ~isempty(IdChannel{comb(o,1)}) && ~isempty(IdChannel{comb(o,2)})
                        signal1 = dFF.(IdChannel{comb(o,1)}).(filt_mode).(zero_mode);
                        signal2 = dFF.(IdChannel{comb(o,2)}).(filt_mode).(zero_mode);

                        save_name = ['Overall correlation between ',filt_mode,' ',...
                            zero_mode,' ',IdChannel{comb(o,1)},'-',IdChannel{comb(o,2)},...
                            ' with ',IdChannel{comb(o,1)},' lag'];

                        [lag2plot,Ovrl_corr.([IdChannel{comb(o,1)},'_',IdChannel{comb(o,2)}])] = ovrl_corr_calculation...
                            (signal1,signal2,time_dwn,Max_lag,corr_type,show_plot,save_plot,PATH2SAVE,save_name);
                    end
                end
            end
            %% Trials analysis
            for trial_analysis = 1
                moving_corr = [];
                %% Trials analysis
                % Define data mode to be used:
                filt_mode = 'unfilt'; % unfilt or filt. Filt referes to the highpass.
                zero_mode = 'uncorrected'; % uncorrected or corrected. Corrected referes to the 8th percentile correction.
                % Define time index for the trial period
                temp_t = time_dwn - time_dwn(1);
                dummie = 1:length(temp_t);
                dummie = dummie(temp_t >= abs(TRANGE(1)));
                idx_Init = dummie(1);
                dummie = 1:length(temp_t);
                dummie = dummie(temp_t >= abs(TRANGE(2)));
                idx_End = dummie(1);

                n = idx_Init + idx_End; % Length of each trial
                t_trials = temp_t(1:n) - temp_t(idx_Init);

                % Compute the aligned results of the trials
                Trial_data = [];
                Behavior_code = fieldnames(Choice_epochs);

                cinnamon_choice = 0;
                if ~isempty(Choice_epochs.cinnamon)
                    cinnamon_choice = 1;
                end
                if (~isempty(intersect(str2double(mice_list(m).name(end)),Animals_cinnamon_correct)) && cinnamon_choice) || ...
                        (isempty(intersect(str2double(mice_list(m).name(end)),Animals_cinnamon_correct)) && cinnamon_choice == 0)
                    Trial_correct = 1;
                    Outcome_name = 'Correct';
                else
                    Trial_correct = 0;
                    Outcome_name = 'Incorrect';
                end
                % Raw Trial data
                for o = 1:length(Behavior_code)
                    if ~isempty(Choice_epochs.(Behavior_code{o}))
                        tmp_name = [Behavior_code{o},'Onset'];
                        Trial_data.(Outcome_name).(tmp_name).idx = ones(size(Choice_epochs.(Behavior_code{o}),1),n)*nan;
                        Trial_data.(Outcome_name).(tmp_name).corr = ones(size(Choice_epochs.(Behavior_code{o}),1),n)*nan;
                        for i = 1:length(IdChannel)
                            if ~isempty(IdChannel{i})
                                Trial_data.(Outcome_name).(tmp_name).dFF.(IdChannel{i}).raw = ones(size(Choice_epochs.(Behavior_code{o}),1),n)*nan;
                                Trial_data.(Outcome_name).(tmp_name).dFF.(IdChannel{i}).baseline_corrected = ones(size(Choice_epochs.(Behavior_code{o}),1),n)*nan;
                            end
                        end
                        for u = 1:size(Choice_epochs.(Behavior_code{o}),1)
                            ix = find(abs(time_dwn-Choice_epochs.(Behavior_code{o})(u,1)) == ...
                                min(abs(time_dwn-Choice_epochs.(Behavior_code{o})(u,1))));
                            if ix+idx_End <= length(time_dwn) && ix-(idx_Init-1) >= 1
                                tmp = ix - (idx_Init-1):ix + idx_End;
                                Trial_data.(Outcome_name).(tmp_name).idx(u,:) = tmp;
                                for i = 1:length(IdChannel)
                                    if ~isempty(IdChannel{i})
                                        Trial_data.(Outcome_name).(tmp_name).dFF.(IdChannel{i}).raw(u,:) = ...
                                            dFF.(IdChannel{i}).(filt_mode).(zero_mode)(tmp);%dFF.(IdChannel{i}).filt(tmp);
                                    end
                                end
                                if ~isempty(moving_corr)
                                    Trial_data.(Outcome_name).(tmp_name).corr(u,:) = moving_corr(tmp);
                                end
                            elseif ix+idx_End > length(time_dwn)
                                tmp = ix - (idx_Init-1):length(time_dwn);
                                Trial_data.(Outcome_name).(tmp_name).idx(u,1:length(tmp)) = tmp;
                                for i = 1:length(IdChannel)
                                    if ~isempty(IdChannel{i})
                                        Trial_data.(Outcome_name).(tmp_name).dFF.(IdChannel{i}).raw(u,1:length(tmp)) = ...
                                            dFF.(IdChannel{i}).(filt_mode).(zero_mode)(tmp);
                                    end
                                end
                                if ~isempty(moving_corr)
                                    Trial_data.(Outcome_name).(tmp_name).corr(u,1:length(tmp)) = moving_corr(tmp);
                                end
                            elseif ix-(idx_Init-1) < 1
                                tmp = 1:(ix + idx_End);
                                Trial_data.(Outcome_name).(tmp_name).idx(u,(n-length(tmp)+1):end) = tmp;
                                for i = 1:length(IdChannel)
                                    if ~isempty(IdChannel{i})
                                        Trial_data.(Outcome_name).(tmp_name).dFF.(IdChannel{i}).raw(u,(n-length(tmp)+1):end) = ...
                                            dFF.(IdChannel{i}).(filt_mode).(zero_mode)(tmp);
                                    end
                                end
                                if ~isempty(moving_corr)
                                    Trial_data.(Outcome_name).(tmp_name).corr(u,(n-length(tmp)+1):end) = moving_corr(tmp);
                                end
                            end
                        end
                        for i = 1:length(IdChannel)
                            if ~isempty(IdChannel{i})
                                %                 Trial_data.TrialOnset.dFF.(IdChannel{i}).baseline_corrected = ...
                                %                     Trial_data.TrialOnset.dFF.(IdChannel{i}).raw - nanmean...
                                %                     (Trial_data.TrialOnset.dFF.(IdChannel{i}).raw...
                                %                     (:,t_trials >= BASELINE_PER(1) & t_trials <= BASELINE_PER(2)),2);
                                Trial_data.(Outcome_name).(tmp_name).dFF.(IdChannel{i}).baseline_corrected = ...
                                    Trial_data.(Outcome_name).(tmp_name).dFF.(IdChannel{i}).raw - nanmedian...
                                    (Trial_data.(Outcome_name).(tmp_name).dFF.(IdChannel{i}).raw...
                                    (:,t_trials >= BASELINE_PER.Onset(1) & ...
                                    t_trials <= BASELINE_PER.Onset(2)),2);
                                Trial_data.(Outcome_name).(tmp_name).dFF.(IdChannel{i}).zscored = ...
                                    (Trial_data.(Outcome_name).(tmp_name).dFF.(IdChannel{i}).raw - mean...
                                    (Trial_data.(Outcome_name).(tmp_name).dFF.(IdChannel{i}).raw...
                                    (:,t_trials >= BASELINE_PER.Onset(1) & t_trials <= BASELINE_PER.Onset(2)),2,'omitnan'))./std(...
                                    Trial_data.(Outcome_name).(tmp_name).dFF.(IdChannel{i}).raw...
                                    (:,t_trials >= BASELINE_PER.Onset(1) & t_trials <= BASELINE_PER.Onset(2)),1,2,'omitnan');
                            end
                        end

                        tmp_name = [Behavior_code{o},'Offset'];
                        Trial_data.(Outcome_name).(tmp_name).idx = ones(size(Choice_epochs.(Behavior_code{o}),1),n)*nan;
                        Trial_data.(Outcome_name).(tmp_name).corr = ones(size(Choice_epochs.(Behavior_code{o}),1),n)*nan;
                        for i = 1:length(IdChannel)
                            if ~isempty(IdChannel{i})
                                Trial_data.(Outcome_name).(tmp_name).dFF.(IdChannel{i}).raw = ones(size(Choice_epochs.(Behavior_code{o}),1),n)*nan;
                                Trial_data.(Outcome_name).(tmp_name).dFF.(IdChannel{i}).baseline_corrected = ones(size(Choice_epochs.(Behavior_code{o}),1),n)*nan;
                            end
                        end
                        for u = 1:size(Choice_epochs.(Behavior_code{o}),1)
                            ix = find(abs(time_dwn-Choice_epochs.(Behavior_code{o})(u,2)) == ...
                                min(abs(time_dwn-Choice_epochs.(Behavior_code{o})(u,2))));
                            if ix+idx_End <= length(time_dwn) && ix-(idx_Init-1) >= 1
                                tmp = ix - (idx_Init-1):ix + idx_End;
                                Trial_data.(Outcome_name).(tmp_name).idx(u,:) = tmp;
                                for i = 1:length(IdChannel)
                                    if ~isempty(IdChannel{i})
                                        Trial_data.(Outcome_name).(tmp_name).dFF.(IdChannel{i}).raw(u,:) = ...
                                            dFF.(IdChannel{i}).(filt_mode).(zero_mode)(tmp);%dFF.(IdChannel{i}).filt(tmp);
                                    end
                                end
                                if ~isempty(moving_corr)
                                    Trial_data.(Outcome_name).(tmp_name).corr(u,:) = moving_corr(tmp);
                                end
                            elseif ix+idx_End > length(time_dwn)
                                tmp = ix - (idx_Init-1):length(time_dwn);
                                Trial_data.(Outcome_name).(tmp_name).idx(u,1:length(tmp)) = tmp;
                                for i = 1:length(IdChannel)
                                    if ~isempty(IdChannel{i})
                                        Trial_data.(Outcome_name).(tmp_name).dFF.(IdChannel{i}).raw(u,1:length(tmp)) = ...
                                            dFF.(IdChannel{i}).(filt_mode).(zero_mode)(tmp);
                                    end
                                end
                                if ~isempty(moving_corr)
                                    Trial_data.(Outcome_name).(tmp_name).corr(u,1:length(tmp)) = moving_corr(tmp);
                                end
                            elseif ix-(idx_Init-1) < 1
                                tmp = 1:(ix + idx_End);
                                Trial_data.(Outcome_name).(tmp_name).idx(u,(n-length(tmp)+1):end) = tmp;
                                for i = 1:length(IdChannel)
                                    if ~isempty(IdChannel{i})
                                        Trial_data.(Outcome_name).(tmp_name).dFF.(IdChannel{i}).raw(u,(n-length(tmp)+1):end) = ...
                                            dFF.(IdChannel{i}).(filt_mode).(zero_mode)(tmp);
                                    end
                                end
                                if ~isempty(moving_corr)
                                    Trial_data.(Outcome_name).(tmp_name).corr(u,(n-length(tmp)+1):end) = moving_corr(tmp);
                                end
                            end
                        end
                        for i = 1:length(IdChannel)
                            if ~isempty(IdChannel{i})
                                %                 Trial_data.TrialOnset.dFF.(IdChannel{i}).baseline_corrected = ...
                                %                     Trial_data.TrialOnset.dFF.(IdChannel{i}).raw - nanmean...
                                %                     (Trial_data.TrialOnset.dFF.(IdChannel{i}).raw...
                                %                     (:,t_trials >= BASELINE_PER(1) & t_trials <= BASELINE_PER(2)),2);
                                Trial_data.(Outcome_name).(tmp_name).dFF.(IdChannel{i}).baseline_corrected = ...
                                    Trial_data.(Outcome_name).(tmp_name).dFF.(IdChannel{i}).raw - nanmedian...
                                    (Trial_data.(Outcome_name).(tmp_name).dFF.(IdChannel{i}).raw...
                                    (:,t_trials >= BASELINE_PER.Offset(1) & ...
                                    t_trials <= BASELINE_PER.Offset(2)),2);
                                Trial_data.(Outcome_name).(tmp_name).dFF.(IdChannel{i}).zscored = ...
                                    (Trial_data.(Outcome_name).(tmp_name).dFF.(IdChannel{i}).raw - mean...
                                    (Trial_data.(Outcome_name).(tmp_name).dFF.(IdChannel{i}).raw...
                                    (:,t_trials >= BASELINE_PER.Onset(1) & t_trials <= BASELINE_PER.Offset(2)),2,'omitnan'))./std(...
                                    Trial_data.(Outcome_name).(tmp_name).dFF.(IdChannel{i}).raw...
                                    (:,t_trials >= BASELINE_PER.Onset(1) & t_trials <= BASELINE_PER.Offset(2)),1,2,'omitnan');
                            end
                        end
                    else
                        tmp_name = [Behavior_code{o},'Onset'];
                        Trial_data.(Outcome_name).(tmp_name) = [];
                        tmp_name = [Behavior_code{o},'Offset'];
                        Trial_data.(Outcome_name).(tmp_name) = [];
                    end
                end
                % Z-score Trial data
                for o = 1:length(Behavior_code)
                    if ~isempty(Choice_epochs.(Behavior_code{o}))
                        tmp_name = [Behavior_code{o},'Onset'];
                        zscored_Trial_data.(Outcome_name).(tmp_name).idx = ones(size(Choice_epochs.(Behavior_code{o}),1),n)*nan;
                        zscored_Trial_data.(Outcome_name).(tmp_name).corr = ones(size(Choice_epochs.(Behavior_code{o}),1),n)*nan;
                        for i = 1:length(IdChannel)
                            if ~isempty(IdChannel{i})
                                zscored_Trial_data.(Outcome_name).(tmp_name).dFF.(IdChannel{i}).raw = ones(size(Choice_epochs.(Behavior_code{o}),1),n)*nan;
                                zscored_Trial_data.(Outcome_name).(tmp_name).dFF.(IdChannel{i}).baseline_corrected = ones(size(Choice_epochs.(Behavior_code{o}),1),n)*nan;
                            end
                        end
                        for u = 1:size(Choice_epochs.(Behavior_code{o}),1)
                            ix = find(abs(time_dwn-Choice_epochs.(Behavior_code{o})(u,1)) == ...
                                min(abs(time_dwn-Choice_epochs.(Behavior_code{o})(u,1))));
                            if ix+idx_End <= length(time_dwn) && ix-(idx_Init-1) >= 1
                                tmp = ix - (idx_Init-1):ix + idx_End;
                                zscored_Trial_data.(Outcome_name).(tmp_name).idx(u,:) = tmp;
                                for i = 1:length(IdChannel)
                                    if ~isempty(IdChannel{i})
                                        zscored_Trial_data.(Outcome_name).(tmp_name).dFF.(IdChannel{i}).raw(u,:) = ...
                                            zscored_dFF.(IdChannel{i}).(filt_mode).(zero_mode)(tmp);%dFF.(IdChannel{i}).filt(tmp);
                                    end
                                end
                                if ~isempty(moving_corr)
                                    zscored_Trial_data.(Outcome_name).(tmp_name).corr(u,:) = moving_corr(tmp);
                                end
                            elseif ix+idx_End > length(time_dwn)
                                tmp = ix - (idx_Init-1):length(time_dwn);
                                zscored_Trial_data.(Outcome_name).(tmp_name).idx(u,1:length(tmp)) = tmp;
                                for i = 1:length(IdChannel)
                                    if ~isempty(IdChannel{i})
                                        zscored_Trial_data.(Outcome_name).(tmp_name).dFF.(IdChannel{i}).raw(u,1:length(tmp)) = ...
                                            zscored_dFF.(IdChannel{i}).(filt_mode).(zero_mode)(tmp);
                                    end
                                end
                                if ~isempty(moving_corr)
                                    zscored_Trial_data.(Outcome_name).(tmp_name).corr(u,1:length(tmp)) = moving_corr(tmp);
                                end
                            elseif ix-(idx_Init-1) < 1
                                tmp = 1:(ix + idx_End);
                                zscored_Trial_data.(Outcome_name).(tmp_name).idx(u,(n-length(tmp)+1):end) = tmp;
                                for i = 1:length(IdChannel)
                                    if ~isempty(IdChannel{i})
                                        zscored_Trial_data.(Outcome_name).(tmp_name).dFF.(IdChannel{i}).raw(u,(n-length(tmp)+1):end) = ...
                                            zscored_dFF.(IdChannel{i}).(filt_mode).(zero_mode)(tmp);
                                    end
                                end
                                if ~isempty(moving_corr)
                                    zscored_Trial_data.(Outcome_name).(tmp_name).corr(u,(n-length(tmp)+1):end) = moving_corr(tmp);
                                end
                            end
                        end
                        for i = 1:length(IdChannel)
                            if ~isempty(IdChannel{i})
                                %                 zscored_Trial_data.TrialOnset.dFF.(IdChannel{i}).baseline_corrected = ...
                                %                     zscored_Trial_data.TrialOnset.dFF.(IdChannel{i}).raw - nanmean...
                                %                     (zscored_Trial_data.TrialOnset.dFF.(IdChannel{i}).raw...
                                %                     (:,t_trials >= BASELINE_PER(1) & t_trials <= BASELINE_PER(2)),2);
                                zscored_Trial_data.(Outcome_name).(tmp_name).dFF.(IdChannel{i}).baseline_corrected = ...
                                    zscored_Trial_data.(Outcome_name).(tmp_name).dFF.(IdChannel{i}).raw - nanmedian...
                                    (zscored_Trial_data.(Outcome_name).(tmp_name).dFF.(IdChannel{i}).raw...
                                    (:,t_trials >= BASELINE_PER.Onset(1) & ...
                                    t_trials <= BASELINE_PER.Onset(2)),2);
                                zscored_Trial_data.(Outcome_name).(tmp_name).dFF.(IdChannel{i}).zscored = ...
                                    (zscored_Trial_data.(Outcome_name).(tmp_name).dFF.(IdChannel{i}).raw - mean...
                                    (zscored_Trial_data.(Outcome_name).(tmp_name).dFF.(IdChannel{i}).raw...
                                    (:,t_trials >= BASELINE_PER.Onset(1) & t_trials <= BASELINE_PER.Onset(2)),2,'omitnan'))./std(...
                                    zscored_Trial_data.(Outcome_name).(tmp_name).dFF.(IdChannel{i}).raw...
                                    (:,t_trials >= BASELINE_PER.Onset(1) & t_trials <= BASELINE_PER.Onset(2)),1,2,'omitnan');
                            end
                        end

                        tmp_name = [Behavior_code{o},'Offset'];
                        zscored_Trial_data.(Outcome_name).(tmp_name).idx = ones(size(Choice_epochs.(Behavior_code{o}),1),n)*nan;
                        zscored_Trial_data.(Outcome_name).(tmp_name).corr = ones(size(Choice_epochs.(Behavior_code{o}),1),n)*nan;
                        for i = 1:length(IdChannel)
                            if ~isempty(IdChannel{i})
                                zscored_Trial_data.(Outcome_name).(tmp_name).dFF.(IdChannel{i}).raw = ones(size(Choice_epochs.(Behavior_code{o}),1),n)*nan;
                                zscored_Trial_data.(Outcome_name).(tmp_name).dFF.(IdChannel{i}).baseline_corrected = ones(size(Choice_epochs.(Behavior_code{o}),1),n)*nan;
                            end
                        end
                        for u = 1:size(Choice_epochs.(Behavior_code{o}),1)
                            ix = find(abs(time_dwn-Choice_epochs.(Behavior_code{o})(u,2)) == ...
                                min(abs(time_dwn-Choice_epochs.(Behavior_code{o})(u,2))));
                            if ix+idx_End <= length(time_dwn) && ix-(idx_Init-1) >= 1
                                tmp = ix - (idx_Init-1):ix + idx_End;
                                zscored_Trial_data.(Outcome_name).(tmp_name).idx(u,:) = tmp;
                                for i = 1:length(IdChannel)
                                    if ~isempty(IdChannel{i})
                                        zscored_Trial_data.(Outcome_name).(tmp_name).dFF.(IdChannel{i}).raw(u,:) = ...
                                            zscored_dFF.(IdChannel{i}).(filt_mode).(zero_mode)(tmp);%dFF.(IdChannel{i}).filt(tmp);
                                    end
                                end
                                if ~isempty(moving_corr)
                                    zscored_Trial_data.(Outcome_name).(tmp_name).corr(u,:) = moving_corr(tmp);
                                end
                            elseif ix+idx_End > length(time_dwn)
                                tmp = ix - (idx_Init-1):length(time_dwn);
                                zscored_Trial_data.(Outcome_name).(tmp_name).idx(u,1:length(tmp)) = tmp;
                                for i = 1:length(IdChannel)
                                    if ~isempty(IdChannel{i})
                                        zscored_Trial_data.(Outcome_name).(tmp_name).dFF.(IdChannel{i}).raw(u,1:length(tmp)) = ...
                                            zscored_dFF.(IdChannel{i}).(filt_mode).(zero_mode)(tmp);
                                    end
                                end
                                if ~isempty(moving_corr)
                                    zscored_Trial_data.(Outcome_name).(tmp_name).corr(u,1:length(tmp)) = moving_corr(tmp);
                                end
                            elseif ix-(idx_Init-1) < 1
                                tmp = 1:(ix + idx_End);
                                zscored_Trial_data.(Outcome_name).(tmp_name).idx(u,(n-length(tmp)+1):end) = tmp;
                                for i = 1:length(IdChannel)
                                    if ~isempty(IdChannel{i})
                                        zscored_Trial_data.(Outcome_name).(tmp_name).dFF.(IdChannel{i}).raw(u,(n-length(tmp)+1):end) = ...
                                            zscored_dFF.(IdChannel{i}).(filt_mode).(zero_mode)(tmp);
                                    end
                                end
                                if ~isempty(moving_corr)
                                    zscored_Trial_data.(Outcome_name).(tmp_name).corr(u,(n-length(tmp)+1):end) = moving_corr(tmp);
                                end
                            end
                        end
                        for i = 1:length(IdChannel)
                            if ~isempty(IdChannel{i})
                                %                 zscored_Trial_data.TrialOnset.dFF.(IdChannel{i}).baseline_corrected = ...
                                %                     zscored_Trial_data.TrialOnset.dFF.(IdChannel{i}).raw - nanmean...
                                %                     (zscored_Trial_data.TrialOnset.dFF.(IdChannel{i}).raw...
                                %                     (:,t_trials >= BASELINE_PER(1) & t_trials <= BASELINE_PER(2)),2);
                                zscored_Trial_data.(Outcome_name).(tmp_name).dFF.(IdChannel{i}).baseline_corrected = ...
                                    zscored_Trial_data.(Outcome_name).(tmp_name).dFF.(IdChannel{i}).raw - nanmedian...
                                    (zscored_Trial_data.(Outcome_name).(tmp_name).dFF.(IdChannel{i}).raw...
                                    (:,t_trials >= BASELINE_PER.Offset(1) & ...
                                    t_trials <= BASELINE_PER.Offset(2)),2);
                                zscored_Trial_data.(Outcome_name).(tmp_name).dFF.(IdChannel{i}).zscored = ...
                                    (zscored_Trial_data.(Outcome_name).(tmp_name).dFF.(IdChannel{i}).raw - mean...
                                    (zscored_Trial_data.(Outcome_name).(tmp_name).dFF.(IdChannel{i}).raw...
                                    (:,t_trials >= BASELINE_PER.Onset(1) & t_trials <= BASELINE_PER.Offset(2)),2,'omitnan'))./std(...
                                    zscored_Trial_data.(Outcome_name).(tmp_name).dFF.(IdChannel{i}).raw...
                                    (:,t_trials >= BASELINE_PER.Onset(1) & t_trials <= BASELINE_PER.Offset(2)),1,2,'omitnan');
                            end
                        end
                    else
                        tmp_name = [Behavior_code{o},'Onset'];
                        zscored_Trial_data.(Outcome_name).(tmp_name) = [];
                        tmp_name = [Behavior_code{o},'Offset'];
                        zscored_Trial_data.(Outcome_name).(tmp_name) = [];
                    end
                end
            end

            %% Add trial to Stim_data:
            for create_Stim_data_structure = 1
                if s == 1
                    % Raw Stim_data
                    Stim_data.All = [];
                    Stim_data.Correct = [];
                    Stim_data.Incorrect = [];

                    Outcome = fieldnames(Trial_data);
                    Outcome = Outcome{1};
                    Stim_data.All = Trial_data.(Outcome);
                    Stim_data.(Outcome) = Trial_data.(Outcome);
                    % Zscored Stim_data
                    zscored_Stim_data.All = [];
                    zscored_Stim_data.Correct = [];
                    zscored_Stim_data.Incorrect = [];

                    Outcome = fieldnames(zscored_Trial_data);
                    Outcome = Outcome{1};
                    zscored_Stim_data.All = zscored_Trial_data.(Outcome);
                    zscored_Stim_data.(Outcome) = zscored_Trial_data.(Outcome);
                else
                    % Raw Stim_data
                    trial_mode = {'raw','baseline_corrected','zscored'};
                    Outcome = fieldnames(Trial_data);
                    Outcome = Outcome{1};
                    if isempty(Stim_data.(Outcome))
                        Stim_data.(Outcome) = Trial_data.(Outcome);
                        fld1 = fieldnames(Stim_data.All);
                        for i = 1:length(fld1)
                            if ~isempty(Stim_data.All.(fld1{i})) && ~isempty(Trial_data.(Outcome).(fld1{i}))
                                fld2 = fieldnames(Stim_data.All.(fld1{i}));
                                for ii = 1:length(fld2)
                                    if isa(Stim_data.All.(fld1{i}).(fld2{ii}),'double')
                                        [Stim_data.All.(fld1{i}).(fld2{ii})] = [Stim_data.All.(fld1{i}).(fld2{ii});...
                                            Trial_data.(Outcome).(fld1{i}).(fld2{ii})];
                                    else
                                        for iii = 1:length(IdChannel)
                                            for iv = 1:length(trial_mode)
                                                [Stim_data.All.(fld1{i}).(fld2{ii}).(IdChannel{iii}).(trial_mode{iv})] = ...
                                                    [Stim_data.All.(fld1{i}).(fld2{ii}).(IdChannel{iii}).(trial_mode{iv});...
                                                    Trial_data.(Outcome).(fld1{i}).(fld2{ii}).(IdChannel{iii}).(trial_mode{iv})];
                                            end
                                        end
                                    end
                                end
                            elseif ~isempty(Trial_data.(Outcome).(fld1{i}))
                                Stim_data.All.(fld1{i}) = Trial_data.(Outcome).(fld1{i});
                            end
                        end
                    else
                        fld1 = fieldnames(Stim_data.All);
                        for i = 1:length(fld1)
                            if ~isempty(Stim_data.All.(fld1{i})) && ~isempty(Trial_data.(Outcome).(fld1{i}))
                                fld2 = fieldnames(Stim_data.All.(fld1{i}));
                                for ii = 1:length(fld2)
                                    if isa(Stim_data.All.(fld1{i}).(fld2{ii}),'double')
                                        [Stim_data.All.(fld1{i}).(fld2{ii})] = [Stim_data.All.(fld1{i}).(fld2{ii});...
                                            Trial_data.(Outcome).(fld1{i}).(fld2{ii})];
                                    else
                                        for iii = 1:length(IdChannel)
                                            for iv = 1:length(trial_mode)
                                                [Stim_data.All.(fld1{i}).(fld2{ii}).(IdChannel{iii}).(trial_mode{iv})] = ...
                                                    [Stim_data.All.(fld1{i}).(fld2{ii}).(IdChannel{iii}).(trial_mode{iv});...
                                                    Trial_data.(Outcome).(fld1{i}).(fld2{ii}).(IdChannel{iii}).(trial_mode{iv})];
                                            end
                                        end
                                    end
                                end
                            elseif ~isempty(Trial_data.(Outcome).(fld1{i}))
                                Stim_data.All.(fld1{i}) = Trial_data.(Outcome).(fld1{i});
                            end
                        end
                        fld1 = fieldnames(Stim_data.(Outcome));
                        for i = 1:length(fld1)
                            if ~isempty(Stim_data.(Outcome).(fld1{i})) && ~isempty(Trial_data.(Outcome).(fld1{i}))
                                fld2 = fieldnames(Stim_data.(Outcome).(fld1{i}));
                                for ii = 1:length(fld2)
                                    if isa(Stim_data.(Outcome).(fld1{i}).(fld2{ii}),'double')
                                        [Stim_data.(Outcome).(fld1{i}).(fld2{ii})] = [Stim_data.(Outcome).(fld1{i}).(fld2{ii});...
                                            Trial_data.(Outcome).(fld1{i}).(fld2{ii})];
                                    else
                                        for iii = 1:length(IdChannel)
                                            for iv = 1:length(trial_mode)
                                                [Stim_data.(Outcome).(fld1{i}).(fld2{ii}).(IdChannel{iii}).(trial_mode{iv})] = ...
                                                    [Stim_data.(Outcome).(fld1{i}).(fld2{ii}).(IdChannel{iii}).(trial_mode{iv});...
                                                    Trial_data.(Outcome).(fld1{i}).(fld2{ii}).(IdChannel{iii}).(trial_mode{iv})];
                                            end
                                        end
                                    end
                                end
                            elseif ~isempty(Trial_data.(Outcome).(fld1{i}))
                                Stim_data.(Outcome).(fld1{i}) = Trial_data.(Outcome).(fld1{i});
                            end
                        end
                    end
                    % Zscored Stim_data
                    Outcome = fieldnames(zscored_Trial_data);
                    Outcome = Outcome{1};
                    if isempty(zscored_Stim_data.(Outcome))
                        zscored_Stim_data.(Outcome) = zscored_Trial_data.(Outcome);
                        fld1 = fieldnames(zscored_Stim_data.All);
                        for i = 1:length(fld1)
                            if ~isempty(zscored_Stim_data.All.(fld1{i})) && ~isempty(zscored_Trial_data.(Outcome).(fld1{i}))
                                fld2 = fieldnames(zscored_Stim_data.All.(fld1{i}));
                                for ii = 1:length(fld2)
                                    if isa(zscored_Stim_data.All.(fld1{i}).(fld2{ii}),'double')
                                        [zscored_Stim_data.All.(fld1{i}).(fld2{ii})] = [zscored_Stim_data.All.(fld1{i}).(fld2{ii});...
                                            zscored_Trial_data.(Outcome).(fld1{i}).(fld2{ii})];
                                    else
                                        for iii = 1:length(IdChannel)
                                            for iv = 1:length(trial_mode)
                                                [zscored_Stim_data.All.(fld1{i}).(fld2{ii}).(IdChannel{iii}).(trial_mode{iv})] = ...
                                                    [zscored_Stim_data.All.(fld1{i}).(fld2{ii}).(IdChannel{iii}).(trial_mode{iv});...
                                                    zscored_Trial_data.(Outcome).(fld1{i}).(fld2{ii}).(IdChannel{iii}).(trial_mode{iv})];
                                            end
                                        end
                                    end
                                end
                            elseif ~isempty(zscored_Trial_data.(Outcome).(fld1{i}))
                                zscored_Stim_data.All.(fld1{i}) = zscored_Trial_data.(Outcome).(fld1{i});
                            end
                        end

                    else
                        fld1 = fieldnames(zscored_Stim_data.All);
                        for i = 1:length(fld1)
                            if ~isempty(zscored_Stim_data.All.(fld1{i})) && ~isempty(zscored_Trial_data.(Outcome).(fld1{i}))
                                fld2 = fieldnames(zscored_Stim_data.All.(fld1{i}));
                                for ii = 1:length(fld2)
                                    if isa(zscored_Stim_data.All.(fld1{i}).(fld2{ii}),'double')
                                        [zscored_Stim_data.All.(fld1{i}).(fld2{ii})] = [zscored_Stim_data.All.(fld1{i}).(fld2{ii});...
                                            zscored_Trial_data.(Outcome).(fld1{i}).(fld2{ii})];
                                    else
                                        for iii = 1:length(IdChannel)
                                            for iv = 1:length(trial_mode)
                                                [zscored_Stim_data.All.(fld1{i}).(fld2{ii}).(IdChannel{iii}).(trial_mode{iv})] = ...
                                                    [zscored_Stim_data.All.(fld1{i}).(fld2{ii}).(IdChannel{iii}).(trial_mode{iv});...
                                                    zscored_Trial_data.(Outcome).(fld1{i}).(fld2{ii}).(IdChannel{iii}).(trial_mode{iv})];
                                            end
                                        end
                                    end
                                end
                            elseif ~isempty(zscored_Trial_data.(Outcome).(fld1{i}))
                                zscored_Stim_data.All.(fld1{i}) = zscored_Trial_data.(Outcome).(fld1{i});
                            end
                        end

                        fld1 = fieldnames(zscored_Stim_data.(Outcome));
                        for i = 1:length(fld1)
                            if ~isempty(zscored_Stim_data.(Outcome).(fld1{i})) && ~isempty(zscored_Trial_data.(Outcome).(fld1{i}))
                                fld2 = fieldnames(zscored_Stim_data.(Outcome).(fld1{i}));
                                for ii = 1:length(fld2)
                                    if isa(zscored_Stim_data.(Outcome).(fld1{i}).(fld2{ii}),'double')
                                        [zscored_Stim_data.(Outcome).(fld1{i}).(fld2{ii})] = [zscored_Stim_data.(Outcome).(fld1{i}).(fld2{ii});...
                                            zscored_Trial_data.(Outcome).(fld1{i}).(fld2{ii})];
                                    else
                                        for iii = 1:length(IdChannel)
                                            for iv = 1:length(trial_mode)
                                                [zscored_Stim_data.(Outcome).(fld1{i}).(fld2{ii}).(IdChannel{iii}).(trial_mode{iv})] = ...
                                                    [zscored_Stim_data.(Outcome).(fld1{i}).(fld2{ii}).(IdChannel{iii}).(trial_mode{iv});...
                                                    zscored_Trial_data.(Outcome).(fld1{i}).(fld2{ii}).(IdChannel{iii}).(trial_mode{iv})];
                                            end
                                        end
                                    end
                                end
                            elseif ~isempty(zscored_Trial_data.(Outcome).(fld1{i}))
                                zscored_Stim_data.(Outcome).(fld1{i}) = zscored_Trial_data.(Outcome).(fld1{i});
                            end
                        end
                    end
                end
            end
            %% Compute measurements for each of the behavioral states:
            filt_mode = 'unfilt'; % unfilt or filt. Filt referes to the highpass.
            zero_mode = 'uncorrected'; % uncorrected or corrected. Corrected referes to the 8th percentile correction.

            threshold_criterion = 'mean';% 'median'
            n_std = 1;

            % Raw dFF
            for o = 1:length(IdChannel)
                if ~isempty(IdChannel{o})
                    tmp = dFF.(IdChannel{o}).(filt_mode).(zero_mode);
                    [~,~,~,test_prom] = findpeaks(tmp);
                    if strcmp(threshold_criterion,'mean')
                        threshold = mean(test_prom)+(n_std*std(test_prom));
                    else
                        threshold = median(test_prom)+(n_std*mad(test_prom));
                    end
                    [tmp_peaks,tmp_locs,tmp_w,tmp_prom] = findpeaks(tmp,'MinPeakProminence',threshold,'Annotate','extents');

                    Trial_Overall_AUC.(Outcome_name).(IdChannel{o}) = trapz(tmp)./(time_dwn(end)-time_dwn(1));
                    Trial_Overall_meandFF.(Outcome_name).(IdChannel{o}) = mean(tmp,'omitnan');

                    for behav = 1:length(Behavior_code)
                        if ~isempty(Choice_epochs.(Behavior_code{behav}))
                            Trial_Measurements.(Outcome_name).norm_AUC.(IdChannel{o}).(Behavior_code{behav}) = ...
                                trapz(tmp(time_dwn>= Choice_epochs.(Behavior_code{behav})(1) & ...
                                time_dwn < Choice_epochs.(Behavior_code{behav})(2)))./(diff(Choice_epochs.(Behavior_code{behav})));
                            Trial_Measurements.(Outcome_name).mean_dFF.(IdChannel{o}).(Behavior_code{behav}) = ...
                                mean(tmp(time_dwn>= Choice_epochs.(Behavior_code{behav})(1) & ...
                                time_dwn < Choice_epochs.(Behavior_code{behav})(2)),'omitnan');

                            Trial_Measurements.(Outcome_name).peaks.(IdChannel{o}).(Behavior_code{behav}) = ...
                                tmp_peaks(time_dwn(tmp_locs)>= Choice_epochs.(Behavior_code{behav})(1) & ...
                                time_dwn(tmp_locs)< Choice_epochs.(Behavior_code{behav})(2));
                            Trial_Measurements.(Outcome_name).peak_width.(IdChannel{o}).(Behavior_code{behav}) = ...
                                tmp_w(time_dwn(tmp_locs)>= Choice_epochs.(Behavior_code{behav})(1) & ...
                                time_dwn(tmp_locs)< Choice_epochs.(Behavior_code{behav})(2));
                            Trial_Measurements.(Outcome_name).peak_prom.(IdChannel{o}).(Behavior_code{behav}) = ...
                                tmp_prom(time_dwn(tmp_locs)>= Choice_epochs.(Behavior_code{behav})(1) & ...
                                time_dwn(tmp_locs)< Choice_epochs.(Behavior_code{behav})(2));

                            Trial_Measurements.(Outcome_name).peaks_rate.(IdChannel{o}).(Behavior_code{behav}) = ...
                                length(tmp_peaks(time_dwn(tmp_locs)>= Choice_epochs.(Behavior_code{behav})(1) & ...
                                time_dwn(tmp_locs)< Choice_epochs.(Behavior_code{behav})(2)))/...
                                (diff(Choice_epochs.(Behavior_code{behav})));
                            Trial_Measurements.(Outcome_name).mean_peaks_val.(IdChannel{o}).(Behavior_code{behav}) = ...
                                mean(tmp_peaks(time_dwn(tmp_locs)>= Choice_epochs.(Behavior_code{behav})(1) & ...
                                time_dwn(tmp_locs)< Choice_epochs.(Behavior_code{behav})(2)),'omitnan');
                            Trial_Measurements.(Outcome_name).mean_peaks_width.(IdChannel{o}).(Behavior_code{behav}) = ...
                                mean(tmp_w(time_dwn(tmp_locs)>= Choice_epochs.(Behavior_code{behav})(1) & ...
                                time_dwn(tmp_locs)< Choice_epochs.(Behavior_code{behav})(2)),'omitnan');
                            Trial_Measurements.(Outcome_name).mean_peaks_prom.(IdChannel{o}).(Behavior_code{behav}) = ...
                                mean(tmp_prom(time_dwn(tmp_locs)>= Choice_epochs.(Behavior_code{behav})(1) & ...
                                time_dwn(tmp_locs)< Choice_epochs.(Behavior_code{behav})(2)),'omitnan');
                        else
                            Trial_Measurements.(Outcome_name).norm_AUC.(IdChannel{o}).(Behavior_code{behav}) = [];
                            Trial_Measurements.(Outcome_name).mean_dFF.(IdChannel{o}).(Behavior_code{behav}) = [];

                            Trial_Measurements.(Outcome_name).peaks.(IdChannel{o}).(Behavior_code{behav}) = [];
                            Trial_Measurements.(Outcome_name).peak_width.(IdChannel{o}).(Behavior_code{behav}) = [];
                            Trial_Measurements.(Outcome_name).peak_prom.(IdChannel{o}).(Behavior_code{behav}) = [];

                            Trial_Measurements.(Outcome_name).peaks_rate.(IdChannel{o}).(Behavior_code{behav}) = [];
                            Trial_Measurements.(Outcome_name).mean_peaks_val.(IdChannel{o}).(Behavior_code{behav}) = [];
                            Trial_Measurements.(Outcome_name).mean_peaks_width.(IdChannel{o}).(Behavior_code{behav}) = [];
                            Trial_Measurements.(Outcome_name).mean_peaks_prom.(IdChannel{o}).(Behavior_code{behav}) = [];
                        end
                    end
                end
            end
            % Zscored dFF
            for o = 1:length(IdChannel)
                if ~isempty(IdChannel{o})
                    tmp = zscored_dFF.(IdChannel{o}).(filt_mode).(zero_mode);
                    [~,~,~,test_prom] = findpeaks(tmp);
                    if strcmp(threshold_criterion,'mean')
                        threshold = mean(test_prom)+(n_std*std(test_prom));
                    else
                        threshold = median(test_prom)+(n_std*mad(test_prom));
                    end
                    [tmp_peaks,tmp_locs,tmp_w,tmp_prom] = findpeaks(tmp,'MinPeakProminence',threshold,'Annotate','extents');

                    zscored_Trial_Overall_AUC.(Outcome_name).(IdChannel{o}) = trapz(tmp)./(time_dwn(end)-time_dwn(1));
                    zscored_Trial_Overall_meandFF.(Outcome_name).(IdChannel{o}) = mean(tmp,'omitnan');

                    for behav = 1:length(Behavior_code)
                        if ~isempty(Choice_epochs.(Behavior_code{behav}))
                            zscored_Trial_Measurements.(Outcome_name).norm_AUC.(IdChannel{o}).(Behavior_code{behav}) = ...
                                trapz(tmp(time_dwn>= Choice_epochs.(Behavior_code{behav})(1) & ...
                                time_dwn < Choice_epochs.(Behavior_code{behav})(2)))./(diff(Choice_epochs.(Behavior_code{behav})));
                            zscored_Trial_Measurements.(Outcome_name).mean_dFF.(IdChannel{o}).(Behavior_code{behav}) = ...
                                mean(tmp(time_dwn>= Choice_epochs.(Behavior_code{behav})(1) & ...
                                time_dwn < Choice_epochs.(Behavior_code{behav})(2)),'omitnan');

                            zscored_Trial_Measurements.(Outcome_name).peaks.(IdChannel{o}).(Behavior_code{behav}) = ...
                                tmp_peaks(time_dwn(tmp_locs)>= Choice_epochs.(Behavior_code{behav})(1) & ...
                                time_dwn(tmp_locs)< Choice_epochs.(Behavior_code{behav})(2));
                            zscored_Trial_Measurements.(Outcome_name).peak_width.(IdChannel{o}).(Behavior_code{behav}) = ...
                                tmp_w(time_dwn(tmp_locs)>= Choice_epochs.(Behavior_code{behav})(1) & ...
                                time_dwn(tmp_locs)< Choice_epochs.(Behavior_code{behav})(2));
                            zscored_Trial_Measurements.(Outcome_name).peak_prom.(IdChannel{o}).(Behavior_code{behav}) = ...
                                tmp_prom(time_dwn(tmp_locs)>= Choice_epochs.(Behavior_code{behav})(1) & ...
                                time_dwn(tmp_locs)< Choice_epochs.(Behavior_code{behav})(2));

                            zscored_Trial_Measurements.(Outcome_name).peaks_rate.(IdChannel{o}).(Behavior_code{behav}) = ...
                                length(tmp_peaks(time_dwn(tmp_locs)>= Choice_epochs.(Behavior_code{behav})(1) & ...
                                time_dwn(tmp_locs)< Choice_epochs.(Behavior_code{behav})(2)))/...
                                (diff(Choice_epochs.(Behavior_code{behav})));
                            zscored_Trial_Measurements.(Outcome_name).mean_peaks_val.(IdChannel{o}).(Behavior_code{behav}) = ...
                                mean(tmp_peaks(time_dwn(tmp_locs)>= Choice_epochs.(Behavior_code{behav})(1) & ...
                                time_dwn(tmp_locs)< Choice_epochs.(Behavior_code{behav})(2)),'omitnan');
                            zscored_Trial_Measurements.(Outcome_name).mean_peaks_width.(IdChannel{o}).(Behavior_code{behav}) = ...
                                mean(tmp_w(time_dwn(tmp_locs)>= Choice_epochs.(Behavior_code{behav})(1) & ...
                                time_dwn(tmp_locs)< Choice_epochs.(Behavior_code{behav})(2)),'omitnan');
                            zscored_Trial_Measurements.(Outcome_name).mean_peaks_prom.(IdChannel{o}).(Behavior_code{behav}) = ...
                                mean(tmp_prom(time_dwn(tmp_locs)>= Choice_epochs.(Behavior_code{behav})(1) & ...
                                time_dwn(tmp_locs)< Choice_epochs.(Behavior_code{behav})(2)),'omitnan');
                        else
                            zscored_Trial_Measurements.(Outcome_name).norm_AUC.(IdChannel{o}).(Behavior_code{behav}) = [];
                            zscored_Trial_Measurements.(Outcome_name).mean_dFF.(IdChannel{o}).(Behavior_code{behav}) = [];

                            zscored_Trial_Measurements.(Outcome_name).peaks.(IdChannel{o}).(Behavior_code{behav}) = [];
                            zscored_Trial_Measurements.(Outcome_name).peak_width.(IdChannel{o}).(Behavior_code{behav}) = [];
                            zscored_Trial_Measurements.(Outcome_name).peak_prom.(IdChannel{o}).(Behavior_code{behav}) = [];

                            zscored_Trial_Measurements.(Outcome_name).peaks_rate.(IdChannel{o}).(Behavior_code{behav}) = [];
                            zscored_Trial_Measurements.(Outcome_name).mean_peaks_val.(IdChannel{o}).(Behavior_code{behav}) = [];
                            zscored_Trial_Measurements.(Outcome_name).mean_peaks_width.(IdChannel{o}).(Behavior_code{behav}) = [];
                            zscored_Trial_Measurements.(Outcome_name).mean_peaks_prom.(IdChannel{o}).(Behavior_code{behav}) = [];
                        end
                    end
                end
            end
            %% Add Trial_Measurements to Measurements:
            for create_Measurements_structure = 1
                % Raw dFF
                if s == 1
                    Measurements.Correct = [];
                    Measurements.Incorrect = [];

                    fld1 = fieldnames(Trial_Measurements.(Outcome_name));
                    for i = 1:length(fld1)
                        for channel = 1:length(IdChannel)
                            fld2 = fieldnames(Trial_Measurements.(Outcome_name).(fld1{i}).(IdChannel{channel}));
                            for ii = 1:length(fld2)
                                if length(Trial_Measurements.(Outcome_name).(fld1{i}).(IdChannel{channel}).(fld2{ii})) == 1
                                    Measurements.(Outcome_name).(fld1{i}).(IdChannel{channel}).(fld2{ii}) = ...
                                        Trial_Measurements.(Outcome_name).(fld1{i}).(IdChannel{channel}).(fld2{ii});
                                else
                                    Measurements.(Outcome_name).(fld1{i}).(IdChannel{channel}).(fld2{ii}) = ...
                                        mean(Trial_Measurements.(Outcome_name).(fld1{i}).(IdChannel{channel}).(fld2{ii}),'omitnan');
                                end
                            end
                        end
                    end

                    Overall_AUC.Correct = [];
                    Overall_meandFF.Correct = [];
                    Overall_AUC.Incorrect = [];
                    Overall_meandFF.Incorrect = [];

                    Overall_AUC.(Outcome_name) = Trial_Overall_AUC.(Outcome_name);
                    Overall_meandFF.(Outcome_name) = Trial_Overall_meandFF.(Outcome_name);
                else
                    Outcome = fieldnames(Trial_data);
                    Outcome = Outcome{1};
                    if isempty(Measurements.(Outcome_name))
                        fld1 = fieldnames(Trial_Measurements.(Outcome_name));
                        for i = 1:length(fld1)
                            for channel = 1:length(IdChannel)
                                fld2 = fieldnames(Trial_Measurements.(Outcome_name).(fld1{i}).(IdChannel{channel}));
                                for ii = 1:length(fld2)
                                    if length(Trial_Measurements.(Outcome_name).(fld1{i}).(IdChannel{channel}).(fld2{ii})) == 1
                                        Measurements.(Outcome_name).(fld1{i}).(IdChannel{channel}).(fld2{ii}) = ...
                                            Trial_Measurements.(Outcome_name).(fld1{i}).(IdChannel{channel}).(fld2{ii});
                                    else
                                        Measurements.(Outcome_name).(fld1{i}).(IdChannel{channel}).(fld2{ii}) = ...
                                            mean(Trial_Measurements.(Outcome_name).(fld1{i}).(IdChannel{channel}).(fld2{ii}),'omitnan');
                                    end
                                end
                            end
                        end
                    else
                        fld1 = fieldnames(Trial_Measurements.(Outcome_name));
                        for i = 1:length(fld1)
                            if ~isempty(Measurements.(Outcome_name).(fld1{i})) && ~isempty(Trial_Measurements.(Outcome_name).(fld1{i}))
                                % fld2 = fieldnames(Measurements.(Outcome_name).(fld1{i}));
                                for channel = 1:length(IdChannel)
                                    fld2 = fieldnames(Measurements.(Outcome_name).(fld1{i}).(IdChannel{channel}));
                                    for ii = 1:length(fld2)
                                        if length(Trial_Measurements.(Outcome_name).(fld1{i}).(IdChannel{channel}).(fld2{ii})) == 1
                                            [Measurements.(Outcome_name).(fld1{i}).(IdChannel{channel}).(fld2{ii})] = ...
                                                [Measurements.(Outcome_name).(fld1{i}).(IdChannel{channel}).(fld2{ii});...
                                                Trial_Measurements.(Outcome_name).(fld1{i}).(IdChannel{channel}).(fld2{ii})];
                                        else
                                            [Measurements.(Outcome_name).(fld1{i}).(IdChannel{channel}).(fld2{ii})] = ...
                                                [Measurements.(Outcome_name).(fld1{i}).(IdChannel{channel}).(fld2{ii});...
                                                mean(Trial_Measurements.(Outcome_name).(fld1{i}).(IdChannel{channel}).(fld2{ii}),'omitnan')];
                                        end
                                    end
                                end
                            elseif ~isempty(Trial_Measurements.(Outcome_name).(fld1{i}))
                                for channel = 1:length(IdChannel)
                                    fld2 = fieldnames(Trial_Measurements.(Outcome_name).(fld1{i}).(IdChannel{channel}));
                                    for ii = 1:length(fld2)
                                        if length(Trial_Measurements.(Outcome_name).(fld1{i}).(IdChannel{channel}).(fld2{ii})) == 1
                                            Measurements.(Outcome_name).(fld1{i}).(IdChannel{channel}).(fld2{ii}) = ...
                                                Trial_Measurements.(Outcome_name).(fld1{i}).(IdChannel{channel}).(fld2{ii});
                                        else
                                            Measurements.(Outcome_name).(fld1{i}).(IdChannel{channel}).(fld2{ii}) = ...
                                                mean(Trial_Measurements.(Outcome_name).(fld1{i}).(IdChannel{channel}).(fld2{ii}),'omitnan');
                                        end
                                    end
                                end
                            end
                        end
                    end

                    if isempty(Overall_AUC.(Outcome_name))
                        Overall_AUC.(Outcome_name) = Trial_Overall_AUC.(Outcome_name);
                        Overall_meandFF.(Outcome_name) = Trial_Overall_meandFF.(Outcome_name);
                    else
                        for channel = 1:length(IdChannel)
                            [Overall_AUC.(Outcome_name).(IdChannel{channel})] = [Overall_AUC.(Outcome_name).(IdChannel{channel});...
                               Trial_Overall_AUC.(Outcome_name).(IdChannel{channel})];
                            [Overall_meandFF.(Outcome_name).(IdChannel{channel})] = [Overall_meandFF.(Outcome_name).(IdChannel{channel});...
                               Trial_Overall_meandFF.(Outcome_name).(IdChannel{channel})];
                        end
                    end
                end

                % Zscored dFF
                if s == 1
                    zscored_Measurements.Correct = [];
                    zscored_Measurements.Incorrect = [];

                    fld1 = fieldnames(zscored_Trial_Measurements.(Outcome_name));
                    for i = 1:length(fld1)
                        for channel = 1:length(IdChannel)
                            fld2 = fieldnames(zscored_Trial_Measurements.(Outcome_name).(fld1{i}).(IdChannel{channel}));
                            for ii = 1:length(fld2)
                                if length(zscored_Trial_Measurements.(Outcome_name).(fld1{i}).(IdChannel{channel}).(fld2{ii})) == 1
                                    zscored_Measurements.(Outcome_name).(fld1{i}).(IdChannel{channel}).(fld2{ii}) = ...
                                        zscored_Trial_Measurements.(Outcome_name).(fld1{i}).(IdChannel{channel}).(fld2{ii});
                                else
                                    zscored_Measurements.(Outcome_name).(fld1{i}).(IdChannel{channel}).(fld2{ii}) = ...
                                        mean(zscored_Trial_Measurements.(Outcome_name).(fld1{i}).(IdChannel{channel}).(fld2{ii}),'omitnan');
                                end
                            end
                        end
                    end

                    zscored_Overall_AUC.Correct = [];
                    zscored_Overall_meandFF.Correct = [];
                    zscored_Overall_AUC.Incorrect = [];
                    zscored_Overall_meandFF.Incorrect = [];

                    zscored_Overall_AUC.(Outcome_name) = zscored_Trial_Overall_AUC.(Outcome_name);
                    zscored_Overall_meandFF.(Outcome_name) = zscored_Trial_Overall_meandFF.(Outcome_name);
                else
                    Outcome = fieldnames(zscored_Trial_data);
                    Outcome = Outcome{1};
                    if isempty(zscored_Measurements.(Outcome_name))
                        fld1 = fieldnames(zscored_Trial_Measurements.(Outcome_name));
                        for i = 1:length(fld1)
                            for channel = 1:length(IdChannel)
                                fld2 = fieldnames(zscored_Trial_Measurements.(Outcome_name).(fld1{i}).(IdChannel{channel}));
                                for ii = 1:length(fld2)
                                    if length(zscored_Trial_Measurements.(Outcome_name).(fld1{i}).(IdChannel{channel}).(fld2{ii})) == 1
                                        zscored_Measurements.(Outcome_name).(fld1{i}).(IdChannel{channel}).(fld2{ii}) = ...
                                            zscored_Trial_Measurements.(Outcome_name).(fld1{i}).(IdChannel{channel}).(fld2{ii});
                                    else
                                        zscored_Measurements.(Outcome_name).(fld1{i}).(IdChannel{channel}).(fld2{ii}) = ...
                                            mean(zscored_Trial_Measurements.(Outcome_name).(fld1{i}).(IdChannel{channel}).(fld2{ii}),'omitnan');
                                    end
                                end
                            end
                        end
                    else
                        fld1 = fieldnames(zscored_Trial_Measurements.(Outcome_name));
                        for i = 1:length(fld1)
                            if ~isempty(zscored_Measurements.(Outcome_name).(fld1{i})) && ~isempty(zscored_Trial_Measurements.(Outcome_name).(fld1{i}))
                                % fld2 = fieldnames(Measurements.(Outcome_name).(fld1{i}));
                                for channel = 1:length(IdChannel)
                                    fld2 = fieldnames(zscored_Measurements.(Outcome_name).(fld1{i}).(IdChannel{channel}));
                                    for ii = 1:length(fld2)
                                        if length(zscored_Trial_Measurements.(Outcome_name).(fld1{i}).(IdChannel{channel}).(fld2{ii})) == 1
                                            [zscored_Measurements.(Outcome_name).(fld1{i}).(IdChannel{channel}).(fld2{ii})] = ...
                                                [zscored_Measurements.(Outcome_name).(fld1{i}).(IdChannel{channel}).(fld2{ii});...
                                                zscored_Trial_Measurements.(Outcome_name).(fld1{i}).(IdChannel{channel}).(fld2{ii})];
                                        else
                                            [zscored_Measurements.(Outcome_name).(fld1{i}).(IdChannel{channel}).(fld2{ii})] = ...
                                                [zscored_Measurements.(Outcome_name).(fld1{i}).(IdChannel{channel}).(fld2{ii});...
                                                mean(zscored_Trial_Measurements.(Outcome_name).(fld1{i}).(IdChannel{channel}).(fld2{ii}),'omitnan')];
                                        end
                                    end
                                end
                            elseif ~isempty(zscored_Trial_Measurements.(Outcome_name).(fld1{i}))
                                for channel = 1:length(IdChannel)
                                    fld2 = fieldnames(zscored_Trial_Measurements.(Outcome_name).(fld1{i}).(IdChannel{channel}));
                                    for ii = 1:length(fld2)
                                        if length(zscored_Trial_Measurements.(Outcome_name).(fld1{i}).(IdChannel{channel}).(fld2{ii})) == 1
                                            zscored_Measurements.(Outcome_name).(fld1{i}).(IdChannel{channel}).(fld2{ii}) = ...
                                                zscored_Trial_Measurements.(Outcome_name).(fld1{i}).(IdChannel{channel}).(fld2{ii});
                                        else
                                            zscored_Measurements.(Outcome_name).(fld1{i}).(IdChannel{channel}).(fld2{ii}) = ...
                                                mean(zscored_Trial_Measurements.(Outcome_name).(fld1{i}).(IdChannel{channel}).(fld2{ii}),'omitnan');
                                        end
                                    end
                                end
                            end
                        end
                    end

                    if isempty(zscored_Overall_AUC.(Outcome_name))
                        zscored_Overall_AUC.(Outcome_name) = zscored_Trial_Overall_AUC.(Outcome_name);
                        zscored_Overall_meandFF.(Outcome_name) = zscored_Trial_Overall_meandFF.(Outcome_name);
                    else
                        for channel = 1:length(IdChannel)
                            [zscored_Overall_AUC.(Outcome_name).(IdChannel{channel})] = [zscored_Overall_AUC.(Outcome_name).(IdChannel{channel});...
                               zscored_Trial_Overall_AUC.(Outcome_name).(IdChannel{channel})];
                            [zscored_Overall_meandFF.(Outcome_name).(IdChannel{channel})] = [zscored_Overall_meandFF.(Outcome_name).(IdChannel{channel});...
                               zscored_Trial_Overall_meandFF.(Outcome_name).(IdChannel{channel})];
                        end
                    end
                end
            end

            %% Save
            if done == 0 || overwrite == 1
                Params.MICEId = mice_list(m).name;
                Params.min2remove = min2remove; % Minutes to remove from the beginning of the recording
                Params.TRANGE = TRANGE; % window size [start time relative to epoc onset, window duration]
                Params.BASELINE_PER = BASELINE_PER; % baseline period within our window
                Params.Limits_Selected_Epochs = Limits_Selected_Epochs; % Duration criteria for selected Choice epochs (in seconds)
                Params.Epoch2analyze = Epoch2analyze; %Time in seconds to get measurements aligned to behavior
                % Params.concatenation_threshold = concatenation_threshold;
                % Params.duration_threshold = duration_threshold;

                % save([PATH2SAVE,'Trial_analysis.mat'],'time_dwn','dFF','lag2plot','Ovrl_corr',...
                %     'moving_corr','Stim_data','Ovrl_Measurements','Measurements','t_trials',...
                %     'Behav_time','Choice_score','Choice_epochs','shock_evoked_response','Params')
                save([PATH2SAVE,'Trial_analysis.mat'],'time_dwn','dFF','zscored_dFF','lag2plot','Ovrl_corr',...
                    'moving_corr','Trial_data','zscored_Trial_data','t_trials','Trial_Measurements','zscored_Trial_Measurements','Behav_time',...
                    'Trial_Overall_AUC','zscored_Trial_Overall_AUC','Trial_Overall_meandFF','zscored_Trial_Overall_meandFF','Choice_score','Choice_epochs','Params')
            end
        end
        clear Trial_data zscored_Trial_data Trial_Measurements zscored_Trial_Measurements Trial_Overall_AUC zscored_Trial_Overall_AUC Trial_Overall_meandFF zscored_Trial_Overall_meandFF
    end  
    
    %% Calculate AUC and absolute peak within 3 seconds after alignment
    
    window2analyze = [0 5];
    fld1 = fieldnames(Stim_data);
    for i = 1:length(fld1)
        if ~isempty(Stim_data.(fld1{i}))
            fld2 = fieldnames(Stim_data.(fld1{i}));
            for ii = 1:length(fld2)
                if ~isempty(Stim_data.(fld1{i}).(fld2{ii}))
                    for iii = 1:length(IdChannel)
                        if ~isempty(IdChannel{iii})
                            fld3 = fieldnames(Stim_data.(fld1{i}).(fld2{ii}).dFF.(IdChannel{iii}));
                            for iv = 1:length(fld3)
                                % For raw
                                tmp = Stim_data.(fld1{i}).(fld2{ii}).dFF.(IdChannel{iii}).(fld3{iv})...
                                    (:,t_trials >= window2analyze(1) & t_trials < window2analyze(2));
                                alignment_measurements.(fld1{i}).(fld2{ii}).(IdChannel{iii}).(fld3{iv}).AUC = ...
                                    trapz(tmp,2);

                                [~,idx] = max(abs(tmp),[],2);
                                alignment_measurements.(fld1{i}).(fld2{ii}).(IdChannel{iii}).(fld3{iv}).peak = tmp(idx);

                                alignment_zscored_measurements.(fld1{i}).(fld2{ii}).(IdChannel{iii}).(fld3{iv}).AUC = ...
                                    trapz(tmp,2);
                                % For zscored
                                tmp = zscored_Stim_data.(fld1{i}).(fld2{ii}).dFF.(IdChannel{iii}).(fld3{iv})...
                                    (:,t_trials >= window2analyze(1) & t_trials < window2analyze(2));
                                alignment_zscored_measurements.(fld1{i}).(fld2{ii}).(IdChannel{iii}).(fld3{iv}).AUC = ...
                                    trapz(tmp,2);

                                [~,idx] = max(abs(tmp),[],2);
                                alignment_zscored_measurements.(fld1{i}).(fld2{ii}).(IdChannel{iii}).(fld3{iv}).peak = tmp(idx);
                            end
                        end
                    end
                else
                    alignment_measurements.(fld1{i}).(fld2{ii}) = [];
                    alignment_zscored_measurements.(fld1{i}).(fld2{ii}) = [];
                end
            end
        end
    end

    %% Plot the results
    PATH2SAVEALLTRIALS = path2save_mouse;
    % Raw
    if ~isempty(Stim_data)
        trial_mode = {'raw','baseline_corrected','zscored'}; % raw or baseline_corrected.
        Outcome = fieldnames(Stim_data);
        for outcome = 2:length(Outcome)
            if ~isempty(Stim_data.(Outcome{outcome}))
                Condition = fieldnames(Stim_data.(Outcome{outcome}));%{'LeverExtension','DipperOn','Reward'};
                avg_mode = 1; % 1, use mean and SEM; 2, use median and MAD
                order_crit = 1;
                for cond = 1:length(Condition)
                    if ~isempty(Stim_data.(Outcome{outcome}).(Condition{cond}))
                        for t_mode = 1:length(trial_mode)
                            plot_trials_data_Maryam(Stim_data.(Outcome{outcome}),t_trials,Outcome{outcome},Condition{cond},IdChannel,trial_mode{t_mode},avg_mode,color2plotchannels,order_crit,limits2plot,LimitXaxis,Color_scale,show_plot,save_plot,PATH2SAVEALLTRIALS)
                        end
                    end
                end
                close all
            end
        end
    end

    %Zscored
    if ~isempty(zscored_Stim_data)
        trial_mode = {'raw','baseline_corrected','zscored'}; % raw or baseline_corrected.
        Outcome = fieldnames(zscored_Stim_data);
        for outcome = 2:length(Outcome)
            if ~isempty(zscored_Stim_data.(Outcome{outcome}))
                Condition = fieldnames(zscored_Stim_data.(Outcome{outcome}));%{'LeverExtension','DipperOn','Reward'};
                avg_mode = 1; % 1, use mean and SEM; 2, use median and MAD
                order_crit = 0; % 0, order of the trials; 1, order based on average intensity from 0 to 10s
                for cond = 1:length(Condition)
                    if ~isempty(zscored_Stim_data.(Outcome{outcome}).(Condition{cond}))
                        for t_mode = 1:length(trial_mode)
                            plot_z_scored_trials_data_Maryam(zscored_Stim_data.(Outcome{outcome}),t_trials,Outcome{outcome},Condition{cond},IdChannel,trial_mode{t_mode},avg_mode,color2plotchannels,order_crit,limits2plot,LimitXaxis,Color_scale,show_plot,save_plot,PATH2SAVEALLTRIALS)
                        end
                    end
                end
                close all
            end
        end
    end
    
    %% All trials Heatmaps 
    raw
    trial_mode = {'raw','baseline_corrected','zscored'};
    if ~isempty(Stim_data)
        Outcome = fieldnames(Stim_data);
        if ~isempty(Stim_data.(Outcome{1}))
            Condition = fieldnames(Stim_data.(Outcome{1}));
        else
            Condition = fieldnames(Stim_data.(Outcome{2}));
        end
        avg_mode = 1; % 1, use mean and SEM; 2, use median and MAD
        for cond = 1:length(Condition)
            for t_mode = 1:length(trial_mode)
                %%% All trials in order
                figure
                for channel = 1:length(IdChannel)
                    tmp = Stim_data.All.(Condition{cond}).dFF.(IdChannel{channel}).(trial_mode{t_mode});
                    subplot(length(IdChannel),1,channel)
                    imagesc(t_trials,1,tmp)
                    colormap('jet')
                    c = colorbar;
                    c.Label.String = 'Intensity (A.U.)';
                    if ~isempty(Color_scale)
                        c_limits = Color_scale;
                        caxis(c_limits)
                    end
                    hold on
                    if ~isempty(LimitXaxis)
                        xlim([LimitXaxis(1) LimitXaxis(end)])
                    else
                        xlim([t_trials(1) t_trials(end)])
                    end
                    xlabel('Time (s)')
                    ylabel('Trial')
                    title([IdChannel{channel},' aligned to ',Condition{cond}],'Interpreter', 'none')
                    box off
                    ylimits = get(gca,'YLim');
                    plot([0 0],ylimits,'k','LineWidth',2)
                end
                if save_plot == 1
                    saveas(gcf,[PATH2SAVEALLTRIALS,'All trials in order Heatmaps of dFF signals centered on ',Condition{cond},' ',trial_mode{t_mode},'.jpg'])
                    saveas(gcf,[PATH2SAVEALLTRIALS,'All trials in order Heatmaps of dFF signals centered on ',Condition{cond},' ',trial_mode{t_mode},'.fig'])
                end
                %%% All trials separated by outcome
                count = 0;
                figure
                for channel = 1:length(IdChannel)
                    tmp = [];
                    for outcome = 2:length(Outcome)
                        if ~isempty(Stim_data.(Outcome{outcome}))
                            if ~isempty(Stim_data.(Outcome{outcome}).(Condition{cond}))
                                count = count + 1;
                                [tmp] = [tmp;Stim_data.(Outcome{outcome}).(Condition{cond}).dFF.(IdChannel{channel}).(trial_mode{t_mode})];
                                if count == 1
                                    trial_limit = size(tmp,1);
                                end
                            end
                        end
                    end
                    subplot(length(IdChannel),1,channel)
                    imagesc(t_trials,1,tmp)
                    colormap('jet')
                    c = colorbar;
                    c.Label.String = 'Intensity (A.U.)';
                    if ~isempty(Color_scale)
                        c_limits = Color_scale;
                        caxis(c_limits)
                    end
                    hold on
                    yline(trial_limit+0.5,'w','LineWidth',2)
                    if ~isempty(LimitXaxis)
                        xlim([LimitXaxis(1) LimitXaxis(end)])
                    else
                        xlim([t_trials(1) t_trials(end)])
                    end
                    xlabel('Time (s)')
                    ylabel('Trial')
                    title([IdChannel{channel},' aligned to ',Condition{cond}],'Interpreter', 'none')
                    box off
                    ylimits = get(gca,'YLim');
                    plot([0 0],ylimits,'k','LineWidth',2)
                end
                if save_plot == 1
                    saveas(gcf,[PATH2SAVEALLTRIALS,'All trials Heatmaps of dFF signals centered on ',Condition{cond},' ',trial_mode{t_mode},'.jpg'])
                    saveas(gcf,[PATH2SAVEALLTRIALS,'All trials Heatmaps of dFF signals centered on ',Condition{cond},' ',trial_mode{t_mode},'.fig'])
                end
            end
        end
    end
    % Zscore
    if ~isempty(zscored_Stim_data)
        Outcome = fieldnames(Stim_data);
        if ~isempty(zscored_Stim_data.(Outcome{1}))
            Condition = fieldnames(zscored_Stim_data.(Outcome{1}));
        else
            Condition = fieldnames(zscored_Stim_data.(Outcome{2}));
        end
        avg_mode = 1; % 1, use mean and SEM; 2, use median and MAD
        for cond = 1:length(Condition)
            for t_mode = 1:length(trial_mode)
                %%% All trials in order
                figure
                for channel = 1:length(IdChannel)
                    tmp = zscored_Stim_data.All.(Condition{cond}).dFF.(IdChannel{channel}).(trial_mode{t_mode});
                    subplot(length(IdChannel),1,channel)
                    imagesc(t_trials,1,tmp)
                    colormap('jet')
                    c = colorbar;
                    c.Label.String = 'Intensity (A.U.)';
                    if ~isempty(Color_scale)
                        c_limits = Color_scale;
                        caxis(c_limits)
                    end
                    hold on
                    if ~isempty(LimitXaxis)
                        xlim([LimitXaxis(1) LimitXaxis(end)])
                    else
                        xlim([t_trials(1) t_trials(end)])
                    end
                    xlabel('Time (s)')
                    ylabel('Trial')
                    title([IdChannel{channel},' aligned to ',Condition{cond}],'Interpreter', 'none')
                    box off
                    ylimits = get(gca,'YLim');
                    plot([0 0],ylimits,'k','LineWidth',2)
                end
                if save_plot == 1
                    saveas(gcf,[PATH2SAVEALLTRIALS,'All trials in order Heatmaps of Zscored dFF signals centered on ',Condition{cond},' ',trial_mode{t_mode},'.jpg'])
                    saveas(gcf,[PATH2SAVEALLTRIALS,'All trials in order Heatmaps of Zscored dFF signals centered on ',Condition{cond},' ',trial_mode{t_mode},'.fig'])
                end
                %%% All trials separated by outcome
                count = 0;
                figure
                for channel = 1:length(IdChannel)
                    tmp = [];
                    for outcome = 2:length(Outcome)
                        if ~isempty(zscored_Stim_data.(Outcome{outcome}))
                            if ~isempty(zscored_Stim_data.(Outcome{outcome}).(Condition{cond}))
                                count = count + 1;
                                [tmp] = [tmp;zscored_Stim_data.(Outcome{outcome}).(Condition{cond}).dFF.(IdChannel{channel}).(trial_mode{t_mode})];
                                if count == 1
                                    trial_limit = size(tmp,1);
                                end
                            end
                        end
                    end
                    subplot(length(IdChannel),1,channel)
                    imagesc(t_trials,1,tmp)
                    colormap('jet')
                    c = colorbar;
                    c.Label.String = 'Intensity (A.U.)';
                    if ~isempty(Color_scale)
                        c_limits = Color_scale;
                        caxis(c_limits)
                    end
                    hold on
                    yline(trial_limit+0.5,'w','LineWidth',2)
                    if ~isempty(LimitXaxis)
                        xlim([LimitXaxis(1) LimitXaxis(end)])
                    else
                        xlim([t_trials(1) t_trials(end)])
                    end
                    xlabel('Time (s)')
                    ylabel('Trial')
                    title([IdChannel{channel},' aligned to ',Condition{cond}],'Interpreter', 'none')
                    box off
                    ylimits = get(gca,'YLim');
                    plot([0 0],ylimits,'k','LineWidth',2)
                end
                if save_plot == 1
                    saveas(gcf,[PATH2SAVEALLTRIALS,'All trials Heatmaps of Zscored dFF signals centered on ',Condition{cond},' ',trial_mode{t_mode},'.jpg'])
                    saveas(gcf,[PATH2SAVEALLTRIALS,'All trials Heatmaps of Zscored dFF signals centered on ',Condition{cond},' ',trial_mode{t_mode},'.fig'])
                end
            end
        end
    end

    % %
    % figure
    % for o = 1:length(IdChannel)
    %     if ~isempty(IdChannel{o})
    %         subplot(length(IdChannel),1,o)
    %         imagesc(t_trials,1,Stim_data.(Condition).dFF.(IdChannel{o}).(trial_mode)(order,:))
    %         colormap('jet')
    %         c = colorbar;
    %         c.Label.String = 'Intensity (A.U.)';
    %         if ~isempty(Color_scale)
    %             c_limits = Color_scale;
    %             caxis(c_limits)
    %         end
    %         hold on
    %         xlim([t_trials(1) t_trials(end)])
    %         xlabel('Time (s)')
    %         ylabel('Animal #')
    %         title([IdChannel{o},' ',trial_mode,' aligned to ',Condition],'Interpreter', 'none')
    %         box off
    %         ylimits = get(gca,'YLim');
    %         plot([0 0],ylimits,'k','LineWidth',2)
    %         %     xline(0,'-k');
    %     end
    % end

    %% Save combined trials
    if done == 0 || overwrite == 1
        Params.MICEId = mice_list(m).name;
        Params.min2remove = min2remove; % Minutes to remove from the beginning of the recording
        Params.TRANGE = TRANGE; % window size [start time relative to epoc onset, window duration]
        Params.BASELINE_PER = BASELINE_PER; % baseline period within our window
        Params.Limits_Selected_Epochs = Limits_Selected_Epochs; % Duration criteria for selected Choice epochs (in seconds)
        Params.Epoch2analyze = Epoch2analyze; %Time in seconds to get measurements aligned to behavior
        % Params.concatenation_threshold = concatenation_threshold;
        % Params.duration_threshold = duration_threshold;

        % save([PATH2SAVE,'Trial_analysis.mat'],'time_dwn','dFF','lag2plot','Ovrl_corr',...
        %     'moving_corr','Stim_data','Ovrl_Measurements','Measurements','t_trials',...
        %     'Behav_time','Choice_score','Choice_epochs','shock_evoked_response','Params')
        
        save([path2save_mouse,'Session_analysis.mat'],'Stim_data','zscored_Stim_data','t_trials','alignment_measurements',...
                'alignment_zscored_measurements','Measurements','zscored_Measurements',...
                'Overall_AUC','zscored_Overall_AUC','Overall_meandFF','zscored_Overall_meandFF','Params')
    end
    clear Stim_data alignment_measurements alignment_zscored_measurements Measurements zscored_Stim_data zscored_Measurements Overall_AUC zscored_Overall_AUC Overall_meandFF zscored_Overall_meandFF

