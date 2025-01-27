% Importing TTL data into MatlabRewardOnset
close all; clear variables; clc;
set(0,'defaultfigurecolor',[1 1 1])


%% Define paths and data to analyze
path2data = 'C:\Users\Anacker1\Desktop\Fiber-Data\Fibers\MH\Y MAZE\'; %Location of the data
path2savefolder = 'C:\Users\Anacker1\Desktop\Fiber-Data\Results\'; %Path to save

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

Animals_Right_correct = [1 2 3 4];
if contains(path2data,'Reversal')
    Animals_Right_correct = setdiff(1:Nmice,Animals_Right_correct);
end

FPrig = 1; % 1 for new rig, 2 for old rig (Pedro's data)
Behav_data = 'Ymaze';%'Digging';'Barnes'
IdChannel = {'dCA1GCamp6f','vCA1GCamp6f'}; %{'Psychlight','jRGECO1b','GCamp6f'};{'GCamp6f',[], 'Psychlight'}in some animals/ some trials are in oposite direction
color2plotchannels = {'r','g'};%{'r',[],'g'}
MiceperFile = 1;
min2remove = 0.1; % Minutes to remove from the beginning of the recording
TRANGE = [-5 30]; % window size [start time relative to epoc onset, window duration]
BASELINE_PER.Onset = [-5 0]; % baseline period within our window
BASELINE_PER.Offset = [-5 0]; % baseline period within our window
BASELINE_PER.ShockOnset = [-5 0]; % baseline period within our window
% BASELINE_PER.Reward = [-3 -1]; % baseline period within our window

Limits_Selected_Epochs.short = [2 3]; % Duration criteria for selected Choice epochs (in seconds)
Limits_Selected_Epochs.long = [4 180];
Limits_IEI = 3;
Epoch2analyze = 2; %Time in seconds to get measurements aligned to behavior:Nmice

for i = 1:length(IdChannel)
    if ~isempty(IdChannel{i})
        limits2plot.(IdChannel{i}) = []; %If empty, automatically adjusted to the data.\
    end
end
limits2plot.all = [];
Color_scale = []; % For the heatmaps. If empty, automatically adjusted to the data.

show_plot = 0; % If 0, plots are not displayed##### if you want the plots to pop up : show_plot= 1
save_plot = 1; % If 0, plots are not saved
reanalysis = 0; % If 1, the code runs for sessions already analyze. If 0, session excluded if analysed
overwrite = 1; % If reanalyzing data, 1 if you want to save them and overwrite results and 0 if not

for m =1:Nmice
    if contains(path2data,'\')
    % % Define the path to the mouse and find the folders to analyze:
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

            %% Define behavioral data
            filt_mode = 'unfilt'; % unfilt or filt. Filt referes to the highpass.
            zero_mode = 'uncorrected';
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
                    Behavior_code = {'Start','Decision','Left','Right'};
                    Color_code = {'c','r','g','b'};
                elseif strcmp(Behav_data,'Barnes')
                    Behavior_code = {'Start','Decision','Left','Right'};
                    Color_code = {'c','r','g','b'};
                end

                %%% For Matlab 2020 or later %%%
                [Behav_time,Choice_score,Choice_epochs] = ExtractBehavData_Maryam...
                    (behav_files(behav_file),mice_list(m).name,sessions(s).name,time_dwn(1),Behavior_code);
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
                
                zoom_limits = [time_dwn(1) 30]; % Window to plot zoom of the session (In seconds)
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
            for ooo = 1
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

                Right_choice = 0;
                if ~isempty(Choice_epochs.Right)
                    Right_choice = 1;
                end
                if ~isempty(intersect(str2double(mice_list(m).name(end)),Animals_Right_correct)) && Right_choice
                    Trial_correct = 1;
                else
                    Trial_correct = 0;
                end
                for o = 1:length(Behavior_code)
                    if ~isempty(Choice_epochs.(Behavior_code{o}))
                        if Trial_correct
                            tmp_name = [Behavior_code{o},'Onset'];
                            Trial_data.Correct.(tmp_name).idx = ones(size(Choice_epochs.(Behavior_code{o}),1),n)*nan;
                            Trial_data.Correct.(tmp_name).corr = ones(size(Choice_epochs.(Behavior_code{o}),1),n)*nan;
                            for i = 1:length(IdChannel)
                                if ~isempty(IdChannel{i})
                                    Trial_data.Correct.(tmp_name).dFF.(IdChannel{i}).raw = ones(size(Choice_epochs.(Behavior_code{o}),1),n)*nan;
                                    Trial_data.Correct.(tmp_name).dFF.(IdChannel{i}).baseline_corrected = ones(size(Choice_epochs.(Behavior_code{o}),1),n)*nan;
                                end
                            end
                            for u = 1:size(Choice_epochs.(Behavior_code{o}),1)
                                ix = find(abs(time_dwn-Choice_epochs.(Behavior_code{o})(u,1)) == ...
                                    min(abs(time_dwn-Choice_epochs.(Behavior_code{o})(u,1))));
                                if ix+idx_End <= length(time_dwn) && ix-(idx_Init-1) >= 1
                                    tmp = ix - (idx_Init-1):ix + idx_End;
                                    Trial_data.Correct.(tmp_name).idx(u,:) = tmp;
                                    for i = 1:length(IdChannel)
                                        if ~isempty(IdChannel{i})
                                            Trial_data.Correct.(tmp_name).dFF.(IdChannel{i}).raw(u,:) = ...
                                                dFF.(IdChannel{i}).(filt_mode).(zero_mode)(tmp);%dFF.(IdChannel{i}).filt(tmp);
                                        end
                                    end
                                    if ~isempty(moving_corr)
                                        Trial_data.Correct.(tmp_name).corr(u,:) = moving_corr(tmp);
                                    end
                                elseif ix+idx_End > length(time_dwn)
                                    tmp = ix - (idx_Init-1):length(time_dwn);
                                    Trial_data.Correct.(tmp_name).idx(u,1:length(tmp)) = tmp;
                                    for i = 1:length(IdChannel)
                                        if ~isempty(IdChannel{i})
                                            Trial_data.Correct.(tmp_name).dFF.(IdChannel{i}).raw(u,1:length(tmp)) = ...
                                                dFF.(IdChannel{i}).(filt_mode).(zero_mode)(tmp);
                                        end
                                    end
                                    if ~isempty(moving_corr)
                                        Trial_data.Correct.(tmp_name).corr(u,1:length(tmp)) = moving_corr(tmp);
                                    end
                                elseif ix-(idx_Init-1) < 1
                                    tmp = 1:(ix + idx_End);
                                    Trial_data.Correct.(tmp_name).idx(u,(n-length(tmp)+1):end) = tmp;
                                    for i = 1:length(IdChannel)
                                        if ~isempty(IdChannel{i})
                                            Trial_data.Correct.(tmp_name).dFF.(IdChannel{i}).raw(u,(n-length(tmp)+1):end) = ...
                                                dFF.(IdChannel{i}).(filt_mode).(zero_mode)(tmp);
                                        end
                                    end
                                    if ~isempty(moving_corr)
                                        Trial_data.Correct.(tmp_name).corr(u,(n-length(tmp)+1):end) = moving_corr(tmp);
                                    end
                                end
                            end
                            for i = 1:length(IdChannel)
                                if ~isempty(IdChannel{i})
                                    %                 Trial_data.TrialOnset.dFF.(IdChannel{i}).baseline_corrected = ...
                                    %                     Trial_data.TrialOnset.dFF.(IdChannel{i}).raw - nanmean...
                                    %                     (Trial_data.TrialOnset.dFF.(IdChannel{i}).raw...
                                    %                     (:,t_trials >= BASELINE_PER(1) & t_trials <= BASELINE_PER(2)),2);
                                    Trial_data.Correct.(tmp_name).dFF.(IdChannel{i}).baseline_corrected = ...
                                        Trial_data.Correct.(tmp_name).dFF.(IdChannel{i}).raw - nanmedian...
                                        (Trial_data.Correct.(tmp_name).dFF.(IdChannel{i}).raw...
                                        (:,t_trials >= BASELINE_PER.Onset(1) & ...
                                        t_trials <= BASELINE_PER.Onset(2)),2);
                                    Trial_data.Correct.(tmp_name).dFF.(IdChannel{i}).zscored = ...
                                        (Trial_data.Correct.(tmp_name).dFF.(IdChannel{i}).raw - mean...
                                        (Trial_data.Correct.(tmp_name).dFF.(IdChannel{i}).raw...
                                        (:,t_trials >= BASELINE_PER.Onset(1) & t_trials <= BASELINE_PER.Onset(2)),2,'omitnan'))./std(...
                                        Trial_data.Correct.(tmp_name).dFF.(IdChannel{i}).raw...
                                        (:,t_trials >= BASELINE_PER.Onset(1) & t_trials <= BASELINE_PER.Onset(2)),1,2,'omitnan');
                                end
                            end

                            tmp_name = [Behavior_code{o},'Offset'];
                            Trial_data.Correct.(tmp_name).idx = ones(size(Choice_epochs.(Behavior_code{o}),1),n)*nan;
                            Trial_data.Correct.(tmp_name).corr = ones(size(Choice_epochs.(Behavior_code{o}),1),n)*nan;
                            for i = 1:length(IdChannel)
                                if ~isempty(IdChannel{i})
                                    Trial_data.Correct.(tmp_name).dFF.(IdChannel{i}).raw = ones(size(Choice_epochs.(Behavior_code{o}),1),n)*nan;
                                    Trial_data.Correct.(tmp_name).dFF.(IdChannel{i}).baseline_corrected = ones(size(Choice_epochs.(Behavior_code{o}),1),n)*nan;
                                end
                            end
                            for u = 1:size(Choice_epochs.(Behavior_code{o}),1)
                                ix = find(abs(time_dwn-Choice_epochs.(Behavior_code{o})(u,2)) == ...
                                    min(abs(time_dwn-Choice_epochs.(Behavior_code{o})(u,2))));
                                if ix+idx_End <= length(time_dwn) && ix-(idx_Init-1) >= 1
                                    tmp = ix - (idx_Init-1):ix + idx_End;
                                    Trial_data.Correct.(tmp_name).idx(u,:) = tmp;
                                    for i = 1:length(IdChannel)
                                        if ~isempty(IdChannel{i})
                                            Trial_data.Correct.(tmp_name).dFF.(IdChannel{i}).raw(u,:) = ...
                                                dFF.(IdChannel{i}).(filt_mode).(zero_mode)(tmp);%dFF.(IdChannel{i}).filt(tmp);
                                        end
                                    end
                                    if ~isempty(moving_corr)
                                        Trial_data.Correct.(tmp_name).corr(u,:) = moving_corr(tmp);
                                    end
                                elseif ix+idx_End > length(time_dwn)
                                    tmp = ix - (idx_Init-1):length(time_dwn);
                                    Trial_data.Correct.(tmp_name).idx(u,1:length(tmp)) = tmp;
                                    for i = 1:length(IdChannel)
                                        if ~isempty(IdChannel{i})
                                            Trial_data.Correct.(tmp_name).dFF.(IdChannel{i}).raw(u,1:length(tmp)) = ...
                                                dFF.(IdChannel{i}).(filt_mode).(zero_mode)(tmp);
                                        end
                                    end
                                    if ~isempty(moving_corr)
                                        Trial_data.Correct.(tmp_name).corr(u,1:length(tmp)) = moving_corr(tmp);
                                    end
                                elseif ix-(idx_Init-1) < 1
                                    tmp = 1:(ix + idx_End);
                                    Trial_data.Correct.(tmp_name).idx(u,(n-length(tmp)+1):end) = tmp;
                                    for i = 1:length(IdChannel)
                                        if ~isempty(IdChannel{i})
                                            Trial_data.Correct.(tmp_name).dFF.(IdChannel{i}).raw(u,(n-length(tmp)+1):end) = ...
                                                dFF.(IdChannel{i}).(filt_mode).(zero_mode)(tmp);
                                        end
                                    end
                                    if ~isempty(moving_corr)
                                        Trial_data.Correct.(tmp_name).corr(u,(n-length(tmp)+1):end) = moving_corr(tmp);
                                    end
                                end
                            end
                            for i = 1:length(IdChannel)
                                if ~isempty(IdChannel{i})
                                    %                 Trial_data.TrialOnset.dFF.(IdChannel{i}).baseline_corrected = ...
                                    %                     Trial_data.TrialOnset.dFF.(IdChannel{i}).raw - nanmean...
                                    %                     (Trial_data.TrialOnset.dFF.(IdChannel{i}).raw...
                                    %                     (:,t_trials >= BASELINE_PER(1) & t_trials <= BASELINE_PER(2)),2);
                                    Trial_data.Correct.(tmp_name).dFF.(IdChannel{i}).baseline_corrected = ...
                                        Trial_data.Correct.(tmp_name).dFF.(IdChannel{i}).raw - nanmedian...
                                        (Trial_data.Correct.(tmp_name).dFF.(IdChannel{i}).raw...
                                        (:,t_trials >= BASELINE_PER.Offset(1) & ...
                                        t_trials <= BASELINE_PER.Offset(2)),2);
                                    Trial_data.Correct.(tmp_name).dFF.(IdChannel{i}).zscored = ...
                                        (Trial_data.Correct.(tmp_name).dFF.(IdChannel{i}).raw - mean...
                                        (Trial_data.Correct.(tmp_name).dFF.(IdChannel{i}).raw...
                                        (:,t_trials >= BASELINE_PER.Onset(1) & t_trials <= BASELINE_PER.Offset(2)),2,'omitnan'))./std(...
                                        Trial_data.Correct.(tmp_name).dFF.(IdChannel{i}).raw...
                                        (:,t_trials >= BASELINE_PER.Onset(1) & t_trials <= BASELINE_PER.Offset(2)),1,2,'omitnan');
                                end
                            end
                        else
                            tmp_name = [Behavior_code{o},'Onset'];
                            Trial_data.Incorrect.(tmp_name).idx = ones(size(Choice_epochs.(Behavior_code{o}),1),n)*nan;
                            Trial_data.Incorrect.(tmp_name).corr = ones(size(Choice_epochs.(Behavior_code{o}),1),n)*nan;
                            for i = 1:length(IdChannel)
                                if ~isempty(IdChannel{i})
                                    Trial_data.Incorrect.(tmp_name).dFF.(IdChannel{i}).raw = ones(size(Choice_epochs.(Behavior_code{o}),1),n)*nan;
                                    Trial_data.Incorrect.(tmp_name).dFF.(IdChannel{i}).baseline_corrected = ones(size(Choice_epochs.(Behavior_code{o}),1),n)*nan;
                                end
                            end
                            for u = 1:size(Choice_epochs.(Behavior_code{o}),1)
                                ix = find(abs(time_dwn-Choice_epochs.(Behavior_code{o})(u,1)) == ...
                                    min(abs(time_dwn-Choice_epochs.(Behavior_code{o})(u,1))));
                                if ix+idx_End <= length(time_dwn) && ix-(idx_Init-1) >= 1
                                    tmp = ix - (idx_Init-1):ix + idx_End;
                                    Trial_data.Incorrect.(tmp_name).idx(u,:) = tmp;
                                    for i = 1:length(IdChannel)
                                        if ~isempty(IdChannel{i})
                                            Trial_data.Incorrect.(tmp_name).dFF.(IdChannel{i}).raw(u,:) = ...
                                                dFF.(IdChannel{i}).(filt_mode).(zero_mode)(tmp);%dFF.(IdChannel{i}).filt(tmp);
                                        end
                                    end
                                    if ~isempty(moving_corr)
                                        Trial_data.Incorrect.(tmp_name).corr(u,:) = moving_corr(tmp);
                                    end
                                elseif ix+idx_End > length(time_dwn)
                                    tmp = ix - (idx_Init-1):length(time_dwn);
                                    Trial_data.Incorrect.(tmp_name).idx(u,1:length(tmp)) = tmp;
                                    for i = 1:length(IdChannel)
                                        if ~isempty(IdChannel{i})
                                            Trial_data.Incorrect.(tmp_name).dFF.(IdChannel{i}).raw(u,1:length(tmp)) = ...
                                                dFF.(IdChannel{i}).(filt_mode).(zero_mode)(tmp);
                                        end
                                    end
                                    if ~isempty(moving_corr)
                                        Trial_data.Incorrect.(tmp_name).corr(u,1:length(tmp)) = moving_corr(tmp);
                                    end
                                elseif ix-(idx_Init-1) < 1
                                    tmp = 1:(ix + idx_End);
                                    Trial_data.Incorrect.(tmp_name).idx(u,(n-length(tmp)+1):end) = tmp;
                                    for i = 1:length(IdChannel)
                                        if ~isempty(IdChannel{i})
                                            Trial_data.Incorrect.(tmp_name).dFF.(IdChannel{i}).raw(u,(n-length(tmp)+1):end) = ...
                                                dFF.(IdChannel{i}).(filt_mode).(zero_mode)(tmp);
                                        end
                                    end
                                    if ~isempty(moving_corr)
                                        Trial_data.Incorrect.(tmp_name).corr(u,(n-length(tmp)+1):end) = moving_corr(tmp);
                                    end
                                end
                            end
                            for i = 1:length(IdChannel)
                                if ~isempty(IdChannel{i})
                                    %                 Trial_data.TrialOnset.dFF.(IdChannel{i}).baseline_corrected = ...
                                    %                     Trial_data.TrialOnset.dFF.(IdChannel{i}).raw - nanmean...
                                    %                     (Trial_data.TrialOnset.dFF.(IdChannel{i}).raw...
                                    %                     (:,t_trials >= BASELINE_PER(1) & t_trials <= BASELINE_PER(2)),2);
                                    Trial_data.Incorrect.(tmp_name).dFF.(IdChannel{i}).baseline_corrected = ...
                                        Trial_data.Incorrect.(tmp_name).dFF.(IdChannel{i}).raw - nanmedian...
                                        (Trial_data.Incorrect.(tmp_name).dFF.(IdChannel{i}).raw...
                                        (:,t_trials >= BASELINE_PER.Onset(1) & ...
                                        t_trials <= BASELINE_PER.Onset(2)),2);
                                    Trial_data.Incorrect.(tmp_name).dFF.(IdChannel{i}).zscored = ...
                                        (Trial_data.Incorrect.(tmp_name).dFF.(IdChannel{i}).raw - mean...
                                        (Trial_data.Incorrect.(tmp_name).dFF.(IdChannel{i}).raw...
                                        (:,t_trials >= BASELINE_PER.Onset(1) & t_trials <= BASELINE_PER.Onset(2)),2,'omitnan'))./std(...
                                        Trial_data.Incorrect.(tmp_name).dFF.(IdChannel{i}).raw...
                                        (:,t_trials >= BASELINE_PER.Onset(1) & t_trials <= BASELINE_PER.Onset(2)),1,2,'omitnan');
                                end
                            end

                            tmp_name = [Behavior_code{o},'Offset'];
                            Trial_data.Incorrect.(tmp_name).idx = ones(size(Choice_epochs.(Behavior_code{o}),1),n)*nan;
                            Trial_data.Incorrect.(tmp_name).corr = ones(size(Choice_epochs.(Behavior_code{o}),1),n)*nan;
                            for i = 1:length(IdChannel)
                                if ~isempty(IdChannel{i})
                                    Trial_data.Incorrect.(tmp_name).dFF.(IdChannel{i}).raw = ones(size(Choice_epochs.(Behavior_code{o}),1),n)*nan;
                                    Trial_data.Incorrect.(tmp_name).dFF.(IdChannel{i}).baseline_corrected = ones(size(Choice_epochs.(Behavior_code{o}),1),n)*nan;
                                end
                            end
                            for u = 1:size(Choice_epochs.(Behavior_code{o}),1)
                                ix = find(abs(time_dwn-Choice_epochs.(Behavior_code{o})(u,2)) == ...
                                    min(abs(time_dwn-Choice_epochs.(Behavior_code{o})(u,2))));
                                if ix+idx_End <= length(time_dwn) && ix-(idx_Init-1) >= 1
                                    tmp = ix - (idx_Init-1):ix + idx_End;
                                    Trial_data.Incorrect.(tmp_name).idx(u,:) = tmp;
                                    for i = 1:length(IdChannel)
                                        if ~isempty(IdChannel{i})
                                            Trial_data.Incorrect.(tmp_name).dFF.(IdChannel{i}).raw(u,:) = ...
                                                dFF.(IdChannel{i}).(filt_mode).(zero_mode)(tmp);%dFF.(IdChannel{i}).filt(tmp);
                                        end
                                    end
                                    if ~isempty(moving_corr)
                                        Trial_data.Incorrect.(tmp_name).corr(u,:) = moving_corr(tmp);
                                    end
                                elseif ix+idx_End > length(time_dwn)
                                    tmp = ix - (idx_Init-1):length(time_dwn);
                                    Trial_data.Incorrect.(tmp_name).idx(u,1:length(tmp)) = tmp;
                                    for i = 1:length(IdChannel)
                                        if ~isempty(IdChannel{i})
                                            Trial_data.Incorrect.(tmp_name).dFF.(IdChannel{i}).raw(u,1:length(tmp)) = ...
                                                dFF.(IdChannel{i}).(filt_mode).(zero_mode)(tmp);
                                        end
                                    end
                                    if ~isempty(moving_corr)
                                        Trial_data.Incorrect.(tmp_name).corr(u,1:length(tmp)) = moving_corr(tmp);
                                    end
                                elseif ix-(idx_Init-1) < 1
                                    tmp = 1:(ix + idx_End);
                                    Trial_data.Incorrect.(tmp_name).idx(u,(n-length(tmp)+1):end) = tmp;
                                    for i = 1:length(IdChannel)
                                        if ~isempty(IdChannel{i})
                                            Trial_data.Incorrect.(tmp_name).dFF.(IdChannel{i}).raw(u,(n-length(tmp)+1):end) = ...
                                                dFF.(IdChannel{i}).(filt_mode).(zero_mode)(tmp);
                                        end
                                    end
                                    if ~isempty(moving_corr)
                                        Trial_data.Incorrect.(tmp_name).corr(u,(n-length(tmp)+1):end) = moving_corr(tmp);
                                    end
                                end
                            end
                            for i = 1:length(IdChannel)
                                if ~isempty(IdChannel{i})
                                    %                 Trial_data.TrialOnset.dFF.(IdChannel{i}).baseline_corrected = ...
                                    %                     Trial_data.TrialOnset.dFF.(IdChannel{i}).raw - nanmean...
                                    %                     (Trial_data.TrialOnset.dFF.(IdChannel{i}).raw...
                                    %                     (:,t_trials >= BASELINE_PER(1) & t_trials <= BASELINE_PER(2)),2);
                                    Trial_data.Incorrect.(tmp_name).dFF.(IdChannel{i}).baseline_corrected = ...
                                        Trial_data.Incorrect.(tmp_name).dFF.(IdChannel{i}).raw - nanmedian...
                                        (Trial_data.Incorrect.(tmp_name).dFF.(IdChannel{i}).raw...
                                        (:,t_trials >= BASELINE_PER.Offset(1) & ...
                                        t_trials <= BASELINE_PER.Offset(2)),2);
                                    Trial_data.Incorrect.(tmp_name).dFF.(IdChannel{i}).zscored = ...
                                        (Trial_data.Incorrect.(tmp_name).dFF.(IdChannel{i}).raw - mean...
                                        (Trial_data.Incorrect.(tmp_name).dFF.(IdChannel{i}).raw...
                                        (:,t_trials >= BASELINE_PER.Onset(1) & t_trials <= BASELINE_PER.Offset(2)),2,'omitnan'))./std(...
                                        Trial_data.Incorrect.(tmp_name).dFF.(IdChannel{i}).raw...
                                        (:,t_trials >= BASELINE_PER.Onset(1) & t_trials <= BASELINE_PER.Offset(2)),1,2,'omitnan');
                                end
                            end
                        end
                    else
                        if Trial_correct
                            tmp_name = [Behavior_code{o},'Onset'];
                            Trial_data.Correct.(tmp_name) = [];
                            tmp_name = [Behavior_code{o},'Offset'];
                            Trial_data.Correct.(tmp_name) = [];
                        else
                            tmp_name = [Behavior_code{o},'Onset'];
                            Trial_data.Incorrect.(tmp_name) = [];
                            tmp_name = [Behavior_code{o},'Offset'];
                            Trial_data.Incorrect.(tmp_name) = [];
                        end
                    end
                end
                %% Add trial to Stim_data:
                if s == 1
                    Stim_data = Trial_data;
                else
                    trial_mode = {'raw','baseline_corrected','zscored'};
                    fld1 = fieldnames(Stim_data);
                    for i = 1:length(fld1)
                        fld2 = fieldnames(Stim_data.(fld1{i}));
                        for outcome = 1:length(fld2)
                            if ~isempty(Stim_data.(fld1{i}).(fld2{outcome})) && ~isempty(Trial_data.(fld1{i}).(fld2{outcome}))
                                fld3 = fieldnames(Stim_data.(fld1{i}).(fld2{outcome}));
                                for ii = 1:length(fld3)

                                    if isa(Stim_data.(fld1{i}).(fld2{outcome}).(fld3{ii}),'double')
                                        [Stim_data.(fld1{i}).(fld2{outcome}).(fld3{ii})] = [Stim_data.(fld1{i}).(fld2{outcome}).(fld3{ii});...
                                            Trial_data.(fld1{i}).(fld2{outcome}).(fld3{ii})];
                                    else
                                        for iii = 1:length(IdChannel)
                                            for iv = 1:length(trial_mode)
                                                [Stim_data.(fld1{i}).(fld2{outcome}).(fld3{ii}).(IdChannel{iii}).(trial_mode{iv})] = ...
                                                    [Stim_data.(fld1{i}).(fld2{outcome}).(fld3{ii}).(IdChannel{iii}).(trial_mode{iv});...
                                                    Trial_data.(fld1{i}).(fld2{outcome}).(fld3{ii}).(IdChannel{iii}).(trial_mode{iv})];
                                            end
                                        end
                                    end
                                end
                            elseif ~isempty(Trial_data.(fld1{i}).(fld2{outcome}))
                                Stim_data.(fld1{i}).(fld2{outcome}) = Trial_data.(fld1{i}).(fld2{outcome});
                            end
                        end
                    end
                end

                %% Compute measurements for each of the behavioral states:
                filt_mode = 'unfilt'; % unfilt or filt. Filt referes to the highpass.
                zero_mode = 'uncorrected'; % uncorrected or corrected. Corrected referes to the 8th percentile correction.
                
                threshold_criterion = 'mean';% 'median'
                n_std = 1;

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

                        for behav = 1:length(Behavior_code)
                            if ~isempty(Choice_epochs.(Behavior_code{behav}))
                                Trial_Measurements.norm_AUC.(IdChannel{o}).(Behavior_code{behav}) = ...
                                    trapz(tmp(time_dwn>= Choice_epochs.(Behavior_code{behav})(1) & ...
                                    time_dwn < Choice_epochs.(Behavior_code{behav})(2)))./(diff(Choice_epochs.(Behavior_code{behav})));
                                Trial_Measurements.mean_dFF.(IdChannel{o}).(Behavior_code{behav}) = ...
                                    mean(tmp(time_dwn>= Choice_epochs.(Behavior_code{behav})(1) & ...
                                    time_dwn < Choice_epochs.(Behavior_code{behav})(2)));

                                Trial_Measurements.peaks.(IdChannel{o}).(Behavior_code{behav}) = ...
                                    tmp_peaks(time_dwn(tmp_locs)>= Choice_epochs.(Behavior_code{behav})(1) & ...
                                    time_dwn(tmp_locs)< Choice_epochs.(Behavior_code{behav})(2));
                                Trial_Measurements.peak_width.(IdChannel{o}).(Behavior_code{behav}) = ...
                                    tmp_w(time_dwn(tmp_locs)>= Choice_epochs.(Behavior_code{behav})(1) & ...
                                    time_dwn(tmp_locs)< Choice_epochs.(Behavior_code{behav})(2));
                                Trial_Measurements.peak_prom.(IdChannel{o}).(Behavior_code{behav}) = ...
                                    tmp_prom(time_dwn(tmp_locs)>= Choice_epochs.(Behavior_code{behav})(1) & ...
                                    time_dwn(tmp_locs)< Choice_epochs.(Behavior_code{behav})(2));

                                Trial_Measurements.peaks_rate.(IdChannel{o}).(Behavior_code{behav}) = ...
                                    length(tmp_peaks(time_dwn(tmp_locs)>= Choice_epochs.(Behavior_code{behav})(1) & ...
                                    time_dwn(tmp_locs)< Choice_epochs.(Behavior_code{behav})(2)))/...
                                    (diff(Choice_epochs.(Behavior_code{behav})));
                                Trial_Measurements.mean_peaks_val.(IdChannel{o}).(Behavior_code{behav}) = ...
                                    mean(tmp_peaks(time_dwn(tmp_locs)>= Choice_epochs.(Behavior_code{behav})(1) & ...
                                    time_dwn(tmp_locs)< Choice_epochs.(Behavior_code{behav})(2)));
                                Trial_Measurements.mean_peaks_width.(IdChannel{o}).(Behavior_code{behav}) = ...
                                    mean(tmp_w(time_dwn(tmp_locs)>= Choice_epochs.(Behavior_code{behav})(1) & ...
                                    time_dwn(tmp_locs)< Choice_epochs.(Behavior_code{behav})(2)));
                                Trial_Measurements.mean_peaks_prom.(IdChannel{o}).(Behavior_code{behav}) = ...
                                    length(tmp_prom(time_dwn(tmp_locs)>= Choice_epochs.(Behavior_code{behav})(1) & ...
                                    time_dwn(tmp_locs)< Choice_epochs.(Behavior_code{behav})(2)));
                            else
                                Trial_Measurements.norm_AUC.(IdChannel{o}).(Behavior_code{behav}) = [];
                                Trial_Measurements.mean_dFF.(IdChannel{o}).(Behavior_code{behav}) = [];

                                Trial_Measurements.peaks.(IdChannel{o}).(Behavior_code{behav}) = [];
                                Trial_Measurements.peak_width.(IdChannel{o}).(Behavior_code{behav}) = [];
                                Trial_Measurements.peak_prom.(IdChannel{o}).(Behavior_code{behav}) = [];

                                Trial_Measurements.peaks_rate.(IdChannel{o}).(Behavior_code{behav}) = [];
                                Trial_Measurements.mean_peaks_val.(IdChannel{o}).(Behavior_code{behav}) = [];
                                Trial_Measurements.mean_peaks_width.(IdChannel{o}).(Behavior_code{behav}) = [];
                                Trial_Measurements.mean_peaks_prom.(IdChannel{o}).(Behavior_code{behav}) = [];
                            end
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
                save([PATH2SAVE,'Trial_analysis.mat'],'time_dwn','dFF','lag2plot','Ovrl_cincorrectorr',...
                    'moving_corr','Trial_data','t_trials','Trial_Measurements','Behav_time','Choice_score','Choice_epochs','Params')
            end
        end
    end  
    % Plot the results
    PATH2SAVEALLTRIALS = path2save_mouse;
    if ~isempty(Stim_data)
        trial_mode = {'raw','baseline_corrected','zscored'}; % raw or baseline_corrected.
        Condition = fieldnames(Stim_data);%{'LeverExtension','DipperOn','Reward'};
        avg_mode = 1; % 1, use mean and SEM; 2, use median and MAD
        order_crit = 1;
        for cond = 1:length(Condition)
            if ~isempty(Stim_data.(Condition{cond}))
                for t_mode = 1:length(trial_mode)
                    plot_trials_data_Lauren(Stim_data,t_trials,Condition{cond},IdChannel,trial_mode{t_mode},avg_mode,color2plotchannels,order_crit,limits2plot,Color_scale,show_plot,save_plot,PATH2SAVEALLTRIALS)
                end
            end
        end
        close all
    end

    % %% Overall Choice vs non Choice activity
    % Choice_idx = [];
    % if ~isempty(~isnan(Choice_epochs))
    %     Choice_idx = [];
    %     for o = 1:size(Choice_epochs,1)
    %         ix1 = find(abs(time_dwn-Choice_epochs(o,1)) == min(abs(time_dwn-Choice_epochs(o,1))));
    %         ix2 = find(abs(time_dwn-Choice_epochs(o,2)) == min(abs(time_dwn-Choice_epochs(o,2))));
    %
    %         [Choice_idx] = [Choice_idx ix1:ix2];
    %     end
    %     Choice_idx = Choice_idx';
    % end
    % if ~isempty(Choice_idx)
    %     dur_Choice = time_dwn(length(Choice_idx));
    % else
    %     dur_Choice = 0;
    % end
    %
    % ix1 = find(abs(time_dwn-Behav_time(end)) == min(abs(time_dwn-Behav_time(end))));
    % Behav_time2 = time_dwn(1:ix1);
    % non_Choice_idx = setdiff(1:ix1,Choice_idx);
    % if ~isempty(non_Choice_idx)
    %     dur_non_Choice = time_dwn(length(non_Choice_idx));
    % else
    %     dur_non_Choice = 0;
    % end
    %
    % % %% Compute measurements of behaviorally aligned epochs
    % if ~isempty(Trial_data)
    %     idx_epoch = find(t_trials >= -1 & t_trials <= Epoch2analyze);
    %     alignment = fieldnames(Stim_data);
    %     trial_mode = {'raw','baseline_corrected','zscored'};
    %     for fld = 1:length(alignment)
    %         for o = 1:length(IdChannel)
    %             if ~isempty(IdChannel{o})
    %                 for oo = 1:length(trial_mode)
    %                     tmp_time = t_trials(idx_epoch);
    %                     tmp = Stim_data.(alignment{fld}).dFF.(IdChannel{o}).(trial_mode{oo})(:,idx_epoch);
    %                     Measurements.AUC.(alignment{fld}).(IdChannel{o}).(trial_mode{oo}) = nan(size(tmp,1),1);
    %                     for u = 1:size(tmp,1)
    %                         Measurements.AUC.(alignment{fld}).(IdChannel{o}).(trial_mode{oo})(u) = trapz(tmp(u,~isnan(tmp(u,:))),2);
    %                     end
    %                     Measurements.meanDFF.(alignment{fld}).(IdChannel{o}).(trial_mode{oo}) = mean(tmp,1,'omitnan');
    %                     [~,loc] = max(abs(tmp),[],2);
    %                     % Linear indexing
    %                     I = (1:size(tmp,1)).';
    %                     J = reshape(loc,[],1);
    %                     k = sub2ind(size(tmp),I,J);
    %                     Measurements.peak.(alignment{fld}).(IdChannel{o}).(trial_mode{oo}) = tmp(k);
    %                 end
    %             end
    %         end
    %     end
    % else
    %     Measurements = [];
    % end
    %
    % %% Compute overall measurements during Choice vs non-Choice epochs
    % filt_mode = 'raw'; % unfilt or filt. Filt referes to the highpass.
    % zero_mode = 'uncorrected';
    %
    % tmp = time_dwn-time_dwn(1);
    % idx_60 = find(tmp >= 60);
    % idx_60 = idx_60(1);
    % idx_per_minute = 1:idx_60:length(tmp);
    % if tmp(length(tmp)-idx_per_minute(end)) >= 15
    %     [idx_per_minute] = [idx_per_minute length(tmp)];
    % else
    %     idx_per_minute(end) = length(tmp);
    % end
    %
    % for o = 1:length(IdChannel)
    %     if ~isempty(IdChannel{o})
    %         tmp = dFF.(IdChannel{o}).raw.uncorrected;
    %         %                         [tmp_peaks,tmp_locs,~,tmp_prom] = findpeaks(tmp,'MinPeakProminence',0.2,'Annotate','extents');
    %         %                         [~,~,tmp_w,~] = findpeaks(tmp,time_dwn,'MinPeakProminence',0.2,'Annotate','extents');
    %
    %         Ovrl_Measurements.norm_AUC.(IdChannel{o}).full.all = trapz(tmp)/(time_dwn(end)-time_dwn(1));
    %         Ovrl_Measurements.mean_dFF.(IdChannel{o}).full.all = mean(tmp);
    %
    %         [tmp_peaks,tmp_locs,tmp_w,tmp_prom] = findpeaks(tmp,'MinPeakProminence',0.2,'Annotate','extents');
    %         Ovrl_Measurements.peaks_rate.(IdChannel{o}).full.all = length(tmp_peaks)/(time_dwn(end)-time_dwn(1));
    %         Ovrl_Measurements.peaks.(IdChannel{o}).full.all = tmp_peaks;
    %         Ovrl_Measurements.peak_width.(IdChannel{o}).full.all = tmp_w;
    %         Ovrl_Measurements.peak_prom.(IdChannel{o}).full.all = tmp_prom;
    %
    %         Ovrl_Measurements.mean_peaks_val.(IdChannel{o}).full.all = mean(Ovrl_Measurements.peaks.(IdChannel{o}).full.all);
    %         Ovrl_Measurements.mean_peaks_width.(IdChannel{o}).full.all = mean(Ovrl_Measurements.peak_width.(IdChannel{o}).full.all);
    %         Ovrl_Measurements.mean_peaks_prom.(IdChannel{o}).full.all = mean(Ovrl_Measurements.peak_prom.(IdChannel{o}).full.all);
    %
    %         Ovrl_Measurements.norm_AUC.(IdChannel{o}).full.perminute = ones(length(idx_per_minute)-1,1)*nan;
    %         Ovrl_Measurements.mean_dFF.(IdChannel{o}).full.perminute = ones(length(idx_per_minute)-1,1)*nan;
    %         Ovrl_Measurements.peaks_rate.(IdChannel{o}).full.perminute = ones(length(idx_per_minute)-1,1)*nan;
    %         Ovrl_Measurements.peaks.(IdChannel{o}).full.perminute = cell(length(idx_per_minute)-1,1);
    %         Ovrl_Measurements.peak_width.(IdChannel{o}).full.perminute = cell(length(idx_per_minute)-1,1);
    %         Ovrl_Measurements.peak_prom.(IdChannel{o}).full.perminute = cell(length(idx_per_minute)-1,1);
    %         Ovrl_Measurements.mean_peaks_val.(IdChannel{o}).full.perminute = ones(length(idx_per_minute)-1,1)*nan;
    %         Ovrl_Measurements.mean_peaks_width.(IdChannel{o}).full.perminute = ones(length(idx_per_minute)-1,1)*nan;
    %         Ovrl_Measurements.mean_peaks_prom.(IdChannel{o}).full.perminute = ones(length(idx_per_minute)-1,1)*nan;
    %         for i = 1:length(idx_per_minute)-1
    %             min_idx = idx_per_minute(i):(idx_per_minute(i+1)-1);
    %             Ovrl_Measurements.norm_AUC.(IdChannel{o}).full.perminute(i) = trapz(tmp(min_idx))/(time_dwn(length(min_idx))-time_dwn(1));
    %             Ovrl_Measurements.mean_dFF.(IdChannel{o}).full.perminute(i) = mean(tmp(min_idx));
    %
    %             [~,idx] = intersect(tmp_locs,min_idx);
    %             Ovrl_Measurements.peaks_rate.(IdChannel{o}).full.perminute(i) = length(idx)/(time_dwn(length(min_idx))-time_dwn(1));
    %             Ovrl_Measurements.peaks.(IdChannel{o}).full.perminute{i} = tmp_peaks(idx);
    %             Ovrl_Measurements.peak_width.(IdChannel{o}).full.perminute{i} = tmp_w(idx);
    %             Ovrl_Measurements.peak_prom.(IdChannel{o}).full.perminute{i} = tmp_prom(idx);
    %
    %             Ovrl_Measurements.mean_peaks_val.(IdChannel{o}).full.perminute(i) = mean(Ovrl_Measurements.peaks.(IdChannel{o}).full.perminute{i});
    %             Ovrl_Measurements.mean_peaks_width.(IdChannel{o}).full.perminute(i) = mean(Ovrl_Measurements.peak_width.(IdChannel{o}).full.perminute{i});
    %             Ovrl_Measurements.mean_peaks_prom.(IdChannel{o}).full.perminute(i) = mean(Ovrl_Measurements.peak_prom.(IdChannel{o}).full.perminute{i});
    %
    %         end
    %
    %
    %         if ~isempty(Choice_idx)
    %             Ovrl_Measurements.norm_AUC.(IdChannel{o}).Choice.all = trapz(tmp(Choice_idx))/(time_dwn(length(Choice_idx))-time_dwn(1));
    %             Ovrl_Measurements.mean_dFF.(IdChannel{o}).Choice.all = mean(tmp(Choice_idx));
    %
    %             [~,idx] = intersect(tmp_locs,Choice_idx);
    %             Ovrl_Measurements.peaks_rate.(IdChannel{o}).Choice.all = length(idx)/(time_dwn(length(Choice_idx))-time_dwn(1));
    %             Ovrl_Measurements.peaks.(IdChannel{o}).Choice.all = tmp_peaks(idx);
    %             Ovrl_Measurements.peak_width.(IdChannel{o}).Choice.all = tmp_w(idx);
    %             Ovrl_Measurements.peak_prom.(IdChannel{o}).Choice.all = tmp_prom(idx);
    %
    %             Ovrl_Measurements.mean_peaks_val.(IdChannel{o}).Choice.all = mean(Ovrl_Measurements.peaks.(IdChannel{o}).Choice.all);
    %             Ovrl_Measurements.mean_peaks_width.(IdChannel{o}).Choice.all = mean(Ovrl_Measurements.peak_width.(IdChannel{o}).Choice.all);
    %             Ovrl_Measurements.mean_peaks_prom.(IdChannel{o}).Choice.all = mean(Ovrl_Measurements.peak_prom.(IdChannel{o}).Choice.all);
    %
    %             Ovrl_Measurements.norm_AUC.(IdChannel{o}).Choice.perminute = ones(length(idx_per_minute)-1,1)*nan;
    %             Ovrl_Measurements.mean_dFF.(IdChannel{o}).Choice.perminute = ones(length(idx_per_minute)-1,1)*nan;
    %             Ovrl_Measurements.peaks_rate.(IdChannel{o}).Choice.perminute = ones(length(idx_per_minute)-1,1)*nan;
    %             Ovrl_Measurements.peaks.(IdChannel{o}).Choice.perminute = cell(length(idx_per_minute)-1,1);
    %             Ovrl_Measurements.peak_width.(IdChannel{o}).Choice.perminute = cell(length(idx_per_minute)-1,1);
    %             Ovrl_Measurements.peak_prom.(IdChannel{o}).Choice.perminute = cell(length(idx_per_minute)-1,1);
    %             Ovrl_Measurements.mean_peaks_val.(IdChannel{o}).Choice.perminute = ones(length(idx_per_minute)-1,1)*nan;
    %             Ovrl_Measurements.mean_peaks_width.(IdChannel{o}).Choice.perminute = ones(length(idx_per_minute)-1,1)*nan;
    %             Ovrl_Measurements.mean_peaks_prom.(IdChannel{o}).Choice.perminute = ones(length(idx_per_minute)-1,1)*nan;
    %             for i = 1:length(idx_per_minute)-1
    %                 Choice_idx2 = intersect(Choice_idx,idx_per_minute(i):(idx_per_minute(i+1)-1));
    %                 if ~isempty(Choice_idx2)
    %                     Ovrl_Measurements.norm_AUC.(IdChannel{o}).Choice.perminute(i) = trapz(tmp(Choice_idx2))/(time_dwn(length(Choice_idx2))-time_dwn(1));
    %                     Ovrl_Measurements.mean_dFF.(IdChannel{o}).Choice.perminute(i) = mean(tmp(Choice_idx2));
    %
    %                     [~,idx] = intersect(tmp_locs,Choice_idx2);
    %                     Ovrl_Measurements.peaks_rate.(IdChannel{o}).Choice.perminute(i) = length(idx)/(time_dwn(length(Choice_idx2))-time_dwn(1));
    %                     Ovrl_Measurements.peaks.(IdChannel{o}).Choice.perminute{i} = tmp_peaks(idx);
    %                     Ovrl_Measurements.peak_width.(IdChannel{o}).Choice.perminute{i} = tmp_w(idx);
    %                     Ovrl_Measurements.peak_prom.(IdChannel{o}).Choice.perminute{i} = tmp_prom(idx);
    %
    %                     Ovrl_Measurements.mean_peaks_val.(IdChannel{o}).Choice.perminute(i) = mean(Ovrl_Measurements.peaks.(IdChannel{o}).Choice.perminute{i});
    %                     Ovrl_Measurements.mean_peaks_width.(IdChannel{o}).Choice.perminute(i) = mean(Ovrl_Measurements.peak_width.(IdChannel{o}).Choice.perminute{i});
    %                     Ovrl_Measurements.mean_peaks_prom.(IdChannel{o}).Choice.perminute(i) = mean(Ovrl_Measurements.peak_prom.(IdChannel{o}).Choice.perminute{i});
    %                 end
    %             end
    %         else
    %             Ovrl_Measurements.norm_AUC.(IdChannel{o}).Choice.all = [];
    %             Ovrl_Measurements.mean_dFF.(IdChannel{o}).Choice.all = [];
    %
    %             Ovrl_Measurements.peaks_rate.(IdChannel{o}).Choice.all = [];
    %             Ovrl_Measurements.peaks.(IdChannel{o}).Choice.all = [];
    %             Ovrl_Measurements.peak_width.(IdChannel{o}).Choice.all = [];
    %             Ovrl_Measurements.peak_prom.(IdChannel{o}).Choice.all = [];
    %
    %             Ovrl_Measurements.mean_peaks_val.(IdChannel{o}).Choice.all = [];
    %             Ovrl_Measurements.mean_peaks_width.(IdChannel{o}).Choice.all = [];
    %             Ovrl_Measurements.mean_peaks_prom.(IdChannel{o}).Choice.all = [];
    %
    %             Ovrl_Measurements.norm_AUC.(IdChannel{o}).Choice.perminute = [];
    %             Ovrl_Measurements.mean_dFF.(IdChannel{o}).Choice.perminute = [];
    %             Ovrl_Measurements.peaks_rate.(IdChannel{o}).Choice.perminute = [];
    %             Ovrl_Measurements.peaks.(IdChannel{o}).Choice.perminute = [];
    %             Ovrl_Measurements.peak_width.(IdChannel{o}).Choice.perminute = [];
    %             Ovrl_Measurements.peak_prom.(IdChannel{o}).Choice.perminute = [];
    %             Ovrl_Measurements.mean_peaks_val.(IdChannel{o}).Choice.perminute = [];
    %             Ovrl_Measurements.mean_peaks_width.(IdChannel{o}).Choice.perminute = [];
    %             Ovrl_Measurements.mean_peaks_prom.(IdChannel{o}).Choice.perminute = [];
    %         end
    %         if ~isempty(non_Choice_idx)
    %             Ovrl_Measurements.norm_AUC.(IdChannel{o}).non_Choice.all = trapz(tmp(non_Choice_idx))/(time_dwn(length(non_Choice_idx))-time_dwn(1));
    %             Ovrl_Measurements.mean_dFF.(IdChannel{o}).non_Choice.all = mean(tmp(non_Choice_idx));
    %
    %             [~,idx] = intersect(tmp_locs,non_Choice_idx);
    %             Ovrl_Measurements.peaks_rate.(IdChannel{o}).non_Choice.all = length(idx)/...
    %                 (time_dwn(length(non_Choice_idx))-time_dwn(1));
    %             Ovrl_Measurements.peaks.(IdChannel{o}).non_Choice.all = tmp_peaks(idx);
    %             Ovrl_Measurements.peak_width.(IdChannel{o}).non_Choice.all = tmp_w(idx);
    %             Ovrl_Measurements.peak_prom.(IdChannel{o}).non_Choice.all = tmp_prom(idx);
    %
    %             Ovrl_Measurements.mean_peaks_val.(IdChannel{o}).non_Choice.all = mean(Ovrl_Measurements.peaks.(IdChannel{o}).non_Choice.all);
    %             Ovrl_Measurements.mean_peaks_width.(IdChannel{o}).non_Choice.all = mean(Ovrl_Measurements.peak_width.(IdChannel{o}).non_Choice.all);
    %             Ovrl_Measurements.mean_peaks_prom.(IdChannel{o}).non_Choice.all = mean(Ovrl_Measurements.peak_prom.(IdChannel{o}).non_Choice.all);
    %
    %             Ovrl_Measurements.norm_AUC.(IdChannel{o}).non_Choice.perminute = ones(length(idx_per_minute)-1,1)*nan;
    %             Ovrl_Measurements.mean_dFF.(IdChannel{o}).non_Choice.perminute = ones(length(idx_per_minute)-1,1)*nan;
    %             Ovrl_Measurements.peaks_rate.(IdChannel{o}).non_Choice.perminute = ones(length(idx_per_minute)-1,1)*nan;
    %             Ovrl_Measurements.peaks.(IdChannel{o}).non_Choice.perminute = cell(length(idx_per_minute)-1,1);
    %             Ovrl_Measurements.peak_width.(IdChannel{o}).non_Choice.perminute = cell(length(idx_per_minute)-1,1);
    %             Ovrl_Measurements.peak_prom.(IdChannel{o}).non_Choice.perminute = cell(length(idx_per_minute)-1,1);
    %             Ovrl_Measurements.mean_peaks_val.(IdChannel{o}).non_Choice.perminute = ones(length(idx_per_minute)-1,1)*nan;
    %             Ovrl_Measurements.mean_peaks_width.(IdChannel{o}).non_Choice.perminute = ones(length(idx_per_minute)-1,1)*nan;
    %             Ovrl_Measurements.mean_peaks_prom.(IdChannel{o}).non_Choice.perminute = ones(length(idx_per_minute)-1,1)*nan;
    %             for i = 1:length(idx_per_minute)-1
    %                 non_Choice_idx2 = intersect(non_Choice_idx,idx_per_minute(i):(idx_per_minute(i+1)-1));
    %                 if ~isempty(non_Choice_idx2)
    %                     Ovrl_Measurements.norm_AUC.(IdChannel{o}).non_Choice.perminute(i) = trapz(tmp(non_Choice_idx2))/(time_dwn(length(non_Choice_idx2))-time_dwn(1));
    %                     Ovrl_Measurements.mean_dFF.(IdChannel{o}).non_Choice.perminute(i) = mean(tmp(non_Choice_idx2));
    %
    %                     [~,idx] = intersect(tmp_locs,non_Choice_idx2);
    %                     Ovrl_Measurements.peaks_rate.(IdChannel{o}).non_Choice.perminute(i) = length(idx)/(time_dwn(length(non_Choice_idx2))-time_dwn(1));
    %                     Ovrl_Measurements.peaks.(IdChannel{o}).non_Choice.perminute{i} = tmp_peaks(idx);
    %                     Ovrl_Measurements.peak_width.(IdChannel{o}).non_Choice.perminute{i} = tmp_w(idx);
    %                     Ovrl_Measurements.peak_prom.(IdChannel{o}).non_Choice.perminute{i} = tmp_prom(idx);
    %
    %                     Ovrl_Measurements.mean_peaks_val.(IdChannel{o}).non_Choice.perminute(i) = mean(Ovrl_Measurements.peaks.(IdChannel{o}).non_Choice.perminute{i});
    %                     Ovrl_Measurements.mean_peaks_width.(IdChannel{o}).non_Choice.perminute(i) = mean(Ovrl_Measurements.peak_width.(IdChannel{o}).non_Choice.perminute{i});
    %                     Ovrl_Measurements.mean_peaks_prom.(IdChannel{o}).non_Choice.perminute(i) = mean(Ovrl_Measurements.peak_prom.(IdChannel{o}).non_Choice.perminute{i});
    %                 end
    %             end
    %         else
    %             Ovrl_Measurements.norm_AUC.(IdChannel{o}).non_Choice.all = [];
    %             Ovrl_Measurements.mean_dFF.(IdChannel{o}).non_Choice.all = [];
    %
    %             Ovrl_Measurements.peaks_rate.(IdChannel{o}).non_Choice.all = [];
    %             Ovrl_Measurements.peaks.(IdChannel{o}).non_Choice.all = [];
    %             Ovrl_Measurements.peak_width.(IdChannel{o}).non_Choice.all = [];
    %             Ovrl_Measurements.peak_prom.(IdChannel{o}).non_Choice.all = [];
    %
    %             Ovrl_Measurements.mean_peaks_val.(IdChannel{o}).non_Choice.all = [];
    %             Ovrl_Measurements.mean_peaks_width.(IdChannel{o}).non_Choice.all = [];
    %             Ovrl_Measurements.mean_peaks_prom.(IdChannel{o}).non_Choice.all = [];
    %
    %             Ovrl_Measurements.norm_AUC.(IdChannel{o}).non_Choice.perminute = [];
    %             Ovrl_Measurements.mean_dFF.(IdChannel{o}).non_Choice.perminute = [];
    %             Ovrl_Measurements.peaks_rate.(IdChannel{o}).non_Choice.perminute = [];
    %             Ovrl_Measurements.peaks.(IdChannel{o}).non_Choice.perminute = [];
    %             Ovrl_Measurements.peak_width.(IdChannel{o}).non_Choice.perminute = [];
    %             Ovrl_Measurements.peak_prom.(IdChannel{o}).non_Choice.perminute = [];
    %             Ovrl_Measurements.mean_peaks_val.(IdChannel{o}).non_Choice.perminute = [];
    %             Ovrl_Measurements.mean_peaks_width.(IdChannel{o}).non_Choice.perminute = [];
    %             Ovrl_Measurements.mean_peaks_prom.(IdChannel{o}).non_Choice.perminute = [];
    %
    %         end
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
        save([path2save_mouse,'Session_analysis.mat'],'Stim_data','t_trials','Params')
    end
    clear Stim_data Ovrl_Measurements Measurements shock_evoked_response
end