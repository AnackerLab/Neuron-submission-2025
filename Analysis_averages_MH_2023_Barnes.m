close all; clear variables; clc;
set(0,'defaultfigurecolor',[1 1 1])

%% Define paths and data to analyze
path2data =  'C:\Users\Anacker1\Desktop\Fiber-Data\Fibers\MH\BARNES MAZE\DAY1\'; %Location of the data
path2savefolder =  'C:\Users\Anacker1\Desktop\Fiber-Data\Results\Barnes double check\channel corrected\Day 1\'; %Path to save

if exist([path2savefolder,'All_mice_average'],'dir') == 0
    mkdir([path2savefolder,'All_mice_average'])
end
path2savefolder = [path2savefolder,'All_mice_average','\'];

if exist([path2savefolder,'figures'],'dir') == 0
    mkdir([path2savefolder,'figures'])
end
path2save_figs = [path2savefolder,'figures','\'];


%% Define the list of mice to analyzed
% Insert it manually
% mice_list = {'KM_M3075','KM_MBLK','KM_M2727'};

% Determine it by the folders included in the path2data
mice_list = dir(path2data);
for o = length(mice_list):-1:1
    if mice_list(o).isdir == 0
        mice_list(o) = [];
    else
        if strcmp(mice_list(o).name,'.') == 1 || strcmp(mice_list(o).name,'..') == 1 || strcmp(mice_list(o).name,'All_mice_average') == 1 
            mice_list(o) = [];
        end
    end
end
Nmice = length(mice_list);

IdChannel = {'vCA1GCamp6f','dCA1GCamp6f'}; %%%%{'GCamp6f',[],'Psychlight'}vCA1GCamp6f dCA1GCamp6f
color2plotchannels = {'g','r'};

for i = 1:length(IdChannel)
    if ~isempty(IdChannel{i})
        limits2plot.(IdChannel{i}) = []; %If empty, automatically adjusted to the data.\
    end
end
limits2plot.all = [];
Color_scale = []; % For the heatmaps. If empty, automatically adjusted to the data.
LimitXaxis = [-3 3]; % Window to visualize in the plot. If empty, full length of TRANGE

show_plot = 1; % If 0, plots are not displayed
save_plot = 1; % If 0, plots are not saved
reanalysis = 1; % If 1, the code runs for sessions already analyze. If 0, session excluded if analysed
overwrite = 0; % If reanalyzing data, 1 if you want to save them and overwrite results and 0 if not

%% Define variables:
count = 0;
for m = 3
    if count == 0
        path2mouse = [path2data,mice_list(m).name,'\'];%%%%:Nmice

        load([path2mouse,'\Session_analysis.mat']) %%Session_analysis
        if ~isempty(Stim_data)
            count = count+1;
            Outcome = fieldnames(Stim_data);
            Condition = fieldnames(Stim_data.(Outcome{1}));
            normalization = fieldnames(Stim_data.(Outcome{1}).(Condition{1}).dFF.(IdChannel{1}));
            Trial_length = length(t_trials);
            measurements_fields = {'AUC','peak'};
            
            for outcome = 2:length(Outcome)
                if ~isempty(Stim_data.(Outcome{outcome}))
                    Conditions2analyze = fieldnames(Stim_data.(Outcome{outcome}));
                    for cond = 1:length(Conditions2analyze)
                        for channel = 1:length(IdChannel)
                            if ~isempty(IdChannel{channel})
                                for norm = 1:length(normalization)
                                    all_mice_trials.(Outcome{outcome}).(Conditions2analyze{cond}).(IdChannel{channel})...
                                        .(normalization{norm}) = nan(Trial_length,Nmice);
                                    all_mice_zscored_trials.(Outcome{outcome}).(Conditions2analyze{cond}).(IdChannel{channel})...
                                        .(normalization{norm}) = nan(Trial_length,Nmice);
                                    for fld = 1:length(measurements_fields)
                                        mean_alignment_measurements_per_mouse.(Outcome{outcome}).(Conditions2analyze{cond}).(IdChannel{channel}).(normalization{norm}).(measurements_fields{fld}) = nan(Nmice,1);
                                        mean_zscored_alignment_measurements_per_mouse.(Outcome{outcome}).(Conditions2analyze{cond}).(IdChannel{channel}).(normalization{norm}).(measurements_fields{fld}) = nan(Nmice,1);
                                    end
                                end
                            end
                        end
                    end
                end
            end
            
            fld = fieldnames(Measurements);
            for ii = 1:length(fld)
                if ~isempty(Measurements.(fld{ii}))
                    fld2 = fieldnames(Measurements.(fld{ii}));
                    for iii = 1:length(fld2)
                        for i = 1:length(IdChannel)
                            if ~isempty(IdChannel{i})
                                fld3 = fieldnames(Measurements.(fld{ii}).(fld2{iii}).(IdChannel{i}));
                                for iv = 1:length(fld3)
                                    pooled_all_mice_Measurements.(fld{ii}).(fld2{iii}).(IdChannel{i}).(fld3{iv}) = [];
                                    mean_Measurements_per_mouse.(fld{ii}).(fld2{iii}).(IdChannel{i}).(fld3{iv}) = nan(Nmice,1);

                                    pooled_all_mice_zscored_Measurements.(fld{ii}).(fld2{iii}).(IdChannel{i}).(fld3{iv}) = [];
                                    mean_zscored_Measurements_per_mouse.(fld{ii}).(fld2{iii}).(IdChannel{i}).(fld3{iv}) = nan(Nmice,1);
                                end
                            end
                        end
                    end
                end
            end

            fld = fieldnames(Overall_AUC);
            for ii = 1:length(fld)
                for i = 1:length(IdChannel)
                    mean_Overall_AUC_per_mouse.(fld{ii}).(IdChannel{i}) = nan(Nmice,1);
                    mean_zscored_Overall_AUC_per_mouse.(fld{ii}).(IdChannel{i}) = nan(Nmice,1);
                    mean_Overall_meandFF_per_mouse.(fld{ii}).(IdChannel{i}) = nan(Nmice,1);
                    mean_zscored_Overall_meandFF_per_mouse.(fld{ii}).(IdChannel{i}) = nan(Nmice,1);
                end
            end

            break
        end
    end
end

%% Get the data
for m = 1:Nmice
    %     if any(m == [2 4])
    % Define the path to the mouse and find the folders to analyze:
    path2mouse = [path2data,mice_list(m).name,'\'];

    load([path2mouse,'\Session_analysis.mat'])%%

    if ~isempty(Stim_data)
        Outcome = fieldnames(Stim_data);
        for outcome = 1:length(Outcome)
            if ~isempty(Stim_data.(Outcome{outcome}))
                Conditions2analyze = fieldnames(Stim_data.(Outcome{outcome}));
                for cond = 1:length(Conditions2analyze)
                    for channel = 1:length(IdChannel)
                        if ~isempty(IdChannel{channel})
                            for norm = 1:length(normalization)
                                if ~isempty(Stim_data.(Outcome{outcome}).(Conditions2analyze{cond}))
                                    all_mice_trials.(Outcome{outcome}).(Conditions2analyze{cond}).(IdChannel{channel})...
                                        .(normalization{norm})(:,m) = mean...
                                        (Stim_data.(Outcome{outcome}).(Conditions2analyze{cond})...
                                        .dFF.(IdChannel{channel}).(normalization{norm}),1,'omitnan')';
                                    all_mice_zscored_trials.(Outcome{outcome}).(Conditions2analyze{cond}).(IdChannel{channel})...
                                        .(normalization{norm})(:,m) = mean...
                                        (zscored_Stim_data.(Outcome{outcome}).(Conditions2analyze{cond})...
                                        .dFF.(IdChannel{channel}).(normalization{norm}),1,'omitnan')';
                                    if ~isempty(alignment_measurements.(Outcome{outcome}).(Conditions2analyze{cond}))
                                        for fld = 1:length(measurements_fields)
                                            mean_alignment_measurements_per_mouse.(Outcome{outcome}).(Conditions2analyze{cond}).(IdChannel{channel})...
                                                .(normalization{norm}).(measurements_fields{fld})(m) = mean...
                                                (alignment_measurements.(Outcome{outcome}).(Conditions2analyze{cond})...
                                                .(IdChannel{channel}).(normalization{norm}).(measurements_fields{fld}),1,'omitnan')';
                                            mean_zscored_alignment_measurements_per_mouse.(Outcome{outcome}).(Conditions2analyze{cond}).(IdChannel{channel})...
                                                .(normalization{norm}).(measurements_fields{fld})(m) = mean...
                                                (alignment_zscored_measurements.(Outcome{outcome}).(Conditions2analyze{cond})...
                                                .(IdChannel{channel}).(normalization{norm}).(measurements_fields{fld}),1,'omitnan')';
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
   
    if ~isempty(Measurements)
        fld = fieldnames(Measurements);
        for ii = 1:length(fld)
            if ~isempty(Measurements.(fld{ii}))
                fld2 = fieldnames(Measurements.(fld{ii}));
                for iii = 1:length(fld2)
                    for i = 1:length(IdChannel)
                        if ~isempty(IdChannel{i})
                            fld3 = fieldnames(Measurements.(fld{ii}).(fld2{iii}).(IdChannel{i}));
                            for iv = 1:length(fld3)
                                [pooled_all_mice_Measurements.(fld{ii}).(fld2{iii}).(IdChannel{i}).(fld3{iv})] = ...
                                    [pooled_all_mice_Measurements.(fld{ii}).(fld2{iii}).(IdChannel{i}).(fld3{iv});...
                                    Measurements.(fld{ii}).(fld2{iii}).(IdChannel{i}).(fld3{iv})];
                                mean_Measurements_per_mouse.(fld{ii}).(fld2{iii}).(IdChannel{i}).(fld3{iv})(m) = ...
                                    mean(Measurements.(fld{ii}).(fld2{iii}).(IdChannel{i}).(fld3{iv}),'omitnan');

                                [pooled_all_mice_zscored_Measurements.(fld{ii}).(fld2{iii}).(IdChannel{i}).(fld3{iv})] = ...
                                    [pooled_all_mice_zscored_Measurements.(fld{ii}).(fld2{iii}).(IdChannel{i}).(fld3{iv});...
                                    zscored_Measurements.(fld{ii}).(fld2{iii}).(IdChannel{i}).(fld3{iv})];
                                mean_zscored_Measurements_per_mouse.(fld{ii}).(fld2{iii}).(IdChannel{i}).(fld3{iv})(m) = ...
                                    mean(zscored_Measurements.(fld{ii}).(fld2{iii}).(IdChannel{i}).(fld3{iv}),'omitnan');
                            end
                        end
                    end
                end
            end
        end
    end

    fld = fieldnames(Overall_AUC);
    for ii = 1:length(fld)
        if ~isempty(Overall_AUC.(fld{ii}))
            for i = 1:length(IdChannel)
                mean_Overall_AUC_per_mouse.(fld{ii}).(IdChannel{i})(m) = mean(Overall_AUC.(fld{ii}).(IdChannel{i}),'omitnan');
                mean_zscored_Overall_AUC_per_mouse.(fld{ii}).(IdChannel{i})(m) = mean(zscored_Overall_AUC.(fld{ii}).(IdChannel{i}),'omitnan');
                mean_Overall_meandFF_per_mouse.(fld{ii}).(IdChannel{i})(m) = mean(Overall_meandFF.(fld{ii}).(IdChannel{i}),'omitnan');
                mean_zscored_Overall_meandFF_per_mouse.(fld{ii}).(IdChannel{i})(m) = mean(zscored_Overall_meandFF.(fld{ii}).(IdChannel{i}),'omitnan');
            end
        end
    end
end

%% Save the results for all sessions for each mouse
%     save([path2save_mouse,'All_sessions_results.mat'],'trials','pressing_delay','ITI')


%% Plot the results
% trial_mode = {'raw','baseline_corrected'}; % raw or baseline_corrected.
% Condition = fieldnames(all_mice_trials);%{'LeverExtension','DipperOn','Reward'};
% avg_mode = 1; % 1, use mean and SEM; 2, use median and MAD
% for cond = 1:length(Condition)
%     for t_mode = 2%1:length(trial_mode)
%         session_names = fieldnames(all_mice_trials.(Conditions2analyze{cond}).dFF.(IdChannel{channel})...
%             .(normalization{norm}));
%         for s = 1:length(session_names)
%             plot_group_data_Ashlea(all_mice_trials,t_trials,Condition{cond},session_names{s},IdChannel,trial_mode{t_mode},avg_mode,limits2plot,Color_scale,show_plot,save_plot,path2save_figs)
%         end
%     end
% end

trial_mode = {'raw','baseline_corrected','zscored'}; % raw or baseline_corrected.
avg_mode = 1; % 1, use mean and SEM; 2, use median and MAD
order_crit = 0;
Outcome = {'Correct'};
for outcome = 1:length(Outcome)
    if ~isempty(all_mice_trials.(Outcome{outcome}))
        Condition = fieldnames(all_mice_trials.(Outcome{outcome}));%{'LeverExtension','DipperOn','Reward'};
        for cond = 1:length(Condition)
            for t_mode = 1:length(trial_mode)
                plot_group_data_Maryam(all_mice_trials.(Outcome{outcome}),t_trials,Outcome{outcome},Condition{cond},IdChannel,trial_mode{t_mode},avg_mode,color2plotchannels,order_crit,limits2plot,Color_scale,show_plot,save_plot,path2save_figs)
                plot_zscored_group_data_Maryam(all_mice_zscored_trials.(Outcome{outcome}),t_trials,Outcome{outcome},Condition{cond},IdChannel,trial_mode{t_mode},avg_mode,color2plotchannels,order_crit,limits2plot,Color_scale,show_plot,save_plot,path2save_figs)
            end
        end
    end
end

save([path2savefolder,'All_mice_results.mat'],'all_mice_trials','all_mice_zscored_trials','t_trials',...
    'mean_Measurements_per_mouse','mean_zscored_Measurements_per_mouse','pooled_all_mice_Measurements',...
    'pooled_all_mice_zscored_Measurements','mean_Overall_AUC_per_mouse','mean_Overall_meandFF_per_mouse',...
    'mean_zscored_Overall_AUC_per_mouse','mean_zscored_Overall_meandFF_per_mouse')
