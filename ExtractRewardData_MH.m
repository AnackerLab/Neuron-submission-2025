function [Reward_times] = ExtractRewardData_Maryam...
    (files,mouseId,sessionId,substracted_onset,Behavior_code)

if  ~isempty(files)
    if contains(files(1).name,'\')
        T = readtable([files(1).folder,'\',files(1).name]);
    else
        T = readtable([files(1).folder,'/',files(1).name]);
    end

    %% Find the right column for this trial:
    ColumnNames = T.Properties.VariableNames;
    TrialColumn = strcmp(sessionId,ColumnNames);
    T = T(:,TrialColumn);
    T = table2array(T);
    T(T == 0) = [];
    T = T';
    Reward_times = T;
    Reward_times(isnan(Reward_times)) = [];

   
else
    Reward_times = [];
end
end