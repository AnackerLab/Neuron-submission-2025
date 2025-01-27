function [timestamps,Choice_score,Choice_epochs] = ExtractBehavData_Maryam...
    (files,mouseId,sessionId,substracted_onset,Behavior_code,Behav_data)

if  ~isempty(files)
    if contains(files(1).name,'\')
        T = readtable([files(1).folder,'\',files(1).name]);
    else
        T = readtable([files(1).folder,'/',files(1).name]);
    end
    timestamps = table2array(T(1:end,1));

    %% Find the right column for this trial:
    ColumnNames = T.Properties.VariableNames;
    TrialColumn = strcmp(sessionId,ColumnNames);
    T = T(:,TrialColumn);
    T = table2array(T);
    if T(1) == 0
        if strcmp(Behav_data,'Digging')
            idx = find(T==1);
            idx = idx(1);
            T(1:idx-1) = 7;
            [Behavior_code] = [Behavior_code,'PreTrial'];
        else
            idx = find(T==1);
            idx = idx(1);
            T(1:idx-1) = 8;
            [Behavior_code] = [Behavior_code,'PreTrial'];
        end
    end
    T(T == 0) = [];
    %% Comment if not working
    if ~isa(timestamps,'double')
        for o = length(timestamps):-1:1
            if ~any(regexp(timestamps{o},'[0-9]'))
                timestamps(o) = [];
            end
        end
        timestamps = str2double(timestamps);
    end
    %% Adjusting the length of timestamps
    idx = find(isnan(timestamps));
    if ~isempty(idx)
        timestamps = timestamps(1:idx(1)-1);
    end
    if length(timestamps) > length(T)
        timestamps(length(T)+1:end) = [];
    end
    
    idx = round(substracted_onset);
    T(1:idx) = [];
    timestamps(1:idx) = [];
    
    %% Get epochs for each behavior
    % Code for behavior:
    % Behaviour1 = 1; Behaviour2 = 2; Behaviour3 = 3; Behaviour4 = 4; Behaviour5 = 5;
    
    for o = 1:length(Behavior_code)
        BehavIdx.(Behavior_code{o}) = find(T == o);
        if ~isempty(BehavIdx.(Behavior_code{o}))
            A = diff(BehavIdx.(Behavior_code{o}));
            B = find(A > 1);
            C = [0;B;length(BehavIdx.(Behavior_code{o}))];
            if ~isempty(C)
                Epochs.(Behavior_code{o}) = ones(length(C)-1,2)*nan;
                for i = 1:length(C)-1
                    if BehavIdx.(Behavior_code{o})(C(i+1)) <= length(timestamps)
                        Epochs.(Behavior_code{o})(i,:) = [timestamps(BehavIdx.(Behavior_code{o})(C(i)+1)) ...
                            timestamps(BehavIdx.(Behavior_code{o})(C(i+1)))];
                    else
                        Epochs.(Behavior_code{o})(i,:) = [timestamps(BehavIdx.(Behavior_code{o})(C(i)+1)) timestamps(end)];
                    end
                end
            else
                Epochs.(Behavior_code{o}) = [];
            end
        else
            Epochs.(Behavior_code{o}) = [];
        end
    end
    Choice_score = T;
    Choice_epochs = Epochs;
else
    Choice_score = [];
    Choice_epochs = [];
end
end
