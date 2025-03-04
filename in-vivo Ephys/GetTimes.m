function [Statistics]= GetTimes(Unit,Variables)
    MotifsName={'Jelly','Empty','Rearing','TurnLeft','TurnRight','Stopping','Walking','Trotting','Running','Chow'};
    JellyIndex = find(strcmp(Unit.ConditionNames, 'Jelly'));
    ChowIndex = find(strcmp(Unit.ConditionNames, 'Chow'));
    LaserIndex = find(strcmp(Unit.ConditionNames, 'Laser'));
    % Extract necessary variables from the Unit structure
    UnitSpecificTimestamps = Unit.UnitSpecificTimestamps;  % Neural activity timestamps
    Range=[0,max(UnitSpecificTimestamps)];
    DLC_codes = Unit.DLC_codes;  % Event types (1, 2, 3, ..., 10)
    DLC_times = Unit.DLC_times;  % Event times corresponding to DLC_codes
    ConditionTypes=Unit.TTLTimestamps_PiezoType;
    % Additional event types (Piezo and Laser)
    Piezo_times = Unit.TTLTimestamps_Piezo;  % Timestamps for Piezo events
    Laser_times = Unit.TTLTimestamps_Laser;  % Timestamps for Laser events
    % Define the unique event types
    event_types = unique(DLC_codes);
    % deal with missing event types
    full_range = 1:10;
% Create an array with NaNs where values are missing
filled_array = NaN(1, length(full_range));
% Assign values from event_types where they match the full range
for i = 1:length(full_range)
    if ismember(full_range(i), event_types)
        filled_array(i) = full_range(i);
    end
end
 event_types=  filled_array; 
      % Loop through shared event types and get statistics 
for i = 1:10
    
        MotifName=MotifsName(i);
    if ismember(i, event_types)
        event_type = event_types(i);  
        UnitSpecificTimestampstemp=UnitSpecificTimestamps;
   %% Jelly video events
if event_type==1  
% Extract the event times for the current event type
        event_times_for_type = DLC_times(DLC_codes == event_type);  
        DLC_Jelly = DLC_times(DLC_codes == event_type);  
        % find the number of Jelly
% trim the unit timestamps to match the condition
        UnitSpecificTimestampstemp=UnitSpecificTimestampstemp(Unit.ConditionNumber==JellyIndex);  
        [Statistics(i)]=StatisticalTests(UnitSpecificTimestampstemp,event_times_for_type,char(MotifsName(event_type)),Range,Variables);

        % Apply DTW between UnitSpecificTimestamps and event_times_for_type
    %% Chow video events
% figure;scatter(event_times_for_type,ones(size(event_times_for_type)),'r|');hold on
% scatter(UnitSpecificTimestampstemp,2*ones(size(UnitSpecificTimestampstemp)),'k|');     
elseif event_type==10  
% Extract the event times for the current event type
        event_times_for_type = DLC_times(DLC_codes == event_type);  
        DLC_Chow = DLC_times(DLC_codes == event_type);  
% find the number of Chow
        UnitSpecificTimestampstemp=UnitSpecificTimestampstemp(Unit.ConditionNumber==ChowIndex);
[Statistics(i)]=StatisticalTests(UnitSpecificTimestampstemp,event_times_for_type,char(MotifsName(event_type)),Range,Variables);
%% All other video events
else             
  % Extract the event times for the current event type
        event_times_for_type = DLC_times(DLC_codes == event_type);  
[Statistics(i)]=StatisticalTests(UnitSpecificTimestampstemp,event_times_for_type,char(MotifsName(event_type)),Range,Variables);
end
    else % event type doesnt exist in this dataset
        try
EventTimes=[];
[Statistics(i)]=StatisticalTests(UnitSpecificTimestampstemp,EventTimes,...
    char(MotifsName(event_type)),Range,Variables);
        catch
        end
    end

end  
%% Jelly Piezo events
              UnitSpecificTimestampstemp=UnitSpecificTimestamps(Unit.ConditionNumber==JellyIndex);
Piezo_timesJelly=Piezo_times(Piezo_times<UnitSpecificTimestampstemp(end)&Piezo_times>UnitSpecificTimestampstemp(1));
    [Statistics(11)]=StatisticalTests(UnitSpecificTimestampstemp,Piezo_timesJelly,'Piezo_Jelly',Range,Variables); 
%%  Chow Piezo events
try
UnitSpecificTimestampstemp=UnitSpecificTimestamps(Unit.ConditionNumber==ChowIndex);
Piezo_timesChow=Piezo_times(Piezo_times<UnitSpecificTimestampstemp(end)&Piezo_times>UnitSpecificTimestampstemp(1));
[Statistics(12)]=StatisticalTests(UnitSpecificTimestampstemp,Piezo_timesChow,'Piezo_Chow',Range,Variables);
catch;end
%% Laser events
UnitSpecificTimestampstemp=UnitSpecificTimestamps(Unit.ConditionNumber==LaserIndex);
try
[Statistics(13)]=StatisticalTests(UnitSpecificTimestampstemp,Laser_times,'Laser',Range,Variables);
catch
    
end
%%
% % Correlation of Jelly and Chow detection with DLC and Piezo
% %% Jelly
% binSize = 1; % 1 second bins
% minTime = min([DLC_Jelly, Piezo_timesJelly]);
% maxTime = max([DLC_Jelly, Piezo_timesJelly]);
% timeEdges = minTime:binSize:maxTime;

% % Convert timestamps to binary event vectors
% DLC_binary = histcounts(DLC_Jelly, timeEdges) > 0; 
% Piezo_binary = histcounts(Piezo_timesJelly, timeEdges) > 0;
% % Calculate correlation on binary event vectors
% [Correlation.spearmanRhoJelly, Correlation.spearmanPJelly] = corr(DLC_binary(:), Piezo_binary(:), 'Type', 'Spearman', 'Rows', 'complete');
% %% Chow
% minTime = min([DLC_Chow; Piezo_timesChow]);
% maxTime = max([DLC_Chow; Piezo_timesChow]);
% timeEdges = minTime:binSize:maxTime;
% % Convert timestamps to binary event vectors
% DLC_binary = histcounts(DLC_Chow, timeEdges) > 0;
% Piezo_binary = histcounts(Piezo_timesChow, timeEdges) > 0;
% % Calculate correlation on binary event vectors
% [Correlation.spearmanRhoChow, Correlation.spearmanPChow] = corr(DLC_binary(:), Piezo_binary(:), 'Type', 'Spearman', 'Rows', 'complete');
end

