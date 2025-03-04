 function [Condition,Index,Units_in_Condition_Titles,Variables,...
     Unit2Screen, Units_in_Condition] = CollectData(Variables...
     ,Index,Units_in_Condition_Titles, Unit2Screen,...
     Units_in_Condition,ConditionOrder) 
 % uses the function "GetMotif.m"
% Timestamps are in microseconds
% Waveforms are in AD - need to cnvert to uV.
%% Load unit data:
Condition.ConditionOrder=ConditionOrder;
Condition.ConditionName=Variables.Type(ConditionOrder);
Condition.UnitName=Variables.UnitName;
Condition.RealFileName=char(Variables.RealFileNames{Condition.ConditionOrder, 1});
% Make the path
Condition.DataPath= [Variables.ComputerDir,'\',Variables.Date...
    ,'\',Variables.MouseName,'\',Condition.RealFileName];
% get cell numbers from sorted file
% Initialize a variable to store the status message
warningMessage = '';
try
[Condition.TimeStamps_Allcells,~,Condition.CellNumbers, Condition.Features, Condition.Samples, Condition.DataHeader] = ...
Nlx2MatSpike([Condition.DataPath,'\TT',num2str(Variables.TetrodeNumber),'_s.ntt'], [1 1 1 1 1], 1, 1, [] ); %clustering file 
Condition.RawCell=Condition.TimeStamps_Allcells(Condition.CellNumbers==Variables.UnitNumber);
catch ME  % Catch the error or warning
    % Store the message in a variable
    warningMessage = ME.message;
try Units_in_Condition(Unit2Screen, Index) = {Condition.warningMessage}; catch end
Units_in_Condition_Titles(Unit2Screen, Index) = {[char(Condition.ConditionName),'.warningMessage']};Index=Index+1;
end
% only keep data relevant to the selected unit
Condition.UnitFeatures=Condition.Features(:,find(Condition.CellNumbers==Variables.UnitNumber));
Condition.UnitSamples=Condition.Samples(:,:,find(Condition.CellNumbers==Variables.UnitNumber));
% get continuous data
[Condition.Ch1.TimeStamps,Condition.Ch1.ChannelNumbers, Condition.Ch1.SampleFrequencies, Condition.Ch1.NumberOfValidSamples, Condition.Ch1.Samples, Condition.Ch1.Header] = Nlx2MatCSC([Condition.DataPath,'\CSC',num2str(4*(Variables.TetrodeNumber-1)+1),'.ncs'], [1 1 1 1 1], 1, 1, []); 
[Condition.Ch2.TimeStamps,Condition.Ch2.ChannelNumbers, Condition.Ch2.SampleFrequencies, Condition.Ch2.NumberOfValidSamples, Condition.Ch2.Samples, Condition.Ch2.Header] = Nlx2MatCSC([Condition.DataPath,'\CSC',num2str(4*(Variables.TetrodeNumber-1)+1),'.ncs'], [1 1 1 1 1], 1, 1, []); 
[Condition.Ch3.TimeStamps,Condition.Ch3.ChannelNumbers, Condition.Ch3.SampleFrequencies, Condition.Ch3.NumberOfValidSamples, Condition.Ch3.Samples, Condition.Ch3.Header] = Nlx2MatCSC([Condition.DataPath,'\CSC',num2str(4*(Variables.TetrodeNumber-1)+1),'.ncs'], [1 1 1 1 1], 1, 1, []); 
[Condition.Ch4.TimeStamps,Condition.Ch4.ChannelNumbers, Condition.Ch4.SampleFrequencies, Condition.Ch4.NumberOfValidSamples, Condition.Ch4.Samples, Condition.Ch4.Header] = Nlx2MatCSC([Condition.DataPath,'\CSC',num2str(4*(Variables.TetrodeNumber-1)+1),'.ncs'], [1 1 1 1 1], 1, 1, []); 
Ch1= reshape(Condition.Ch1.Samples, 1, []); % This will turn it into a 1xN vector
AllCSC=nan(4,length(Ch1));
AllCSC(1,:)=reshape(Condition.Ch1.Samples, 1, []);
AllCSC(2,:)=reshape(Condition.Ch2.Samples, 1, []);
AllCSC(3,:)=reshape(Condition.Ch3.Samples, 1, []);
AllCSC(4,:)=reshape(Condition.Ch4.Samples, 1, []);
Condition.AllCSC=AllCSC;
%% get information from the header:
% get ADBitVolts
lineWithADBitVolts = '';
% Loop through each line in DataHeader
for i = 1:length(Condition.DataHeader)
    currentLine = Condition.DataHeader{i};  % Get the current line
    % Check if the current line contains the string 'ADBitVolts'
    if contains(currentLine, 'ADBitVolts')
        lineWithADBitVolts = currentLine;
        break;  % Exit the loop once the line is found
    end
end
% Check if the line was found
if ~isempty(lineWithADBitVolts)
    % Extract the 4 numerical values from the line
    numericValues = regexp(lineWithADBitVolts, '\d+\.\d+', 'match');
        % Convert the extracted strings to numbers
    Condition.ADBitVolts = str2double(numericValues);
end
% get File created info 
% Initialize variables to store the date and time string
lineWithTimeCreated = '';
% Loop through each line in Jelly.DataHeader
for i = 1:length(Condition.DataHeader)
    currentLine = Condition.DataHeader{i};  % Get the current line
    % Check if the current line contains the string 'TimeCreated'
    if contains(currentLine, 'TimeCreated')
        lineWithTimeCreated = currentLine;
        break;  % Exit the loop once the line is found
    end
end
% Check if the line was found
if ~isempty(lineWithTimeCreated)
    % Extract the date and time from the line
    dateTimeStr = regexp(lineWithTimeCreated, '\d{4}/\d{2}/\d{2} \d{2}:\d{2}:\d{2}', 'match', 'once');
    % Convert the date and time string to MATLAB datetime format
    Condition.TimeCreated = datetime(dateTimeStr, 'InputFormat', 'yyyy/MM/dd HH:mm:ss');
end
%% get information about events:
[Condition.TimeStamps_events, ~,Condition.TTLs, ~, Condition.EventStrings, Condition.EventHeader] = ...
    Nlx2MatEV([Condition.DataPath,'\Events.nev'], [1 1 1 1 1], 1, 1, []); %event file
Condition.RecordingStart=Condition.TimeStamps_events((ismember(Condition.EventStrings,'Starting Recording')));
Condition.RecordingEnd=Condition.TimeStamps_events((ismember(Condition.EventStrings,'Stopping Recording')));
%% Figure out the TTLs
TTLtypesfound=unique(Condition.TTLs);
TTLtypesfound(TTLtypesfound==0)=[];
TTLNumbers=Condition.TTLs;
TTLTimestamps=Condition.TimeStamps_events;
%% classify the TTL types in file
   % % plot ttls over time
%     figure
%  for i=1:length(TTLtypesfound)
%   tempTTLtimestamps=(Condition.TimeStamps_events(Condition.TTLs==TTLtypesfound(i))-Condition.TimeStamps_events(1))/1000000;
%   scatter(tempTTLtimestamps,TTLtypesfound(i)*ones(size( tempTTLtimestamps))); 
%   hold on  
%  end
%  xlim([90 140])
 if ~isempty(TTLtypesfound)
       Condition.PiezoTTLType=[];
        Condition.VideoTTLType=[];
        Condition.LaserTTLType=[];
 % try to find the video TTL (it is given every five seconds)
for i=1:length(TTLtypesfound)
tempTTLtimestamps=(Condition.TimeStamps_events(Condition.TTLs==TTLtypesfound(i))-Condition.TimeStamps_events(1))/1000000;
TimestampDiffrence=[tempTTLtimestamps(2:end),tempTTLtimestamps(end)]-tempTTLtimestamps(1:end); 
tempSTDEV=std((TimestampDiffrence(1:end-1)));
if round(median(TimestampDiffrence(1:end-1)))==5 % if its 5 sec apart
    Condition.VideoTTLType=TTLtypesfound(i);
TTLtypesfound(TTLtypesfound==TTLtypesfound(i))=[];
        break
end
end
if isempty(Condition.VideoTTLType) %  video TTL not found, set to deafult number 20
Condition.VideoTTLType=20;
end
 % try to find the Laser TTL (it is given in a constant frequency of 2 hz)
for i=1:length(TTLtypesfound)
tempTTLtimestamps=(Condition.TimeStamps_events(Condition.TTLs==TTLtypesfound(i))-Condition.TimeStamps_events(1))/1000000;
TimestampDiffrence=[tempTTLtimestamps(2:end),tempTTLtimestamps(end)]-tempTTLtimestamps(1:end); 
if round(median(TimestampDiffrence(1:end-1)),2)==0.5 % if its 0.5 sec apart (2hz)
    Condition.LaserTTLType=TTLtypesfound(i);
    TTLtypesfound(TTLtypesfound==TTLtypesfound(i))=[];
    if Condition.ConditionOrder==3
Variables.TTLs{Condition.ConditionOrder}=Condition.LaserTTLType;
end
        break
end
end
if isempty(Condition.LaserTTLType) % no Laser TTL found, set to deafult number 20
Condition.LaserTTLType=20;
end
% the remaining TTLtype is Piezo (if exists)
if ~isempty(TTLtypesfound)
Condition.PiezoTTLType= TTLtypesfound  ;
if Condition.ConditionOrder<3
Variables.TTLs{Condition.ConditionOrder}=Condition.PiezoTTLType;
end
else
 Condition.PiezoTTLType=20;  
end
else
        Condition.PiezoTTLType=20;
        Condition.VideoTTLType=20;
        Condition.LaserTTLType=20;
 end
 
Condition.RawPiezo = Condition.TimeStamps_events(Condition.TTLs==Condition.PiezoTTLType);
% plot to confirm it makes sense
% TTLtypes = unique(Condition.TTLs);
% TTLtypes(TTLtypes==0)=[];
% TTLTimes = cell(length(TTLtypes), 1);  % Initialize as a cell array to store each type's times
% figure
% for Types = 1:length(TTLtypes)
%     % Get the indices of the events matching the current TTL type
%     matching_indices = Condition.TTLs == TTLtypes(Types);
%         % Store the timestamps for the current TTL type in the cell array
%     TTLTimes{Types} = Condition.TimeStamps_events(matching_indices);
% scatter(   TTLTimes{Types},Types*ones(size( TTLTimes{Types})));hold on
% end
% legend(num2str(TTLtypes(1)),num2str(TTLtypes(2)))

try Units_in_Condition(Unit2Screen, Index) = {Condition.VideoTTLType}; catch end
Units_in_Condition_Titles(Unit2Screen, Index) = {[char(Condition.ConditionName),'.VideoTTLType']};Index=Index+1;
try Units_in_Condition(Unit2Screen, Index) = {Condition.PiezoTTLType}; catch end
Units_in_Condition_Titles(Unit2Screen, Index) = {[char(Condition.ConditionName),'.PiezoTTLType']};Index=Index+1;
try Units_in_Condition(Unit2Screen, Index) = {Condition.LaserTTLType}; catch end
Units_in_Condition_Titles(Unit2Screen, Index) = {[char(Condition.ConditionName),'.LaserTTLType']};Index=Index+1;

%%
% read the log file to find the time this recording was created (doesnt
% matter which condition the cration date is the same
LogFile = fopen([Condition.DataPath,'\CheetahLogFile.txt'],'r');
LogFile=textscan(LogFile, '%s', 'Delimiter', '\n');
LogFile=LogFile{1, 1};  
%% find time created of Log file
lineWithTimeOpened = '';
% Loop through each line in Jelly.DataHeader
for i = 1:length(LogFile)
    currentLine = LogFile{i};  % Get the current line
    % Check if the current line contains the string 'TimeCreated'
    if contains(currentLine, 'Time Opened')
        lineWithTimeOpened = currentLine;
        break;  % Exit the loop once the line is found
    end
end
% Check if the line was found
if ~isempty(lineWithTimeOpened)
    % Extract the date and time from the line
    dateTimeStr = regexp(lineWithTimeOpened, '\d{4}/\d{2}/\d{2} \d{2}:\d{2}:\d{2}', 'match', 'once');
    % Convert the date and time string to MATLAB datetime format
    ExperimentCreated=datetime(dateTimeStr, 'InputFormat', 'yyyy/MM/dd HH:mm:ss');
end 
%% find time created of channel file
lineWithTimeCreated = '';
% Loop through each line in Jelly.DataHeader
for i = 1:length(Condition.Ch1.Header)
    currentLine = Condition.Ch1.Header{i};  % Get the current line
    % Check if the current line contains the string 'TimeCreated'
    if contains(currentLine, 'TimeCreated')
        lineWithTimeCreated = currentLine;
        break;  % Exit the loop once the line is found
    end
end
% Check if the line was found
if ~isempty(lineWithTimeCreated)
    % Extract the date and time from the line
    dateTimeStr = regexp(lineWithTimeCreated, '\d{4}/\d{2}/\d{2} \d{2}:\d{2}:\d{2}', 'match', 'once');
    % Convert the date and time string to MATLAB datetime format
    ChannelTimeCreated=datetime(dateTimeStr, 'InputFormat', 'yyyy/MM/dd HH:mm:ss');
end  
% Convert the creation times to microseconds since the start of the experiment
exp_start_microseconds = posixtime(ExperimentCreated) * 1e6;
channel_created_microseconds = posixtime(ChannelTimeCreated) * 1e6;
% Calculate the offset in microseconds between each channel and the experiment start
offset_microseconds = channel_created_microseconds - exp_start_microseconds;
Condition.FirstTimestamp=offset_microseconds;
Condition.TimeStamps_AllcellsNorm=Condition.TimeStamps_Allcells-Condition.FirstTimestamp; 
%% only keep data relevant to the selected unit
Condition.UnitTimeStamps_Unit=Condition.TimeStamps_Allcells(find(Condition.CellNumbers==Variables.UnitNumber)); 
Condition.TimeStamps_unit_Norm=Condition.UnitTimeStamps_Unit-Condition.FirstTimestamp; 
% screen for relevent events by event type
% only keep data relevant to the TTL type
%% Collect the data into the list
try Units_in_Condition(Unit2Screen, Index) = {Condition.FirstTimestamp}; catch end
Units_in_Condition_Titles(Unit2Screen, Index) = {[char(Condition.ConditionName),'.FirstTimestamp']};Index=Index+1;
try Units_in_Condition(Unit2Screen, Index) = {Condition.TimeStamps_unit_Norm}; catch end
Units_in_Condition_Titles(Unit2Screen, Index) = {[char(Condition.ConditionName),'.TimeStamps_unit_Norm']};Index=Index+1;
try Units_in_Condition(Unit2Screen, Index) = {Condition.UnitSamples}; catch end
Units_in_Condition_Titles(Unit2Screen, Index) = {[char(Condition.ConditionName),'.UnitSamples']};Index=Index+1;
if Condition.ConditionOrder<3
Condition.TimeStamps_Piezo=Condition.TimeStamps_events(:,find(Condition.TTLs==Condition.PiezoTTLType))';
Condition.TimeStamps_PizeoNorm=Condition.TimeStamps_Piezo-Condition.FirstTimestamp; %specific to the feeding events
try Units_in_Condition(Unit2Screen, Index) = {Condition.TimeStamps_PizeoNorm}; catch end
Units_in_Condition_Titles(Unit2Screen, Index) = {[char(Condition.ConditionName),'.TimeStamps_PizeoNorm']};Index=Index+1;
else
Condition.TimeStamps_Laser=Condition.TimeStamps_events(Condition.TTLs==Condition.LaserTTLType)';
Condition.TimeStamps_LaserNorm=Condition.TimeStamps_Laser-Condition.FirstTimestamp;
try Units_in_Condition(Unit2Screen, Index) = {Condition.TimeStamps_LaserNorm}; catch end
Units_in_Condition_Titles(Unit2Screen, Index) = {[char(Condition.ConditionName),'.TimeStamps_LaserNorm']};Index=Index+1;
end
% Enforce feeding detection criteria: At least 2 activations & max 6s gap
validFeedingEvents = [];
for i = 1:length(Condition.TimeStamps_Piezo) - 1
    % Calculate time difference between successive piezo activations
    time_diff = Condition.TimeStamps_Piezo(i+1) - Condition.TimeStamps_Piezo(i);

    % If activation occurs twice & within 6s, keep it
    if time_diff <= 6
        validFeedingEvents = [validFeedingEvents; Condition.TimeStamps_Piezo(i)];
    end
end

% Update feeding timestamps with valid ones
Condition.ValidTimeStamps_Piezo = validFeedingEvents;

if ConditionOrder<3
%% get information about video for feeding files:
[Condition] = GetMotif(Condition,Variables); %(EphysObj,Plot(true/false),Threshold,MinimalBoutTime(seconds)),20,4
try Units_in_Condition(Unit2Screen, Index) = {Condition.TimestampMotifs}; catch end
Units_in_Condition_Titles(Unit2Screen, Index) = {[char(Condition.ConditionName),'.TimestampMotifs']};Index=Index+1;
try Units_in_Condition(Unit2Screen, Index) = {Condition.MotifTTL}; catch end
Units_in_Condition_Titles(Unit2Screen, Index) = {[char(Condition.ConditionName),'.MotifTTL']};Index=Index+1;
end
end
