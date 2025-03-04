function [Unit]=ConnectDataFiles(Jelly,Chow,Laser,Variables)
ConditionList={Chow,Jelly,Laser};
ConditionNames={'Chow','Jelly','Laser'};
for ConditionIndex=1:3
Condition=ConditionList{1,ConditionIndex};    
% get File created info 
% Initialize variables to store the date and time string
%% find tetrode file created
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
    TetrodesTimeCreated(ConditionIndex)=datetime(dateTimeStr, 'InputFormat', 'yyyy/MM/dd HH:mm:ss');
end   
TetrodeFirst(ConditionIndex)=Condition.TimeStamps_events(1);
TetrodeLast(ConditionIndex)=Condition.TimeStamps_events(end);
%% find time created of Event file
lineWithTimeCreated = '';
% Loop through each line in Jelly.DataHeader
for i = 1:length(Condition.EventHeader)
    currentLine = Condition.EventHeader{i};  % Get the current line
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
    EventsTimeCreated(ConditionIndex)=datetime(dateTimeStr, 'InputFormat', 'yyyy/MM/dd HH:mm:ss');
end     
EventFirst(ConditionIndex)=Condition.Ch1.TimeStamps(1);
EventLast(ConditionIndex)=Condition.Ch1.TimeStamps(end);
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
    ChannelTimeCreated(ConditionIndex)=datetime(dateTimeStr, 'InputFormat', 'yyyy/MM/dd HH:mm:ss');
end  
ChannelFirst(ConditionIndex)=Condition.TimeStamps_events(1);
ChannelLast(ConditionIndex)=Condition.TimeStamps_events(end);
end
% read the log file to find the time this recording was created (doesnt
% matter which condition the cration date is the same
LogFile = fopen([Jelly.DataPath,'\CheetahLogFile.txt'],'r');
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
%% find out how to normalise the data:
% Convert the creation times to microseconds since the start of the experiment
exp_start_microseconds = posixtime(ExperimentCreated) * 1e6;
channel_created_microseconds = posixtime(ChannelTimeCreated) * 1e6;
% Calculate the offset in microseconds between each channel and the experiment start
offset_microseconds = channel_created_microseconds - exp_start_microseconds;
[~,SortOrder]=sort(offset_microseconds,'ascend');
Variables.SortOrder=SortOrder;
ConditionList=ConditionList(SortOrder);
ConditionNames=ConditionNames(SortOrder);
Variables.ConditionNames=ConditionNames;
offset_microseconds=offset_microseconds(SortOrder);
% figure;scatter(offset_microseconds,ones(size(offset_microseconds)));
% Connect the data in the correct way
UnitSpecificTimestamps=[];ConditionNumber=[];
Channel_1=[];Channel_2=[];Channel_3=[];Channel_4=[];
TTLTimestamps_Video=[];TTLTimestamps_Piezo=[];TTLTimestamps_Laser=[];TTLTimestamps_PiezoType=[];
DLC_codes=[];DLC_times=[];
for ConditionIndex=1:3
Condition=ConditionList{1, ConditionIndex};
offset=offset_microseconds(1);
% TTLTimestamps_VideoTemp=Condition.TimeStamps_events(Condition.TTLs==Condition.VideoTTLType)-offset;
if ~strcmp(Condition.ConditionName, 'Laser')
TempUnit=Condition.RawCell-offset;
TTLTimestamps_PiezoTemp=Condition.RawPiezo-offset;
TempTimestampMotifs=Condition.alignedtimestampsVideoMSec'-offset+ Variables.Factor*1e6;
% % plot to check time elingmnet
% figure; scatter(TempUnit,1*ones(size(TempUnit)));hold on
% scatter(TTLTimestamps_PiezoTemp,2*ones(size(TTLTimestamps_PiezoTemp)));hold on
% scatter(TempTimestampMotifs,3*ones(size(TempTimestampMotifs)));hold off
% legend('Unit','Piezo','Video')
if isempty(TempTimestampMotifs)
TempTimestampMotifs=Jelly.TimestampMotifs'* 1e6+Condition.Ch1.TimeStamps(1)-offset+ Variables.Factor*1e6;
end
try
Tempfiltered_MotifCodes=Condition.filtered_MotifCodes';   
if   strcmp(Condition.ConditionName, 'Chow')
Tempfiltered_MotifCodes(Tempfiltered_MotifCodes==1)=10;     
end     
catch
Tempfiltered_MotifCodes=[] ;   
end
if strcmp(Condition.ConditionName, 'Chow')
tempTTLTimestamps_PiezoType=2;
else
tempTTLTimestamps_PiezoType=1;
end 
TTLTimestamps_LaserTemp=[];
else
TempUnit=Condition.TimeStamps_Allcells(Condition.CellNumbers==Variables.UnitNumber)-offset;
TTLTimestamps_LaserTemp=Condition.TimeStamps_events(Condition.TTLs==Condition.LaserTTLType)-offset;
TTLTimestamps_PiezoTemp=[];
tempTTLTimestamps_PiezoType=[];
end
%% connect data
UnitSpecificTimestamps=[UnitSpecificTimestamps,TempUnit];
ConditionNumber=[ConditionNumber,ConditionIndex*ones(size(TempUnit))];
Channel_1=[Channel_1,Condition.Ch1.TimeStamps-offset];
Channel_2=[Channel_2,Condition.Ch2.TimeStamps-offset];
Channel_3=[Channel_3,Condition.Ch3.TimeStamps-offset];
Channel_4=[Channel_4,Condition.Ch4.TimeStamps-offset];
% TTLTimestamps_Video=[TTLTimestamps_Video,TTLTimestamps_VideoTemp];
if ~strcmp(Condition.ConditionName, 'Laser')
TTLTimestamps_Piezo=[TTLTimestamps_Piezo,TTLTimestamps_PiezoTemp];
TTLTimestamps_PiezoType=[TTLTimestamps_PiezoType,tempTTLTimestamps_PiezoType*ones(size(TTLTimestamps_PiezoTemp))];   
try
DLC_codes=[DLC_codes,Tempfiltered_MotifCodes];
catch
    DLC_codes=[DLC_codes,Tempfiltered_MotifCodes'];
end
try
DLC_times=[DLC_times,TempTimestampMotifs];
catch
DLC_times=[DLC_times,TempTimestampMotifs'];
end
else
TTLTimestamps_Laser=[TTLTimestamps_Laser,TTLTimestamps_LaserTemp];
end
end



FirstTime=Channel_1(1);
UnitSpecificTimestamps=UnitSpecificTimestamps-FirstTime;
Channel_1=Channel_1-FirstTime;
Channel_2=Channel_2-FirstTime;
Channel_3=Channel_2-FirstTime;
Channel_4=Channel_2-FirstTime;
TTLTimestamps_Video=TTLTimestamps_Video-FirstTime;
TTLTimestamps_Piezo=TTLTimestamps_Piezo-FirstTime;
TTLTimestamps_Laser=TTLTimestamps_Laser-FirstTime;
DLC_times=DLC_times-FirstTime;
%% save data to structure
Unit.UnitSpecificTimestamps=UnitSpecificTimestamps;
Unit.Channel_1=Channel_1;
Unit.TTLTimestamps_Video=TTLTimestamps_Video;
Unit.TTLTimestamps_Piezo=TTLTimestamps_Piezo;
Unit.TTLTimestamps_Laser=TTLTimestamps_Laser;
Unit.DLC_codes=DLC_codes;
Unit.DLC_times=DLC_times;
Unit.ConditionNumber=ConditionNumber;
Unit.TTLTimestamps_PiezoType=TTLTimestamps_PiezoType;
Unit.TTLTimestamps_PiezoType=TTLTimestamps_PiezoType;
Unit.ConditionNames=ConditionNames;
%% plot
if Variables.DisplayPlot
fig = figure('units', 'normalized', 'outerposition', [0 0 1 1]);  % Create a full-screen figure with 6 subplots% Randomly select 5 LED-on and 5 LED-off frames for display
scatter(UnitSpecificTimestamps(ConditionNumber==1),ones(size(UnitSpecificTimestamps(ConditionNumber==1))),'r|'); hold on
scatter(UnitSpecificTimestamps(ConditionNumber==2),ones(size(UnitSpecificTimestamps(ConditionNumber==2))),'g|'); hold on
scatter(UnitSpecificTimestamps(ConditionNumber==3),ones(size(UnitSpecificTimestamps(ConditionNumber==3))),'b|'); hold on
scatter(TTLTimestamps_Video,2*ones(size(TTLTimestamps_Video)),'|m'); hold on
scatter(TTLTimestamps_Piezo,2*ones(size(TTLTimestamps_Piezo)),'|c'); hold on
scatter(TTLTimestamps_Laser,2*ones(size(TTLTimestamps_Laser)),'|y'); hold on
scatter(DLC_times,DLC_codes+2,'|k')
% legend(char(ConditionNames(1)),char(ConditionNames(2)),char(ConditionNames(3))...
%     ,'Video','Piezo','Laser','DLC');
ylim([0 12]);
yticks(1:12);
YaxisName={'Unit FR','TTLs','Jelly','Empty','Rearing','TurnLeft','TurnRight','Stopping','Walking','Trotting','Running','Chow'};
set(gca, 'YTickLabel', YaxisName,'FontSize', 12, 'FontWeight', 'bold');
% Concatenate the path with the desired file name
fileName = [num2str(Variables.TetrodeNumber),num2str(Variables.UnitNumber),'overview.png'];
% Define the full save path
fullSavePath = fullfile(Variables.UnitGeneralPath, 'figures', fileName);
% Check if the 'figures' folder exists, and create it if it doesn't
figuresFolderPath = fullfile(Variables.UnitGeneralPath, 'figures');
if ~exist(figuresFolderPath, 'dir')
    mkdir(figuresFolderPath);
end
% Save the figure
sgtitle('Events Overview');
   print(fig, '-painters', '-dpdf', fullfile(Variables.UnitGeneralPath, ...
    'figures\',[num2str(Variables.TetrodeNumber), ...
    num2str(Variables.UnitNumber), 'EventsOverview.pdf']));
%fig = figure('units', 'normalized', 'outerposition', [0 0 1 1]);  % Create a full-screen figure with 6 subplots% Randomly select 5 LED-on and 5 LED-off frames for display


%% Get example images from the video files
fig = figure('units', 'normalized', 'outerposition', [0 0 1 1]);  % Create a full-screen figure with 6 subplots% Randomly select 5 LED-on and 5 LED-off frames for display
CodeNames={'Food','Empty','Rearing','TurnLeft','TurnRight','Stopping','Walking','Trotting','Running'};
frameCodes = Jelly.filtered_MotifCodes;  % Assuming this is your array with values from 1 to 10
% Initialize an array to hold the selected frame indices
selectedFrames = zeros(1, 9);  % One frame for each code from 1 to 9
% Loop through each unique number from 1 to 10 and randomly select a frame
for code = 1:9
    % Find indices in frameCodes that match the current code
    indices = find(frameCodes == code);   
        if ~isempty(indices)
        % Randomly select one index from these matching indices
        selectedFrames(code) = indices(randi(length(indices)));
    else
        % If no matching frame is found, skip this code
        disp(['Code ' num2str(code) ' not found in Unit.DLC_codes']);
    end
end
videoObj = VideoReader([Variables.VideoPath,'\',Jelly.Movie_AVI_Cropped]);
% Loop through each selected frame, read and display it
for i = 1:9
    try
       frameNumber = selectedFrames(i);
    videoObj.CurrentTime = (frameNumber - 1) / videoObj.FrameRate;  % Set to the specific frame time
    frame = readFrame(videoObj);  % Read the frame
    % Display frame in a subplot
    subplot(2, 5, i);  % Arrange subplots in a 2x5 grid
    imshow(frame);
    title([num2str(i), CodeNames(i)]);
        catch
        continue
    end
end
sgtitle('ExampleFrames Jelly');
   print(fig, '-painters', '-dpdf', fullfile(Variables.UnitGeneralPath, ...
    'figures\',[num2str(Variables.TetrodeNumber), ...
    num2str(Variables.UnitNumber), 'ExampleFrames_Jelly.pdf']));%%
%%
fig = figure('units', 'normalized', 'outerposition', [0 0 1 1]);  % Create a full-screen figure with 6 subplots% Randomly select 5 LED-on and 5 LED-off frames for display
frameCodes = Chow.filtered_MotifCodes;  % Assuming this is your array with values from 1 to 10
% Initialize an array to hold the selected frame indices
selectedFrames = zeros(1, 9);  % One frame for each code from 1 to 9
% Loop through each unique number from 1 to 10 and randomly select a frame
for code = 1:9
    % Find indices in frameCodes that match the current code
    indices = find(frameCodes == code);   
        if ~isempty(indices)
        % Randomly select one index from these matching indices
        selectedFrames(code) = indices(randi(length(indices)));
    else
        % If no matching frame is found, skip this code
        disp(['Code ' num2str(code) ' not found in Unit.DLC_codes']);
    end
end
videoObj = VideoReader([Variables.VideoPath,'\',Chow.Movie_AVI_Cropped]);
% Loop through each selected frame, read and display it
for i = 1:9
    try
       frameNumber = selectedFrames(i);
    videoObj.CurrentTime = (frameNumber - 1) / videoObj.FrameRate;  % Set to the specific frame time
    frame = readFrame(videoObj);  % Read the frame
    % Display frame in a subplot
    subplot(2, 5, i);  % Arrange subplots in a 2x5 grid
    imshow(frame);
    title([num2str(i), CodeNames(i)]);
    catch
        continue
    end
end
sgtitle('ExampleFrames Chow');
   print(fig, '-painters', '-dpdf', fullfile(Variables.UnitGeneralPath, ...
    'figures\',[num2str(Variables.TetrodeNumber), ...
    num2str(Variables.UnitNumber), 'ExampleFrames_Chow.pdf']));%%

end
end   
  





