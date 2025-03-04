function [Condition] = VideoSync(Condition, Variables)
%% this code uses the LED signal in the arduino to figure out the video start time relative to the NLX recording
% Close all existing figure windows
close all;
% Variables
% Load the video
VideoFolder = fullfile(Variables.VideoPath, Condition.Movie_AVI);
video = VideoReader(VideoFolder);
% Number of frames to display and initialize array to store frame numbers
numFrames =length(Condition.MotifTTL); % Total number of frames
frameRate = video.FrameRate;
% Generate timestamps for each frame
timestampsVideoMSec = 1000000*(0:numFrames-1) / frameRate; % Array of timestamps in seconds
% Define the name of the csv file saved
originalFilename = Condition.Movie_AVI;
% Use fileparts to split the filename into parts
[filepath, name, ext] = fileparts(originalFilename);
% Check if the extension is '.avi'
if strcmp(ext, '.avi')
    % Construct the new filename with '_Piezo.csv' suffix
    newFilename = fullfile(filepath, [name, '_Piezo.csv']);
    newVarname=fullfile(filepath, [name, '_FrameStatus.mat']);
    newOffsetName=fullfile(filepath, [name, '_Offset.mat']);
else
    error('The file does not have a .avi extension.');
end
SaveFileName=[Variables.VideoPath,'\',newFilename];
SaveVarName=[Variables.VideoPath,'\',newVarname];
SaveOffsetName=[Variables.VideoPath,'\',newOffsetName];
%% if this was done for one of the units in this experiment, you can skip for the rest
for GetFrames=1:1
if  exist(SaveVarName, 'file')==2 || exist(SaveFileName, 'file')==2
        try
        load(SaveVarName, 'FrameStatus'); 
    catch
FrameStatus = readmatrix(SaveFileName);   
FrameStatus = logical(FrameStatus);
    end
else % skip this if file already exists

% Display the first frame for ROI selection
video.CurrentTime = 0;  % Reset video to the start time if needed
firstFrame = readFrame(video);
if size(firstFrame, 3) == 3
    firstFrame = rgb2gray(firstFrame);  % Convert to grayscale if needed
end
% Display the first frame and let user draw the ROI
imshow(firstFrame);
title('Draw a rectangle around the area of the Arduino LED');
roi = drawrectangle;  % Allows you to draw an ROI interactively
wait(roi);  % Wait until the rectangle is drawn and confirmed
% Get the position of the ROI as [x, y, width, height]
roiPosition = round(roi.Position);
frameNumbers = zeros(1, numFrames);  % Array to hold frame numbers
% Initialize storage for LED on and off frame data
ledOnData = struct();
ledOffData = struct();
% Collect all frames for computing pixel-wise baseline statistics
allFrames = [];
for idx = 1:numFrames
    if hasFrame(video)
        frame = readFrame(video);
        croppedFrame = imcrop(frame, roiPosition);
        if size(croppedFrame, 3) == 3
            croppedFrame = rgb2gray(croppedFrame);  % Convert to grayscale if needed
        end
        allFrames = cat(3, allFrames, croppedFrame);
%         figure; heatmap(croppedFrame);title(num2str(idx))
    else
        break;
    end
end
% Compute the pixel-wise mean and standard deviation across all frames
% Flatten all pixel values across all frames into a single vector
allPixels = double(allFrames(:));  % Convert to double for accuracy in calculations
% Calculate the overall mean and standard deviation of all pixel values
overallMean = mean(allPixels);
overallStd = std(allPixels);
% define threshold value
maxValues = squeeze(max(max(allFrames, [], 1), [], 2));
MeanmaxValues=mean(maxValues);
stdmaxValues=std(double(maxValues));
outlierThreshold = 3;
zScores = (maxValues - MeanmaxValues) / stdmaxValues;
ledOnFrames = find(abs(zScores) > outlierThreshold)+1;
%  figure;plot(maxValues);hold on; scatter(ledOnFrames,ones(size(ledOnFrames)),'r')
% % Convert LED-on frame indices to timestamps
% Get all frame indices (1 to total frames)
allFrameIndices = 1:numFrames;
% Find LED-off frames by removing LED-on frames from the total frames
% ledOnFrameSet = cell2mat(ledOnIndices); % Convert LED-on indices to array
ledOffFrames = setdiff(allFrameIndices, ledOnFrames)';
%% show a frame
video.FrameRate; %sets the video position to the start of the desiredFrameNumber frame.
if Variables.DisplayPlot
fig = figure('units', 'normalized', 'outerposition', [0 0 1 1]);  % Create a full-screen figure with 6 subplots% Randomly select 5 LED-on and 5 LED-off frames for display
sampleLedOnFrames = randsample(ledOnFrames, 5);
sampleLedOffFrames = randsample(ledOffFrames, 5);
for i = 1:5
    % Display LED-on sample frames
    video.CurrentTime = sampleLedOnFrames(i) / frameRate;
    frame = readFrame(video);
    croppedFrame = imcrop(frame, roiPosition);  % Crop to ROI
    subplot(2, 5, i);
    imshow(croppedFrame);
    title(sprintf('LED On: Frame %d', sampleLedOnFrames(i)));
    
    % Display LED-off sample frames
    video.CurrentTime = sampleLedOffFrames(i) / frameRate;
    frame = readFrame(video);
    croppedFrame = imcrop(frame, roiPosition);  % Crop to ROI
    subplot(2, 5, i + 5);
    imshow(croppedFrame);
    title(sprintf('LED Off: Frame %d', sampleLedOffFrames(i)));
end
sgtitle(char(Condition.ConditionName));
   print(fig, '-painters', '-dpdf', fullfile(Variables.UnitGeneralPath, ...
    'figures\',[char(Condition.ConditionName), num2str(Variables.TetrodeNumber), ...
    num2str(Variables.UnitNumber), 'DetectionOfLED.pdf']));
end
%% export list of frames as CSV
FrameStatus=false(length(allFrameIndices),1);
validLedOnFrames = ledOnFrames(ledOnFrames <= numFrames);
FrameStatus(validLedOnFrames) = true;
% Save FrameStatus as a CSV file
writematrix(FrameStatus, SaveFileName);
save(SaveVarName, 'FrameStatus');
end
end
Condition.FrameStatus=FrameStatus;
for AlignSignal=1:1
%% align the signals
% get timestamps in seconds of both signals 
if exist(SaveOffsetName, 'file')==2 
load(SaveOffsetName, 'totalAlignmentOffset');
else
    try
Piezo=double(Condition.RawPiezo-Condition.RawPiezo(1))/1000000;
LEDPiezo=Condition.TimestampMotifs(FrameStatus)';
ttlTimestamps = Condition.RawPiezo;  % Ensure this matches the variable in your workspace
LEDPiezoOriginals=find(FrameStatus);
%% Figure out overlap with piezo TTLs
% Step 1: Calculate the initial alignment offset to roughly align the video timestamps
initialOffset = Condition.Ch1.TimeStamps(1);
% Step 2: Align the LED timestamps initially
initialAlignedLEDTimestamps = timestampsVideoMSec(LEDPiezoOriginals) + initialOffset;
% Step 3: Perform cross-correlation using only the small LED and TTL datasets
% Define the min and max time range based on data in aligned LED and TTL timestamps
minTime = min([initialAlignedLEDTimestamps(1), ttlTimestamps(1)]);
maxTime = max([initialAlignedLEDTimestamps(end), ttlTimestamps(end)]);
% Generate binary signals for cross-correlation based on LED and TTL timestamps
ledSignal = ismember(round(initialAlignedLEDTimestamps), round(minTime:maxTime));  % LED signal
ttlSignal = ismember(round(ttlTimestamps), round(minTime:maxTime));                % TTL signal
% Perform cross-correlation to determine the refined alignment offset
try
    [correlation, lag] = xcorr(ttlSignal, ledSignal); 
[~, maxIndex] = max(correlation);
refinedLag = lag(maxIndex);  % Optimal lag for refined alignment
catch
    refinedLag=0;
end
% Step 4: Calculate the total alignment offset as a single number
totalAlignmentOffset = initialOffset + refinedLag;  % Single number for full alignment
save(SaveOffsetName, 'totalAlignmentOffset');
    catch
    initialOffset = Condition.Ch1.TimeStamps(1);
    refinedLag=0;
    totalAlignmentOffset = initialOffset + refinedLag;  % Single number for full alignment

    end
end
% Step 5: Apply this single offset to all video timestamps
alignedtimestampsVideoMSec = timestampsVideoMSec + totalAlignmentOffset;
Condition.alignedtimestampsVideoMSec=double(alignedtimestampsVideoMSec)';
% % Plot the aligned timestamps to check code
% Step 6 - plot
% % Define the interval used for cross-correlation (e.g., a 20-second gap)
% intervalStart = 60 * 1000;  % Start after the first minute (60 seconds, in ms)
% intervalDuration = 20 * 1000;  % Duration of 20 seconds (in ms)
% intervalEnd = intervalStart + intervalDuration;
% figure;
% hold on;
% % Plot all aligned video timestamps
% plot((finalAlignedVideoTimestamps(FrameStatus)-Condition.Ch1.TimeStamps  (1)/1000000), ones(size(finalAlignedVideoTimestamps)), 'bo', 'DisplayName', 'Aligned Video Timestamps'); 
% % Plot all TTL timestamps from the unit activity signal
% plot((ttlTimestamps-Condition.Ch1.TimeStamps  (1)/1000000), ones(size(ttlTimestamps)), 'rx', 'DisplayName', 'TTL Timestamps in Unit Activity');
% % Add plot labels and legend
% xlabel('Aligned Timestamps (ms)');
% ylabel('Signal Presence');
% title('Final Alignment of Video and Unit Activity Timestamps with Selected Time Gap');
% legend('show');
% hold off;
end
end