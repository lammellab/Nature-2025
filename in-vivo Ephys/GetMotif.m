 function [Condition] = GetMotif(Condition,Variables) 
%% find video relevant files
% Define the real-world dimensions of the box (in cm)
real_width = 30;  %  cm (X-axis)%%%%% real size 30x17 cm!!! not estimated, real
real_height = 17; %  cm (Y-axis)
% Define thresholds
Threshold_Body2Head=[0 3]; 
Threshold_Plate2Head=[0 5]; 
Threshold_Velocity.Stopping=[0,1];
Threshold_Velocity.Walking=[1,5];
Threshold_Velocity.Trotting=[5,10];
Threshold_Velocity.Running=[10,50];
turn_threshold = 15;% deg
VideoDirectory=dir(Variables.VideoPath);VideoDirectory=extractfield(VideoDirectory,'name')';VideoDirectory=VideoDirectory(3:end);
% Get the specific condition name for this order and remove spaces
ConditionName = strrep(Variables.ConditionName{1, Condition.ConditionOrder}, ' ', ''); % Remove spaces
% Loop through the VideoDirectory
for i = 1:length(VideoDirectory)
    currentFile = VideoDirectory{i}; % Get the current file name
    % 1. Files containing "DLC", ConditionName, and ".csv"
    if contains(lower(currentFile), 'dlc') && contains(lower(currentFile), lower(ConditionName)) && contains(lower(currentFile), '.csv')
        Condition.DLC_CSV = currentFile;
    % 2. Files containing ConditionName and ".csv", but not "DLC"
    elseif contains(lower(currentFile), lower(ConditionName)) && contains(lower(currentFile), '.csv') && ~contains(lower(currentFile), 'dlc')
        Condition.biobserveCSV = currentFile;
    % 3. Files containing ConditionName and ".avi" (case insensitive for file extension)
elseif contains(lower(currentFile), lower(ConditionName)) && contains(lower(currentFile), '.avi') && contains(lower(currentFile), 'cropped')
        Condition.Movie_AVI_Cropped = currentFile;
        % 4. Files containing ConditionName and ".avi" (case insensitive for file extension)
    elseif contains(lower(currentFile), lower(ConditionName)) && contains(lower(currentFile), '.avi')
        Condition.Movie_AVI = currentFile;
      %5. Condition.Video.VideoMp4=VideoDirectory(FindRelevantFiles,1).name;
 elseif contains(lower(currentFile), lower(ConditionName)) && contains(lower(currentFile), 'mp4')
        Condition.Movie_Mp4 = currentFile;
    end
end
% look at the track from bioobserve
try
Condition.Track = biobserve_import([Variables.VideoPath,'\',Condition.biobserveCSV]);
catch
end
% other variables
timestamp_interval = 5; %Sec - this is the interval of the TTL signal from the bioobserve to NLX
% Read the entire CSV file
DLCInfo = readtable([Variables.VideoPath,'\',Condition.DLC_CSV], 'ReadVariableNames', false);
% Assign the processed headers to the table
VariableNames = {'bodyparts', 'NoseX', 'NoseY', 'NoseLikelihood', ...
                 'ImplantX', 'ImplantY', 'ImplantLikelihood', ...
                 'BodyX', 'BodyY', 'BodyLikelihood', ...
                 'TailX', 'TailY', 'TailLikelihood', ...
                 'TLBoxX', 'TLBoxY', 'TLBoxLikelihood', ...
                 'BRBoxX', 'BRBoxY', 'BRBoxLikelihood'};

% Assign the VariableNames as the column headers for the table
DLCInfo.Properties.VariableNames = VariableNames;
Condition.DLCInfo=DLCInfo;
InfoArray=table2array(DLCInfo);
% read the video file
VideoInfo=VideoReader([Variables.VideoPath,'\',Condition.Movie_AVI_Cropped]);Condition.VideoInfo=VideoInfo;
Condition.VideoFrameTimes=linspace(0,Condition.VideoInfo.Duration,Condition.VideoInfo.NumFrames); % X values for frame basd analysis
Condition.FrameEdgesSec=linspace(0,Condition.VideoFrameTimes(end)*2-mean([Condition.VideoFrameTimes(end),Condition.VideoFrameTimes(end-1)]),Condition.VideoInfo.NumFrames+1);
Condition.FrameEdgesMin=Condition.FrameEdgesSec/60;
% Identify the X and Y position columns
position_columns = contains(VariableNames, {'X', 'Y'});
% Get likelihood columns (assuming they end with 'Likelihood')
likelihood_columns = contains(VariableNames, 'Likelihood');
% Loop through each body part to replace low-likelihood positions with NaN
for i = 1:sum(contains(VariableNames, 'Likelihood'))
    % Get the corresponding position columns for this body part
    x_column = find(position_columns, 1, 'first') + (i-1)*2;  % X column
    y_column = x_column + 1;      % Y column
    
    % Get the likelihood column for this body part
    likelihood_column = find(likelihood_columns, 1, 'first') + i - 1;
    
    % Identify frames with low likelihood
    low_likelihood_frames = DLCInfo{:, likelihood_column} < 0.85;% threshold for including data
    
    % Replace X and Y positions with NaN where likelihood is low
    DLCInfo{low_likelihood_frames, [x_column, y_column]} = NaN;
end
DLCInfoSmooth=DLCInfo;
% Smooth the X and Y data to account for natural mouse movement
% DLCInfoSmooth{:, contains(VariableNames, {'X', 'Y'})} = smoothdata(DLCInfo{:, contains(VariableNames, {'X', 'Y'})}, 'movmean', 5);
Condition.DLCInfoSmooth=DLCInfoSmooth;
% define the arena
% Extract the X and Y positions
implant_X = DLCInfoSmooth.ImplantX;
implant_Y = DLCInfoSmooth.ImplantY;
body_X = DLCInfoSmooth.BodyX;
body_Y = DLCInfoSmooth.BodyY;
% Calculate the minimum and maximum X and Y to define the edges of movement
min_X = min([implant_X; body_X], [], 'omitnan'); % Left edge
max_X = max([implant_X; body_X], [], 'omitnan'); % Right edge
min_Y = min([implant_Y; body_Y], [], 'omitnan'); % Bottom edge
max_Y = max([implant_Y; body_Y], [], 'omitnan'); % Top edge
% Define the four corners of the rectangular movement area
Condition.Arena.bottom_left = [min_X, max_Y];
Condition.Arena.bottom_right = [max_X, max_Y];
Condition.Arena.top_left = [min_X, min_Y];
Condition.Arena.top_right = [max_X, min_Y];
% Find out the locations of the boxes
TLBox_X = DLCInfoSmooth.TLBoxX;
TLBox_Y = DLCInfoSmooth.TLBoxY;
BRBox_X = DLCInfoSmooth.BRBoxX;
BRBox_Y = DLCInfoSmooth.BRBoxY;
% Calculate the average locations of the Empty Box and Food Box
Condition.Center_TLBox = [median(TLBox_X, 'omitnan'),median(TLBox_Y, 'omitnan')];
Condition.Center_BRBox = [median(BRBox_X, 'omitnan'),median(BRBox_Y, 'omitnan')];
% Extract the X and Y coordinates of the Empty Box and Food Box
TLBox_coords = Condition.Center_TLBox;
BRBox_coords = Condition.Center_BRBox;
% Calculate the Euclidean distance between Empty Box and Food Box
Condition.PlateDistance = sqrt((BRBox_coords(1) - TLBox_coords(1))^2 + (BRBox_coords(2) - TLBox_coords(2))^2);
v = VideoReader([Variables.VideoPath,'\',Condition.Movie_AVI_Cropped]);
% Show the frame
fig=figure;
imshow(read(v, 500));
hold on; 
% plot(implant_X, implant_Y, '-b', 'DisplayName', 'Implant'); % Blue circles 
% plot(body_X, body_Y, '-r', 'DisplayName', 'Implant'); % Blue circles 
scatter(Condition.Arena.top_left(1),Condition.Arena.top_left(2)); hold on%bl
scatter(Condition.Arena.top_right(1),Condition.Arena.top_right(2)); hold on%br
scatter(Condition.Arena.bottom_left(1),Condition.Arena.bottom_left(2)); hold on%tl
scatter(Condition.Arena.bottom_right(1),Condition.Arena.bottom_right(2)); hold on%tr
scatter(Condition.Center_TLBox(1),Condition.Center_TLBox(2)); hold on
scatter(Condition.Center_BRBox(1) ,Condition.Center_BRBox(2)); hold on
legend('ArenaTL','ArenaTR','ArenaBL','ArenaBR','PlateTL', 'PlateBR')
fileName = [Variables.MouseName,' ',Variables.Date,' ',num2str(Variables.TetrodeNumber),...
num2str(Variables.UnitNumber),'_',char(Condition.ConditionName),' Plate Placement.pdf'];
fullSavePath = fullfile(Variables.UnitGeneralPath, 'figures', fileName);
%% Manual selection of plate location (Optional step - in case DLC did not do a good job) 
for DefinePlate=1:1
% %% enter plate location manually if its off
% Define the path for saving/loading plate coordinates
%%
saveFilePath = fullfile(Variables.VideoPath, 'plate_coords.mat');
% Check if the file with saved coordinates exists
if isfile(saveFilePath)
    % Load saved coordinates
    load(saveFilePath, 'BRBox_coords', 'TLBox_coords');
    disp('Loaded previously saved plate coordinates.');
else
    % Manual selection of BR plate
    title('Select BR plate');
    [x, y] = ginput(1); % Select Bottom-Right plate
    BRBox_coords = [x, y];

    % Manual selection of TL plate
    title('Select TL plate');
    [x, y] = ginput(1); % Select Top-Left plate
    TLBox_coords = [x, y];

    % Save the coordinates to a file for future use
    save(saveFilePath, 'BRBox_coords', 'TLBox_coords');
    disp(['Plate coordinates saved to: ', saveFilePath]);
end
end 
saveas(fig, fullSavePath);  
for FindMotifs=1:1
%% Calculate Parameters:
ParametersArray=nan(height(Condition.DLCInfo),4);
ParametersArray(:,1)=sqrt((body_X-implant_X).^2 + (body_Y-implant_Y).^2);
ParametersArray(:,2)=sqrt((implant_X-Condition.Center_BRBox(1)).^2 + (implant_Y-Condition.Center_BRBox(2)).^2);
ParametersArray(:,3)=sqrt((implant_X-Condition.Center_TLBox(1)).^2 + (implant_Y-Condition.Center_TLBox(2)).^2);
%% DISTANCE of body between frames:
Futurebody_X=[body_X(2:end);0];
Futurebody_Y=[body_Y(2:end);0];
ParametersArray(:,4)=sqrt(([body_X(2:end);0]-body_X).^2 + ([body_Y(2:end);0]-body_Y).^2);
ParametersArray(end,4)=0;
% Convert to real world distances using the scaling factor: 
% Define the top-left and bottom-right coordinates of the box in pixels
TLBox_coords = [min(body_X), min(body_Y)];  % Top-left corner in pixels
BRBox_coords = [max(body_X), max(body_Y)]; % Bottom-right corner in pixels
% Calculate the pixel dimensions of the box
pixel_width = BRBox_coords(1) - TLBox_coords(1);   % Width in pixels (X-axis)
pixel_height = BRBox_coords(2) - TLBox_coords(2);  % Height in pixels (Y-axis)
% Calculate the diagonal distance in pixels (Pythagorean theorem)
diagonal_pixels = sqrt(pixel_width^2 + pixel_height^2);
% Calculate the diagonal distance in real-world size (cm)
diagonal_real = sqrt(real_width^2 + real_height^2);  % Diagonal in cm
% Calculate the scaling factor (cm per pixel)
scaling_factor = diagonal_real / diagonal_pixels;
% convert the distances form pixel to cm:
ParametersArray=ParametersArray*scaling_factor;
% convert the distance moved between frames to cm/sec
ParametersArray(:,5)=VideoInfo.FrameRate*ParametersArray(:,4);
%% ANGLE OF HEAD AND BODY:
% make a vector from the head and body position
num_frames = length(body_X); % Total number of frames
angles = zeros(1, num_frames - 1); % Preallocate space for angles between frames
for i = 1:num_frames-1
      % Calculate the vector for the current frame using implant and body coordinates
    v1 = [implant_X(i) - body_X(i), implant_Y(i) - body_Y(i)];
    
    % Calculate the vector for the next frame using implant and body coordinates
    v2 = [implant_X(i+1) - body_X(i+1), implant_Y(i+1) - body_Y(i+1)];
    
    % Normalize the vectors to make them unit vectors (for accurate angle calculation)
    v1 = v1 / norm(v1);
    v2 = v2 / norm(v2);
    
    % Calculate the signed angle using atan2 (2D cross product)
    cross_product = v1(1) * v2(2) - v1(2) * v2(1); % 2D cross product
    dot_product = dot(v1, v2); % Dot product of the two vectors
    
    % Compute the angle using atan2 to get signed angles
    angle_rad = atan2(cross_product, dot_product);
       
        % Convert the angle to degrees
    angles(i) = rad2deg(real(angle_rad)); % Use real() to discard small imaginary parts
end
angles=[angles';0];
ParametersArray(:,6)=angles;
% work on the angle data to detect real turning events and not just move head
% Smooth the angle data using a moving average
window_size = 5; % Define the size of the moving window
smoothed_angles = movmean(angles, window_size);
% Detect sustained turns for left and right separately - logic:
% Turn Threshold: We define a threshold for the angle to determine if a turn is happening. 
% For example, if the angle between two consecutive frames exceeds a certain value 
%(e.g., 2 degrees), it is considered a "turn."
% Positive angles (> turn_threshold) indicate a left turn.
% Negative angles (< -turn_threshold) indicate a right turn.
% Sustained Turn: A real turn is considered "sustained" if the mouse continues 
% turning in the same direction for multiple consecutive frames 
%(in our case, at least 3 frames). We want to avoid detecting brief, 
%momentary head movements as actual turns.
% Define a threshold for real turns (e.g., 2 degrees)
% Identify frames where the mouse is turning (left or right)
turns_right = smoothed_angles > turn_threshold;  % Positive angles -> left turn
turns_left = smoothed_angles < -turn_threshold; % Negative angles -> right turn
% Initialize the event array (0 for no turn, -1 for left, 1 for right)
event_array = zeros(1, length(smoothed_angles));
min_turn_duration = 2; % Minimum number of consecutive frames for a sustained turn
% Mark left turns (-1) in the event array
count = 0;
for i = 1:length(turns_left)
    if turns_left(i)
        if count == 0  % Mark the start of a turn
            start_frame = i;
        end
        count = count + 1;
    else
        if count >= min_turn_duration
            event_array(start_frame:i-1) = -1; % Mark the left turn in the array
        end
        count = 0;  % Reset count if turn ends
    end
end
% Catch the last left event if it ends at the last frame
if count >= min_turn_duration
    event_array(start_frame:end) = -1;
end
% Mark right turns (1) in the event array
count = 0;
for i = 1:length(turns_right)
    if turns_right(i)
        if count == 0  % Mark the start of a turn
            start_frame = i;
        end
        count = count + 1;
    else
        if count >= min_turn_duration
            event_array(start_frame:i-1) = 1; % Mark the right turn in the array
        end
        count = 0;  % Reset count if turn ends
    end
end
% Catch the last right event if it ends at the last frame
if count >= min_turn_duration
    event_array(start_frame:end) = 1;
end
ParametersArray(:,7)=event_array';
ParametersArrayNames={'Head2Body', 'Head2BR','Head2TL','BodyFromLastFrame','Velocity_cm/sec','head/body angle','Turn'};
Motifs=nan(height(ParametersArray),1);
MotifsName={'Food','Empty','Rearing','TurnLeft','TurnRight','Stopping','Walking','Trotting','Running'};
for FrameNumber=1:height(ParametersArray)
% check if we are close to any plate, if they are in a small range of a
% plate but not close enough to deterine feeding, ignore that frame
   Distance2BR= ParametersArray(FrameNumber,2);
   Distance2TL= ParametersArray(FrameNumber,3);
   if Distance2BR>Threshold_Plate2Head(2)+2 && Distance2BR>Threshold_Plate2Head(2)+2
SafeDistance=true;
   else
SafeDistance=false;
   end
if ParametersArray(FrameNumber,1)<Threshold_Body2Head(2)%'Rearing'
    if SafeDistance
MotifNames(FrameNumber)=MotifsName(3);
MotifCodes(FrameNumber)=3; 
    else
MotifNames(FrameNumber)={'Unclassified'};
MotifCodes(FrameNumber)=nan; 
    end
elseif ParametersArray(FrameNumber,2)<Threshold_Plate2Head(2)%Close2BR
    % mouse is close to BR plate. 
    if Variables.FoodBR %plate contains food
MotifNames(FrameNumber)=MotifsName(1);
MotifCodes(FrameNumber)=1;
    else
MotifNames(FrameNumber)=MotifsName(2);
MotifCodes(FrameNumber)=2;  
    end   
elseif ParametersArray(FrameNumber,3)<Threshold_Plate2Head(2)%Close2TL
        % mouse is close to TL plate. 
    if ~Variables.FoodBR %plate contains food
MotifNames(FrameNumber)=MotifsName(1);
MotifCodes(FrameNumber)=1;
    else %plate is empty 
MotifNames(FrameNumber)=MotifsName(2);
MotifCodes(FrameNumber)=2;   
    end  
    
elseif ParametersArray(FrameNumber,7)<0 %TurnLeft
MotifNames(FrameNumber)=MotifsName(4);
MotifCodes(FrameNumber)=4;      
elseif ParametersArray(FrameNumber,7)>0 %TurnRight
MotifNames(FrameNumber)=MotifsName(5);
MotifCodes(FrameNumber)=5; 
elseif ParametersArray(FrameNumber,5)<Threshold_Velocity.Stopping(2) %Stopping
if SafeDistance
MotifNames(FrameNumber)=MotifsName(6);
MotifCodes(FrameNumber)=6;      
    else
MotifNames(FrameNumber)={'Unclassified'};
MotifCodes(FrameNumber)=nan; 
    end
elseif Threshold_Velocity.Walking(1)<ParametersArray(FrameNumber,5)...
        && ParametersArray(FrameNumber,5)<Threshold_Velocity.Walking(2) %Walking
    if SafeDistance
MotifNames(FrameNumber)=MotifsName(7);
MotifCodes(FrameNumber)=7;
    else
MotifNames(FrameNumber)={'Unclassified'};
MotifCodes(FrameNumber)=nan; 
    end
elseif Threshold_Velocity.Trotting(1)<ParametersArray(FrameNumber,5)...
        && ParametersArray(FrameNumber,5)<Threshold_Velocity.Trotting(2)%Trotting
if SafeDistance
MotifNames(FrameNumber)=MotifsName(8);
MotifCodes(FrameNumber)=8;  
    else
MotifNames(FrameNumber)={'Unclassified'};
MotifCodes(FrameNumber)=nan; 
    end
elseif Threshold_Velocity.Running(1)<ParametersArray(FrameNumber,5)...
        && ParametersArray(FrameNumber,5)<Threshold_Velocity.Running(2) %Running
    if SafeDistance
MotifNames(FrameNumber)=MotifsName(9);
MotifCodes(FrameNumber)=9; 
    else
MotifNames(FrameNumber)={'Unclassified'};
MotifCodes(FrameNumber)=nan; 
    end
end
end
MotifCodes=MotifCodes';
% % % % plot to test code
% figure
% Frame=find(MotifCodes==1); 
% randomFrame = randsample(Frame, 1);
% % Load the video file
% videoFile = [Variables.VideoPath,'\',Condition.Movie_AVI]; % Replace with your video file name
% vidObj = VideoReader(videoFile);
% % Specify the frame number you want to display
% frameNumber = randomFrame; % Change this to the frame number you want
% % Calculate the time corresponding to the desired frame
% frameTime = (frameNumber - 1) / VideoInfo.FrameRate;
% % Set the current time of the video to the frame time
% VideoInfo.CurrentTime = frameTime;
% % Read the specific frame
% frame = readFrame(VideoInfo);
% % Display the frame
% imshow(frame);
% title('Should be eating');
% %%
%% Initialize a filtered_data array, starting with the original data
filtered_MotifCodes = MotifCodes;
% Loop over the data, skipping the first and last elements to avoid boundary issues
for i = 2:length(MotifCodes)-1
    % Check if the current value is an error (different from both neighbors)
    if MotifCodes(i) ~= MotifCodes(i-1) && MotifCodes(i) ~= MotifCodes(i+1)
        % If it's an error, replace it with the previous value (or next, as they are the same)
        filtered_MotifCodes(i) = MotifCodes(i-1);
    end
end
% Handle special case where two consecutive values are different from neighbors
for i = 2:length(MotifCodes)-2
    if MotifCodes(i) ~= MotifCodes(i-1) && MotifCodes(i) == MotifCodes(i+2) && MotifCodes(i+1) ~= MotifCodes(i)
        filtered_MotifCodes(i+1) = MotifCodes(i);
    end
end
end
% % Plot the data
% figure;
% % Frame rate and corresponding time axis
% frame_rate = 15; % Frames per second
% num_frames = length(filtered_MotifCodes); % Number of frames in the data
% time = (0:num_frames-1) / frame_rate; % Time axis in seconds
% % Create a color map for each unique code in the data
% unique_values = unique(filtered_MotifCodes); % Unique values in the data
% colors = lines(length(unique_values)); % Generate distinct colors
% % Plot the data with different colors for each code
% hold on;
% for i = 1:length(unique_values)
%     % Find indices of the current unique value
%     indices = (filtered_MotifCodes == unique_values(i));
%     
%     % Plot those indices with the corresponding color
%     plot(time(indices), filtered_MotifCodes(indices), '|', 'Color', colors(i,:), 'MarkerSize', 12);
% end
% % Replace the y-ticks with the MotifsName
% yticks(unique_values); % Set y-ticks to be the unique data values
% yticklabels(MotifsName(unique_values)); % Replace y-tick labels with corresponding behavior labels
% % Customize the plot
% xlabel('Time (seconds)');
% ylabel('Code');
% title('Data Over Time with Different Codes');
% legend(arrayfun(@num2str, unique_values, 'UniformOutput', false));
% grid on;
% hold off;
MotifStart=[filtered_MotifCodes(1);filtered_MotifCodes(1:end-1)]-filtered_MotifCodes;
MotifEnd=filtered_MotifCodes-[filtered_MotifCodes(2:end);filtered_MotifCodes(end)];
MotifTTL=nan(size(MotifStart));MotifTTL(1)=filtered_MotifCodes(1);MotifTTL(end)=filtered_MotifCodes(end)+10;
MotifTTL(MotifStart ~= 0)=filtered_MotifCodes((MotifStart ~= 0));% show start of event by number
MotifTTL(MotifEnd ~= 0)=filtered_MotifCodes((MotifEnd ~= 0))+10;% show end of event by number+10
Condition.MotifTTL=MotifTTL;
Condition.filtered_MotifCodes=filtered_MotifCodes;
Condition.MotifsName=MotifsName';
Condition.TimestampMotifs=linspace(0,Condition.VideoInfo.Duration,length(MotifTTL))';
 end