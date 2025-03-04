function [Variables,Units_in_Condition] = CreateListOfFiles(Variables)
% Define file path
filePath = [Variables.ComputerDir,'\ListOfFilesToAnalyse.xlsx'];
% Check if the file exists
if isfile(filePath)
    % Read the Excel file
[~, ~, DataFromExcel] = xlsread(filePath, 'All Units', '', 'basic');
else
    disp('File not found');
    return;  % Exit if the file is not found
end
% Remove first row (titles)
DataFromExcel = DataFromExcel(3:end, :);
% Initialize an empty cell array for Units_in_Condition
Units_in_Condition = {};
% Loop through the data
for Unit2Screen = 1:size(DataFromExcel, 1)
    % Get relevant information from the current row
    Date = DataFromExcel{Unit2Screen, 2};           % Date (Column 2)
    MouseName = DataFromExcel{Unit2Screen, 3};      % Mouse_Name (Column 3)
    JellyFile = DataFromExcel{Unit2Screen, 4};      % Jelly file (Column 4)
    ChowFile = DataFromExcel{Unit2Screen, 5};       % Chow file (Column 5)
    LaserFile = DataFromExcel{Unit2Screen, 6};      % Laser file (Column 6)
    
    % Ensure that JellyFile, ChowFile, and LaserFile are numeric and not NaN
    if isnumeric(JellyFile) && ~isnan(JellyFile) && ...
       isnumeric(ChowFile) && ~isnan(ChowFile) && ...
       isnumeric(LaserFile) && ~isnan(LaserFile)
    
        % Get the values from columns 7 onward to handle any number of tetrodes
        numCols = size(DataFromExcel, 2) - 6; % Number of tetrode columns
        for j = 1:numCols
            % Check if the value is numeric or string in columns beyond 6
            UnitValue = DataFromExcel{Unit2Screen, j+6};
            
            if ischar(UnitValue)  % If it's a string (e.g., '1,2,3')
                UnitNumbers = str2num(UnitValue);  % Convert string to numbers
            elseif isnumeric(UnitValue) && ~isnan(UnitValue)  % If it's already numeric
                UnitNumbers = UnitValue;  % Treat it as a single unit number
            else
                continue;  % If it's neither, skip to the next iteration
            end
            
            % Loop over each unit number
            for unit = UnitNumbers
                % Skip if the unit is NaN
                if isnan(unit)
                    continue;
                end
                
                % Generate filenames for Jelly, Chow, and Laser
                formattedString_Jelly = sprintf('%s_%s_Tetrode_%d_Unit_%d_File_%02d', Date, MouseName, j, unit, JellyFile);
                formattedString_Chow = sprintf('%s_%s_Tetrode_%d_Unit_%d_File_%02d', Date, MouseName, j, unit, ChowFile);
                formattedString_Laser = sprintf('%s_%s_Tetrode_%d_Unit_%d_File_%02d', Date, MouseName, j, unit, LaserFile);
                
                % Store these formatted strings in the same row (three columns)
                Units_in_Condition{end+1, 1} = formattedString_Jelly; % Column 1 for Jelly
                Units_in_Condition{end, 2} = formattedString_Chow;    % Column 2 for Chow
                Units_in_Condition{end, 3} = formattedString_Laser;   % Column 3 for Laser
            end
        end
    end
end

end

    

