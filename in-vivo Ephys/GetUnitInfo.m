function [Index,Units_in_Condition_Titles,Variables,  Unit2Screen, Units_in_Condition]=GetUnitInfo( Variables, Unit2Screen, Units_in_Condition,Units_in_Condition_Titles) 
%% extract available variables
%  Unit Name
Variables.UnitName = Units_in_Condition{Unit2Screen, 1};
Variables.UnitName=Variables.UnitName(1:end-8);
disp(Variables.UnitName)
% Use regexp to extract the parts of the string
tokens = regexp(Variables.UnitName, '(?<Date>\d{4}-\d{2}-\d{2})_(?<MouseName>[\w-]+)_Tetrode_(?<TetrodeNumber>\d+)_Unit_(?<UnitNumber>\d+)', 'names');
% Check if tokens are returned
if ~isempty(tokens)
    % Store extracted values into Variables structure
    Variables.Date = tokens.Date;
    Variables.TetrodeNumber = str2double(tokens.TetrodeNumber);
    Variables.UnitNumber = str2double(tokens.UnitNumber);
    Variables.MouseName = tokens.MouseName;
else
    disp('No tokens extracted. Please check the format of UnitName.');
end
% Get the diet type
REGlist={'13-3', '13-9','11-7','18-9','19-9','18-5','21-5','20-7'};
for u=1:length(REGlist)
    if Variables.MouseName(end-3:end)== REGlist{u}
        Variables.DietType='REG';
        break
    else
       Variables.DietType='HFD'; 
    end
end
Variables.UnitGeneralPath=[Variables.ComputerDir,'\',Variables.Date,'\',Variables.MouseName,'\'];
Variables.ChannelNumbers=[4*(Variables.TetrodeNumber-1)+1,...
4*(Variables.TetrodeNumber-1)+2,4*(Variables.TetrodeNumber-1)+3,4*(Variables.TetrodeNumber-1)+4];
%% Get the information from all conditions for a specific unit
Variables.FileNumbers=nan(1,3);
for ConditionNumber=1:3
TempString=Units_in_Condition{Unit2Screen, ConditionNumber};
Variables.FileNumbers(ConditionNumber)=str2num(TempString(end-1:end));
end
Variables.Type = {'Jelly','Chow','Laser'};  % Trim any whitespace
Variables.TTLType={'Piezo','Piezo','Laser'};
Variables.EventID={10,10,11};
Variables.TTLs={4,4,1};
% determine the location of the food plate
Variables.FoodBR=datetime(Variables.Date)<datetime('2021-01-01');
% Path names
Variables.FolderPath= [Variables.ComputerDir,'\',Variables.Date...
    ,'\',Variables.MouseName];
% Get the path to video files for the Jelly and Chow conditions
Variables.VideoPath= [Variables.ComputerDir,'\',Variables.Date...
    ,'\',Variables.MouseName,'\Video'];
% define path to data
% Get the exact file names
RealFileNames=dir(Variables.FolderPath);
RealFileNames=extractfield(RealFileNames,'name')';
Variables.RealFileNames=RealFileNames(Variables.FileNumbers+2);
Variables.FoodConsumed=nan(1,3);
for ConditionNumber=1:3
% Use regexp to find the food consumed pattern (e.g., '0.541-gr')
tokensReal = regexp(Variables.RealFileNames{ConditionNumber, 1}, '_(?<FoodConsumed>\d+(\.\d+)?)-gr', 'names');
% Check if tokens are returned
if ~isempty(tokensReal)
    % Store the extracted value into Variables.FoodConsumed
    Variables.FoodConsumed(ConditionNumber) = str2num(tokensReal.FoodConsumed);  
else
    % If no value found, store NaN
    Variables.FoodConsumed(ConditionNumber) = NaN;
end
end
% determine the referance to condition
ConditionList={'Jelly1','Jelly 1', 'Chow', '5ms_2hz','Laser_Single-5ms_2hz','Jelly - Exposure','Jelly','jelly','chow','Laser'};
for ConditionNumber=1:3
for u=1:length(ConditionList)
if contains(Variables.RealFileNames{ConditionNumber, 1},ConditionList{u})
  Variables.ConditionName(ConditionNumber)={ConditionList{u}};
  break
end
end
end
Index=1;
Units_in_Condition_Titles(Unit2Screen,Index)={'Jelly_File'};Index=Index+1;
Units_in_Condition_Titles(Unit2Screen,Index)={'Chow_File'};Index=Index+1;
Units_in_Condition_Titles(Unit2Screen,Index)={'Laser_File'};Index=Index+1;
try Units_in_Condition(Unit2Screen, Index) = {Variables.TetrodeNumber}; catch end
Units_in_Condition_Titles(Unit2Screen, Index) = {'Variables.TetrodeNumber'};Index=Index+1;
try Units_in_Condition(Unit2Screen, Index) = {Variables.UnitNumber  }; catch end
Units_in_Condition_Titles(Unit2Screen, Index) = {'Variables.UnitNumber'};Index=Index+1;
try Units_in_Condition(Unit2Screen, Index) = {Variables.MouseName}; catch end
Units_in_Condition_Titles(Unit2Screen, Index) = {'Variables.MouseName'};Index=Index+1;
try Units_in_Condition(Unit2Screen, Index) = {Variables.DietType}; catch end
Units_in_Condition_Titles(Unit2Screen, Index) = {'Variables.DietType'};Index=Index+1;
try Units_in_Condition(Unit2Screen, Index) = {Variables.Date}; catch end
Units_in_Condition_Titles(Unit2Screen, Index) = {'Variables.Date'};Index=Index+1;
end


