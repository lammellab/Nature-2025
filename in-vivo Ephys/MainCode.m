clear; close all;clc;
%% Neta Gazit Shimoni Updated 11/02/24
for ChooseVariables=1:1
%% this code uses the following functions: 
% CreateListOfFiles.m % GetUnitInfo.m % CollectData.m % GetMotif.m % VideoSync.m
% EventsLaser.m % ConnectDataFiles.m % GetTimes.m % StatisticalTests.m % UpdateListInfo.m
%% this code uses the following libraries: 
% NeuralynxMatlabImportExport_v6.0.0
%% VARIABLES RELATED TO ANALYSIS
Variables.MinimalboutIntervalSec=6; % interval between bouts
Variables.BoutLengthSecondsLimit=1; % minimal bout length
Variables.TimeLimit=15; % time in minutes-limits the timeframe of the analysis from start to the specified time
Variables.DisplayPlot=false;
Variables.RunLaser=true;
Variables.RunFood=true;
Variables.DisplaySummaryPlot=false;
Variables.WindowForBehavior=5;%sec
Variables.MinimalBoatInterval=6;%sec
Variables.Factor=0; 
Variables.ComputerDir='E:';
% add NLX libraries
addpath('C:\Users\netas\Documents\MATLAB\Obesity\EphysNature\NeuralynxMatlabImportExport_v6.0.0')
end
% Create list of units to analyse :
[Variables,Units_in_Condition] = CreateListOfFiles(Variables);  % Correct way to call the method
count=1; Units_in_Condition_Titles={};
combinedStatisticsTable = struct([]);
% combinedCorrelation = struct([]);
UnitIndex=1;
for UnitNumber=[1]%1:length(Units_in_Condition)
    close all
clearvars -except UnitIndex UnitNumber combinedStatisticsTable Units_in_Condition Units_in_Condition_Titles Variables
     try
        Variables.UnitIndex=UnitIndex;
        Variables.Unit=UnitNumber;
%% collect basic data   
[Index,Units_in_Condition_Titles,Variables,  UnitNumber, Units_in_Condition]...
    =GetUnitInfo( Variables, UnitNumber, Units_in_Condition,Units_in_Condition_Titles) ;
%% Load data
% % create Jelly object
[Jelly,Index,Units_in_Condition_Titles,Variables, UnitNumber,...
    Units_in_Condition] = CollectData(Variables,Index,...
    Units_in_Condition_Titles, UnitNumber, Units_in_Condition,1) ;
[Jelly] = VideoSync(Jelly, Variables);
% % create Chow object
[Chow,Index,Units_in_Condition_Titles,Variables,  UnitNumber,...
    Units_in_Condition] = CollectData(Variables,Index,...
    Units_in_Condition_Titles, UnitNumber, Units_in_Condition,2) ;
[Chow] = VideoSync(Chow, Variables);
% % create Laser object
[Laser,Index,Units_in_Condition_Titles,Variables,  UnitNumber,...
    Units_in_Condition] = CollectData(Variables,Index,...
    Units_in_Condition_Titles, UnitNumber, Units_in_Condition,3) ;
try
Laser.Xmin=-5;Laser.Xmax=20; % range to plot around the Laser pulse in ms
Laser.TimeXmin=50.2;Laser.TimeXmax=50.4; % range to plot around the Laser pulse in ms
[Laser] = EventsLaser(Laser,Variables); 
[Laser,Units_in_Condition,Units_in_Condition_Titles,Index]=UpdateListInfo...
    (Laser,Units_in_Condition, UnitNumber,Index,Units_in_Condition_Titles);
catch;end
%% Run the analysis and get the relevant information of each unit and event type
%  [EventTimes]=DefineEvents(Jelly,Variables,EventType)
 [Unit]=ConnectDataFiles(Jelly,Chow,Laser,Variables);    
[StatisticsTable]= GetTimes(Unit,Variables);    % Append currentStruct to combinedStruct
for t=1:1   
    try
StatisticsTable(13).Unit=StatisticsTable(12).Unit;
StatisticsTable(13).UnitIndex=StatisticsTable(12).UnitIndex;
StatisticsTable(13).DietType=StatisticsTable(12).DietType;
StatisticsTable(13).MouseName=StatisticsTable(12).MouseName;
StatisticsTable(13).Date=StatisticsTable(12).Date;
StatisticsTable(13).TetrodeNumber=StatisticsTable(12).TetrodeNumber;
StatisticsTable(13).UnitNumber=StatisticsTable(12).UnitNumber;
StatisticsTable(13).Condition={'Laser'};
StatisticsTable(13).Factor=StatisticsTable(12).Factor;
StatisticsTable(13).Bouts=nan;
StatisticsTable(13).TotalFiringRateBoutHz=nan;
StatisticsTable(13).TotalBoutTime=nan;
StatisticsTable(13).TotalNonBoutTime=nan;
StatisticsTable(13).TotalFiringRateBaselineHz=nan;
StatisticsTable(13).TotalFiringRateBoutZ=nan;
StatisticsTable(13).TotalFiringRateBaselineZ=nan;
StatisticsTable(13).p_value=nan;
StatisticsTable(13).Decision=Units_in_Condition{UnitNumber, 40} ;
StatisticsTable(13).pre_bout_rates=nan;
StatisticsTable(13).in_bout_rates=nan;
StatisticsTable(13).pre_bout_ratesMedian=nan;
StatisticsTable(13).in_bout_ratesMedian=nan;
StatisticsTable(13).pre_bout_ratesMean=Units_in_Condition{UnitNumber, 46}  ;
StatisticsTable(13).in_bout_ratesMean=Units_in_Condition{UnitNumber, 44}  ;
StatisticsTable(13).Type=nan;
StatisticsTable(13).unit_mean=nan;
StatisticsTable(13).unit_sem=nan;
StatisticsTable(13).Latency_ms=Units_in_Condition{UnitNumber, 42} ;
StatisticsTable(13).SignificanceMarker=Units_in_Condition{UnitNumber, 41} ;
if Units_in_Condition{UnitNumber, 40}
    for m=1:13
  StatisticsTable(m).Tagged=true;
    end
else
        for m=1:13
  StatisticsTable(m).Tagged=false;
  ends
end
end
end
combinedStatisticsTable = [combinedStatisticsTable, StatisticsTable];
% combinedCorrelation = [combinedCorrelation, Correlation];
UnitIndex=UnitIndex+1;
% close all

end
        catch
        continue
        disp('Unit did not run')
end
save('Units_in_Condition', 'Units_in_Condition', '-v7.3');
save('combinedStatisticsTable', 'combinedStatisticsTable', '-v7.3');
end