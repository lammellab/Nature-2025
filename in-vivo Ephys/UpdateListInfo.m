function [Condition,Units_in_Condition,Units_in_Condition_Titles, FinalIndex]=UpdateListInfo(Condition,Units_in_Condition, Unit2Screen,Index,Units_in_Condition_Titles)
% (Condition,Unit2Screen, Units_in_Condition,Units_in_Condition_Titles,Index);
if Condition.ConditionOrder==3 
%% Info about the unit%% Info about the optotagging
% Decision
try Units_in_Condition(Unit2Screen, Index) = {Condition.Tagged}; catch end
Units_in_Condition_Titles(Unit2Screen, Index) = {[char(Condition.ConditionName),'.Tagged']};Index=Index+1;
try Units_in_Condition(Unit2Screen, Index) = {Condition.SignificanceMarker}; catch end
Units_in_Condition_Titles(Unit2Screen, Index) = {[char(Condition.ConditionName),'.SignificanceMarker']};Index=Index+1;
% Stats
% Information about the tagging
try Units_in_Condition(Unit2Screen, Index) = {Condition.Latency_ms}; catch end
Units_in_Condition_Titles(Unit2Screen, Index) = {[char(Condition.ConditionName),'.Latency_ms']};Index=Index+1;
try Units_in_Condition(Unit2Screen, Index) = {Condition.MaxHeight_uV}; catch end
Units_in_Condition_Titles(Unit2Screen, Index) = {[char(Condition.ConditionName),'.MaxHeight_uV']};Index=Index+1;
try Units_in_Condition(Unit2Screen, Index) = {Condition.AverageMaxValueHz}; catch end
Units_in_Condition_Titles(Unit2Screen, Index) = {[char(Condition.ConditionName),'.AverageMaxValueHz']};Index=Index+1;
try Units_in_Condition(Unit2Screen, Index) = {Condition.MedianValueHz}; catch end
Units_in_Condition_Titles(Unit2Screen, Index) = {[char(Condition.ConditionName),'.AverageMedianValueHz']};Index=Index+1;
try Units_in_Condition(Unit2Screen, Index) = {Condition.BaselineValueHz}; catch end
Units_in_Condition_Titles(Unit2Screen, Index) = {[char(Condition.ConditionName),'.BaselineValueHz']};Index=Index+1;
%% Plots info
try Units_in_Condition(Unit2Screen, Index) = {Condition.avg_spike_rate}; catch end
Units_in_Condition_Titles(Unit2Screen, Index) = {[char(Condition.ConditionName),'.avg_spike_rate']};Index=Index+1;
try Units_in_Condition(Unit2Screen, Index) = {Condition.avg_spike_rate_z}; catch end
Units_in_Condition_Titles(Unit2Screen, Index) = {[char(Condition.ConditionName),'.avg_spike_rate_z']};Index=Index+1;
try Units_in_Condition(Unit2Screen, Index) = {Condition.std_spike_rate}; catch end
Units_in_Condition_Titles(Unit2Screen, Index) = {[char(Condition.ConditionName),'.std_spike_rate']};Index=Index+1;
try Units_in_Condition(Unit2Screen, Index) = {round(Condition.time_vector)'}; catch end
Units_in_Condition_Titles(Unit2Screen, Index) = {[char(Condition.ConditionName),'.time_vector']};Index=Index+1;
try Units_in_Condition(Unit2Screen, Index) = {[Condition.TimeXmin, Condition.TimeXmax]}; catch end
Units_in_Condition_Titles(Unit2Screen, Index) = {[char(Condition.ConditionName),'.Times2Diplay']};Index=Index+1;
try Units_in_Condition(Unit2Screen, Index) = {Condition.Heatmap}; catch end
Units_in_Condition_Titles(Unit2Screen, Index) = {[char(Condition.ConditionName),'.Heatmap']};Index=Index+1;
end
if Condition.ConditionOrder<3 
try Units_in_Condition(Unit2Screen, Index) = {Condition.MaxChannel}; catch end
Units_in_Condition_Titles(Unit2Screen, Index) = {[char(Condition.ConditionName),'.MaxChannel']};Index=Index+1;     

% % Information about the behavioral responses -adapt 
try Units_in_Condition(Unit2Screen, Index) = {Condition.mean_z_scored_spikes}; catch end
Units_in_Condition_Titles(Unit2Screen, Index) = {[char(Condition.ConditionName),'.mean_z_scored_spikes']};Index=Index+1;     
try Units_in_Condition(Unit2Screen, Index) = {Condition.sem_z_scored_spikes}; catch end
Units_in_Condition_Titles(Unit2Screen, Index) = {[char(Condition.ConditionName),'.sem_z_scored_spikes']};Index=Index+1;     
try Units_in_Condition(Unit2Screen, Index) = {Condition.Type}; catch end
Units_in_Condition_Titles(Unit2Screen, Index) = {[char(Condition.ConditionName),'.Type']};Index=Index+1;     
try Units_in_Condition(Unit2Screen, Index) = {Condition.Significant}; catch end
Units_in_Condition_Titles(Unit2Screen, Index) = {[char(Condition.ConditionName),'.Significant']};Index=Index+1;     
try Units_in_Condition(Unit2Screen, Index) = {Condition.SignificanceMarker}; catch end
Units_in_Condition_Titles(Unit2Screen, Index) = {[char(Condition.ConditionName),'.SignificanceMarker']};Index=Index+1;    
try Units_in_Condition(Unit2Screen, Index) = {Condition.Pvalue}; catch end
Units_in_Condition_Titles(Unit2Screen, Index) = {[char(Condition.ConditionName),'.Pvalue']};Index=Index+1;    
try Units_in_Condition(Unit2Screen, Index) = {Condition.statsW}; catch end
Units_in_Condition_Titles(Unit2Screen, Index) = {[char(Condition.ConditionName),'.statsW']};Index=Index+1; try Units_in_Condition(Unit2Screen, Index) = {Condition.MotifsName}; catch end
Units_in_Condition_Titles(Unit2Screen, Index) = {char(Condition.MotifsName)};Index=Index+1;
try Units_in_Condition(Unit2Screen, Index) = {Condition.TimestampMotifs }; catch end
Units_in_Condition_Titles(Unit2Screen, Index) = {[char(Condition.ConditionName),'.TimestampMotifsSN ']};Index=Index+1;
try Units_in_Condition(Unit2Screen, Index) = { Condition.MotifTTL}; catch end
Units_in_Condition_Titles(Unit2Screen, Index) = {[char(Condition.ConditionName),'.MotifTTL ']};Index=Index+1;
%% Plots info
%event related plots
try Units_in_Condition(Unit2Screen, Index) = {Condition.avg_spike_rate}; catch end
Units_in_Condition_Titles(Unit2Screen, Index) = {[char(Condition.ConditionName),'.avg_spike_rate']};Index=Index+1;
try Units_in_Condition(Unit2Screen, Index) = {Condition.time_bins}; catch end
Units_in_Condition_Titles(Unit2Screen, Index) = {[char(Condition.ConditionName),'.time_bins']};Index=Index+1;
try Units_in_Condition(Unit2Screen, Index) = {Condition.HeatmapZ}; catch end
Units_in_Condition_Titles(Unit2Screen, Index) = {[char(Condition.ConditionName),'.Heatmap']};Index=Index+1;
% time plot
try Units_in_Condition(Unit2Screen, Index) = {Condition.firing_rateHz }; catch end
Units_in_Condition_Titles(Unit2Screen, Index) = {[char(Condition.ConditionName),'.firing_rateHz']};Index=Index+1;
try Units_in_Condition(Unit2Screen, Index) = {Condition.z_scored_firing_rate }; catch end
Units_in_Condition_Titles(Unit2Screen, Index) = {[char(Condition.ConditionName),'.z_scored_firing_rate']};Index=Index+1;
try Units_in_Condition(Unit2Screen, Index) = {Condition.time_bins_minutes }; catch end
Units_in_Condition_Titles(Unit2Screen, Index) = {[char(Condition.ConditionName),'.time_bins_minutes']};Index=Index+1;
end
FinalIndex=Index;
end






