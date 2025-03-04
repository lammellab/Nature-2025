 function [Condition] = EventsLaser(Condition,Variables) 
% Load the unit and event data
Variables.DisplayPlot=true;
EventTimeStamps=double(Condition.TimeStamps_LaserNorm)'/1000000; %r
UnitTimeStamps_cellsSeconds=double(Condition.TimeStamps_unit_Norm)'/1000000;% c
% figure
% scatter(EventTimeStamps,ones(size(EventTimeStamps)));
% hold on
% scatter(UnitTimeStamps_cellsSeconds',2*ones(size(UnitTimeStamps_cellsSeconds')));
% xlim([150 155]);
% Load the unit and event data% Load the unit and event data
if Variables.DisplayPlot; fig = figure('units', 'normalized', 'outerposition', [0 0 1 1]); end % Create a full-screen figure with 6 subplots
%% First subplot: Raster plot
for RasterPlot=1:1
    % Define bin size in seconds (you can change this value)
bin_size = 0.001; % Time bin size in seconds (e.g., 1 ms bins)
% Initialize arrays for storing spike rates, raster plot data, and heatmap data
all_spike_rates = []; raster_spike_times = {}; raster_event_indices = []; heatmap_data = [];
% Iterate through each laser event, except the last one
for EventNumber = 1:length(EventTimeStamps)
    try
    % Calculate the time to the next event
    event_Starttime = EventTimeStamps(EventNumber);
    next_event_time = EventTimeStamps(EventNumber+1);
    event_interval = next_event_time - event_Starttime;    
    % Define pre-event and post-event windows
    pre_event_window = 0.3 * event_interval;
    post_event_window = event_interval - pre_event_window;    
    % Number of bins in the time window (with smaller bins)
    num_bins = ceil((pre_event_window + post_event_window) / bin_size);  
    % Initialize spike counts for this event
    spike_counts = zeros(num_bins, 1);   
    % Select unit spikes that fall within this window
    time_window_start = event_Starttime - pre_event_window;
    time_window_end = event_Starttime + post_event_window;
    spikes_in_window = UnitTimeStamps_cellsSeconds(UnitTimeStamps_cellsSeconds>= ...
    time_window_start & UnitTimeStamps_cellsSeconds <= time_window_end) - event_Starttime;
        % Bin the spikes and update the spike count
    for Simulated = 1:num_bins
        bin_start = (Simulated - 1) * bin_size - pre_event_window;
        bin_end = Simulated * bin_size - pre_event_window;
        spike_counts(Simulated) = sum(spikes_in_window >= bin_start...
            & spikes_in_window < bin_end);
    end   
    % Normalize spike counts to Hz
    spike_rate_per_bin = spike_counts / bin_size; % This keeps the rate in Hz
    if size(all_spike_rates,1)==size(spike_rate_per_bin,1)
    all_spike_rates = [all_spike_rates, spike_rate_per_bin]; %#ok<AGROW>
    elseif EventNumber>1
    % Calculate the difference in length
difference = size(all_spike_rates,1) - size(spike_rate_per_bin,1);
    % Append zeros to the shorter vector (vector2)
spike_rate_per_bin = [spike_rate_per_bin, zeros(1, difference)];
    all_spike_rates = [all_spike_rates, spike_rate_per_bin]; %#ok<AGROW>
    end
   
    % Store the spike times for the raster plot
    raster_spike_times{EventNumber} = spikes_in_window * 1000; % Convert to ms
    raster_event_indices = [raster_event_indices; repmat(EventNumber, length(spikes_in_window), 1)]; %#ok<AGROW>
    
    % Store spike rate for the heatmap
    heatmap_data = [heatmap_data; spike_rate_per_bin']; %#ok<AGROW>
    catch
        continue
end
end
% Invert the order of events for the heatmap (event 1 at bottom)
heatmap_data = flipud(heatmap_data);
% Compute the average and standard deviation of spike counts across events
Condition.all_spike_rates=all_spike_rates;
avg_spike_rate = mean(all_spike_rates, 2);
std_spike_rate = std(all_spike_rates, 0, 2);
Condition.avg_spike_rate=avg_spike_rate;Condition.std_spike_rate=std_spike_rate;
% Generate the time vector for plotting (converted to milliseconds)
time_vector = linspace(-pre_event_window, post_event_window, num_bins) * 1000; % Convert to ms
Condition.time_vector=time_vector;
[Condition.AverageMaxValueHz,Latency_ms]=max(avg_spike_rate);
Condition.Latency_ms=time_vector(Latency_ms);
BaselineTimes=find(-10< time_vector & time_vector <0);
LaserRelevantTimes=find(0< time_vector & time_vector <10);
for EventNumber=1:size(all_spike_rates,2)
SpikesInLatency(EventNumber)=sum(all_spike_rates(BaselineTimes,EventNumber)/1000);%number of spikes
SpikesBaseline(EventNumber)=sum(all_spike_rates(LaserRelevantTimes,EventNumber)/1000);

end
Condition.SpikesInLatency=SpikesInLatency';
Condition.SpikesBaseline=SpikesBaseline';
% test data for normal distribution
[NormalSpikesInLatency, ~] = lillietest(Condition.SpikesInLatency);
[NormalSpikesBaseline, ~] = lillietest(Condition.SpikesBaseline);
if NormalSpikesInLatency==1 || NormalSpikesBaseline==1
    %if not normally distributed use non-parametric test
    Condition.NormalDistribution=false;
end
% Perform the Mann-Whitney U test (Wilcoxon rank-sum test) to determine if spikes after laser are diffrent then baseline
[Condition.WilcoxonRankSum_Pvalue, Significant_WilcoxonRankSum, WilcoxonRankSum_Stats] = ...
    ranksum(Condition.SpikesInLatency, Condition.SpikesBaseline);
Condition.WilcoxonRankSum_zval=WilcoxonRankSum_Stats.zval;
Condition.WilcoxonRankSum_ranksum=WilcoxonRankSum_Stats.ranksum;
if Significant_WilcoxonRankSum==1 
    if sum(avg_spike_rate(LaserRelevantTimes))>sum(avg_spike_rate(BaselineTimes))
Condition.WilcoxonTagged=true;
    end
else
Condition.WilcoxonTagged=false;
end 
Condition.MedianValueHz=median(avg_spike_rate);
% get the baseline firing rate 10 sec before the laser pulse
Condition.BaselineValueHz=mean(Condition.avg_spike_rate(-10 <Condition.time_vector & Condition.time_vector<0));
% Calculate the average intensity per bin for the new line plot
avg_intensity_per_bin = mean(heatmap_data, 1);
    if Variables.DisplayPlot
subplot(3,2,1);
hold on;
for EventNumber = 1:length(raster_spike_times)
    plot(raster_spike_times{EventNumber}, raster_event_indices(raster_event_indices == EventNumber), 'k.', 'MarkerSize', 8);
end
xlim([Condition.Xmin Condition.Xmax]); % Show -5 to 20 ms
xlabel('Time (ms)');
ylabel('Event (Laser Pulses)');
title('Raster Plot of Neural Activity Around Laser Events');
% Add vertical red line at time = 0 for the laser event
plot([0 0], [0 length(EventTimeStamps)], 'r--', 'LineWidth', 2);
hold off;
    end
end
%% Second subplot: Z-score of Firing Rate histogram
for FiringRateHistogram=1:1
if Variables.DisplayPlot
subplot(3,2,3);
hold on;
% Calculate z-scores for average spike rate
avg_spike_rate_z = (avg_spike_rate - mean(avg_spike_rate)) / std(avg_spike_rate);
std_spike_rate_z = std_spike_rate / std(avg_spike_rate); % Normalize the standard deviation as well
Condition.avg_spike_rate_z=avg_spike_rate_z;
Condition.std_spike_rate_z=std_spike_rate_z;
% Plot shaded area for standard deviation in z-score
fill([time_vector, fliplr(time_vector)], ...
     [avg_spike_rate_z + std_spike_rate_z; flipud(avg_spike_rate_z - std_spike_rate_z)], ...
     [0.9 0.9 0.9], 'EdgeColor', 'none'); % Shaded area
% Plot average spike rate in z-score
plot(time_vector, avg_spike_rate_z, 'b', 'LineWidth', 2);
% Add vertical red line at time = 0 for the laser event
plot([0 0], ylim, 'r--', 'LineWidth', 2);
% Customize the plot
xlim([Condition.Xmin Condition.Xmax]); % Set x-axis limits to show -5 to 20 ms
xlabel('Time (ms)');
ylabel('Firing Rate (Z-score)');
title('Average Firing Rate (Z-score) Around Laser Events');
legend('Standard Deviation', 'Average Firing Rate', 'Laser Event');
hold off;
    end
end
%% Third subplot: Average Heatmap
for HeatMapPlot=1:1
% limit avg_spike_rate data to the selected frame [-5 20]
avg_spike_rateInRange=avg_spike_rate(find(round(time_vector)==Condition.Xmin):find(round(time_vector)==Condition.Xmax));
avg_spike_rateInRangeZscored=(avg_spike_rateInRange - mean(avg_spike_rateInRange)) ./ std(avg_spike_rateInRange);
Condition.TimeRageHistogramX=round(time_vector(find(round(time_vector)==Condition.Xmin):find(round(time_vector)==Condition.Xmax)));
Condition.Heatmap=avg_spike_rateInRangeZscored;
if Variables.DisplayPlot
subplot(3,2,5);
Sub4=heatmap(avg_spike_rateInRangeZscored');
Sub4.XDisplayLabels = string(Condition.TimeRageHistogramX);
colormap('jet');
Sub4.ColorLimits = [0 7];
xlabel('Time (ms)');
title('Intensity Averaged Heatmap');
% Add vertical red line at time = 0 for the laser event
% plot([0 0], ylim, 'r--', 'LineWidth', 2);
% Save the figure as a PDF
    end
end
% determine if a cell is opto-tagged (use the bootstrapping strategy from Cerniauskas 2019)
%% Fotrh subplot: Bootstrapping
% try
for DetectEvents=1:1
    [TimeStamps_events, EventIDs,EventTTLs, EventFeatures, EventventStrings, Header]...
        = Nlx2MatEV([Condition.DataPath,'\',  'Events.nev'], [1 1 1 1 1], 1, 1, []);
    try
    [TimeStamps_Cells, ~, CellNumbers, ~, ~, ~] = Nlx2MatSpike([Condition.DataPath, '\TT',num2str(Variables.TetrodeNumber),'_AutoSort.ntt'], [1 1 1 1 1], 1, 1, []);     
    catch
    [TimeStamps_Cells, ~, CellNumbers, ~, ~, ~] = Nlx2MatSpike([Condition.DataPath, '\TT',num2str(Variables.TetrodeNumber),'_s.ntt'], [1 1 1 1 1], 1, 1, []);     
    end
    
    TimeStamps_Unit=TimeStamps_Cells(CellNumbers==Variables.UnitNumber);
    CellNumbers=CellNumbers(CellNumbers==Variables.UnitNumber);
    % Zero out timestamps based on condition
if TimeStamps_Unit(1) > TimeStamps_events(1)
    TimeStamps_events_zeroed_s = (TimeStamps_events-TimeStamps_events(1))/1000000;
    TimeStamps_Unit_zeroed_s = (TimeStamps_Unit-TimeStamps_events(1))/1000000;
else
    TimeStamps_events_zeroed_s = (TimeStamps_events-TimeStamps_Unit(1))/1000000;
    TimeStamps_Unit_zeroed_s = (TimeStamps_Unit-TimeStamps_Unit(1))/1000000;
end
cell = cat(1,TimeStamps_Unit_zeroed_s,CellNumbers); 
NonUnitLocations = find(CellNumbers ~=Variables.UnitNumber); % "~=1" looks at cell#1, "~=2" looks at cell#2, etc.
cell(:,NonUnitLocations) = []; % clear out cells that do not belong to this unit
laser = cat(1,TimeStamps_events_zeroed_s,EventTTLs);
laser(laser == 0) = NaN;
eventcount = find(laser(2,:) ==Variables.TTLs{1, Condition.ConditionOrder}); %counts the number o laser pulses
NotEvents=find(laser(2,:)~=Variables.TTLs{1, Condition.ConditionOrder});
laser(2,NotEvents) = NaN; %exclude other events
for EventNumber = 1:length(eventcount) %numbers the laser pulses from 1 to last
    laser(2,eventcount(EventNumber)) =  EventNumber;
end
%% Information about time-locked events
% find the number of spikes in each 2ms timebin for all the events
%creat array of timespikes around laser pulses 
% find number of spikes for each event in 1 ms bins
EventBins=[-100/1000:1/1000:100/1000];
% create variable to store the data
EventArray=zeros(length(eventcount),length(EventBins));
PlotArray=zeros(length(eventcount),length(EventBins));
LaserSpikeLocations=false(1,length(TimeStamps_Unit_zeroed_s));
for EventNumber=1:length(eventcount)
%normalise spike times to event onset
tempSpikeTimes=TimeStamps_Unit_zeroed_s-TimeStamps_events_zeroed_s(eventcount(EventNumber));
% store spikes detected in each bin
for Bin=1:length(EventBins)-1
    SpikesFound=find(EventBins(Bin)<=tempSpikeTimes  &  tempSpikeTimes<=EventBins(Bin+1));
SpikeLocations=(EventBins(Bin)<=tempSpikeTimes  &  tempSpikeTimes<=EventBins(Bin+1));
    LaserSpikeLocations(SpikeLocations==1)=true;
if isempty(SpikesFound)
    EventArray(EventNumber,Bin)=0;
    PlotArray(EventNumber,Bin)=0;
else
    EventArray(EventNumber,Bin)=length(SpikesFound);
    PlotArray(EventNumber,Bin)=EventNumber*length(SpikesFound);
end    
end  
end
% sum the values of spikes in each bin
AverageEventArray(1,:)=sum(EventArray,1);
%find the highest response bin time in sec
for EventBinsIndex = 1:length(EventBins) %finds number of values in 0-2ms, 1-3ms, 2-4ms, etc. interval
    count_real(EventBinsIndex) = sum([AverageEventArray(EventBinsIndex) AverageEventArray(EventBinsIndex)]);
end
RangeForTagged=find((0<=EventBins & EventBins<=10/1000));
[highestresponse,highestresponseIndex] = max(count_real(RangeForTagged)); %M = highest number of spikes in interval, I = index of M
highestresponseIndex=RangeForTagged(highestresponseIndex);
highestresponseTimesMS=1000*mean([EventBins(highestresponseIndex) EventBins(highestresponseIndex+1)]);
TotalSpikesInRange = sum(AverageEventArray);
% Enforce opto-tagging based on latency (<8 ms)
if highestresponse > percentile999 && highestresponseTimesMS < 8
    Condition.OptoTagged = true;
    Condition.OptoTag_Latency = highestresponseTimesMS; % Store latency
else
    Condition.OptoTagged = false;
end

% Plot each row of Y against the same X-axis
% figure; hold on;
% for i = 1:size(PlotArray, 1)  % Loop through the 40 rows of Y
%     scatter(EventBins, PlotArray(i, :));  % Scatter plot for each row
% end
% hold off;
%% bootstrapping
bootstrap = zeros(1,10000);
for bootstrapNumber = 1:10000 % do this 10000 times
    %generate simulated data
% Initialize the random data array with zeros
simulatedData = zeros(1, length(AverageEventArray));
% Randomly choose positions to place spikes, allowing repeats
randomIndices = randi([1 length(AverageEventArray)], 1, TotalSpikesInRange);  % Choose positions with replacement
% Add spikes to the selected positions
for i = 1:TotalSpikesInRange
    simulatedData(randomIndices(i)) = simulatedData(randomIndices(i)) + 1;  % Increment the spike count at each index
end
% Run the same bin analysis that we did for the real data
for EventBinsIndex = 1:length(EventBins) %finds number of values in 0-2ms, 1-3ms, 2-4ms, etc. interval
    count_shuffle(EventBinsIndex) = sum([simulatedData(EventBinsIndex) simulatedData(EventBinsIndex)]);
end
bootstrap(bootstrapNumber) = max(count_shuffle);
end
bootstrap = transpose(bootstrap);
pd = fitdist(bootstrap,'Normal');
percentile95 = icdf(pd,0.95);
percentile99 = icdf(pd,0.99);
percentile999 = icdf(pd,0.999);
if highestresponse>percentile95
Condition.Tagged=true;
% Find the percentile
    if highestresponse > percentile999
        Condition.BootstrappingPvalue = 0.001;
        Condition.SignificanceMarker = '***';  % p < 0.001
    elseif highestresponse > percentile99
        Condition.BootstrappingPvalue = 0.01;
        Condition.SignificanceMarker = '**';   % p < 0.01
    else
        Condition.BootstrappingPvalue = 0.05;
        Condition.SignificanceMarker = '*';    % p < 0.05
    end
else
    Condition.Tagged = false;
    Condition.BootstrappingPvalue = cdf(pd, highestresponse) * 100;
    Condition.SignificanceMarker = 'ns';       % p > 0.05
end

if Variables.DisplayPlot
subplot(3,2,2);
x_values = 0:0.1:40;
y = pdf(pd,x_values);
yyaxis left
histogram(bootstrap)
hold on
yyaxis right
plot(x_values,y,'LineWidth',2)
xline(highestresponse, '-r', 'Observed response');
xline(percentile95, '-g', '95%');
title(['Opto-tagging analysis (' Condition.SignificanceMarker ')'])
if highestresponse>percentile95
xlim([0 round(highestresponse*1.2)])
else
xlim([0 round(percentile95*1.2)])
end
end
end
% catch
% end
%% Fifth subplot: Example trace
for ExampleTrace=1:1
% Show example traces
TimeStampsEvents_ms = Condition.TimeStamps_LaserNorm' / 1000;  % Laser event times in ms
TimeStampsUnits_ms = Condition.TimeStamps_unit_Norm' / 1000;  % Unit activity times in ms
spikes_after_laser = [];  % To store spikes 0–10 ms after the laser
for EventNumber = 1:length(TimeStampsEvents_ms)
    laser_time = TimeStampsEvents_ms(EventNumber);
        % Find spikes that occurred 0–10 ms after the laser event
    spikes_in_window = find((TimeStampsUnits_ms...
        >= laser_time & TimeStampsUnits_ms <= laser_time + 100));
        % Collect the spike times
    spikes_after_laser = [spikes_after_laser; spikes_in_window]; %#ok<AGROW>
end
WaveformsAfterLaser=Condition.Samples(:,:,spikes_after_laser);
WaveformsLaserAverage=mean(WaveformsAfterLaser,3,'omitnan');
Xaxis=(1:128);
YaxisLaser=reshape(WaveformsLaserAverage, [], 1);
% Condition.correlation_Laser_Baseline = corr(YaxisBaseline, YaxisLaser);
YaxisLaser=median(Condition.ADBitVolts)*YaxisLaser*1000000;
YaxisLaser=YaxisLaser-YaxisLaser(1);
if Variables.DisplayPlot
subplot(3,2,4); plot(Xaxis,YaxisLaser,'r'); 
xlim([0 128]);
% Add axis titles
xlabel('Channel');
ylabel('Voltage (uV)');
title('Unit average waveform');
xticks([16 48 80 112]);               % Set custom tick positions
xticklabels({'Ch1', 'Ch2', 'Ch3', 'Ch4'});  % Set custom labels
end
[Condition.MaxHeight_uV,ChannelMax]=max(YaxisLaser);
 end
%% Six subplot: Time example
for TimeExample=1:1
ChannelTimestamps=(double(Condition.Ch1.TimeStamps)-Condition.FirstTimestamp)/1000000;
LaserTimestamps=Condition.TimeStamps_LaserNorm/1000000;
NonUnitLocations=Condition.TimeStamps_unit_Norm'  /1000000;
% Preallocate the new array to store the results (15289 * 512)
new_length = length(ChannelTimestamps) * 512;
NewTimestamps = zeros(1, new_length);
% Fill the new array by interpolating between consecutive timestamps
index = 1;
for EventNumber = 1:length(ChannelTimestamps) - 1
    % Generate 512 values between each consecutive pair of timestamps
    ticks = linspace(ChannelTimestamps(EventNumber), ChannelTimestamps(EventNumber+1), 512);
    NewTimestamps(index:index+511) = ticks;  % Add 512 values (including start and end)
    index = index + 512;
end
% For the last value, we simply copy it 512 times
NewTimestamps(end-511:end) = ChannelTimestamps(end);
if 0 <= ChannelMax && ChannelMax <= 32 %ch1
Condition.Ch1.NewTimestamps=NewTimestamps';
Condition.Ch1.uV=reshape(Condition.Ch1.Samples, [], 1);
Condition.Ch1.uV=median(Condition.ADBitVolts)*Condition.Ch1.uV;
Condition.MaxChannel=1;
elseif 32 <= ChannelMax && ChannelMax <= 64 %ch2
Condition.Ch2.NewTimestamps=NewTimestamps';
Condition.Ch2.uV=reshape(Condition.Ch2.Samples, [], 1);
Condition.Ch2.uV=median(Condition.ADBitVolts)*Condition.Ch2.uV;
Condition.MaxChannel=2;
elseif 64 <= ChannelMax && ChannelMax <= 96 %ch3
Condition.Ch3.NewTimestamps=NewTimestamps';
Condition.Ch3.uV=reshape(Condition.Ch3.Samples, [], 1);
Condition.Ch3.uV=median(Condition.ADBitVolts)*Condition.Ch3.uV*10000000;
Condition.MaxChannel=3;
elseif 96 <= ChannelMax && ChannelMax <= 128 %ch4
Condition.Ch4.NewTimestamps=NewTimestamps';
Condition.Ch4.uV=reshape(Condition.Ch4.Samples, [], 1);
Condition.Ch4.uV=median(Condition.ADBitVolts)*Condition.Ch4.uV;
Condition.MaxChannel=4;
end
if Variables.DisplayPlot
    if Condition.MaxChannel==1
subplot(3,2,6); plot(NewTimestamps',Condition.Ch1.uV)
    elseif Condition.MaxChannel==2
subplot(3,2,6); plot(NewTimestamps',Condition.Ch2.uV)
    elseif Condition.MaxChannel==3
subplot(3,2,6); plot(NewTimestamps',Condition.Ch3.uV)
    elseif Condition.MaxChannel==4
subplot(3,2,6); plot(NewTimestamps',Condition.Ch4.uV)
    end
xlim([Condition.TimeXmin Condition.TimeXmax]);  
% find laser events collected in the specified timerange
TTLEvents=LaserTimestamps((Condition.TimeXmin <= LaserTimestamps...
    & LaserTimestamps <= Condition.TimeXmax));
for m=1:length(TTLEvents)  
xline(TTLEvents(m),'r')
end
% find units collected in the specified timerange
UnitStamps=UnitTimeStamps_cellsSeconds((Condition.TimeXmin <= UnitTimeStamps_cellsSeconds...
    & UnitTimeStamps_cellsSeconds <= Condition.TimeXmax));
for m=1:length(UnitStamps)  
xline(UnitStamps(m),'g')
end

legend({['Ch',num2str(Condition.MaxChannel)], 'Laser'});
end
end
if Variables.DisplayPlot
    sgtitle(char(Condition.ConditionName));
   print(fig, '-painters', '-dpdf', fullfile(Variables.UnitGeneralPath, ...
    'figures\',[char(Condition.ConditionName), num2str(Variables.TetrodeNumber), ...
    num2str(Variables.UnitNumber), 'Opto-tagging_Figure.pdf']));

close all;end
close all
 end
 

 