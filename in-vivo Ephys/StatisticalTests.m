function [Stats] = StatisticalTests(UnitTimes, EventTimes, EventName, Range, Variables)
Condition2Plot='Jelly';
if strcmp(EventName , Condition2Plot); Variables.DisplayPlot=true;end
    % Initial setup and metadata storage
    Stats.Unit = Variables.Unit;
    Stats.UnitIndex = Variables.UnitIndex;
    Stats.DietType = Variables.DietType;
    Stats.MouseName = Variables.MouseName;
    Stats.Date = Variables.Date;
    Stats.TetrodeNumber = Variables.TetrodeNumber;
    Stats.UnitNumber = Variables.UnitNumber;
    Stats.Condition = EventName;
    Stats.Factor=Variables.Factor;
    try
    if isempty(EventTimes)
        % Set all stats to NaN if EventTimes is empty
        Stats.Bouts = []; % Initialize as empty if no bouts are found
        Stats.EventMatrix = NaN;
        Stats.zscoreRaster = NaN;
        Stats.TotalFiringRateBoutHz = NaN;
        Stats.TotalFiringRateBaselineHz = NaN;
        Stats.TotalFiringRateBoutZ = NaN;
        Stats.TotalFiringRateBaselineZ = NaN;
        Stats.TotalBoutTime = NaN;
        Stats.TotalNonBoutTime = NaN;
        Stats.MeanBaseline = NaN;
        Stats.MeanBout = NaN;
        Stats.MedianBaseline = NaN;
        Stats.MedianBout = NaN;
        Stats.Type = NaN;
        Stats.CohensD = NaN;
        Stats.HodgesLehmann = NaN;
        Stats.p_value = NaN;
        Stats.Decision = false;
        Stats.unit_mean=NaN;
        Stats.unit_sem=NaN;
        disp('EventTimes is empty. All statistics set to NaN.');
    else
        %% Bout Analysis Parameters
    if strcmp(Stats.Condition, 'Laser')
    min_bout_size = 1; % Minimum number of timestamps per bout
    max_within_bout_interval = 20 * 1e3; % 20 milliseconds in microseconds
    analysis_window = 10 * 1e3; % 10 milliseconds for analysis window
    elseif strcmp(Stats.Condition, 'TurnLeft')||strcmp(Stats.Condition, 'TurnRight')
    min_bout_size = 1; % Minimum number of timestamps per bout
    max_within_bout_interval = 3 * 1e6; % 3 seconds in microseconds
    analysis_window = 3 * 1e6; % 3 seconds in microseconds for analysis window        
    else
    min_bout_size = 7; % Minimum number of timestamps per bout
    max_within_bout_interval = 3 * 1e6; % 3 seconds in microseconds
    analysis_window = 3 * 1e6; % 3 seconds in microseconds for analysis window        
    end
    % Initialize bout information storage
    Stats.Bouts = [];
    total_bout_duration = 0;
    total_bout_spikes = 0;
    % Arrays to store 3-second pre-bout and in-bout firing rates for statistical comparison
    pre_bout_rates = [];
    in_bout_rates = [];
    pre_bout_medians = [];
    in_bout_medians = [];
    % Total session firing rate for normalization
    session_duration = (Range(end) - Range(1)) / 1e6; % in seconds
    session_total_spikes = length(UnitTimes);
    session_firing_rate = session_total_spikes / session_duration;
    %% Bout Analysis
    current_bout = [EventTimes(1)]; % Start with the first event time
    for i = 2:length(EventTimes)
        interval = EventTimes(i) - EventTimes(i - 1);
        if interval <= max_within_bout_interval
            % Continue current bout
            current_bout = [current_bout, EventTimes(i)];
        else
            % If the bout meets the minimum size, process it
            if length(current_bout) >= min_bout_size
                bout_start = current_bout(1);
                bout_end = current_bout(end);
                bout_duration = (bout_end - bout_start) / 1e6; % Duration in seconds
                % Calculate event frequency in the bout
                event_frequency_bout = length(current_bout) / bout_duration;
                % Calculate unit spikes in the bout
                bout_unit_times = UnitTimes(UnitTimes >= bout_start & UnitTimes <= bout_end);
                unit_event_frequency_bout = length(bout_unit_times) / bout_duration;
                normalized_unit_frequency_bout = unit_event_frequency_bout / session_firing_rate;
                % Calculate 3-second firing rate before and during the bout
                pre_bout_start = max(bout_start - analysis_window, Range(1)); % Start time for 3 seconds before bout
                pre_bout_unit_times = UnitTimes(UnitTimes >= pre_bout_start & UnitTimes < bout_start);
                pre_bout_duration = (bout_start - pre_bout_start) / 1e6; % Pre-bout duration in seconds
                pre_bout_rate = length(pre_bout_unit_times) / pre_bout_duration;
                pre_bout_median = median(pre_bout_unit_times); % Median firing time before bout
                % Check for bout duration constraints (3 seconds) and calculate in-bout firing rate
                in_bout_unit_times = UnitTimes(UnitTimes >= bout_start & UnitTimes < bout_start + analysis_window);
                in_bout_rate = length(in_bout_unit_times) / 3; % Since analysis window is fixed to 3 seconds
                in_bout_median = median(in_bout_unit_times); % Median firing time during bout
                % Store pre-bout and in-bout rates for statistical comparison
                pre_bout_rates = [pre_bout_rates; pre_bout_rate];
                in_bout_rates = [in_bout_rates; in_bout_rate];
                pre_bout_medians = [pre_bout_medians; pre_bout_median];
                in_bout_medians = [in_bout_medians; in_bout_median];
                % Store bout details
                 Location=length(Stats.Bouts)+1;
                Stats.Bouts(Location).StartTime = bout_start;
                Stats.Bouts(Location).EndTime = bout_end;
                Stats.Bouts(Location).LengthSec = bout_duration;
                Stats.Bouts(Location).EventFrequencyHz = event_frequency_bout;
                Stats.Bouts(Location).UnitEventFrequencyHz = unit_event_frequency_bout;
                Stats.Bouts(Location).NormalizedUnitFrequency = normalized_unit_frequency_bout;
                Stats.Bouts(Location).PreBoutRateHz = pre_bout_rate;
                Stats.Bouts(Location).InBoutRateHz = in_bout_rate;
                Stats.Bouts(Location).PreBoutMedian = pre_bout_median;
                Stats.Bouts(Location).InBoutMedian = in_bout_median;
                % Accumulate totals for later calculations
                total_bout_duration = total_bout_duration + bout_duration;
                total_bout_spikes = total_bout_spikes + length(bout_unit_times);
            end
            % Start a new bout with the current event
            current_bout = [EventTimes(i)];
        end
    end

    % Handle the last bout
    if length(current_bout) >= min_bout_size
        bout_start = current_bout(1);
        bout_end = current_bout(end);
        bout_duration = (bout_end - bout_start) / 1e6; % Duration in seconds

        event_frequency_bout = length(current_bout) / bout_duration;

        bout_unit_times = UnitTimes(UnitTimes >= bout_start & UnitTimes <= bout_end);
        unit_event_frequency_bout = length(bout_unit_times) / bout_duration;
        normalized_unit_frequency_bout = unit_event_frequency_bout / session_firing_rate;

        pre_bout_start = max(bout_start - analysis_window, Range(1));
        pre_bout_unit_times = UnitTimes(UnitTimes >= pre_bout_start & UnitTimes < bout_start);
        pre_bout_duration = (bout_start - pre_bout_start) / 1e6;
        pre_bout_rate = length(pre_bout_unit_times) / pre_bout_duration;
        pre_bout_median = median(pre_bout_unit_times);

        in_bout_unit_times = UnitTimes(UnitTimes >= bout_start & UnitTimes < bout_start + analysis_window);
        in_bout_rate = length(in_bout_unit_times) / 3;
        in_bout_median = median(in_bout_unit_times);

        pre_bout_rates = [pre_bout_rates; pre_bout_rate];
        in_bout_rates = [in_bout_rates; in_bout_rate];
        pre_bout_medians = [pre_bout_medians; pre_bout_median];
        in_bout_medians = [in_bout_medians; in_bout_median];
Location=length(Stats.Bouts)+1;
        Stats.Bouts(Location).StartTime = bout_start;
        Stats.Bouts(Location).EndTime = bout_end;
        Stats.Bouts(Location).LengthSec = bout_duration;
        Stats.Bouts(Location).EventFrequencyHz = event_frequency_bout;
        Stats.Bouts(Location).UnitEventFrequencyHz = unit_event_frequency_bout;
        Stats.Bouts(Location).NormalizedUnitFrequency = normalized_unit_frequency_bout;
        Stats.Bouts(Location).PreBoutRateHz = pre_bout_rate;
        Stats.Bouts(Location).InBoutRateHz = in_bout_rate;
        Stats.Bouts(Location).PreBoutMedian = pre_bout_median;
        Stats.Bouts(Location).InBoutMedian = in_bout_median;
        total_bout_duration = total_bout_duration + bout_duration;
        total_bout_spikes = total_bout_spikes + length(bout_unit_times);
    end
        %% Summary Calculations for Overall Bout Metrics
    % Calculate TotalFiringRateBoutHz for bout periods
    Stats.TotalFiringRateBoutHz = total_bout_spikes / total_bout_duration;
    Stats.TotalBoutTime = total_bout_duration;

    % Calculate TotalNonBoutTime and baseline firing rate during non-bout periods
    Stats.TotalNonBoutTime = session_duration - total_bout_duration;
    non_bout_spikes = session_total_spikes - total_bout_spikes; % Spikes outside of bouts
    Stats.TotalFiringRateBaselineHz = non_bout_spikes / Stats.TotalNonBoutTime;
    % Calculate session-wide mean and std firing rate across entire session (bout + non-bout)
    session_firing_rate = session_total_spikes / session_duration; % Mean firing rate over entire session
    % Compute standard deviation of firing rate across the session
    % by dividing the session into bins (e.g., 1s bins) to capture variability
    bin_width = 1e6; % 1 second in microseconds
    session_bins = histcounts(UnitTimes, Range(1):bin_width:Range(end));
    session_firing_rate_std = std(session_bins / (bin_width / 1e6)); % Std over 1-second bins
    % Calculate z-scored firing rates for bouts and baseline
    Stats.TotalFiringRateBoutZ = (Stats.TotalFiringRateBoutHz - session_firing_rate) / session_firing_rate_std;
    Stats.TotalFiringRateBaselineZ = (Stats.TotalFiringRateBaselineHz - session_firing_rate) / session_firing_rate_std;
    % Statistical analysis: compare pre-bout and in-bout rates using Wilcoxon signed-rank test
    if ~isempty(pre_bout_rates) && ~isempty(in_bout_rates)
[Stats.p_value, Stats.Decision,SignedRank] = signrank(pre_bout_rates, in_bout_rates);

% Compute differences
differences = pre_bout_rates - in_bout_rates;
% Remove zero differences
nonzero_differences = differences(differences ~= 0);
% Compute ranks of absolute differences
[~, ~, ranks] = unique(abs(nonzero_differences));  % Ranks of absolute differences
% Assign signs to the ranks
signed_ranks = ranks .* sign(nonzero_differences);
% Compute W (sum of positive ranks)
Stats.zval = sum(signed_ranks(signed_ranks > 0)); % this is actually "W"
if isfield(SignedRank, 'signedrank')
Stats.signedrank=SignedRank.signedrank;
else
        Stats.signedrank=nan;     
end
        Stats.pre_bout_rates = pre_bout_rates;
        Stats.in_bout_rates = in_bout_rates;
        Stats.pre_bout_ratesMean = mean(pre_bout_rates);
        Stats.in_bout_ratesMean = mean(in_bout_rates);
        Stats.pre_bout_ratesMedian = median(pre_bout_rates);
        Stats.in_bout_ratesMedian = median(in_bout_rates);
        if Stats.in_bout_ratesMean> Stats.pre_bout_ratesMean
             Stats.Type=1;
        elseif Stats.in_bout_ratesMean< Stats.pre_bout_ratesMean
             Stats.Type=2;    
        else
            Stats.Type=NaN;    
        end
    else
        Stats.p_value = NaN;
        Stats.Decision = NaN;
        Stats.pre_bout_rates = NaN;
        Stats.in_bout_rates = NaN;
        Stats.pre_bout_ratesMean = NaN;
        Stats.in_bout_ratesMean = NaN;
        Stats.pre_bout_ratesMedian = NaN;
        Stats.in_bout_ratesMedian = NaN;
        Stats.Type=NaN;    

    end
    try

if Variables.DisplayPlot; fig = figure('units', 'normalized', 'outerposition', [0 0 1 1]); end % Create a full-screen figure with 6 subplots

% Define time range for Laser and other conditions
if strcmp(Stats.Condition, 'Laser')
    pre_bout_time = -3 * 1e3;    % 3 milliseconds in microseconds
    post_bout_time = 50 * 1e3;   % 50 milliseconds in microseconds
    time_range = [-0.003, 0.05]; % Time range in seconds for Laser
else
    pre_bout_time = -5 * 1e6;    % 3 seconds in microseconds
    post_bout_time = 20 * 1e6;   % 10 seconds in microseconds
    time_range = [-5, 20];       % Time range in seconds for other conditions
end

AllboutEventTimes = [];
AllboutUnitTimes = [];

% 1. Bottom subplot - Raster plot of Events and Unit Events
if Variables.DisplayPlot; subplot(2, 1, 2);hold on;end
% limit to maximum 50 events

% if length(Stats.Bouts) > 50
%     % Randomly select 50 unique indices from Stats.Bouts
%     LenthLimitIndices = randperm(length(Stats.Bouts), 50);
%     % Use the selected indices to get 50 random elements
%     Stats.BoutsAnalysed = LenthLimitIndices;
% else
%     LenthLimitIndices=1:length(Stats.Bouts);
% end
for boutIdx = 1:length(Stats.Bouts)
    boutStart = Stats.Bouts(boutIdx).StartTime;
    boutEnd = Stats.Bouts(boutIdx).EndTime;
    
    % Align times relative to the start of the current bout
    alignedEventTimes = EventTimes - boutStart;
    alignedUnitTimes = UnitTimes - boutStart;
    
    % Filter to show events within the specified pre and post-bout times
    boutEventTimes = alignedEventTimes(alignedEventTimes >= pre_bout_time & alignedEventTimes <= post_bout_time);
    boutUnitTimes = alignedUnitTimes(alignedUnitTimes >= pre_bout_time & alignedUnitTimes <= post_bout_time);
    AllboutEventTimes = [AllboutEventTimes, boutEventTimes];
    AllboutUnitTimes = [AllboutUnitTimes, boutUnitTimes];
if Variables.DisplayPlot    
    % Plot events (red) and unit events (black) for this bout
scatter(boutUnitTimes / 1e6, boutIdx * ones(size(boutUnitTimes)),100, '|','MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'MarkerEdgeAlpha', 0.8);
hold on
scatter(boutEventTimes / 1e6, boutIdx * ones(size(boutEventTimes)),100, '|', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r', 'MarkerEdgeAlpha', 0.1);
xlim(time_range);
xline(0,'g--');
xlabel('Time (s) from Bout Start');
ylabel('Bout Index');

if Stats.p_value<=0.05 && Stats.p_value>0.01 
    title([Stats.Condition,num2str(Stats.p_value), '*']);  % Set title to the condition name
elseif Stats.p_value<=0.01 && Stats.p_value>0.001 
        title([Stats.Condition,num2str(Stats.p_value), '**']);  % Set title to the condition name
elseif Stats.p_value<=0.001
        title([Stats.Condition,num2str(Stats.p_value), '***']);  % Set title to the condition name
else
        title([Stats.Condition,num2str(Stats.p_value)]);  % Set title to the condition name
end
end
end
hold off;

%% 2. Top subplot - Histogram of Raster Plot (Event and Unit Frequencies with SEM)
if ~strcmp(Stats.Condition, 'Laser')
    if Variables.DisplayPlot ;subplot(2, 1, 1);end

% Define binning parameters for histogram
bin_width = 0.5 * 1e6;  % 100 ms bin width in microseconds
bins = pre_bout_time:bin_width:post_bout_time;
bin_centers = (bins(1:end-1) + bins(2:end)) / 2;  % Bin centers in microseconds

% Initialize matrices to hold individual bout histograms
event_hist_matrix = zeros(length(Stats.Bouts), length(bin_centers));
unit_hist_matrix = zeros(length(Stats.Bouts), length(bin_centers));

% Calculate histogram for each bout and store in matrix
for boutIdx = 1:length(Stats.Bouts)
    boutStart = Stats.Bouts(boutIdx).StartTime;
    
    % Align times relative to the start of the current bout
    alignedEventTimes = EventTimes - boutStart;
    alignedUnitTimes = UnitTimes - boutStart;
    
    % Calculate histograms for this bout
    event_hist_matrix(boutIdx, :) = histcounts(alignedEventTimes, bins);
    unit_hist_matrix(boutIdx, :) = histcounts(alignedUnitTimes, bins);
end
% Calculate mean and SEM across bouts
event_mean = mean(event_hist_matrix) / (bin_width / 1e6);  % Convert to Hz
event_mean=zscore(event_mean);
event_mean=event_mean-mean(event_mean(1:10));
event_sem = std(event_hist_matrix) / sqrt(length(Stats.Bouts)) / (bin_width / 1e6);
Stats.unit_mean = mean(unit_hist_matrix) / (bin_width / 1e6);    % Convert to Hz
Stats.unit_sem = std(unit_hist_matrix) / sqrt(length(Stats.Bouts)) / (bin_width / 1e6);
% unit_mean = Stats.unit_mean/mean(Stats.unit_mean(:,1:5));
unit_mean = zscore(Stats.unit_mean);
unit_mean = unit_mean-unit_mean(1);
unit_sem = zscore(Stats.unit_sem);
% Plot unit frequency on the left y-axis
if Variables.DisplayPlot
% yyaxis left
fill([bin_centers, fliplr(bin_centers)] / 1e6, ...
    [unit_mean + unit_sem, fliplr(unit_mean - unit_sem)], 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
hold on;
plot(bin_centers / 1e6, unit_mean, 'k-', 'LineWidth', 1.5);  % Black line for unit frequency
ylabel('Unit Frequency (Zscore)');
% Plot event frequency on the right y-axis
% yyaxis right
fill([bin_centers, fliplr(bin_centers)] / 1e6, ...
    [event_mean + event_sem, fliplr(event_mean - event_sem)], 'b', 'FaceAlpha', 0.1, 'EdgeAlpha', 0.1,'EdgeColor', 'none');
bin_centers = smooth(bin_centers, 10); % Moving average with a window size of 5
plot(bin_centers / 1e6, event_mean, 'LineWidth', 0.5);  % Red line for event frequency
ylabel('Event Frequency (Hz)');
% yLimits = ylim;
% % Calculate the range of y-axis and increase it by 20%
% yRange = yLimits(2) - yLimits(1);
% newYRange = yRange * 1.2;
% % Calculate the new y-axis limits, centered around the original midpoint
% yMidpoint = mean(yLimits);
% newYLimits = [yMidpoint - newYRange/2, yMidpoint + newYRange/2];
% % Apply the new y-axis limits
% ylim(newYLimits);
ylim([-5 10])
% Set labels and limits
xlim(time_range);
xlabel('Time (s) from Bout Start');
title([num2str(Variables.TetrodeNumber),num2str(Variables.UnitNumber)]);
% legend({'Unit Frequency SEM', 'Unit Frequency Mean', 'Event Frequency SEM', 'Event Frequency Mean'}, 'Location', 'best');
hold off;
fileName = [Variables.MouseName,'',Variables.Date,'',num2str(Variables.TetrodeNumber),num2str(Variables.UnitNumber),'_',Stats.Condition,'.pdf'];
fullSavePath = fullfile(Variables.UnitGeneralPath, 'figures', fileName);
saveas(fig, fullSavePath);  % gcf gets the current figure
end
end
    catch
Stats.unit_mean = NaN;
Stats.unit_sem =  NaN;
    end
    end
    catch
    end

end
