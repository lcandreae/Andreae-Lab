%% Code was used for the analysis of calcium transients in acute slices treated with NMDA/glycine in different hippocampal layers (SO, SR, SLM)
% Data from imaging is stored in Excel files
% Excel file contains raw fluorescence trace (mean values within a dendritic ROI)
% Excel file also contains background fluorescence trace (mean values within a small square region of the background that exhibits no changes in fluorescence over time)



%% dF/F computation & peak amplitude measurement

% Import data
calcium_data = xlsread(file); % file = file name

% Extract fluorescence within dendritic ROI
raw_fluorescence = calcium_data(:,2); % Dendritic fluorescence data stored in second column within Excel file

% Plot original data
figure;
plot(raw_fluorescence);
title('Raw Fluorescence Data');
xlabel('Time (frames)');
ylabel('Fluorescence');

% Background subtraction to correct for background noise
background = calcium_data(:,3);  % extracts background fluorescence trace stored in third column within Excel file
avgbackground = mean(background);
corrected_fluorescence = raw_fluorescence - avgbackground;

% Plot background-subtracted fluorescence data
figure;
plot(corrected_fluorescence);
title('Background Subtracted Fluorescence');
xlabel('Time (frames)');
ylabel('Fluorescence');

% Detect baseline fluorescence
x_coordinates = round(ginput(2));   % Select two points (via 2 clicks) that border a segment of the fluorescence trace clearly representing baseline fluorescence
baselinepre = mean(corrected_fluorescence(x_coordinates(1,1) : x_coordinates(2,1)));   % mean baseline fluorescence before (pre) dF/F

% dF/F computation
norm_data = (corrected_fluorescence - baselinepre) / baselinepre;

% Plot normalized data
figure;
plot(norm_data);
title('Normalized dF/F Data');
xlabel('Time (frames)');
ylabel('dF/F');

% Determine peak amplitude
baselinepost = mean(norm_data(x_coordinates(1,1) : x_coordinates(2,1)));  % Baseline after (post) dF/F computation
min_value = min(norm_data);
max_value = max(norm_data);
peak_amp = max_value - min_value;  % peak amplitude is determined as the difference between the min and max value within the dF/F trace




%% Code to determine the duration of the calcium transient
% Step 1: determine threshold for event detection (values above the treshold are accepted as part of transient)

% Calculate standard deviation of the baseline noise
baseline = norm_data(x_coordinates(1,1) : x_coordinates(2,1));
std_baseline = std(baseline);  % Standard deviation of the baseline trace

% Define the transient threshold as 'n-times' the standard deviation of the baseline (we used '6' as our constant - may need to adjust based on how noisy your data is)
threshold = 6 * std_baseline;




%% Determine the duration of the transient
% Step 2: find values above the threshold and determine duration of transient

% Find indices where the normalized data exceeds the threshold
above_threshold = norm_data > threshold;

% Visualize normalized data and threshold line for calcium transient
figure;
plot(norm_data, 'k');
hold on;
yline(threshold, 'r', 'LineWidth', 2);
title('Transient Detection with Threshold');
xlabel('Time (frames)');
ylabel('dF/F');

% Identify the first and last values above the threshold
start_idx = find(above_threshold, 1, 'first');  % First value above threshold (start of transient)
end_idx = find(above_threshold, 1, 'last');  % Last value above threshold (end of transient)

% Plot start & end points of the transient onto the dF/F trace for visualization
plot(start_idx, norm_data(start_idx), 'ko');
plot(end_idx, norm_data(end_idx), 'ko');

% Count number of values above the threshold -> this equals the duration of the transient (in frames)
duration_frames = sum(above_threshold);

% Convert from frames to seconds
frame_rate = n;  % Adjust 'n' based on your frame rate (frames per second)
duration = duration_frames / frame_rate;

% Display the duration
disp(['Transient duration in frames: ', num2str(duration_frames)]);
disp(['Transient duration in seconds: ', num2str(duration)]);



