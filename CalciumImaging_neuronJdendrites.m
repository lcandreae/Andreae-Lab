%% Code used for analysis of local calcium transients along reconstructed dendrites
% Code runs with TIFF files in which frames showing widespread fluorescence increases across entire neurons (global calcium transients) were removed, restricting the analysis to the detection of local calcium events

% Short explanation to reconstruct dendrites in NeuronJ and import into MATLAB:
% 1) Open ImageJ
% 2) Open calcium imaging time-series -> create an averaged/maximum projection image -> Save as 8-bit (Image –> type –> 8-bit)
% 3) Open NeuronJ (Plugins -> NeuronJ)
% 4) Load saved image into NeuronJ
% 4) ‘Add tracings’ button to recreate dendrites
% 5) Export tracings as a ‘tab-delimited text file’ (first option)
% 6) Copy & paste data in text file into an excel file to ease importing


%% File names required for code to run


% Name of calcium imaging time-series file
tiffile = 'Example.tif'; % name of TIFF stack

% Name of text file containing x & y coordinates of reconstructed dendrites (exported in NeuronJ)
txtfile = 'Example.txt';

% Name of excel file containing x & y coordinates of reconstructed dendrites
xlfile = 'Example.xls';



%% Import calcium imaging time-series


nFrames = n; % adjust 'n' according to number of frames in calcium imaging time-series

% Load the TIFF image stack
imageStack = imread(tiffile, 'Index', 1); % Load the first frame
for i = 2:nFrames
    frame = imread(tiffile, 'Index', i);
    imageStack = cat(3, imageStack, frame); % Concatenate frames along the third dimension
end



%% Visualize the first frame with dendrites overlayed

% Open image of first frame
frame1 = imageStack(:, :, 1);
figure, imagesc(frame1), colormap(gray), axis image


% Overlay reconstructed dendrites onto dendrites within image
j = 1;
for i = 1:999
    B = importdata(txtfile, '\t', j);
    b1 = size(B);
       if b1(1) > 1
            break
       end
    line(i).data = B.data;
    hold on
    plot(line(i).data(:, 1), line(i).data(:, 2), 'y');
    % Add number labels to the lines
    label_x = mean(line(i).data(:, 1)); % Adjust label position as needed
    label_y = mean(line(i).data(:, 2)); % Adjust label position as needed
    text(label_x, label_y, num2str(i), 'Color', 'y', 'HorizontalAlignment', 'center');   
    ad = size(B.data);
    j = j + ad(1) + 1;
end




%% Read in reconstruction excel file data concatenates all dendrite points

% Code reads in excel data & concatenates all dendrites
[num,txt,all] = xlsread(xlfile); 
x_num = num(:,1);
y_num = num(:,2);
x_num = x_num(~isnan(x_num));
y_num = y_num(~isnan(y_num));

dendritex = x_num;
dendritey = y_num;



%% Create raster plot of raw fluorescence traces for coordinates along all dendrites -> dendrites are concatenated unless each dendrite is run individually


ptsddrite = length(dendritex); % obtain number of datapoints along dendrites

ddrite_rawtrace = zeros(ptsddrite,nFrames);


% Extract fluorescence intensity values from the image stack at the pixel coordinates of the reconstructed dendrites
for i = 1:nFrames
    
    currimage = double(imageStack(:, :, i));
    
    for j = 1:ptsddrite
      ddrite_rawtrace(j,i)  = currimage( dendritey(j) , dendritex(j) );
    end
    
end

figure, imagesc(ddrite_rawtrace)




%% dF/F computation

prctile_constant = 10;  % Value used as percentile cutoff for baseline detection


prct = prctile(ddrite_rawtrace',  prctile_constant  ,1);   % finds 10th percentile value of each trace -> values below this cutoff for each trace will be considered as 'baseline noise values'


% find Baseline (F0) for each trace
for i = 1 : ptsddrite
    idx = find(ddrite_rawtrace(i,:) <= prct(i));          % finds the indices where the values within the trace are <= to the 10th percentile value of that trace
    baselines(1,i) = mean(ddrite_rawtrace(i,idx));        % calculates the mean of the values found in the previous step, which represent the baseline (F0) for the individual trace
end


% Compute dF/F (relative change in fluorescence for each time point in each trace based on the calculated baselines)
dF = (ddrite_rawtrace' - repmat(baselines,nFrames,1)) ./ repmat(baselines,nFrames,1);    
dF = dF';

figure, imagesc(dF)





%% Smooth dF/F traces

window_constant = 10; % smoothing factor

for i = 1:ptsddrite
    curr_ddrite = dF(i,:);
    smooth_dF(i,:) = smooth(curr_ddrite,  window_constant   );
end

figure, imagesc(smooth_dF)

smooth_dendrite = smooth_dF;



%% Threshold determination for event detection in each trace
% The standard deviation of the noise band (lowest 10% values) is multiplied by 30 within each trace to detect the threshold value for event detection

stdMultiply_constant = 30;

ddriteprct = prctile(smooth_dendrite',   prctile_constant   ,1);  % finds 10th percentile value of each trace -> values below this cutoff for each trace will be considered as 'baseline noise values'

% find Threshold for event detection in each trace
for i = 1 : ptsddrite
    idx = find(smooth_dendrite(i,:) <= ddriteprct(i));          % finds the indices where the values within the trace are <= to the 10th percentile value of that trace
    stds(1,i) = std(smooth_dendrite(i,idx));        % calculates the standard deviation of the values found in the previous step, which represent the baseline noise for the individual trace
    
    stdthresholds(1,i) =      stdMultiply_constant      * stds(1,i);   % 30 is used as constant for multiplication with standard deviation of noise band for each trace to determine the cutoff for event detection
end




%% Event detection
% Events contain values in traces that are above the cutoff values (exceed 30 times the
% standard deviation of the individual noise band); events are connected in
% time (across frames) and space (across neighouring coordinates)

% Create a logical matrix to store binary image & double to store non-binary image
binaryImage = false(size(smooth_dendrite));
nonbinaryImage = zeros(size(smooth_dendrite));

for cellIdx = 1 : size(smooth_dendrite,1)
    threshold = stdthresholds(cellIdx);
    
    % Apply the threshold to create the binary image for this cell
    binaryImage(cellIdx, :) = smooth_dendrite(cellIdx, :) > threshold;
    
    % Apply the threshold to create the non-binary image for this cell
    nonbinaryImage(cellIdx, :) = smooth_dendrite(cellIdx, :) .* (smooth_dendrite(cellIdx, :) > threshold);
end

% Display the binary image
figure;
imagesc(binaryImage);

% Display the non-binary image
figure;
imagesc(nonbinaryImage);





%% Exclusion of events smaller than 6 pixels - defined as noise

% Minimum area threshold
minArea = 6;

% Remove small regions from the binary image
binaryImage = bwareaopen(binaryImage, minArea);

% Remove corresponding small regions from the non-binary image
nonbinaryImage(~binaryImage) = 0;

% Display the updated binary image
figure;
imagesc(binaryImage);
title('Size filtered Binary Image');

% Display the updated non-binary image
figure;
imagesc(nonbinaryImage);
title('Size filtered Non-Binary Image');





%% Measurements of properties of events

stats = regionprops(binaryImage, 'Area', 'Centroid', 'BoundingBox');

figure;
imagesc(nonbinaryImage);

% Detects connected components - creates bounding boxes around connected components (events)
for i = 1:numel(stats)
    % Access blob properties
    area = stats(i).Area;
    centroid = stats(i).Centroid;
    boundingBox = stats(i).BoundingBox;

    % Draw bounding box on the original image
    hold on;
    rectangle('Position', boundingBox, 'EdgeColor', 'k', 'LineWidth', 1);
    
    % Calculate label position (center of the bounding box)
    labelX = boundingBox(1) + boundingBox(3) / 2;
    labelY = boundingBox(2) + boundingBox(4) / 2;
    
    % Label the bounding box with component number
    text(labelX, labelY, num2str(i), 'Color', 'y', 'FontSize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    
    hold off;

    % Display properties
    fprintf('Blob %d:\n', i);
    fprintf('Area: %d pixels\n', area);
    fprintf('Centroid: (%.2f, %.2f)\n', centroid(1), centroid(2));
    fprintf('Bounding Box: [%.2f, %.2f, %.2f, %.2f]\n', boundingBox);
end

% Measure properties of connected components in the logical matrix
stats = regionprops(binaryImage, 'PixelIdxList');

% Initialize an array to store the sums of dF/F values in connected components
sumValues = zeros(1, numel(stats));

% Initialize an array to store the peak amplitude of dF/F values in connected components
maxValues = zeros(1, numel(stats));

% Loop through each connected component - obtain data on peak amplitude of events and sum of dF/F values within events
for i = 1:numel(stats)
    pixelIdx = stats(i).PixelIdxList;
    
    valuesInComponent = smooth_dendrite(pixelIdx);

    maxValues(i) = max(valuesInComponent); % Collect peak amplitude dF/F value within each event
    
    % Calculate the sum of dF/F values within the connected component
    sumValues(i) = sum(valuesInComponent);
end



% Obtain data on duration of events
% Measure properties of connected components (blobs) in the binary image
stats = regionprops(binaryImage, 'Area', 'Centroid', 'BoundingBox');

% Initialize arrays to store the widths
widths = zeros(1, numel(stats));

% Loop through each connected component
for i = 1:numel(stats)
    % Access blob properties
    boundingBox = stats(i).BoundingBox;
    
    % Calculate the width (x-axis range)
    width = boundingBox(3);
    
    % Store the width (duration of event)
    widths(i) = width;
end







% Display the sum of dF/F values within each connected component
for i = 1:numel(stats)
    fprintf('Connected Component %d: Sum of dF/F Values = %.2f\n', i, sumValues(i));
end


% Display the widths (duration of events)
fprintf('Widths of all bounding boxes/regions:\n');
disp(widths);


% Display the peak amplitude dF/F value within each event
fprintf('Peak dF/F amplitude of each bounding box/region:\n');
disp(maxValues);



%% DONE
clear all; close all; clc;















% This code detects global calcium events and long local events (>150 s), and tests whether their timing of onset is coupled by computing onset latency 
% The analysis runs on full-frame time-series data.
% It was applied to all neurons showing calcium events lasting longer than 150 s.
% For neurons with two long events, the code was run twice:
%   (1) to test coupling of the longest event with its nearest global event, and
%   (2) again for the second-longest event, after setting the first event to zero 
%       to avoid re-testing the same event.
% Code runs separately from previous code (this code specifically tests global-local event coupling)


%% ------------------------------------------------------------------------
% File names and basic parameters
% -------------------------------------------------------------------------
tiffile = 'Example.tif';   % Calcium imaging TIFF stack
txtfile = 'Example.txt';         % NeuronJ tracing text file
xlfile  = 'Example.xls';         % Excel file with coordinates

nFrames = n;   % Number of frames (adjust)

%% ------------------------------------------------------------------------
% Import calcium imaging time-series
% -------------------------------------------------------------------------
fprintf('Loading image stack: %s (%d frames)...\n', tiffile, nFrames);
imageStack = imread(tiffile, 'Index', 1);
for i = 2:nFrames
    frame = imread(tiffile, 'Index', i);
    imageStack = cat(3, imageStack, frame);
end
fprintf('Image stack loaded (%d × %d × %d).\n', size(imageStack));

%% ========================================================================
%                           GLOBAL EVENT DETECTION
% ========================================================================

%% Visualize first frame with reconstructed dendrites overlayed
frame1 = imageStack(:, :, 1);
figure('Name','First Frame with Dendrites'), imagesc(frame1), colormap(gray), axis image;
title('First Frame with Dendrites Overlay');

% Read and plot dendrite tracings from NeuronJ exported text file
j = 1;
for i = 1:999
    B = importdata(txtfile, '\t', j);
    b1 = size(B);
    if b1(1) > 1
        break
    end
    line(i).data = B.data; %#ok<SAGROW>
    hold on
    plot(line(i).data(:, 1), line(i).data(:, 2), 'y');
    text(mean(line(i).data(:,1)), mean(line(i).data(:,2)), num2str(i), ...
        'Color', 'y', 'HorizontalAlignment', 'center');
    ad = size(B.data);
    j = j + ad(1) + 1;
end
hold off;

%% Read reconstructed dendrite coordinates from Excel
[num, ~, ~] = xlsread(xlfile);
x_num = num(:,1);
y_num = num(:,2);
x_num = x_num(~isnan(x_num));
y_num = y_num(~isnan(y_num));
dendritex = x_num;
dendritey = y_num;

%% Create raster plot of raw fluorescence traces (3x3 pixel averaging)
ptsddrite = length(dendritex);
ddrite_rawtrace = zeros(ptsddrite, nFrames);

windowSize = 1;  % averaging radius (1 -> 3x3 patch)

for fi = 1:nFrames
    currimage = double(imageStack(:, :, fi));
    for pt = 1:ptsddrite
        y = dendritey(pt);
        x = dendritex(pt);
        xMin = max(1, x - windowSize);
        xMax = min(size(currimage, 2), x + windowSize);
        yMin = max(1, y - windowSize);
        yMax = min(size(currimage, 1), y + windowSize);
        localPatch = currimage(yMin:yMax, xMin:xMax);
        ddrite_rawtrace(pt, fi) = mean(localPatch(:));
    end
end

figure('Name','Raw Fluorescence'); imagesc(ddrite_rawtrace);
title('Raw Fluorescence Traces (3×3 Pixel Averaging)');
xlabel('Frame'); ylabel('Dendritic coordinate');

%% dF/F computation
prctile_constant = 10; % percentile cutoff for baseline detection
prct = prctile(ddrite_rawtrace', prctile_constant, 1);

baselines = nan(1, ptsddrite);
for i = 1:ptsddrite
    idx = find(ddrite_rawtrace(i,:) <= prct(i));
    if isempty(idx)
        baselines(i) = mean(ddrite_rawtrace(i,:)); % fallback if no points below percentile
    else
        baselines(i) = mean(ddrite_rawtrace(i, idx));
    end
end

% Avoid dividing by zero baselines
baselines(baselines == 0) = eps;

dF = (ddrite_rawtrace' - repmat(baselines, nFrames, 1)) ./ repmat(baselines, nFrames, 1);
dF = dF';
figure('Name','dF/F Raster'); imagesc(dF);
title('ΔF/F Raster Plot (Unsmoothened)'); xlabel('Frame'); ylabel('Dendritic coordinate');

% For GLOBAL events we do NOT smooth (use raw dF)
smooth_dendrite = dF;

%% Threshold determination for global event detection
stdMultiply_constant = 10; % multiplier for baseline STD to set threshold
ddriteprct = prctile(smooth_dendrite', prctile_constant, 1);

stds = zeros(1, ptsddrite);
stdthresholds = zeros(1, ptsddrite);
for i = 1:ptsddrite
    idx = find(smooth_dendrite(i,:) <= ddriteprct(i));
    if isempty(idx)
        stds(i) = std(smooth_dendrite(i,:));
    else
        stds(i) = std(smooth_dendrite(i, idx));
    end
    stdthresholds(i) = stdMultiply_constant * stds(i);
end

%% Event detection (binary and non-binary masks)
binaryImage = false(size(smooth_dendrite));
nonbinaryImage = zeros(size(smooth_dendrite));

for cellIdx = 1:size(smooth_dendrite,1)
    threshold = stdthresholds(cellIdx);
    binaryImage(cellIdx, :) = smooth_dendrite(cellIdx, :) > threshold;
    nonbinaryImage(cellIdx, :) = smooth_dendrite(cellIdx, :) .* (smooth_dendrite(cellIdx, :) > threshold);
end

figure('Name','Binary Global Events'); imagesc(binaryImage);
title('Binary Event Image'); xlabel('Frame'); ylabel('Dendritic coordinate');
figure('Name','Nonbinary Global Events'); imagesc(nonbinaryImage);
title('Non-Binary Event Image'); xlabel('Frame'); ylabel('Dendritic coordinate');

%% Exclude small events (< 6 pixels)
minArea = 6;
binaryImage = bwareaopen(binaryImage, minArea);
nonbinaryImage(~binaryImage) = 0;

figure('Name','Filtered Binary'), imagesc(binaryImage);
title('Filtered Binary Image'); xlabel('Frame'); ylabel('Dendritic coordinate');
figure('Name','Filtered Nonbinary'), imagesc(nonbinaryImage);
title('Filtered Non-Binary Image'); xlabel('Frame'); ylabel('Dendritic coordinate');

%% Compute fraction of active dendritic coordinates per frame (global activity)
sum_active = sum(binaryImage, 1);
normsum_active = sum_active / size(smooth_dendrite,1);
percent_active = max(normsum_active);

figure('Name','Normalized Active Fraction');
plot(normsum_active, 'k');
ylim([0 1]);
yticks(0:0.25:1);
yticklabels({'0','0.25','0.5','0.75','1'});
yline(0.75, '--', 'LineWidth', 1.5, 'Color', [0.5 0.5 0.5]);
title('Proportion of Active Coordinates per Frame');
xlabel('Frame'); ylabel('Proportion active');

figure('Name','Global Summary');
subplot(3,1,1), imagesc(nonbinaryImage), title('Nonbinary Global Events');
subplot(3,1,2), imagesc(binaryImage), title('Binary Global Events');
subplot(3,1,3), plot(normsum_active, 'k');
ylim([0 1]); yticks(0:0.25:1); yline(0.75, '--', 'LineWidth', 1);
title('Normalized Sum Active');

% Keep for later alignment
binaryglobal = binaryImage;
nonbinaryglobal = nonbinaryImage;






%% ========================================================================
%                           LOCAL EVENT DETECTION
% ========================================================================

%% Visualize first frame again with dendrites (local detection)
figure('Name','First Frame with Dendrites (Local)'), imagesc(frame1), colormap(gray), axis image;
title('First Frame with Dendrites Overlay (Local)');

% Re-read NeuronJ tracings for plotting (same code as above)
j = 1;
for i = 1:999
    B = importdata(txtfile, '\t', j);
    b1 = size(B);
    if b1(1) > 1
        break
    end
    line(i).data = B.data; %#ok<SAGROW>
    hold on
    plot(line(i).data(:, 1), line(i).data(:, 2), 'y');
    text(mean(line(i).data(:,1)), mean(line(i).data(:,2)), num2str(i), ...
        'Color', 'y', 'HorizontalAlignment', 'center');
    ad = size(B.data);
    j = j + ad(1) + 1;
end
hold off;

%% Read coordinates (again to ensure local section has them)
[num, ~, ~] = xlsread(xlfile);
x_num = num(:,1);
y_num = num(:,2);
x_num = x_num(~isnan(x_num));
y_num = y_num(~isnan(y_num));
dendritex = x_num;
dendritey = y_num;

%% Create raster plot of raw fluorescence traces (single pixel here)
ptsddrite = length(dendritex);
ddrite_rawtrace = zeros(ptsddrite, nFrames);

for fi = 1:nFrames
    currimage = double(imageStack(:, :, fi));
    for pt = 1:ptsddrite
        ddrite_rawtrace(pt, fi) = currimage(dendritey(pt), dendritex(pt));
    end
end

figure('Name','Raw Fluorescence (Local)'), imagesc(ddrite_rawtrace);
title('Raw Fluorescence Traces (Local)'); xlabel('Frame'); ylabel('Dendritic coordinate');

%% dF/F computation (local)
prctile_constant = 10;
prct = prctile(ddrite_rawtrace', prctile_constant, 1);

baselines = nan(1, ptsddrite);
for i = 1:ptsddrite
    idx = find(ddrite_rawtrace(i,:) <= prct(i));
    if isempty(idx)
        baselines(i) = mean(ddrite_rawtrace(i,:));
    else
        baselines(i) = mean(ddrite_rawtrace(i, idx));
    end
end

baselines(baselines == 0) = eps; % avoid division by zero
dF = (ddrite_rawtrace' - repmat(baselines, nFrames, 1)) ./ repmat(baselines, nFrames, 1);
dF = dF';
figure('Name','dF/F (Local)'), imagesc(dF);
title('ΔF/F Raster Plot (Local)'); xlabel('Frame'); ylabel('Dendritic coordinate');

%% Smooth dF/F traces for local event detection
window_constant = 10; % smoothing factor (frames)
smooth_dF = nan(size(dF));
for i = 1:ptsddrite
    curr_ddrite = dF(i,:);
    smooth_dF(i,:) = smooth(curr_ddrite, window_constant);
end
figure('Name','Smoothed dF/F (Local)'), imagesc(smooth_dF);
title('Smoothed ΔF/F (Local)');

smooth_dendrite = smooth_dF;

%% Threshold determination for local event detection
stdMultiply_constant = 25; % multiplier used previously for local detection
ddriteprct = prctile(smooth_dendrite', prctile_constant, 1);

stds = zeros(1, ptsddrite);
stdthresholds = zeros(1, ptsddrite);
for i = 1:ptsddrite
    idx = find(smooth_dendrite(i,:) <= ddriteprct(i));
    if isempty(idx)
        stds(i) = std(smooth_dendrite(i,:));
    else
        stds(i) = std(smooth_dendrite(i, idx));
    end
    stdthresholds(i) = stdMultiply_constant * stds(i);
end

%% Event detection (local): binary and non-binary
binaryImage = false(size(smooth_dendrite));
nonbinaryImage = zeros(size(smooth_dendrite));

for cellIdx = 1:size(smooth_dendrite,1)
    threshold = stdthresholds(cellIdx);
    binaryImage(cellIdx, :) = smooth_dendrite(cellIdx, :) > threshold;
    nonbinaryImage(cellIdx, :) = smooth_dendrite(cellIdx, :) .* (smooth_dendrite(cellIdx, :) > threshold);
end

figure('Name','Binary Local Events'); imagesc(binaryImage);
title('Binary Local Events'); xlabel('Frame'); ylabel('Dendritic coordinate');
figure('Name','Nonbinary Local Events'); imagesc(nonbinaryImage);
title('Non-Binary Local Events'); xlabel('Frame'); ylabel('Dendritic coordinate');

%% Remove small local events (<6 pixels)
minArea = 6;
binaryImage = bwareaopen(binaryImage, minArea);
nonbinaryImage(~binaryImage) = 0;

figure('Name','Filtered Binary Local'), imagesc(binaryImage);
title('Filtered Binary Image (Local)'); xlabel('Frame'); ylabel('Dendritic coordinate');
figure('Name','Filtered Nonbinary Local'), imagesc(nonbinaryImage);
title('Filtered Non-Binary Image (Local)'); xlabel('Frame'); ylabel('Dendritic coordinate');

%% Measurements of properties of local events
stats = regionprops(binaryImage, 'Area', 'Centroid', 'BoundingBox');
figure('Name','Local Events with Bounding Boxes'), imagesc(nonbinaryImage);
title('Local Events (non-binary)'); hold on;
for i = 1:numel(stats)
    area = stats(i).Area;
    centroid = stats(i).Centroid;
    boundingBox = stats(i).BoundingBox;
    rectangle('Position', boundingBox, 'EdgeColor', 'k', 'LineWidth', 1);
    labelX = boundingBox(1) + boundingBox(3) / 2;
    labelY = boundingBox(2) + boundingBox(4) / 2;
    text(labelX, labelY, num2str(i), 'Color', 'y', 'FontSize', 12, 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    fprintf('Blob %d: Area=%d, Centroid=(%.2f,%.2f), BBox=[%.2f %.2f %.2f %.2f]\n', ...
        i, area, centroid(1), centroid(2), boundingBox);
end
hold off;

% Sum and peak values per connected component
statsPix = regionprops(binaryImage, 'PixelIdxList');
sumValues = zeros(1, numel(statsPix));
maxValues = zeros(1, numel(statsPix));
for i = 1:numel(statsPix)
    pixIdx = statsPix(i).PixelIdxList;
    valuesInComponent = smooth_dendrite(pixIdx);
    maxValues(i) = max(valuesInComponent);
    sumValues(i) = sum(valuesInComponent);
end

% Get widths (duration) from bounding boxes
statsBBox = regionprops(binaryImage, 'BoundingBox');
widths = zeros(1, numel(statsBBox));
for i = 1:numel(statsBBox)
    widths(i) = statsBBox(i).BoundingBox(3); % width in frames
end

% Display statistics
for i = 1:numel(statsPix)
    fprintf('Connected Component %d: Sum dF/F = %.2f, Peak = %.2f, Width = %.2f\n', ...
        i, sumValues(i), maxValues(i), widths(i));
end

% Save local binary/nonbinary for alignment
binarylocal = binaryImage;
nonbinarylocal = nonbinaryImage;








%% ========================================================================
%                                ALIGNMENT
% ========================================================================

%% STEP 0: Extract Longest Local Event (preserve original nonbinary signal)
% Use widths from local detection to select longest event
if isempty(widths)
    error('No local events found; cannot extract longest event for alignment.');
end

[~, idxLongest] = max(widths); % index of longest local event

% Get pixel indices of connected components (from local binary)
statsPixelIdx_local = regionprops(binarylocal, 'PixelIdxList');
if idxLongest > numel(statsPixelIdx_local)
    error('Index mismatch: longest event index exceeds available connected components.');
end
longestEventPixels = statsPixelIdx_local(idxLongest).PixelIdxList;

longestEvent = zeros(size(binarylocal));
longestEvent(longestEventPixels) = nonbinarylocal(longestEventPixels);

figure('Name','Longest Local Event'), imagesc(longestEvent);
title('Isolated Longest Local Event'); xlabel('Frame'); ylabel('Dendritic coordinate');

%% STEP 1: Extract Global Events Around Peaks
% Identify global peaks using normalized sum active (normsum_active computed earlier)

threshold = 0.75;   % threshold used for global peak detection (proportion)
window = 10;        % ± frames around each global peak for event extraction
nFramesVec = length(normsum_active);

% Use findpeaks if available
try
    [~, peakFrames] = findpeaks(normsum_active, 'MinPeakHeight', threshold);
catch
    % Fallback simple peak detection: any frame where normsum_active > threshold
    peakFrames = find(normsum_active > threshold);
    % Optionally reduce to local maxima:
    peakFrames = peakFrames([true, diff(normsum_active(peakFrames)) > 0]); %#ok<AGROW>
end

% Build mask for ±window around each peak
global_event_mask = false(1, nFramesVec);
for i = 1:length(peakFrames)
    f = peakFrames(i);
    minF = max(1, f - window);
    maxF = min(nFramesVec, f + window);
    global_event_mask(minF:maxF) = true;
end

binary_global_extracted = binaryglobal;
binary_global_extracted(:, ~global_event_mask) = 0;

nonbinary_global_extracted = nonbinaryglobal;
nonbinary_global_extracted(:, ~global_event_mask) = 0;

%% STEP 2: Visualize Binary Overlay (Global = Red, Longest Local = Green)
[rows, cols] = size(binarylocal);
RGB_binary = zeros(rows, cols, 3);
RGB_binary(:,:,1) = binary_global_extracted;     % Red = global
RGB_binary(:,:,2) = longestEvent > 0;            % Green = longest local

figure('Name','Binary Overlay Global vs Longest Local'), image(RGB_binary);
axis tight; xlabel('Frame'); ylabel('Dendritic coordinate');
title('Binary Overlay: Global (Red), Longest Local (Green)');

%% STEP 3: Visualize Nonbinary Overlay (Global = White, Local = Orange)
% Normalize to shared max to maintain relative brightness
maxVal = max([max(nonbinary_global_extracted(:)), max(longestEvent(:))]);
if isempty(maxVal) || maxVal == 0
    maxVal = 1; % avoid division by zero
end

nonbinary_global_norm = nonbinary_global_extracted / maxVal;
nonbinary_local_norm  = longestEvent / maxVal;

% Brighten global for visibility
globalBoost = 1.5;
nonbinary_global_norm = min(nonbinary_global_norm * globalBoost, 1);

RGB_nonbinary = zeros(rows, cols, 3);
% Local in orange (R + 0.5 * G)
RGB_nonbinary(:,:,1) = nonbinary_local_norm;
RGB_nonbinary(:,:,2) = 0.5 * nonbinary_local_norm;
% Add global (white) to all channels (use max to overlay)
RGB_nonbinary(:,:,1) = max(RGB_nonbinary(:,:,1), nonbinary_global_norm);
RGB_nonbinary(:,:,2) = max(RGB_nonbinary(:,:,2), nonbinary_global_norm);
RGB_nonbinary(:,:,3) = max(RGB_nonbinary(:,:,3), nonbinary_global_norm);

figure('Name','Nonbinary Overlay (Global + Local)'), image(RGB_nonbinary);
axis tight; xlabel('Frame'); ylabel('Dendritic coordinate');
title('Nonbinary Overlay: Global (White), Longest Local (Orange)');

%% STEP 4: Extract & Visualize Closest Global Event to Longest Local Event
% Determine onset frame of longest local event
statsLocal = regionprops(logical(longestEvent), 'BoundingBox');
if isempty(statsLocal)
    error('No bounding box found for longestEvent.');
end
local_onset_frame = max(1, round(statsLocal(1).BoundingBox(1)));

% If no peaks were found, throw a warning and use the maximum of normsum_active
if isempty(peakFrames)
    warning('No global peaks found via threshold; using frame with maximal normsum_active as global onset.');
    [~, fmax] = max(normsum_active);
    peakFrames = fmax;
end

% Find the global peak closest in time to the local onset
latencies = abs(peakFrames - local_onset_frame);
[~, idxClosest] = min(latencies);

% Create mask for only the closest global event ± window frames
global_event_mask = false(1, nFramesVec);
f = peakFrames(idxClosest);
minF = max(1, f - window);
maxF = min(nFramesVec, f + window);
global_event_mask(minF:maxF) = true;

binary_global_extracted = binaryglobal;
binary_global_extracted(:, ~global_event_mask) = 0;

nonbinary_global_extracted = nonbinaryglobal;
nonbinary_global_extracted(:, ~global_event_mask) = 0;

% Visualize Binary overlay (closest global vs longest local)
RGB_binary = zeros(rows, cols, 3);
RGB_binary(:,:,1) = binary_global_extracted;        % Red = closest global
RGB_binary(:,:,2) = longestEvent > 0;               % Green = longest local

figure('Name','Binary Overlay Closest Global vs Longest Local'), image(RGB_binary);
axis tight; xlabel('Frame'); ylabel('Dendritic coordinate');
title('Binary Overlay: Closest Global Event (Red), Longest Local Event (Green)');

% Visualize Nonbinary overlay (closest global vs longest local)
maxVal = max([max(nonbinary_global_extracted(:)), max(longestEvent(:))]);
if isempty(maxVal) || maxVal == 0
    maxVal = 1;
end

nonbinary_global_norm = nonbinary_global_extracted / maxVal;
nonbinary_local_norm = longestEvent / maxVal;
nonbinary_global_norm = min(nonbinary_global_norm * globalBoost, 1);

RGB_nonbinary = zeros(rows, cols, 3);
RGB_nonbinary(:,:,1) = nonbinary_local_norm;
RGB_nonbinary(:,:,2) = 0.5 * nonbinary_local_norm;
RGB_nonbinary(:,:,1) = max(RGB_nonbinary(:,:,1), nonbinary_global_norm);
RGB_nonbinary(:,:,2) = max(RGB_nonbinary(:,:,2), nonbinary_global_norm);
RGB_nonbinary(:,:,3) = max(RGB_nonbinary(:,:,3), nonbinary_global_norm);

figure('Name','Nonbinary Overlay Closest Global vs Local'), image(RGB_nonbinary);
axis tight; xlabel('Frame'); ylabel('Dendritic coordinate');
title('Nonbinary Overlay: Closest Global Event (White), Longest Local Event (Orange)');

%% Combine global and local into one nonbinary raster (preserve global values)
combinedImage = nonbinary_global_extracted;  % Start with global
localOnlyMask = (nonbinary_global_extracted == 0) & (longestEvent ~= 0);
combinedImage(localOnlyMask) = longestEvent(localOnlyMask);

figure('Name','Combined Global+Local'), imagesc(combinedImage);
colormap(turbo); colorbar;
axis tight; xlabel('Frame'); ylabel('Dendritic coordinate');
title('Overlay: Global Preserved, Local Added (No Overlap Addition)');

%% ========================================================================
%                      BINARY ONSET LATENCY & VISUALIZATION
% ========================================================================

% Find earliest frame where any pixel is active for global and local masks
[~, global_onset_frame] = max(any(binary_global_extracted, 1));
[~, local_onset_frame]  = max(any(longestEvent > 0, 1));

observed_latency = abs(global_onset_frame - local_onset_frame);
fprintf('Latency (onset) in frames (binary): %d\n', observed_latency);

% Visualization: overlay local and nearest global event and annotate latency
[rows, nFramesTotal] = size(binary_global_extracted);

RGB = zeros(rows, nFramesTotal, 3);
RGB(:,:,1) = binary_global_extracted;  % Red: global event
RGB(:,:,2) = longestEvent > 0;         % Green: local event

figure('Name','Latency Overlay'); image(RGB);
xlabel('Frame'); ylabel('Dendritic coordinate');
title(sprintf('Overlay of Longest Local and Nearest Global Event\nLatency = %d frames', observed_latency), 'FontWeight', 'bold');
hold on;

% Vertical lines at onsets
xline(local_onset_frame, '--g', 'LineWidth', 2, 'Label', 'Local Onset', ...
    'LabelOrientation', 'horizontal', 'LabelVerticalAlignment', 'bottom');
xline(global_onset_frame, '--r', 'LineWidth', 2, 'Label', 'Global Onset', ...
    'LabelOrientation', 'horizontal', 'LabelVerticalAlignment', 'bottom');

% Zoom window around local onset
xlim([max(1, local_onset_frame - 15), min(nFramesTotal, local_onset_frame + 30)]);
ylim([1, rows]);

% Add horizontal arrow/line showing latency
yArrow = rows + 2;
plot([local_onset_frame global_onset_frame], [yArrow yArrow], 'k-', 'LineWidth', 2);
plot(local_onset_frame, yArrow, 'kv', 'MarkerFaceColor', 'g', 'MarkerSize', 8);
plot(global_onset_frame, yArrow, 'kv', 'MarkerFaceColor', 'r', 'MarkerSize', 8);

midpoint = (local_onset_frame + global_onset_frame) / 2;
text(midpoint, yArrow + 1, sprintf('%d frames', observed_latency), ...
    'HorizontalAlignment', 'center', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'k');

% Add a second, bold white-on-black label slightly above for contrast
text(midpoint, yArrow + 3, sprintf('%d frames', observed_latency), ...
    'HorizontalAlignment', 'center', 'FontSize', 12, 'FontWeight', 'bold', ...
    'Color', 'w', 'BackgroundColor', 'k', 'Margin', 2);

ylim([1, yArrow + 5]);
legend({'Global (Red)', 'Local (Green)'});
hold off;

%% DONE
fprintf('Processing complete.\n');
























