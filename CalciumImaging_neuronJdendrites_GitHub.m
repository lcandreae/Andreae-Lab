%% Code used for analysis of calcium transients along reconstructed dendrites

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







