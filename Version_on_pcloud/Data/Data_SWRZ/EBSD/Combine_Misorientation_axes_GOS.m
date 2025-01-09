
clear all 
close all
%% Do we need to add MTEC ptah

addpath("C:\Users\User\Documents\mtex-5.11.2");
startup_mtex;

% Read the Excel spreadsheet


% Filter rows where 'Deformed' column is 'Y'

spreadsheet_path = "/Users/cljd/pCloud Drive/WORK-GENERAL/POSTDOC-UCB/BERKELEY-VIBE/Documents/Projects/Kil_SWRZ/EBSD/Rapid_response_day2/Day2_edited.xlsx";
data_table = readtable(spreadsheet_path);
filtered_data = data_table(strcmp(data_table.Deformed, 'Y'), :);

   

% Initialize an empty cell array to store the loaded data
data_combined = [];
data_combined_norm=[];

% Initialize a cell array to store the filenames used
filenames_used = {};
array_length={};
num_sample_points=10000

% Lets define a function for sampling
function sampled_set = sample_quaternions(quaternion_set, num_points)
    num_data_points = length(quaternion_set);
    if num_data_points >= num_points
        % Randomly sample without replacement if we have enough data points
        indices = randperm(num_data_points, num_points);
        sampled_set = quaternion_set(indices);
    else
        % Replicate the data points to reach the desired number
        replication_factor = ceil(num_points / num_data_points);
        replicated_set = repmat(quaternion_set, replication_factor, 1);
        sampled_set = replicated_set(1:num_points);
    end
end


% Loop through each row in the filtered data
for i = 1:height(filtered_data)
    % Get the filename and modify it
    filename = filtered_data.Filename{i};
    % filename = strrep(filename, '\', ''); % Remove backslashes
    
    % Get the grainID
    grainID = filtered_data.grainID(i);
    
    % Add grainID and _misorientation to the filename
    filename = [filename '_grain' num2str(grainID) '_misorientation'];

    % Determine the folder path based on the value in the 'Day' column
    day = filtered_data.Day(i);
     if day == 2
        folder_path = "/Users/cljd/pCloud Drive/WORK-GENERAL/POSTDOC-UCB/BERKELEY-VIBE/Documents/Projects/Kil_SWRZ/EBSD/Rapid_response_day2/h5oina/Mis_axes";
 
     else
        error('Unexpected value in the Day column');
    end


 

    % Construct the full path to the .mat file
    file_path = fullfile(folder_path, [filename, '.mat']);

    % Load the specific variable from the file
    if isfile(file_path)
        data = load(file_path, 'save_mis');
        len_ar=length(data.save_mis)
        array_length=[array_length, len_ar]

        % Lets try normalizing them

        % Append the loaded data to the combined array
        data_combined = [data_combined; data.save_mis];

        % Now lets try normalizing them
        weighted=sample_quaternions(data.save_mis, num_sample_points);
        
        data_combined_norm=[data_combined_norm; weighted]

        % Append the filename to the list of filenames used
        filenames_used{end+1} = file_path;
    else
        fprintf('File not found: %s\n', file_path);
    end
end
filenames_used = filenames_used';
%% Make a two part figure
close all
figure;

cmax=12
mtexFig = mtexFigure('ncols', 2, 'nrows', 1, 'spacing', 'loose', 'margin', 1);
% Create the 2 by 1 subplot (2 columns, 1 row)
nextAxis;
plotAxisDistribution(data_combined, 'contourf');
mtexColorMap LaboTeXColorMap
%colorbar;
title('Combined Axis Distribution (no scaling)');
clim([0, cmax])
% Assuming you have another set of combined data for normalization per grain
% For the sake of example, we will use the same combined_quaternions
% Replace this with your actual normalized data

nextAxis;
plotAxisDistribution(data_combined_norm, 'contourf');
mtexColorMap LaboTeXColorMap
colorbar;
title('Axis Distribution with Normalization per Grain');
clim([0, cmax])

%% Make a two part figure
close all
figure;

cmax=12
mtexFig = mtexFigure('ncols', 2, 'nrows', 1, 'spacing', 'loose', 'margin', 1);
% Create the 2 by 1 subplot (2 columns, 1 row)
nextAxis;
plotAxisDistribution(data_combined, 'contourf');

%colorbar;
title('Combined Axis Distribution (no scaling)');

% Assuming you have another set of combined data for normalization per grain
% For the sake of example, we will use the same combined_quaternions
% Replace this with your actual normalized data

nextAxis;
plotAxisDistribution(data_combined_norm, 'contourf');

colorbar;
title('Axis Distribution with Normalization per Grain');
