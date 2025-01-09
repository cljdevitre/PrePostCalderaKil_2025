clear all 
close all

%% ADD MTEC ptah
startup_mtex;

%% Specify Crystal and Specimen Symmetries _ Load depending on analysis setup
% crystal symmetry
% This is used to ID subgrains
subgrain_angle=0.3; 
% This is used to pull out what counts as a separate grain when plotting each grain separatly
grain_div_angle=3; 
color_angle=3;
len_filter=1 % filter out shorter grain boundaries than this um
results = {};

CS = {... 
  'notIndexed',...
  crystalSymmetry('mmm', [4.8 10 6], 'mineral', 'Forsterite', 'color', [0.53 0.81 0.98]),...
  crystalSymmetry('m-3m', [8.358 8.358 8.358], 'mineral', 'Chromite', 'color', 'light green'),...
  crystalSymmetry('Pbca', [7.048 7.193 18.159], 'mineral', 'Enstatite', 'color', 'light yellow')};

% plotting convention
setMTEXpref('xAxisDirection','east');
setMTEXpref('zAxisDirection','outOfPlane');
%% Loading individual files

day_set='day1'
% path to filesf
if strcmp(day_set, 'day1');
       pname="P:\WORK-GENERAL\POSTDOC-UCB\BERKELEY-VIBE\Documents\Projects\Kil_SWRZ\EBSD\Rapid_response_day1\h5oina";
else
    error('Invalid day_set value');
end
%% 

folder = fullfile(pname, 'mis2mean_unfilt');

% Load Individual EBSD Map  - special cases, MLP_9, MP2_63, MLP_43,
% MLP_39,MLP_49, MP2_67, MP2_59b_analysis100, MP2_58B_hough60_refined 
% Get a list of all .ctf files in the directory
%files = dir(fullfile(pname, '\MLP_9 - EBSD Data.ctf')); % was '*.ctf'
files = dir(fullfile(pname, '\Rapid_response_day1 Specimen 2 919_70 Map Data 71.h5oina')); % was '*.ctf'#'\*.h5oina'
for i = 1:length(files)
    display(files(i))
% We need to clear variables each time thro

if contains(files(i).name, 'ctf')
    filename = ['\' files(i).name];
    filename=strrep(filename, ' - EBSD Data','');
    filename=strrep(filename, '.ctf','');
    fname = fullfile(pname, files(i).name);

    ebsd = EBSD.load(fname,CS,'interface','ctf',...
  'convertEuler2SpatialReferenceFrame');
%ebsd=rotate(ebsd, 180*degree, 'keepXY')

else
    
    filename = ['\' files(i).name];
    filename=strrep(filename, ' Map Data ','_');
    filename=strrep(filename, '.h5oina','');
    fname = fullfile(pname, files(i).name);
    ebsd = EBSD.load(fname,CS,'interface','h5oina');
end
ebsd.unitCell = ebsd.unitCell;


%% Data Filtering
[grains, ebsd.grainId, ebsd.mis2mean]= calcGrains(ebsd,'angle',10*degree);
ebsd(grains(grains.grainSize<2)).phase=0; %removes misindexed pixels 


%%

%Apply 5 neighbour Kuwahara filter
F = KuwaharaFilter; 
F.numNeighbours = 5;
ebsd= smooth(ebsd('indexed'), F, 'fill',grains);




%% Subgrain reconstruction
[grains,ebsd.grainId,ebsd.mis2mean]= calcGrains(ebsd,'angle',[subgrain_angle*degree,grain_div_angle*degree]); % Want this big, so it pulls out each grain separatly. 
iter=[1,5,8];
for Z=1:length(iter)
    grains=smooth(grains, iter(Z));
end
grains=grains(grains.grainSize > max(grains.grainSize/30)); % Only use grains within 20X size of max grains
grains=grains('Fo'); %Construct grains variable that is restricted to only the phase forsterite

remainingGrainIds = [grains.id];
isRemainingGrain = ismember(ebsd.grainId, remainingGrainIds);
ebsd_Fo = ebsd(isRemainingGrain);


%% Mis 2 mean plot
mtexFig = newMtexFigure('layout',[1, 2]); % 1 row, 2 columns
plot(ebsd_Fo,ebsd_Fo.mis2mean.angle./degree);
hold on
title('Mis2Mean Coloring')
plot(grains.innerBoundary, 'lineColor', 'k', 'width', 2, label='subgrain inner')
plot(grains.boundary, 'lineColor', 'k', label='subgrain bound')

clim([0 color_angle]);
%mtexColorbar;
h = mtexColorbar;
h.Label.String = 'Mis2Mean angle';  % Set the label text



% Loop through grains to label them with their GOS value. 
for j = 1:length(grains)
    
    % Calculate the center of the grain
    centerX = mean(grains(j).boundary.x);
    centerY = mean(grains(j).boundary.y);
    
    % Get the grain ID and GOS value
    grainId = grains(j).id;
    GOS = grains(j).GOS ./ degree;
    
    
    % Create a text label for the grain ID and GOS
    text(centerX, centerY, sprintf('ID: %d\nGOS: %.2f', grainId, GOS), ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
         'FontSize', 8, 'BackgroundColor', 'white', 'EdgeColor', 'black', 'Margin', 1);
end




% Lets add the IPF figure --------------
mtexFig.nextAxis;
title('IPF Coloring')

% Lopp through grains to plot an IPF figure to show misorientation angle +
% direction relative to mean grain 
for z = 1:length(grains)
    
    % Get EBSD data for the current grain
    currentGrainIdsz = grains(z).id;  % IDs of pixels in the current grain
     % Combine grain ID with filename
    namez = [num2str(currentGrainIdsz) '_' filename];
    currentGrainEBSDz = ebsd_Fo(ebsd_Fo.grainId == currentGrainIdsz);  % Subset EBSD data for current grain
    oM = ipfHSVKey(currentGrainEBSDz);
    oM.inversePoleFigureDirection = orientation(grains(z).meanOrientation) * oM.whiteCenter;
    oM.colorPostRotation = reflection(yvector);
    oM.maxAngle = 2*degree;
    color=oM.orientation2color(currentGrainEBSDz.orientations);
    hold on
    [~,mP]=plot(currentGrainEBSDz, color);


    
    
end

% Check if the folder 'mis2mean_unfilt' exists in the specified path; if not, create it

if ~exist(folder, 'dir')
    mkdir(folder);
end

savePath = fullfile(folder, [strrep(filename, '.ctf', ''), '_mis2mean.png']);
saveas(gcf, savePath)

close all 

%-------------------------------------------------------------------------
% Start afresh with a new figure! This is to show each grain, and its
% misorientation signature 
%------------------------------------------------------------------------
mtexFig = newMtexFigure('layout',[2, length(grains)]); % 2 row, N columns
figure(mtexFig.parent);  % Ensure we are modifying the correct figure
set(mtexFig.parent, 'Position', [100, 100, 1200, 1200]);  % Adjust these values as needed
hold on


% Loop through each grain again for IPF coloring. 
for k = 1:length(grains)
    x=1; % x value for plotting
    y=k; % column based on grain
    
    % Get EBSD data for the current grain
    currentGrainIds = grains(k).id;  % IDs of pixels in the current grain
     % Combine grain ID with filename
    name = [num2str(currentGrainIds) '_' filename];
    GOS = grains(k).GOS ./ degree;



    currentGrainEBSD = ebsd_Fo(ebsd_Fo.grainId == currentGrainIds);  % Subset EBSD data for current grain
    oM = ipfHSVKey(currentGrainEBSD);
    oM.inversePoleFigureDirection = grains(k).meanOrientation * oM.whiteCenter; % Previously was orientation(grains(k).meanOrientation) * oM.whiteCenter;
    oM.maxAngle = grain_div_angle*degree;
    color=oM.orientation2color(currentGrainEBSD.orientations);

    % Lets add the grain number
    
    

    % This is plotting the pretty IPF map. 
    mtexFig.nextAxis(x, y)
    hold on
    [~,mP]=plot(currentGrainEBSD, color);
    centerX = mean(grains(k).boundary.x);
    centerY = mean(grains(k).boundary.y);
    text(centerX, centerY, sprintf('ID: %d\nGOS: %.2f', currentGrainIds , GOS), ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
         'FontSize', 8, 'BackgroundColor', 'white', 'EdgeColor', 'black', 'Margin', 1);
   


    % Lets calculate grain reference orientation deviation (GROD) 
    % the misorientation between each pixel orientation to the grain mean orientation
    ori=currentGrainEBSD.orientations;
    mis2mean=inv(grains(k).meanOrientation) .* ori;
    GOS=currentGrainEBSD.grainMean(mis2mean.angle, currentGrainIds);
    


     gb_C_m=[grains(k).innerBoundary];

    if isempty(gb_C_m)
        disp('Grain boundary is empty!');
        continue; % Skip the current iteration if there is no boundary data
    end


   

%----------------------- First, lets filter boundaries based on
%misorientation angles ------------------------------------------%
         
% 
    z = subgrain_angle; % Plots misorientation axes that have angles greater than this
    ang = 10; % Plots misorientation axes that have angles less than this
% 
%     % Conditions
    cond1 = gb_C_m.misorientation.angle./degree > z;
    cond2 = gb_C_m.misorientation.angle./degree < ang;
    cond3 = gb_C_m.misorientation.angle./degree > ang;

    cond=z < gb_C_m.misorientation.angle./degree & gb_C_m.misorientation.angle./degree < ang;

    condlength=gb_C_m.segLength>len_filter;
    

        if ~any(cond&condlength)
        disp('No boundaries wit selected angles');
        continue; % Skip the current iteration if there is no boundary data
        end


    plot(gb_C_m(condlength&cond), 'linecolor','k','linewidth',1.5) % This is all boundaries. 
    %plot(grains.innerBoundary, 'lineColor', 'y', label='subgrain inner')
    mP.micronBar.visible = 'off';

    % Title for each subplot
    formattedFilename = strrep(filename, '\', '\\');
formattedFilename = strrep(formattedFilename, '_', '-');

% Set the title with specified font size
title(['Grain ', num2str(grains(k).id), ' from ', formattedFilename], 'FontSize', 10);

    
  % Now lets work out how to work out how much of the boundaries are tilt
  % and twist for each grain 

 %Calculate the direction of the misorientation axis in specimin coordinates. 
o=ebsd(gb_C_m.ebsdId).orientations;
try
    ROT=axis(o(:,1),o(:,2));
    %Calculate grain boundary directions in specimin coordinates
    GBD=gb_C_m.direction;


    % Considering boundary dip - shallowly dipping boundaries are unlikely to be observed
    % Calculate vector perpendicuular to misorientation axis and grain boundary direction in specimin coordinates. This is the normal to a tilt boundary
    tiltnormal=cross(ROT, GBD);
    % Define vertical vector in specimin coordinates
    z=vector3d(0, 0, 1);
    % Condition tilt 1: Only consider tilt boundaries that dip steeper than 15 degree wrt. the plane of section.
    condtilt1=angle(tiltnormal,z, 'antipodal')>15*degree;
    % Condition twist 1 : only consider twist boundaries that dip steeper than 15 degree wrt. the plane of section.
    condtwist1=angle(ROT,z, 'antipodal')>15*degree;
    % Condition twist 2 : Ideal twist boundary has boundary strike perpendicular to rotation axis. Here we allow 15 degrees leniancy.
    condtwist2=angle(ROT, GBD, 'antipodal')>75*degree;
    % If a boundary meets both twist criteria, classify as a twist boundary. Include angular, noise and length critera.
    TwistLog=condtwist1&condtwist2;
    TwistBoundaries=TwistLog&cond1&cond2&condlength;
    % Classify all boundaries that are not twist as tilt boundaries, Include angular, noise and length critera.
    TiltLog=~TwistLog;
    TiltBoundaries=TiltLog&cond1&cond2&condlength&condtilt1;

    plot(gb_C_m(TiltBoundaries), 'linecolor','red','linewidth',2)
    plot(gb_C_m(TwistBoundaries), 'linecolor','y','linewidth',2)
    TiltLength=sum(gb_C_m(TiltBoundaries).segLength);
    TwistLength=sum(gb_C_m(TwistBoundaries).segLength);



    Tot_length=sum(gb_C_m.segLength);
    PercTilt=100*TiltLength/Tot_length;
    PercTwist=100*TwistLength/Tot_length;
    PercUnc=100*(Tot_length-TiltLength-TwistLength)/Tot_length;
    AvRotAxis=mean(gb_C_m.misorientation(cond&condlength));
    prop_length=Tot_length/((grains(k).grainSize)^0.5);

    results = [results; {filename, grains(k).id, grains(k).GOS./degree, grains(k).grainSize, ...
        grains(k).area, grains(k).subBoundaryLength, grains(k).equivalentRadius, prop_length, Tot_length, TiltLength, TwistLength, PercTilt, PercTwist, PercUnc }];


    % Plot data on next axis
    x=2; %update to move onto the next figure
    mtexFig.nextAxis(x, y); % Now move onto the next axis - to see what you get IPF wise
    hold on
    plotAxisDistribution(gb_C_m.misorientation(cond&condlength), 'contourf');

    % Lets save these axes
    folder_mis = fullfile(pname, 'Mis_axes');
    if ~exist(folder_mis, 'dir')
        mkdir(folder_mis);
    end

    %if (grains(k).GOS./degree >0.25)

    savename=sprintf('%s_grain%d_misorientation', filename, grains(k).id);
    save_mis=gb_C_m.misorientation;
    fullFilePath_axis=fullfile(folder_mis, savename);
    save(fullFilePath_axis, 'save_mis')
    %end
    %save('')
catch ME
    warning('An error occurred while computing the rotation axis: %s', ME.message);
    % Continue to the next iteration of the loop or perform any other desired action
    continue; % Assuming this code is within a loop
end





end
% Save the entire figure with all subplots
savePath = fullfile(folder, [filename, '_IPF_All_Grains.png']);
mkdir(folder); % Ensure the folder exists
saveas(gcf, savePath);




close all

varsToKeep = {'files', 'folder', 'pname', 'len_filter', 'CS', 'i', 'subgrain_angle', 'grain_div_angle', 'color_angle', 'results'};
allVars = who;% Get a list of all variables in the workspace


varsToClear = setdiff(allVars, varsToKeep); % Determine which variables to clear

clear(varsToClear{:}); % Clear the selected variables


end
%% 

%%
fullmatrix=fullfile(pname, 'test.xlsx');
headers={'Filename', 'grainID',	'GOS' 'Grain Size (pixels)', 'Grain Size (um2)',...
    'subBoundaryLength (um)', 'equivalentRadius (um)',	'GB Length/Sqrt Size', 'Total GB length',	'Tilt length',	'Twist length',	'Perc Tilt',	'Perc Twist', 'Perc unclassified'};
results_with_head=[headers; results];
writecell(results_with_head, fullmatrix)

% %% Load one as a test
% % Set the folder path
% folder_path = fullfile(folder, 'Mis_axes');
% 
% % List all files in the folder
% file_list = dir(folder_path);
% 
% % Initialize an empty cell array to store the loaded data
% data_combined = [];
% 
% % Loop through each file in the folder
% for i = 1:length(file_list)
%     % Get the current file name
%     current_file = file_list(i).name;
% 
%     % Check if the current file contains 'MLP_43' in its name
%     if contains(current_file, 'MLP_43') && strcmp(current_file(end-3:end), '.mat')
%         % Construct the full path to the .mat file
%         file_path = fullfile(folder_path, current_file);
% 
%         % Load the specific variable from the file
%         data = load(file_path, 'save_mis');
% 
%         % Append the loaded data to the combined array
%         data_combined = [data_combined; data.save_mis];
%     end
% end
% 
% % Plot the combined data
% plotAxisDistribution(data_combined, 'contourf');


%% Now Load them all


% % Set the folder path
% folder_path = [folder '/Mis_axes'];
% 
% % Get a list of all .mat files in the folder
% file_list = dir(fullfile(folder_path, '*.mat'));
% 
% % Initialize a cell array to store the loaded data
% data_combined = [];
% 
% % Loop through each file and load the data
% for i = 1:length(file_list)
%     file_path = fullfile(folder_path, file_list(i).name);
%     data = load(file_path, 'save_mis'); % Load the specific variable from each file
%     data_combined = [data_combined; data.save_mis]; % Append to the combined variable
% end
% 
% % Plot the combined data
% plotAxisDistribution(data_combined, 'contourf');


  