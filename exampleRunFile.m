%% Example run file for the IDIC
% VARIABLES OPTIONS
% -------------------------------------------------------------------------
%   imagesFolder: subdirectory containing the series of images on which to run
%                 FIDIC. Images are read in alphanumeric order, which is
%                 assumed to be consistent with time steps.
%   Ext: the file extension of the input images.  All images must be of the
%        same type.
%   numImage: optional parameter to pass to img2mat limiting the number of
%             images processed from the subdirectory identified in imagesFolder
%   sSize: interrogation window (subset) size for the first iterations.
%          Must be 32,64,96, or 128 pixels and a two column
%          array (one for each dimenision) or scalar (equal for all
%          dimensions).
%   sSizeMin: Minimum subset size for the final iteration.
%             Must be 32,64,96, or 128 pixels scalar (equal for all
%             dimensions).
%   runMode:  string that defines the method of running IDIC. Options:
%             cumulative (time0 -> time1, time0 -> time2, ...)
%             (Allowable inputs: 'c','cum','cumulative')
%             or
%             incremental (time0 -> time1, time1 -> time2, ...)
%             (Allowable inputs: 'i','inc','incremental')
%             or
%             hybrid (q-factor based reference updating)
%             (Allowable inputs: 'h','hyb','hybrid')
%
% OUTPUTS
% -------------------------------------------------------------------------
%   u:  displacement field vector calculated from FIDIC. Format: cell array,
%      which is a 3D vector (components in x,y)  per each time point
%      (units are in pixels)
%         u{time}{1} = displacement in x-direction at t=time of size MxNxP
%         u{time}{2} = displacement in y-direction at t=time of size MxNxP
%         u{time}{3} = displacement magnitude at t=time of size MxNxP
%   cc: typically large (100Mb+) struct containing runtime diagnostics from
%       the q-factor stages
%   dm: final subset spacing in px
%   gridPoints: final measurement point meshgrid
%
% NOTES
% Please cite
% Landauer, A.K., Patel, M., Henann, D.L. et al. Exp Mech (2018).
% https://doi.org/10.1007/s11340-018-0377-4
% -------------------------------------------------------------------------
%% Set up workspace and images

clear; close all; clc;

%set up runtime variables
sSize = [64 64];
sSizeMin = 16;
runMode = 'c'; %use 'i' for incremental mode, 'c' for cumulative, and 'h' for hybrid
ext_in = 'tif'; %Input image format
folder_in = ['.',filesep,'example_images']; %Folder containing the images series
max_def_idx = 'b'; %Specify where the max deformation occurs
                    %use 'center' or 'c' for the center image,
                    %'end' or 'e' for the last image,
                    %'beginning' or 'b' for the first,
                    %or specific with an integer
numImages = 3; %use only first n images in the folder, 'all' for all images
spacing = 1; %spacing between images to using in input stack
              %Convert input images to .mat
%Compute basic noise floor and measurement resultion metrics
[noise_percent,meas_res,CI_disp_mean,no_im] = image_eval(folder_in,ext_in,sSize,sSizeMin);

%Optionally crop images.  Set crop to 'y' or 'yes' to enable cropping,
%or enter a set of coordinates [x_topleft,y_topleft,x_bottomright,y_bottomright].
crop = 'yes';
% crop = [1,1,511,511];
[crop_nw_loc,folder_out,ext_crp] = imageCropping(folder_in,ext_in,numImages,spacing,max_def_idx,crop);

resultsFolder = ['.',filesep,'Results',filesep];


filter_yes_no = 'yes';
filt_opt = {'gaussian',[3,3],0.5};

mat_file_save_yes_no = 'no'; % Whether or not to save each img or keep in RAM.
                      %For large images or a large number of image, save
                      %them to disk ("yes") if you run out of RAM, otherwise
                      %keeping in RAM ("no") is a faster options.

if strcmp(mat_file_save_yes_no(1),'n')
    %no mat files saved
    if ischar(numImages(1))
        [cellIMG,filename] = img2mat(folder_out,mat_file_save_yes_no,...
            ext_crp,filter_yes_no,filt_opt,spacing); %All images in "imagesFolder"
    else
        [cellIMG,filename] = img2mat(folder_out,mat_file_save_yes_no,ext_crp,...
            filter_yes_no,filt_opt,spacing,numImages); %Images 1 to numImages
    end
else
    %with mat files saved
    if ischar(numImages(1))
        [~,filename] = img2mat(folder_out,mat_file_save_yes_no,...
            ext_crp,filter_yes_no,filt_opt,spacing); %All images in "imagesFolder"
    else
        [~,filename] = img2mat(folder_out,mat_file_save_yes_no,ext_crp,...
            filter_yes_no,filt_opt,spacing,numImages); %Images 1 to numImages
    end
end

%% RUNNING DIC

%set up parallel pool - used in the "DIC.m" function, but
pool = gcp('nocreate');
if isempty(pool)
    curCluster = parcluster();
    pool = parpool(curCluster.NumWorkers);
end

% Estimate displacements via IDIC
if strcmp(mat_file_save_yes_no(1),'n')
[u, cc, dm, gridPoints, tSwitch] = ...
    funIDIC(cellIMG, sSize, sSizeMin, runMode); %pass "cellIMG" if mat_file_save is
                                                 %"no" or "filename" if
                                                 %"yes" (first argument)
else
[u, cc, dm, gridPoints, tSwitch] = funIDIC(filename, sSize, sSizeMin, runMode);
end

% Save the results
if exist(resultsFolder,'dir') ~= 7
    mkdir(resultsFolder)
end

%% Post-processing and reporting
if no_im == 0
    %Build the reporting table struct array
    prefilt_str = strcat(filt_opt{1},', ',num2str(filt_opt{2}),', ',num2str(filt_opt{3}));
    reporting_table = struct('cameraNoise',noise_percent,'prefiltering',prefilt_str,...
        'initialSubset',sSize,'finalSubset',sSizeMin,'step',dm,'xcorrType','normalized','interpolant','spline',...
        'numMeasurementPts',numel(u{1}{1}),'totalImages',length(u)+1,...
        'displacementSpatialRes',mean(sSize),'displacementResX',meas_res(1),...
        'displacementResY',meas_res(2),'RunMode',runMode);

    fprintf('\n-----------------------------------------\n');
    fprintf('Run Parameters and Measurement Specifications\n')
    fprintf('-----------------------------------------\n');
    fprintf('Camera Noise \t\t  %0.2g%%\n',noise_percent);
    fprintf('Prefiltering \t\t  %s, %0.2gx%0.2g, %0.2g\n',filt_opt{1},...
        filt_opt{2}(1),filt_opt{2}(2),filt_opt{3});
    fprintf('Initial subset        %ipx by %ipx\n',sSize(1),sSize(2));
    fprintf('Final (square) subset %ipx\n',sSizeMin);
    fprintf('Step         \t\t  %0.2gpx\n',dm);
    fprintf('Run mode     \t\t  %s\n',runMode);
    fprintf('Correlation type      %s\n','normalized');
    fprintf('Interpolation \t\t  Spline\n');
    fprintf('Measurement points \t  %i\n',numel(u{1}{1}));
    fprintf('Total images \t\t  %i\n',length(u)+1);
    fprintf('Displacement\n   Spatial resolution \t ~%0.4gpx \n   ',mean(sSizeMin/2));
    fprintf('Measurement res, x    %0.2gpx\n   ',meas_res(1));
    fprintf('Measurement res, y    %0.2gpx\n',meas_res(2));
    fprintf('-----------------------------------------\n');

    %Save relavent workspace variables
    save(strcat(resultsFolder,'results_qDIC.mat'),'u','cc','dm','gridPoints','reporting_table');
else
    %Save relavent workspace variables
    save(strcat(resultsFolder,'results_qDIC.mat'),'u','cc','dm','gridPoints');
end


%% CLEAN UP
%Clean up the current set of images from the cd
if strcmp(mat_file_save_yes_no(1),'y')
    delete *IDIC_image*.mat
    delete(strcat(folder_out,'*.',ext_crp));
end
