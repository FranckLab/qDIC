function [crop_nw_loc,folder_out,fmt] = imageCropping(folder_in,ext_in,sSize,max_def_idx,crop)
%This function crops the input images to include only the region of
%interest
%
% INPUTS
% -------------------------------------------------------------------------
%   folder_in: folder containing orginal images
%   ext_in: image formate extention
%   folder_in: folder containing orginal images
%   sSize: interrogation window (subset) size
%   max_def_index: string specifying where the max deformation occurs 
%                   use 'center' or 'c' for the center image,
%                   'end' or 'e' for the last image,
%                   'beginning' or 'b' for the first,
%                   or specific with an integer 
%   crop: string to specify whether to crop images, set to 'y' or 'yes' to crop
% OUTPUTS
% -------------------------------------------------------------------------
%   crop_nw_loc: location of the northwest corner of the cropped region
%   folder_out: the location where the images were placed
%
%NOTES
% -------------------------------------------------------------------------
%

%% Setup

%Output variables
fmt = 'tif';
folder_out = strcat(folder_in,filesep,'cropped_images',filesep);
ext_out = strcat('.',fmt);

%Make a new output folder if none exists
if exist(folder_out,'dir') ~= 7
    mkdir(folder_out);
end

prefixes = cell(1,26^3);
alphab = {'a','b','c','d','e','f','g','h','i','j','k','l','m','n','o',...
    'p','q','r','s','t','u','v','w','x','y','z'};

kk = 0;
for hh = 1:length(alphab)
    for ii = 1:length(alphab)
        for jj = 1:length(alphab)
            kk = kk + 1;
            prefixes{kk} = strcat(alphab{hh},alphab{ii},alphab{jj});
        end
    end
end

%% Read in image filenames
files = dir(strcat(folder_in,filesep,'*',ext_in));
l = length(files);

%% Get Cropping Region
%% Get Cropping Region
if strcmp(crop, 'yes')||strcmp(crop, 'y')
    if strcmp(max_def_idx,'center')||strcmp(max_def_idx,'c')
        im_loc = ceil(l/2);
    elseif strcmp(max_def_idx,'end')||strcmp(max_def_idx,'e')
        im_loc = l;
    elseif strcmp(max_def_idx,'beginning')||strcmp(max_def_idx,'b')
        im_loc = 1;
    end
    figure
    imagesc(imread(strcat(folder_in,filesep,files(im_loc).name)))
    title('Click to select cropping region. Define two points: top left and bottom right')
    axis('image'); colormap gray
    [X,Y] = ginput(2);
    X = ceil(X);
    Y = ceil(Y);
    X_ss(1) = X(1); %- mod(X(1),max(sSize))+max(sSize); %place the point such that an
    %interger number of subsets is used
    %Crop agressively.
    X_ss(2) = X(2);% - mod(X(2),max(sSize));
    Y_ss(1) = Y(1);% - mod(Y(1),max(sSize))+max(sSize);
    Y_ss(2) = Y(2);% - mod(Y(2),max(sSize));
    close
    
    crop_nw_loc = [X_ss(1),Y_ss(1)];
    
else
    crop_nw_loc = [1,1];
    X_ss(1) = 1;
    X_ss(2) = size(imread(strcat(folder_in,filesep,files(1).name)),2);
    Y_ss(1) = 1;
    Y_ss(2) = size(imread(strcat(folder_in,filesep,files(1).name)),1);
    
end

%% Crop and write out files

image_idx = 1:l;
% Loop through files
for ii = 1:length(image_idx)
    READ = imread(strcat(folder_in,filesep,files(image_idx(ii)).name));
    try
    READ = rgb2gray(READ);
    catch
    end
    IMG = READ(Y_ss(1):Y_ss(2),X_ss(1):X_ss(2),1); %Cropped size from ginput
    
    dir_filename = strcat(folder_out,prefixes{image_idx(ii)},'_image_number_',...
        num2str(image_idx(ii)),ext_out); %use prefix to ensure proper ordering
    imwrite(IMG,dir_filename,fmt); %Write the file with the specified settings
end

