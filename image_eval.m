function [noise_percent,spatial_res,CI_disp_mean,no_im] = image_eval(Folder,ext,sSize0,sSizeMin)
%This function performs basic noise floor and spetial resolution analyses
%for images used in the DIC.  To use, take several (2+) completely static
%images of the speckle pattern and label these with the keyword "static" in
%the image directory.
%
% VARIABLES OPTIONS
% -------------------------------------------------------------------------
%   Folder: subdirectory containing the series of images on which to run
%           the evalution, there should be 2 or more images with "static"
%           as part of the filename
%   Ext: the file extension of the input images.  All images must be of the
%        same type.
%
% OUTPUTS
% -------------------------------------------------------------------------
%   noise_percent:  percentage of the full-scale range of the image format
%                   greyscale that the sensor noise occupies
%   spatial_res: the spatial resolution that can be expected from the
%                algorithm, based on the noise, speckle pattern, optics, 
%                and other error sources.
%   CI_disp_mean: mean confidence interval on the displacement measured.
%                 This should be centered on zero, unless bias errors exist
%   no_im: flag indicating that no "static"-labeled images were found
%
% NOTES
% Please cite:
% Landauer, A.K., Patel, M., Henann, D.L. et al. Exp Mech (2018). 
% https://doi.org/10.1007/s11340-018-0377-4

% -------------------------------------------------------------------------

%% Retrieve images

%Load all of the file's directory information
files = dir(strcat(Folder,filesep,'*static*',ext));
l = length(files);

%Only procede if evaluation images are present
if l == 0
    disp('Image evaluation skipped')
    %set failure flag
    no_im = 1;
    noise_percent = nan;
    spatial_res = nan;
    CI_disp_mean = nan;
else
    disp('Running image evaluation DIC')
    no_im = 0;
    
    %read in the image sequence
    for ii = 1:l
        READ = imread(strcat(Folder,filesep,files(ii).name));
        try 
            full_images(:,:,ii) = double(rgb2gray(READ));
        catch
            full_images(:,:,ii) = double(READ);
        end
    end
    
    %find the bitdepth of the images
    S = whos('READ');
    if strcmp(S.class,'uint8')
        depth = 256;
    elseif strcmp(S.class,'uint16')
        depth = 65536;
    else
        depth = max(full_images(:));
    end
    clear READ
    
    %% Select eval region
    imagesc(full_images(:,:,1))
    title('Select noise evaluation ROI. Define two points: top left and bottom right')
    axis('image'); colormap gray
    [X,Y] = ginput(2);
    X = ceil(X);
    Y = ceil(Y);
    close
    

    images = full_images(Y(1):Y(2),X(1):X(2),:);
    
    %% Noise level
    %mean_image = mean(image,3);
    std_image = std(images,0,3);
    
    noise_level = std(std_image(:));
    
    noise_percent = noise_level/depth*100;
    
    %% Spatial resolution
    
    %set up parameters
    image_pair = cell(1,2);
    image_pair{1} = images(:,:,1);
    image_pair{2} = images(:,:,2);
    
    u0 = cell(1,2);
    u0{1} = 0;
    u0{2} = 0;

    
    %Do itereative DIC between the identical images
    [u,~,~,~] = IDIC(image_pair,sSize0,sSizeMin,u0);
    
    %Compute spatial resolutions
    spatial_res(1) = std(u{1},0,'all','omitnan');
    spatial_res(2) = std(u{2},0,'all','omitnan');
    spatial_res(3) = std(u{3},0,'all','omitnan');
    
    z = 1.96;
    %Compute the confidence interval on the displacements
    CI_disp(:,:,1) = u{3} - z*spatial_res(3)/sqrt(2);
    CI_disp_mean(1) = mean(CI_disp(:,:,1),'all','omitnan');
    CI_disp(:,:,2) = u{3} + z*spatial_res(3)/sqrt(2);
    CI_disp_mean(2) = mean(CI_disp(:,:,2),'all','omitnan');
    
end
