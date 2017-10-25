function [IMG,filename,filt_opt] = img2mat(folder_in,ext_in,smoothing,s)
%Read images and write them out in .mat
%
% INPUTS
% -------------------------------------------------------------------------
%   folder_in: folder containing orginal images
%   ext_in: image formate extention
%   folder_in: folder containing orginal images
%   smoothing: on/off for Gaussian smoothing prefilter ([3,3],0.5)
%   s: max number of images to use
%
% OUTPUTS
% -------------------------------------------------------------------------
%   IMG: image stack for all current images
%   filename: regex for filename prefix
%   filt_opt: filter options used if the smoothing param was set to 'yes'
%
% NOTES
% -------------------------------------------------------------------------
%

%Load all of the files directory information
files = dir(strcat('.',filesep,folder_in,filesep,'*',ext_in));

%Determine the number of files
if nargin<4
    s = length(files);
end

if strcmp(smoothing,'on')
    filt_opt = {'gaussian',[3,3],0.5};
    
    filter_gauss = fspecial(filt_opt{1},filt_opt{2},filt_opt{3});
    
    % Loop through files, reading in alpha-numeric order
    for ii = 1:s
        READ = imread(strcat(folder_in,filesep,files(ii).name));
        %store the image, and do a small amount of gaussian blurring to
        %improve contrast gradients
        IMG(:,:,ii) = imfilter(double(READ(:,:,1)),filter_gauss,'replicate');
        
        % Option to plot the images
        %         imshow(IMG(:,:,ii))
        %         drawnow
    end
    
else
    filt_opt = {'none',[nan,nan],nan};
    % Loop through files, reading in alpha-numeric order
    for ii = 1:s
        READ = imread(strcat(folder_in,filesep,files(ii).name));
        %store the image, and do a small amount of gaussian blurring to
        %improve contrast gradients
        IMG(:,:,ii) = double(READ(:,:,1));
        
        % Option to plot the images
        %         imshow(IMG(:,:,ii))
        %         drawnow
    end
    
end
%     cellIMG = cell(1);

%
for ii = 1:s
    cellIMG{1} = IMG(:,:,ii); %Make a new variable to hold the current
    %image, needed for "save" to work properly
    filename = strcat(files(ii).name,'IDIC_image_',num2str(ii+999),'.mat');
    save(filename,'cellIMG');
end

filename = '*IDIC_image*';

end
