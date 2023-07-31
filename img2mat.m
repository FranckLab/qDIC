function [cellIMG,filename,filt_opt] = img2mat(folder_in,mat_file_save,ext_in,smoothing,filt_opt,freq_img,s)
%Read images and write them out in .mat
%
% INPUTS
% -------------------------------------------------------------------------
%   folder_in: folder containing orginal images
%   mat_file_save: 'yes' or 'no' to save out imags to .mat
%   ext_in: image formate extention
%   smoothing: on/off for Gaussian smoothing prefilter
%   filt_opt: Gaussian smoothing prefilter options (default: [3,3],0.5)
%   freq_img: frequency to sample the image stack
%   s: max number of images to use
%
% OUTPUTS
% -------------------------------------------------------------------------
%   cellIMG: image stack for all current images
%   filename: regex for filename prefix
%   filt_opt: filter options used if the smoothing param was set to 'yes'
%
% NOTES
% -------------------------------------------------------------------------
%

%Load all of the files directory information
files = dir(strcat(folder_in,filesep,'*',ext_in));

loc = strcat(folder_in,filesep,'*',ext_in);
if isempty(files)
    fprintf('NO IMAGE FILES FOUND UNDER: %s \n',loc)
end

%Determine the number of files
if nargin<7
    s = length(files);
end

if strcmp(smoothing,'yes')
    
    filter_gauss = fspecial(filt_opt{1},filt_opt{2},filt_opt{3});
    
    % Loop through files, reading in alpha-numeric order
    cnt = 0;
    for ii = 1:freq_img:s
        cnt = cnt+1;
        READ = imread(strcat(folder_in,filesep,files(ii).name));
        %store the image, and do a small amount of gaussian blurring to
        %improve contrast gradients
        cellIMG{cnt} = imfilter(double(READ(:,:,1)),filter_gauss,'replicate');
        
        filename = 'images not saved to disk';
        % Option to plot the images
        %         imshow(IMG(:,:,ii))
        %         drawnow
    end
    
else
    cnt = 0;
    filt_opt = {'none',[nan,nan],nan};
    % Loop through files, reading in alpha-numeric order
    for ii = 1:freq_img:s
        cnt = cnt+1;
        READ = imread(strcat(folder_in,filesep,files(ii).name));
        %store the image, and do a small amount of gaussian blurring to
        %improve contrast gradients
        cellIMG{cnt} = double(READ(:,:,1));
        
        filename = 'images not saved to disk';
        % Option to plot the images
        %         imshow(IMG(:,:,ii))
        %         drawnow
    end
    
end
%     cellIMG = cell(1);

% only save out if requested
if strcmp(mat_file_save(1),'y')
    for ii = 1:length(cellIMG)
        IMG{1} = cellIMG{ii}; %Make a new variable to hold the current img
        
        %image, needed for "save" to work properly
        filename = strcat(files(ii).name(1:end-(length(ext_in)+1)),...
            '_IDIC_image_',num2str(ii+9999),'.mat');
        save(filename,'IMG');
    end
    cellIMG = 'images saved to disk';
    filename = '*IDIC_image*';
else
    filename = files(1).name;
end



end
