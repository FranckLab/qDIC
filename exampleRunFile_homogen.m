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
%   incORcum: string that defines the method of running IDIC. Options:
%             cumulative (time0 -> time1, time0 -> time2, ...)
%             (Allowable inputs: 'c','cum','cumulative')
%             or
%             incremental (time0 -> time1, time1 -> time2, ...)
%             (Allowable inputs: 'i','inc','incremental')
%
% OUTPUTS
% -------------------------------------------------------------------------
%   u:  displacement field vector calculated from FIDIC. Format: cell array,
%      which is a 3D vector (components in x,y)  per each time point
%      (units are in pixels)
%         u{time}{1} = displacement in x-direction at t=time of size MxNxP
%         u{time}{2} = displacement in y-direction at t=time of size MxNxP
%         u{time}{3} = displacement magnitude at t=time of size MxNxP
%   cc: peak values of the cross-correlation for each interrogation
%   dm: final subset spacing in px
%
% NOTES
% -------------------------------------------------------------------------
%% Set up workspace and images

%Image cropping to get a region of interest on the specimen
% [crop_nw_loc,folder_out] = imageCropping(folder_in,ext_in,sSize,max_def_idx);

% clear; close all; clc;

function [] = exampleRunFile_homogen(sSizeMin)

if nargin<1
	sSizeMin = 32;
end

sSize = [64 64];
% sSizeMin = 32;
runMode = 'i'; %use 'i' for incremental mode, 'c' for cumulative, and 'h' for hybrid
folder_out = 'full_image_series_comp';
ext_crp = 'tif'; %output image file form, defined in image_cropping.m
% resultsFolder = ['.',filesep,'Results_hybrid_homDef',filesep];
% numImages = 15;

%Convert input images to .mat
[cellIMG,filename,filt_opt] = img2mat(folder_out,ext_crp,'off'); %All images in "imagesFolder"
% [cellIMG,filename] = img2mat(folder_out,ext_crp,'on',numImages); %Images 1 to
%numImages only
filename = '*IDIC_image*.mat';

%% RUNNING DIC

pool = gcp('nocreate');
if isempty(pool)
    curCluster = parcluster();
%     curCluster.NumWorkers = 6;
%     saveProfile(curCluster);
    pool = parpool(curCluster.NumWorkers);
end

% Estimate displacements via IDIC
[u, ~, dm, gridPoints] = funIDIC(filename, sSize, sSizeMin, runMode);

save(['testing_workspace_homogen_9iter2_ssmin',num2str(sSizeMin)])          

%%
clearvars -except sSizeMin 

for ii = 1:3
    if ii == 1
        sSizeMin = 16;
    elseif ii == 2
        sSizeMin = 32;
    else
        sSizeMin = 64;
    end

load(['testing_workspace_homogen_9iter2_ssmin',num2str(sSizeMin)])
numImages = length(u)+1;
% clear cc 

% u_new = u;
% qfactor_U95_test_plotting(u,u_new,cc,dm);

        % use equilibrium smoothing to remove some non-physical motions (see
        % Soonsung Hong, Huck Beng Chew, Kyung-Suk Kim, Cohesive-zone laws
        % for void growth � I. Experimental field projection of crack-tip
        % crazing in glassy polymers, JMPS, 57,8, 2009, Pp 1357-73, %
        % http://dx.doi.org/10.1016/j.jmps.2009.04.003.)

% [u_hybrid,u_imp_hybrid,u_err_hybrid,u_err_2norm_hybrid,...
%     nan_mask_hybrid] = error_mapping_genDef(u,dm);

% save('.././results_hybrid_genDef_16pxRefine')
if runMode(1) == 'i'
    u = inc2cum(u,dm,gridPoints,'spline');
end

bndry_wd = 0;
[u_masked_hyb,u_imp_masked_hyb,u_err_hyb,l2err_hyb,std_err_hyb,mean_err_hyb,mean_err_mag_hyb,imp_strain_hyb,nan_mask]...
    = error_mapping_homogen(u,dm,gridPoints,numImages,bndry_wd);

clear cc
save(['.././results_hybrid_homogen_9iter2_ssmin',num2str(sSizeMin)])

% [u_masked_hyb,u_imp_masked_hyb,u_err_hyb,l2err_hyb,std_err_hyb,mean_err_hyb,imp_strain_hyb,nan_mask]...
%     = error_mapping_genDef2(u,dm,numImages,bndry_wd);
% run .././error_comparison_qFIDIC_FIDIC_zeroDisp
end

% delete *IDIC_image*.mat
% % Save the results
% if exist(resultsFolder,'dir') ~= 7
%     mkdir(resultsFolder)
% end

end

%Build the reporting table struct array
% prefilt_str = strcat(filt_opt{1},', ',num2str(filt_opt{2}),', ',num2str(filt_opt{3}));
% reporting_table = struct('cameraNoise',noise_percent,'prefiltering',prefilt_str,...
%     'subset',sSize,'step',dm,'xcorrType',norm_xcc,'interpolent','spline',...
%     'numMeasurementPts',numel(u{1}{1}),'totalImages',length(u)+1,...
%     'displacementSpatialRes',mean(sSize),'displacementResX',meas_res(1),...
%     'displacementResY',meas_res(2));

%Save relavent workspace variables
% save(strcat(resultsFolder,'resultsFIDIC.mat'),'u','cc','cellIMG','dm','reporting_table');

%% CLEAN UP
%Clean up the current set of images from the cd
% delete *IDIC_image*.mat
% delete(strcat(folder_out,'*.',ext_crp));
