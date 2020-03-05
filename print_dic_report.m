function [] = print_dic_report(reporting_table)
%helper function for qDIC to print out a "DIC report table" of from the 
%reporting_table variable built from the output of image_eval.m, and 
%output when an image series with static images is run
%
%
% VARIABLES OPTIONS
% -------------------------------------------------------------------------
%   reporting_table: the reporting_table variable save with a qDIC run that
%                    has static images
%
% OUTPUTS
% -------------------------------------------------------------------------
%   none

fprintf('\n-----------------------------------------\n');
fprintf('Run Parameters and Measurement Specifications\n')
fprintf('-----------------------------------------\n');
fprintf('Camera Noise \t\t  %0.2g%%\n',reporting_table.cameraNoise);
fprintf('Prefiltering \t\t  %s\n',reporting_table.prefiltering);
fprintf('Image size            %ipx by %ipx\n',reporting_table.imageSize(1),reporting_table.imageSize(2));
fprintf('Initial subset        %ipx by %ipx\n',reporting_table.initialSubset(1),reporting_table.initialSubset(2));
fprintf('Final (square) subset %ipx\n',reporting_table.finalSubset);
fprintf('Step         \t\t  %0.2gpx\n',reporting_table.step);
fprintf('Run mode     \t\t  %s\n',reporting_table.RunMode);
fprintf('Correlation type      %s\n','normalized');
fprintf('Interpolation \t\t  Spline\n');
fprintf('Measurement points \t  %i\n',reporting_table.numMeasurementPts);
fprintf('Total images \t\t  %i\n',reporting_table.totalImages);
fprintf('Displacement\n   Spatial resolution \t ~%0.4gpx \n   ',reporting_table.finalSubset/2);
fprintf('Measurement res, x    %0.2gpx\n   ',reporting_table.displacementResX);
fprintf('Measurement res, y    %0.2gpx\n',reporting_table.displacementResY);
fprintf('-----------------------------------------\n');



