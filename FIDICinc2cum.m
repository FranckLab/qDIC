function [uCum] = FIDICinc2cum(u,tSwitch,dm,m,option)
% Converting incremental displacement to cumulative displacement
%
%   Inputs:
%   -----------------------------------------------------------------------
%     u: It is cumulative displacement field which may or may not have
%     change of reference image.
%
%     tSwitch: The time point position where the reference image for
%     computating cumulative displacment is updated to bad correlations.
%
%     dm:  Spacing between grid points in 'u.'
%
%     m: meshgrid from DIC
%
%     option: Interpolation method
%
%   Outputs:
%   -----------------------------------------------------------------------
%     uCum: Updated cumulative field with same reference image.

% Find first timepoint where reference image is updated
idx = find(tSwitch>0,1,'first');
uCum = u(1:idx-1);

if isempty(idx)
    uCum = u;
else
    
    % Loop through all timepoints after first timepoint where reference image
    % is updated
    for i = idx:length(u)
        
        %uTemp stores displacement in a format to utilize converstion between
        %incremental to cumulative displacment.
        %uTemp{1} stores value of cumulative displacement till the last
        %timepoint where the reference image is updated.
        %uTemp{2}: Cumulative displacement calculated at current timepoint wrt
        %updated reference image
        
        if tSwitch(i) >= 1
            if i == 1
                uTemp{1} = cell(1,2);
                for ii = 1:2
                    uTemp{1}{ii} = zeros(size(u{1}{ii}));
                end
            else
                uTemp{1} = uCum{i-1};
            end
        end
        
        uTemp{2} = u{i};
        [uTemp] = inc2cum(uTemp,dm,m,option);
        uCum{i} = uTemp{2};
    end
end

end

% % close all
% %
% % % for ii = 1:8
% % %     figure,imagesc(isnan(cc{ii}.qfactors_accept{1}(1+trim_size(1):end-trim_size(1),...
% % %     1+trim_size(2):end-trim_size(2)))),colorbar, axis image
% % % figure,histogram(cc{ii}.qfactors(3,:),200)
% % % figure,histogram(cc{ii}.qfactors(4,:),200)
% % % pause
% % % end