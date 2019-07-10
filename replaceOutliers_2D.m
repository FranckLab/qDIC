function [u,alpha_mask,nan_mask,edge_pts] = replaceOutliers_2D(u,cc)
%function that takes in the disp field, and a cc struct with fields for
%flagged bad points and overwrites poor values with information from nearby
%nodes
%
% u = flagOutliers(u,cc,thr,epsilon) removes outliers identified in
% previous steps
%
% INPUTS
% -------------------------------------------------------------------------
%   u: cell containing the input displacement field. (u{1:3} = {u_x, u_y,
%   	u_mag})
%   cc: cross-correlation diagnostic struct
%
% OUTPUTS
% -------------------------------------------------------------------------
%   u: displacement field with outliers removed and interpolated
%   alpha_mask: the mask of points removed, as an aphla-bending format for
%               plotting
%   nan_mask: the mask of points removed, as nans
%   edge_pts: edge points of the image that are nan'ed because they lack
%             support for interpolation
%
% NOTES
% -------------------------------------------------------------------------
%Please cite:
% Landauer, A.K., Patel, M., Henann, D.L. et al. Exp Mech (2018).
% https://doi.org/10.1007/s11340-018-0377-4


inpaint_opt = 0;

%convert to alpha masks for plotting
%zero == pass
alpha_mask.qf1 = zeros(size(cc.qfactors_accept{1}));
alpha_mask.qf2 = zeros(size(cc.qfactors_accept{2}));

for i = 1:length(u)
    alpha_mask.u{i} = zeros(size(cc.u_outlier{i}));
    alpha_mask.u{i}(~isfinite(cc.u_outlier{i})) = 0.4;
end

%0.4 == reject
alpha_mask.qf1(~isfinite(cc.qfactors_accept{1})) = 0.4;
alpha_mask.qf2(~isfinite(cc.qfactors_accept{2})) = 0.4;

%convert to NaN masks for processing
%one == pass
nan_mask.qf1 = ones(size(cc.qfactors_accept{1}));
nan_mask.qf2 = ones(size(cc.qfactors_accept{2}));

for i = 1:length(u)
    nan_mask.u{i} = ones(size(cc.u_outlier{i}));
    nan_mask.u{i}(~isfinite(cc.u_outlier{i})) = nan;
end

%nan == reject
nan_mask.qf1(~isfinite(cc.qfactors_accept{1})) = nan;
nan_mask.qf2(~isfinite(cc.qfactors_accept{2})) = nan;

edgesq = ones(size(nan_mask.qf2));
edge_pts = ones(size(nan_mask.qf2));
edgesq(4:end-3,4:end-3) = 0;

edges1 = ones(size(nan_mask.qf2));
edges1(2:end-1,2:end-1) = 0;
edges1(edges1==1) = nan;
edges1 = edges1+1;

edge_pts(edgesq==1) = nan_mask.qf1(edgesq==1);
edge_pts(edgesq==1) = edge_pts(edgesq==1).*nan_mask.qf1(edgesq==1).*edges1(edgesq==1);


nan_mask.qf1(edgesq==1) = 1;
nan_mask.qf2(edgesq==1) = 1;

for i = 1:length(u)
    % Union nan mask
    union_mask_u{i} = nan_mask.qf1.*nan_mask.qf2.*nan_mask.u{i};

    % apply nan mask
    u{i} = u{i}.*union_mask_u{i};

    % fill in nan'd values
    u{i} = inpaint_nans(double(u{i}),inpaint_opt);

end
