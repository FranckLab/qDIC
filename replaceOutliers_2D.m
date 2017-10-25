function [u,alpha_mask,nan_mask,edge_pts] = replaceOutliers_2D(u,cc)
%function that takes in the disp field, and a cc struct with fields for
%flagged bad points and overwrites poor values with information from nearby
%nodes


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

% nan_mask.U951 = ones(size(cc.U95_accept{1}));
% nan_mask.U952 = ones(size(cc.U95_accept{2}));

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


% nan_mask.qf1(~isfinite(cc.qfactors_accept{1})) = nan;
% nan_mask.qf2(~isfinite(cc.qfactors_accept{2})) = nan;

% nan_mask.U951(~isfinite(cc.U95_accept{1})) = nan;
% nan_mask.U952(~isfinite(cc.U95_accept{2})) = nan;

for i = 1:length(u)
    % Union nan mask
%     union_mask_u{i} = nan_mask.qf1.*nan_mask.qf2.*nan_mask.U951...
%         .*nan_mask.U952.*nan_mask.u{i};
    
    union_mask_u{i} = nan_mask.qf1.*nan_mask.qf2.*nan_mask.u{i};
    
    % apply nan mask
    u{i} = u{i}.*union_mask_u{i};
    
    % fill in nan'd values
    u{i} = inpaint_nans(double(u{i}),inpaint_opt);
    
end



