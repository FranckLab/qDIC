function [u, cc, dm, mFinal, decorrFlag] = IDIC(varargin)
%
% INPUTS
% -------------------------------------------------------------------------
%   I0: cell containing the undeformed, I0{1}, and deformed, I0{2} images
%   sSize: interrogation window (subset) size
%   sSizeMin: interrogation window (subset) minimum size
%   u0: pre-estimated displacement field (typically zeros)
%
% OUTPUTS
% -------------------------------------------------------------------------
%   u: displacement field vector defined at every meshgrid point with
%      spacing dm. Format: cell array, each containing a 3D matrix
%         (components in x,y)
%         u{1} = displacement in x-direction
%         u{2} = displacement in y-direction
%         u{3} = magnitude
%   cc: peak values of the cross-correlation for each interrogation
%   dm: meshgrid spacing (8 by default)
%   m_final: final meshgrid
%   decorr_flag: flag indicating that decorrelation is occuring
%
% NOTES
% -------------------------------------------------------------------------
% To run you may need a compatible C compiler. Please see
% (http://www.mathworks.com/support/compilers/R2014a/index.html)
%
% If used please cite:
% Landauer, A.K., Patel, M., Henann, D.L. et al. Exp Mech (2018). 
% https://doi.org/10.1007/s11340-018-0377-4

% PRESET CONSTANTS
maxIterations = 9; % maximum number of iterations
dm = 8; % desired output mesh spacing
norm_xcc = 'n'; %switch to 'u' for un-normalized xcc - faster but less accurate
convergenceCrit = [0.15, 0.25, 0.075]; % convergence criteria
ccThreshold = 1e-4; % bad cross-correlation threshold (un-normalized)
%[not used in norm version]
sizeThresh = 196;  %Threshold for maximum size of bad correlation regions
percentThresh = 25; %Threshold for max % of total measurement pts failing q-factor testing
stDevThresh = 0.15; %Threshold for max st. dev. of the fitted gaussian to the second peak

cc = cell(1);
cc{1} = struct('A',[],'max',[],'maxIdx',[],'qfactors',[],'q_thresh',[],...
    'qfactors_accept',[]);

[I0, sSize, sSizeMin, sSpacing, padSize, DICPadSize, u] = parseInputs(varargin{:});

% START ITERATING
i = 2; converged01 = 0; SSE = []; I = I0;

t0 = tic;
while ~converged01 && i - 1 < maxIterations
    ti = tic;

    % Check for convergence
    [converged01, SSE(i-1) , sSize(i,:), sSpacing(i,:)] = ...
        checkConvergenceSSD_2D(I,SSE,sSize,sSizeMin,sSpacing,convergenceCrit);
    
    if ~converged01
        
        finalSize = sSize(i,:);
        [I, m] = parseImages(I,sSize(i,:),sSpacing(i,:));
        
        %warp images with displacement guess if a cumulative step
        if i == 2 %only on the first iteration
            if numel(u{1}) == 1 %on the first image the disp guess is zero
                u{1} = zeros(length(m{1}),length(m{2}));
                u{2} = zeros(length(m{1}),length(m{2}));
            end
            u{1} = inpaint_nans(u{1}); %inpaint nans from last time's edge pts
            u{2} = inpaint_nans(u{2});
            I = areaMapping_2D(I,m,u); %otherwise map the images w/ the initial guess
            [I, m] = parseImages(I,sSize(i,:),sSpacing(i,:));
%             u{1} = 0; %reset disp to zero
%             u{2} = 0;
        end
        
        % run cross-correlation to get an estimate of the displacements
        [du, cc{i-1}] = DIC(I,sSize,sSpacing(i,:),DICPadSize,ccThreshold,norm_xcc);

        % add the displacements from previous iteration to current
        [u, ~, cc{i-1}, mFinal] = addDisplacements_2D(u,du,cc{i-1},cc{i-1},m,dm);
        
        % filter the displacements using a predictor filter
        u = filterDisplacements_2D(u,sSize(i,:)/dm);
        cc{i-1} = flagOutliers_2D(u,cc{i-1},1.25,0.1);
        
        % remove questionable displacement points
        [u,cc{i-1}.alpha_mask,cc{i-1}.nan_mask,cc{i-1}.edge_pts] = replaceOutliers_2D(u,cc{i-1});
        
        % mesh and pad images based on new subset size and spacing
        [I, m] = parseImages(I0,sSize(i,:),sSpacing(i,:));
        
        % map areas based on displacment field
        I = areaMapping_2D(I,m,u);
        
        if sum(isnan(u{1}(:)+isnan(u{2}(:)))) > 0 || sum(isnan(I{1}(:)+isnan(I{2}(:)))) > 0
            disp('nan found')
        end
        
        disp(['Elapsed time (iteration ',num2str(i-1),'): ',num2str(toc(ti))]);
        i = i + 1;
    end
    
end

%prepare outputs
decorrFlag = decorrelationCheck(cc,sizeThresh,percentThresh,stDevThresh);
u{1} = cc{end}.edge_pts.*u{1};
u{2} = cc{end}.edge_pts.*u{2};
[u,cc,mFinal] = parseOutputs(u,cc,finalSize,padSize,mFinal);

disp(['Convergence at iteration ',num2str(i-1)]);
disp(['Total time: ',num2str(toc(t0))]);
end



%% ========================================================================
function varargout = parseImages(varargin)
% pads images and creates meshgrid

I{1} = single(varargin{1}{1});
I{2} = single(varargin{1}{2});
sSize = varargin{2};
sSpacing = varargin{3};

prePad = sSize/2;
postPad = sSize/2;

sizeI = size(I{1});
I{1} = padarray(I{1},prePad,'replicate','pre');
I{1} = padarray(I{1},postPad,'replicate','post');

I{2} = padarray(I{2},prePad,'replicate','pre');
I{2} = padarray(I{2},postPad,'replicate','post');

idx = cell(1,2);
for i = 1:2, idx{i} = (1:sSpacing(i):(sizeI(i) + 1)) + sSize(i)/2; end

% [m{1},m{2}] = meshgrid(idx{:});

varargout{    1} = I;
varargout{end+1} = idx;

end

%% ========================================================================
function varargout = parseInputs(varargin)
% parses inputs and pads images so that there is an divisable meshgrid number.

I0{1} = single(varargin{1}{1});
I0{2} = single(varargin{1}{2});

% I0{1} = permute(I0{1},[2 1]);
% I0{2} = permute(I0{2},[2 1]);

sSize = varargin{2};
sSize = [sSize(2), sSize(1)];

sSizeMin = varargin{3};

sSpacing = sSize/2;
u0_ = varargin{4};
u0 = u0_(1:2);

DICPadSize = sSpacing/2;

sizeI0 = size(I0{1});
sizeI = ceil(sizeI0./sSpacing).*sSpacing;
prePad = ceil((sizeI - sizeI0)/2);
postPad = floor((sizeI - sizeI0)/2);

I{1} = padarray(I0{1},prePad,'replicate','pre');
I{1} = padarray(I{1},postPad,'replicate','post');

I{2} = padarray(I0{2},prePad,'replicate','pre');
I{2} = padarray(I{2},postPad,'replicate','post');

varargout{    1} = I;
varargout{end+1} = sSize;
varargout{end+1} = sSizeMin;
varargout{end+1} = sSpacing;
varargout{end+1} = [prePad; postPad];
varargout{end+1} = DICPadSize;
varargout{end+1} = u0;

end


function [u,cc,m] = parseOutputs(u,cc,sSize,padSize,m_)
% parses outputs and unpads the displacment field and cc.

for i = 1:2
    m{i} = m_{i}-padSize(1,i)-sSize(i)/2;
end
u{3} = sqrt(u{1}.^2 + u{2}.^2);

%remove the memory-intesive and generally unneeded parts of the cc struct
for jj = 1:length(cc)
    cc{jj}.A = [];
    cc{jj}.max = [];
end

end

function decorrFlag = decorrelationCheck(cc,sizeThresh,percentThresh,stDevThresh)
%This function takes in a cc structure and checks the q-factor maps for
%decorrelation indications.

%Initilize
decorrFlag = 0;

if nargin < 4
    stDevThresh = 0.13;
elseif nargin < 3
    percentThresh = 15;
elseif nargin < 2
    sizeThresh = 121;
end

trim_size = ceil((cc{end}.sSize/2)./cc{end}.sSpacing)+1;

%get the q-factors maps, exclude boundary pixels since they are
%reliably poor.
qf1 = cc{end}.qfactors_accept{1}(1+trim_size(1):end-trim_size(1),...
    1+trim_size(2):end-trim_size(2));
qf2 = cc{end}.qfactors_accept{2}(1+trim_size(1):end-trim_size(1),...
    1+trim_size(2):end-trim_size(2));

%compute the % of bad correlation pts
num_pts = numel(qf1)+numel(qf2);
num_fail = sum(isnan(qf1(:))+isnan(qf2(:)));

if num_fail/num_pts*100 >= percentThresh
    decorrFlag = 1;
    disp('percent decorr')
end

%binarize the qf maps
qf1_bin = isnan(qf1);
qf2_bin = isnan(qf2);

%get connect regions
qf1_connectivity = bwconncomp(qf1_bin,8);
qf2_connectivity = bwconncomp(qf2_bin,8);

%find the number of px in each connect region
qf1_numPixels = cellfun(@numel,qf1_connectivity.PixelIdxList);
qf2_numPixels = cellfun(@numel,qf2_connectivity.PixelIdxList);

%find the largest area of the connected regions
[~,idx1] = max(qf1_numPixels);
[~,idx2] = max(qf2_numPixels);

%compute the region properties for the largest region, normalized by is
%eccentricity (since its harder to correctly predict far away from known data
%with inpaint nans)
if ~isempty(idx1)
    S1 = regionprops(qf1_connectivity,'Eccentricity','Area');
    ecc(1) = S1(idx1).Eccentricity;
    area(1) = S1(idx1).Area;
    if ecc(1)+0.25 > 1
        knockdown(1) = 1;
    else
        knockdown(1) = ecc(1)+0.25;
    end
    norm_area(1) = area(1)/knockdown(1);
else
    norm_area(1) = 0;
end

if ~isempty(idx2)
    S2 = regionprops(qf2_connectivity,'Eccentricity','Area');
    ecc(2) = S2(idx2).Eccentricity;
    area(2) = S2(idx2).Area;
    if ecc(2)+0.1 > 1
        knockdown(2) = 1;
    else
        knockdown(2) = ecc(2)+0.1;
    end
    norm_area(2) = area(2)*knockdown(2);
else
    norm_area(2) = 0;
end
% qf_holeSize(idx1) = 0;
% [qf_second,idx2] = max([qf1_numPixels,qf2_numPixels]);
% norm_area
if max(norm_area(:)) > sizeThresh
    decorrFlag = 1;
    disp('size decorr')
end

% do the st dev based thresholding
complete_qfactors = cc{end}.qfactors(3:4,:);
for ii = 1:size(complete_qfactors,1)
    
    x = sort(complete_qfactors(ii,:));
    
    %do the parameter estimation
    options = statset('MaxIter',2000, 'MaxFunEvals',4000);
    % options.FunValCheck = 'off';
    
    %     disp('Parameters estimated for single peak Gaussian')
    paramEsts = mle(x,'options',options);%,'optimfun','fmincon'
    paramEsts = [0.5,paramEsts(1),paramEsts(1),paramEsts(2),...
        paramEsts(2)];
    
    qfactor_fit_means(ii) = paramEsts(3);
    qfactor_fit_stds(ii) = paramEsts(5);
    
end

if cc{end}.sSize(1) == 16
    stDevThresh = stDevThresh + 0.01;
elseif cc{end}.sSize(1) == 64
    stDevThresh = stDevThresh - 0.01;
end

if qfactor_fit_stds(1) >= stDevThresh || qfactor_fit_stds(2) >= stDevThresh
    decorrFlag = 1;
    disp('stdev decorr')
end


end