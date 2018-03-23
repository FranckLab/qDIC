function [u, cc] = DIC(varargin)
% [du, cc] = DIC(I,sSize,sSpacing,ccThreshold) estimates
% displacements between two images through digital image
% correlation.
%
% INPUTS
% -------------------------------------------------------------------------
%   I: cell containing the undeformed, I{1}, and deformed, I{2} 2-D images
%   sSize: interrogation window (subset) size
%   sSpacing: interrogation window (subset) spacing.  Determines window
%             overlap factor
%   DICPadSize: padding arround the current window
%   ccThreshold: threshold values that define a bad cross-correlation
%   norm_xcc: flag for normalize cross-correlation cross-correlation
%
% OUTPUTS
% -------------------------------------------------------------------------
%   u: cell containing the displacement field (u{1:2} = {u_x, u_y})
%   cc: cc space maps, q-factors, and diagnostic struct
%
% NOTES
% -------------------------------------------------------------------------
% all functions are self contained
%
% If used please cite:
% Landauer, A.K., Patel, M., Henann, D.L. et al. Exp Mech (2018). 
% https://doi.org/10.1007/s11340-018-0377-4

% Parse inputs and create meshgrid
[I,m,mSize,sSize,sSizeFull,sSpacing,MTF,M,ccThreshold,norm_xcc] = ...
    parseInputs(varargin{:});

% Initialize variables
mSize_ = prod(mSize);
u12 = zeros(mSize_,2);
% cc = zeros(mSize_,1);

cc = struct('A',[],'max',[],'sSpacing',[],'sSize',[],'maxIdx',[],'qfactors',[],'q_thresh',[],...
    'qfactors_accept',[]);

%inject small noise to avoid flat subsets
noise = rand(size(I{1}))/10000;
I{1} = I{1}+noise;
I{2} = I{2}+noise;

A = cell(1,mSize_);

if strcmp(norm_xcc,'n')||strcmp(norm_xcc,'norm')||strcmp(norm_xcc,'normalized')
    
    parfor k = 1:mSize_
        
        %-----------------------------------------------------------------------
        % grab the moving subset from the images
        subst = I{1}(m{1}(k,:),m{2}(k,:));
        B = I{2}(m{1}(k,:),m{2}(k,:));
        %     size_sbst = size(subst);
        
        % multiply by the modular transfer function to alter frequency content
        subst = MTF.*subst; B = MTF.*B;
        
        % run cross-correlation
        A_ = normxcorr2(subst,B);
        A_small = A_(size(B,1)/2:size(B,1)+size(B,1)/2-1, ...
            size(B,2)/2:size(B,2)+size(B,2)/2-1);
        A{k} = imrotate(A_small,180);
       
        % find maximum index of the cross-correlaiton
        [max_val(k), maxIdx{k}] = max(A{k}(:));
        
        % compute pixel resolution displacements
        [u1, u2] = ind2sub(sSize,maxIdx{k});
        
        % gather the 3x3 pixel neighborhood around the peak
        try xCorrPeak = reshape(A{k}(u1 + (-1:1), u2 + (-1:1)),9,1);
            % least squares fitting of the peak to calculate sub-pixel disp
            du12 = lsqPolyFit2(xCorrPeak, M{1}, M{2});
            u12(k,:) = [u1 u2] + du12' - (sSize/2);
            %----------------------------------------------------------
        catch
            u12(k,:) = nan;
        end
    end
    
else
    parfor k = 1:mSize_
        
        %-----------------------------------------------------------------------
        % grab the moving subset from the images
        subst = I{1}(m{1}(k,:),m{2}(k,:));
        B = I{2}(m{1}(k,:),m{2}(k,:));
        %     size_sbst = size(subst);
        
        
        % multiply by the modular transfer function to alter frequency content
        subst = MTF.*subst; B = MTF.*B;
        A{k} = xCorr2(subst,B,sSize); %Custom fft-based correllation - fast,
        %but not normalized. Be careful using this
        %if there is a visable intensity gradient in
        %an image
        
        % find maximum index of the cross-correlaiton
        [max_val(k), maxIdx{k}] = max(A{k}(:));
        
        % compute pixel resolution displacements
        [u1, u2] = ind2sub(sSize,maxIdx{k});
        
        % gather the 3x3 pixel neighborhood around the peak
        try xCorrPeak = reshape(A{k}(u1 + (-1:1), u2 + (-1:1)),9,1);
            % least squares fitting of the peak to find sub-pixel disp
            du12 = lsqPolyFit2(xCorrPeak, M{1}, M{2});
            u12(k,:) = [u1 u2] + du12' - (sSize/2) - 1;
            %--------------------------------------------------------------
        catch
            u12(k,:) = nan;
        end
        
    end
    
    
end

cc.A = A;
cc.maxIdx = maxIdx;
cc.max = max_val;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reshape displacements and discover bad correlation points
if size(sSpacing,1) == 1
    spacingChange = 1;
elseif sSpacing(end,1) == sSpacing(end-1,1)
    spacingChange = 0;
else
    spacingChange = 1;
end

cc.max = reshape(double(cc.max),mSize);
cc.sSpacing = sSpacing(end,:);
cc.sSize = sSize;
[cc, ccMask_] = ...
    removeBadCorrelations(I,cc,ccThreshold,norm_xcc,spacingChange,mSize,mSize_);

for ii = 1:2
    ccMask{ii} = reshape(double(ccMask_(ii,:)),mSize);
end

u{1} = reshape(double(u12(:,2)),mSize);%.*ccMask{1}.*ccMask{2};
u{2} = reshape(double(u12(:,1)),mSize);%.*ccMask{1}.*ccMask{2};

end

%% ========================================================================
function varargout = parseInputs(varargin)
% Parse inputs and create meshgrid

I{1} = varargin{1}{1};
I{2} = varargin{1}{2};
sSizeFull = varargin{2};
sSpacing = varargin{3};
padSize = varargin{4};
ccThreshold = varargin{5};
norm_xcc = varargin{6};
% cc = varargin{7};

% pad images with zeros so that we don't grab any subset outside of the image
% domain. This would produce an error
I{1} = padarray(I{1},padSize,'replicate','both');
I{2} = padarray(I{2},padSize,'replicate','both');
sizeV = size(I{1});

sSize = sSizeFull(end,:);

% Initialize Mesh Variables
idx = cell(1,2);
for i = 1:2
    idx{i} = (1+padSize(i)) : sSpacing(i) : (sizeV(i)-sSize(i)-padSize(i)+1);
end
[m{1},m{2}] = ndgrid(idx{:});


% sSize = [sSize(2) sSize(1)];
mSize = size(m{1});
mSize_ = prod(mSize);

m_ = cell(1,2);
for k = 1:2
    m_{k} = zeros([mSize_,sSize(k)],'uint16');
    repmat_ = repmat((1:sSize(k))-1,mSize_,1);
    m_{k} = bsxfun(@plus, repmat_,m{k}(:));
end

% Initialize quadratic least squares fitting coefficients
[mx, my] = meshgrid((-1:1),(-1:1));
m = [mx(:), my(:)];

M{1} = zeros(size(m,1),6);
for i = 1:size(m,1)
    x = m(i,1); y = m(i,2);
    M{1}(i,:) = [1,x,y,x^2,x*y,y^2];
end

M{2} = M{1}'*M{1};

% Generate Moduluar transfer function (see eq. 3)
[~,~,MTF] = generateMTF(sSize);

%% Parse outputs

varargout{    1} = I;
varargout{end+1} = m_;
varargout{end+1} = mSize;
varargout{end+1} = sSize;
varargout{end+1} = sSizeFull;
varargout{end+1} = sSpacing;
varargout{end+1} = MTF;
varargout{end+1} = M;
varargout{end+1} = ccThreshold;
varargout{end+1} = norm_xcc;
% varargout{end+1} = cc;

end

%% ========================================================================
function A = xCorr2(A,B,sSize)
% performs fft based cross correlation of A and B (see equation 2)

A = fftn(A,sSize);
B = fftn(B,sSize);
B = conj(B);
A = A.*B;
A = ifftn(A);
A = real(A);
A = fftshift(A);

end

%% ========================================================================
function    duvw = lsqPolyFit2(b, M, trMM)
% LeastSqPoly performs a polynomial fit in the least squares sense
% Solves M*x = b,
% trMM = transpose(M)*M
% trMb = transpose(M)*b
%
% If you need to generate the coefficients then uncomment the following
% [mx, my, mz] = meshgrid(-1:1,-1:1,-1:1);
% m = [mx(:), my(:), mz(:)];
%
% for i = 1:size(m,1)
%    x = m(i,1); y = m(i,2); z = m(i,3);
%    M1(i,:) = [1,x,y,z,x^2,x*y,x*z,y^2,y*z,z^2];  %3D case
%               1 2 3 4  5   6   7   8   9   10

%    M1(i,:) = [1,x,y,x^2,x*y,y^2];                %2D Case
%               1 2 3  5   6   8
%               1 2 3  4   5   6
% end
%
% trMM1 = M'*M;

% b = log(b);
trMb = sum(bsxfun(@times, M, b));

x = trMM\trMb'; %solve for unknown coefficients

A = [x(5), 2*x(4);
    2*x(6),  x(5)];

duvw = (A\(-x([2 3])));

end

%% ========================================================================
function varargout = generateMTF(sSize)
% MTF functions taken from
% J. Nogueira, A Lecuona, P. A. Rodriguez, J. A. Alfaro, and A. Acosta.
% Limits on the resolution of correlation PIV iterative methods. Practical
% implementation and design of weighting functions. Exp. Fluids,
% 39(2):314{321, July 2005. doi: 10.1007/s00348-005-1017-1

%% equation 4

if prod(single(sSize == 16)) || prod(single(sSize == 32)) || prod(single(sSize == 64))
    sSize = sSize(1);
    
    x = cell(1,2);
    [x{1}, x{2}] = meshgrid(1:sSize,1:sSize);
    
    nu{1} = 1;
    for i = 1:2
        x{i} = x{i} - sSize/2 - 0.5;
        x{i} = abs(x{i}/sSize);
        nu{1} = nu{1}.*(3*(4*x{i}.^2-4*x{i}+1));
    end
    
    %% equation 5
    [x{1}, x{2}] = meshgrid(1:sSize,1:sSize);
    
    for i = 1:2, x{i} = x{i} - sSize/2 - 0.5; end
    
    r = abs(sqrt(x{1}.^2 + x{2}.^2)/sSize);
    nu{2}  = zeros(size(r));
    nu{2}(r < 0.5) = 24/pi*(4*r(r < 0.5).^2-4*r(r < 0.5)+1);
    
    %% equation 6
    [x{1}, x{2}] = meshgrid(1:sSize,1:sSize);
    
    nu{3} = 1;
    for i = 1:2
        x{i} = x{i} - sSize/2 - 0.5;
        x{i} = (x{i}/sSize);
        
        nu{3} = nu{3}.*(12*abs(x{i}).^2 - 12*abs(x{i}) + 3 + ...
            0.15*cos(4*pi*x{i}) + 0.20*cos(6*pi*x{i}) + ...
            0.10*cos(8*pi*x{i}) + 0.05*cos(10*pi*x{i}));
        
    end
    nu{3}(nu{3} < 0) = 0;
    
else
    
    nu{1} = ones(sSize(1),sSize(2));
    nu{2} = nu{1};
    nu{3} = nu{1};   %==== Based on the usage and the paper this comes from
    %nu should remain 3-dims even for the 2D case
    
end

nu = cellfun(@(x) x/sum(x(:)), nu, 'UniformOutput',0);
nu = cellfun(@sqrt, nu, 'UniformOutput',0);

varargout = nu;

end

%% ========================================================================
function [cc, ccMask] = ...
    removeBadCorrelations(I,cc,ccThreshold,norm_xcc,sizeChange,mSize,mSize_)
% flag bad correlation points, using q-factor analysis

if strcmp(norm_xcc,'n')||strcmp(norm_xcc,'norm')||strcmp(norm_xcc,'normalized')
    
    ccMask = ones(size(cc.max));
    
    for k = 1:mSize_
        
        cc_min = cc.A{k};% - min(cc.A{k}(:));, not needed since this is nxcc
        
        phi = sort(cc_min(:));
        
        P_img = imregionalmax(cc_min);
        
        peaks = unique(cc_min(P_img)); %some can be non-unique,
        %but this means its a bad xcc
        
        %compute several quality metrics, as given in "Xue Z, Particle Image
        % Velocimetry Correlation Signal-to-noise Metrics, Particle Image
        % Pattern Mutual Information and Measurement uncertainty Quantification.
        % MS Thesis, Virginia Tech, 2014.
        try
            ppr = peaks(end)/peaks(end-1);  %peak to second peak ratio (NOT USED)
            %min value = 1 (worst-case)
        catch
            ppr = 100; %a large value, since one peak only is good
        end
        prmsr = ((peaks(end)^2)/(rms(phi(phi<phi(end)/2))^2)); %peak to RMS ratio (NOT USED)
        %min value = 2; (worst case)
        
        %peak to corr. energy ratio
        pce = ((peaks(end)^2)/(1/numel(phi)*(sum(abs(phi(:).^2)))));
        %min value -> 1 (worst case)
        
        [cc_hist,~] = histcounts(cc_min,30);
        entropy = 0;
        for i = 1:30
            p(i) = cc_hist(i)/sum(cc_hist);
            if p(i) == 0
                entropy = entropy+p(i);
            else
                entropy = entropy+p(i)*log(1/p(i));
            end
        end
        
        ppe = 1/entropy; %peak to cc (information) entropy
        %min value -> 0 (worst case)
        
        qfactors_(:,k) = [ppr,prmsr,pce,ppe];
        
    end
    
    for k = 1:size(qfactors_,1)
        qf_ = (qfactors_(k,:)-min(qfactors_(k,:)));
        cc.qfactors(k,:) = qf_/max(qf_);
    end
    
    if sizeChange == 1
        %recompute threshold, only use pce & ppe since these give the best
        %results emprically.
        
        stdev_prefact = 1.00;
        
        for ii = 1:2
            
            [qf_para{ii},single_distro] = bimodal_gauss_fit(cc.qfactors(2+ii,:));
            
            if single_distro == 0%(qf_para{ii}(2) + 2*qf_para{ii}(4)) < (qf_para{ii}(3) - 2*qf_para{ii}(5))
                cc.q_thresh{ii} = qf_para{ii}(3) - stdev_prefact*qf_para{ii}(5);
            elseif single_distro == 1
                cc.q_thresh{ii} = qf_para{ii}(3) - stdev_prefact*qf_para{ii}(5);
            else
                cc.q_thresh{ii} = qf_para{ii}(3) - stdev_prefact*qf_para{ii}(5);
            end 
        end
        q_trim = [cc.q_thresh{1};cc.q_thresh{2}];
    else
        q_trim = [cc.q_thresh{1};cc.q_thresh{2}];
    end
    
    %NaN the qfactor values that are below the threshold
    temp = bsxfun(@le,cc.qfactors(3:4,:),q_trim);
    qfactors_accept = cc.qfactors(3:4,:);
    qfactors_accept(temp) = NaN;
        
    for ii = 1:2
        cc.qfactors_accept{ii} = reshape(double(qfactors_accept(ii,:)),mSize);
        %         cc.U95_accept{ii} = reshape(double(U95_accept(ii,:)),mSize);
    end
    
    ccMask = ones(size(qfactors_accept)) + ...
        zeros(size(qfactors_accept)).*qfactors_accept;
    
else
    minOS = 1;
    for i = 1:2
        zeroIdx = I{i} == 0;
        threshold = mean2(I{i}(~zeroIdx));
        I_ = I{i}.*(I{i} < threshold);
        minOS = minOS*sum(I_(:))/sum(~zeroIdx(:));
    end
    cc.max = cc.max - minOS;
    cc.max = cc.max/(max(I{1}(:))*max(I{2}(:)));
    ccMask = double(cc.max >= ccThreshold);
    
    CC = bwconncomp(~ccMask);
    if length(CC.PixelIdxList) > 0
        [~,idx] = max(cellfun(@numel,CC.PixelIdxList));
        ccMask(CC.PixelIdxList{idx}) = inf;
    end
    ccMask(cc.max == 0) = nan;
    ccMask(~isfinite(ccMask)) = 0;
    cc = cc.max.*ccMask;
    
end

end

%% ========================================================================
function [paramEsts,single_distro] = bimodal_gauss_fit(x)
%This function takes a dataset and fits a bimodal Gaussian distro to it.

x = sort(x);

%set function for bimodal Gaussian
pdf_normmixture = @(x,p,mu1,mu2,sigma1,sigma2) ...
    p*normpdf(x,mu1,sigma1) + (1-p)*normpdf(x,mu2,sigma2);
pdf_single = @(x,mu1,sigma1) ...
    normpdf(x,mu1,sigma1);

%starting params, biased mixture toward "good" values,
%centered at quartiles, equal std dev.
pStart = 0.25;
muStart = quantile(x,[.10 .75]);
sigmaStart(1) = sqrt(var(x(1:round(length(x)/5))));
%- 0.25*diff(quantile(x,[0.01 0.25])).^2);
sigmaStart(2) = sqrt(var(x(round(length(x)/10):round(3*length(x)/4))));
%... - 0.25*diff(quantile(x,[0.25 0.75])).^2);%1:round(length(x)/2)
start = [pStart muStart sigmaStart];

%set lower and upper bounds
lb = [0 -inf -inf 0.00001 0.00001];
ub = [1 inf inf inf inf];

%do the parameter estimation
options = statset('MaxIter',2000, 'MaxFunEvals',4000);
% options.FunValCheck = 'off';
try
    single_distro = 0;
    paramEsts = mle(x, 'pdf',pdf_normmixture, 'start',start, ...
        'lower',lb, 'upper',ub, 'options',options);%,'optimfun','fmincon'
    
    if paramEsts(2)-paramEsts(4) >= paramEsts(3)+paramEsts(5) || ...
            paramEsts(2)+paramEsts(4) <= paramEsts(3)-paramEsts(5)
        
        single_distro = 1;
        %     disp('Parameters estimated for single peak Gaussian')
        paramEsts = mle(x,'options',options);%,'optimfun','fmincon'
        paramEsts = [0.5,paramEsts(1),paramEsts(1),paramEsts(2),...
            paramEsts(2)];
        
    end
    
catch
    single_distro = 1;
%     disp('Parameters estimated for single peak Gaussian')
    paramEsts = mle(x,'options',options);%,'optimfun','fmincon'
    paramEsts = [0.5,paramEsts(1),paramEsts(1),paramEsts(2),...
        paramEsts(2)];
end

% % %show the result
% % figure
% % % [~, bins] = 
% % histogram(x,100);
% % % bins = -2.5:.5:7.5;
% % % h = bar(bins,histc(x,bins)/(length(x)*0.5),'histc');
% % % histogram(x,100)
% % % h.FaceColor = [0.9 0.9 0.9];
% % xgrid = linspace(1.1*min(x),1.1*max(x),200);
% % pdfgrid = pdf_normmixture(xgrid,paramEsts(1),paramEsts(2),paramEsts(3),...
% %     paramEsts(4),paramEsts(5));
% % hold on
% % plot((paramEsts(3) - 2*paramEsts(5)),pdfgrid,'or')
% % plot((paramEsts(2) + 2*paramEsts(4)),pdfgrid,'*r')
% % plot(xgrid,pdfgrid,'-b')
% % hold off
% % xlabel('x')
% % ylabel('Probability Density')

end