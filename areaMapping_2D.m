function I = areaMapping_2D(varargin)
% I = areaMapping(I0,m,u0) symmetrically warps undeformed
% and deformed 2-D images by the displacement field from previous iteration
% using trilinear interpolation.
%
% INPUTS
% -------------------------------------------------------------------------
%   I0: cell containing the undeformed, I0{1}, and deformed, I0{2} images
%   m: meshgrid of the displacement field
%   u0: displacement field vector defined at every meshgrid point with 
%      spacing dm. Format: cell array, each containing a 2D matrix 
%         (components in x,y,z)
%         u0{1} = displacement in x-direction
%         u0{2} = displacement in y-direction
%         u0{3} = magnitude
%
% OUTPUTS
% -------------------------------------------------------------------------
%   I: cell containing the symmetrically warped images of I0{1} and I0{2}
%
% NOTES
% -------------------------------------------------------------------------
% if interp2FastMex and mirt2D-mexinterp are used to perform the interpolation
% since they are faster and less memory demanding than MATLAB's interp3
% function. To run you need a compatible C compiler. Please see
% (http://www.mathworks.com/support/compilers/R2014a/index.html)
%
% For more information please see section 2.2.
% If used please cite:
% Bar-Kochba E., Toyjanova J., Andrews E., Kim K., Franck C. (2014) A fast 
% iterative digital volume correlation algorithm for large deformations. 
% Experimental Mechanics. doi: 10.1007/s11340-014-9874-2

interp_opt = 'spline';
inpaint_opt = 0;

[I0,m0,u0] = parseInputs(varargin{:});

idx = cell(1,2); idx_u0 = cell(1,2);
for i = 1:2
    idx{i} = m0{i}(1):(m0{i}(end) - 1); % get index for valid indices of image
    idx_u0{i} = linspace(1,size(u0{1},i),length(idx{i})+1); % get index for displacement meshgrid
    idx_u0{i} = idx_u0{i}(1:length(idx{i}));
    u0{i} = double(u0{i}*0.5);
end

[m_u0{2}, m_u0{1}] = ndgrid(idx_u0{:});
[m{2}, m{1}] = ndgrid(idx{:});

u = cell(1,2);
for i = 1:2
    u{i} = interp2(u0{i},m_u0{1}, m_u0{2},interp_opt);
    %u{i} = mirt2D_mexinterp(u0{i}, m_u0{1}, m_u0{2}); %Uses cpp for
                                                        %interpolation, minor
                                                        %performance gain at
                                                        %the cost of using a mex 
end
%% Warp images (see eq. 8)
mForward = cellfun(@(x,y) x + y, m, u, 'UniformOutput',false);
mBackward = cellfun(@(x,y) x - y, m, u, 'UniformOutput',false);

% interpolate images based on the deformation
% 
I{1} = inpaint_nans(interp2(I0{1},mForward{1}, mForward{2},interp_opt),inpaint_opt);
I{2} = inpaint_nans(interp2(I0{2},mBackward{1}, mBackward{2},interp_opt),inpaint_opt);

% I{1} = (interp2(I0{1},mForward{1}, mForward{2},'cubic'));
% I{2} = (interp2(I0{2},mBackward{1}, mBackward{2},'cubic'));
%I{1} = mirt2D_mexinterp(I0{1},  mForward{1}, mForward{2}); %mex-interp fcns
%I{2} = mirt2D_mexinterp(I0{2},  mBackward{1}, mBackward{2});
end

%% ========================================================================
function varargout = parseInputs(varargin)
I0 = varargin{1};
m  = varargin{2};
u0 = varargin{3};

% convert I0 to double.  Other datatypes will produce rounding errors
I0 = cellfun(@double, I0, 'UniformOutput', false);

varargout{      1} = I0;
varargout{end + 1} = m;
varargout{end + 1} = u0;

end
