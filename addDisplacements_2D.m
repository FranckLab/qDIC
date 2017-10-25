function [u, du, cc, m] = addDisplacements_2D(u0,du0,cc,cc0,m0,dm)
% u = addDisplacements(u0,du0,cc,cc0,m0,dm) adds displacements from previous
% iteratations to the current iterate's displacement field (du).
%
% INPUTS
% -------------------------------------------------------------------------
%   u0: past displacement field vector defined at every meshgrid point with
%      spacing dm. Format: cell array, each containing a 2D matrix
%         (components in x,y)
%         u0{1} = displacement in x-direction
%         u0{2} = displacement in y-direction
%         u0{3} = magnitude
%   du0: current displacements as
%         du0{1} = displacement in x-direction
%         du0{2} = displacement in y-direction
%   cc: final cross correlation matrix used to define interpolant locations
%   cc0: initial cross correlation matrix used to define interpolant locations
%   m0: mesh grid parameter
%   dm: subset spacing paramter used to set up grids
%
% OUTPUTS
% -------------------------------------------------------------------------
%   u: cell containing the added displacement fields
%   du: cell containing interpolated incremental displacements
%   cc: cross-correlation matrix interpolated to fit with u and du
%   m: final meshgrid
%
% NOTES
% -------------------------------------------------------------------------
%
% If used please cite:
% Bar-Kochba E., Toyjanova J., Andrews E., Kim K., Franck C. (2014) A fast
% iterative digital volume correlation algorithm for large deformations.
% Experimental Mechanics. doi: 10.1007/s11340-014-9874-2

interp_opt = 'spline';
inpaint_opt = 0;

for i = 1:2, du0{i} = inpaint_nans(du0{i},inpaint_opt); end % remove NaNs if present

idx = cell(1,2);
for i = 1:2, idx{i} = m0{i}(1):dm:m0{i}(end); end % construct new meshgrid

[m0_{1}, m0_{2}] = ndgrid(m0{1},m0{2});
[m{1}, m{2}] = ndgrid(idx{1},idx{2});

% sample to desired mesh spacing
du = cell(1,2);
for i = 1:2
    F = griddedInterpolant(m0_{1}, m0_{2}, du0{i}, interp_opt);
    du{i} = F(m{1},m{2});
    
    F_qf = griddedInterpolant(m0_{1}, m0_{2}, cc0.qfactors_accept{i}, 'nearest');
    cc.qfactors_accept{i} = F_qf(m{1},m{2});
    
end

F_cc = griddedInterpolant(m0_{1}, m0_{2}, cc0.max, interp_opt);
cc.max = F_cc(m{1},m{2});

if  sum(cellfun(@numel,u0)) == 2 || mean(u0{1}(:)) == 0, u = du; % on first iteration u = du
else
    u = cellfun(@plus,u0,du,'UniformOutput',0);
    % else u^(k) = sum(u^(k-1)) + du (see eq. 7)
end

end
