function [ucum] = inc2cum(u,dm,m,option)
%Conver incremental FIDIC/FIDVC displacement to cumulative displacements
%   Inputs:
%   -----------------------------------------------------------------------
%   u:  displacement field vector calculated from FIDIC/FIDVC. 
%       Format: cell array, which is a 2D/3D vector (components in x,y,z)  
%       per each time point (units are in voxels)
%         u{time}{1} = displacement in x-direction at t=time of size MxNxP
%         u{time}{2} = displacement in y-direction at t=time of size MxNxP
%         if 3D:
%         u{time}{3} = displacement in z-direction at t=time of size MxNxP
%   dm:  Spacing between grid points in 'u.'
%   m: meshgrid from DIC
%   option: Interpolation method
%   
%   Outputs:
%   u:  Cumulative displacement in the same format as input


%% Parse Inputs

% Change order of dimension 1 and 2 and replace edge nans. 
for i = 1:length(u)
    
    nan_mask{i}{1} = u{i}{1}./u{i}{1};
    nan_mask{i}{2} = u{i}{2}./u{i}{2};
    
    tempU = inpaint_nans(u{i}{2});
    u{i}{2} = inpaint_nans(u{i}{1});
    u{i}{1} = tempU;
end
clear tempU

% Find number of dimensions
nDim = length(size(u{1}{1})); %2D or 3D

%% Start inc to cum

% Create reference and deformed grid at t=0
% sizeU = num2cell(size(u{1}{1}));
% for i = 1:nDim
%     sizeU{i} = 1:dm:dm*sizeU{i};
% end
% [refGrid{1:length(sizeU)}] = ndgrid(sizeU{:}); % reference grid

refGrid = m;
defGrid = refGrid;  %Deformed grid is same as refGrid at t=0;


% Convert inc to cum at all time points
for t = 1:length(u)
    
    % For each dimensions
    for i = 1:nDim
        if t==1
            udef{i} = u{t}{i};
        else
            
            % Find u at deformed grid points
            udef{i} = interpn(refGrid{:},u{t}{i},defGrid{:},option,NaN);
        end
    end
    for i = 1:nDim
        % Update deformed grid points
        defGrid{i} = defGrid{i} + udef{i};        
        % Calculate cumulative u from deformed grid points
        ucum{t}{i} = defGrid{i} - refGrid{i};
    end
end


% Revert back u to original format (Change order of dimension 1 and 2.) 
for i = 1:length(ucum)
    tempU = ucum{i}{2};
    ucum{i}{2} = ucum{i}{1};
    ucum{i}{1} = tempU;
    ucum{i}{3} = sqrt(ucum{i}{1}.^2 + ucum{i}{2}.^2);
    
    ucum{i}{1} = nan_mask{i}{1}.*ucum{i}{1};
    ucum{i}{2} = nan_mask{i}{2}.*ucum{i}{2};
end

end