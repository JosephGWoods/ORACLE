function data_mode = NPointFT(data, dim, varargin)
% NPointFT - DFT over bSSFP RF phase cycle increments to obtain Fourier coefficients (bSSFP modes)
%
% data_mode = NPointFT(data, dim, varargin)
%
% in:
%      data  - complex phase-cycled bSSFP data (arbitrary dimensions, but must include a dimension of phase-cycled data along dim)
%      dim   - dimension along which to perform the DFT (i.e. the phase cycle dimension)
%      varargin - optional parameters as name-value pairs:
%           'phi':      phase cycle increments in radians
%                       (default: phi = linspace(0, 2*pi, size(data, dim) + 1); phi(end) = [])
%           'mode':     mode index (default = [-1, 0, 1])
%           'use_norm': normalize by 1/N (for N phase cycles) (default = true)
%           'use_sum':  sum across phase cycle increments during the Fourier transform (default = true)
%                       If false, does not sum across phase cycle increments, and instead outputs the
%                       Fourier coefficients for each phase cycle increment along a new mode dimension
%                       at the end of the output array.
%
% out:
%      data_mode - complex mode bSSFP data
%                  (if use_sum = true, size(data_mode) = size(data) except along dim, where size(data_mode, dim) = length(mode))
%                  (if use_sum = false, size(data_mode) = [size(data), length(mode)], with mode dimension at the end)
%
if ~exist('data','var') || isempty(data); error('data must be provided!'); end
if ~exist('dim','var')  || isempty(dim);  error('dim must be provided!');  end

% Parse inputs
p = inputParser;
p.addParameter('phi'  , []         , @(x) isnumeric(x) || isempty(x));
p.addParameter('mode' , [-1, 0, 1] , @isvector);
p.addParameter('bnorm', true       , @islogical);
p.addParameter('bsum' , true       , @islogical);
p.parse(varargin{:});
p = p.Results;
phi   = p.phi;
mode  = p.mode;
bnorm = p.bnorm;
bsum  = p.bsum;
if isempty(phi)
    phi = linspace(0, 2*pi, size(data, dim) + 1);
    phi(end) = [];
end


% Get size of input data
size_data = size(data);

% Set up PCycInc dimensions for broadcasting
idx      = ones(1, length(size_data)); % makes a cell array of '1'
idx(dim) = length(phi);
phi      = reshape(phi, idx);

% Setup size of output mode data
if bsum
    % If the Fourier transform sum is to be performed,
    % update size of PCyc dimension to number of modes
    size_data_mode      = size_data;
    size_data_mode(dim) = length(mode);
else
    % If no Fourier transform sum, put mode dimension at the end.
    size_data_mode = [];
    for indMode = 1:length(mode)
        size_data_mode = cat(length(size_data)+1,size_data_mode,size_data);
    end
end

% Indexing array for selecting mode dimension
idx = repmat({':'}, 1, length(size_data_mode)); % makes a cell array of ':'

% Perform DFT along dim
data_mode = zeros(size_data_mode);
if bsum
    for indMode = 1:length(mode)
        idx{dim}          = indMode; % replace ':' with indices along desired direction
        data_mode(idx{:}) = sum( data .* exp(1i .* mode(indMode) .* phi) , dim );
    end
else
    for indMode = 1:length(mode)
        idx{end}          = indMode; % replace ':' with indices along the mode dimension
        data_mode(idx{:}) = data .* exp(1i .* mode(indMode) .* phi);
    end
end

% Normalize by number of phase cycles (mean, rather than sum)
if bnorm
    data_mode = data_mode / size_data(dim);
end
