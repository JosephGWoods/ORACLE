function data_mode = NPointFT(data, dim, PCycInc, mode, bnorm)
% Description: 
% DFT over bSSFP RF phase cycle increments to obtain Fourier coefficients (bSSFP modes)
%
% in:
%      data    - complex phase-cycled bSSFP data [length(PCycInc),1] or [1,length(PCycInc)] or [x,y,z,length(PCycInc)]
%      mode    - mode index (e.g., -1, 0, 1)
%      PCycInc - phase cycle increments (radians)
%      bnorm   - normalize by 1/N (for N phase cycles) (default = true)
%
% out:
%      data_mode - complex mode bSSFP data scalar or [x,y,z,1]
%
if ~exist('data','var') || isempty(data)
    error('data must be provided!')
end
if ~exist('dim','var') || isempty(dim)
    error('dim must be provided!')
end
if ~exist('PCycInc','var') || isempty(PCycInc)
    PCycInc = linspace(0, 2*pi, size(data, dim) + 1);
    PCycInc(end) = [];
end
if ~exist('mode','var') || isempty(mode)
    PCycInc = [-1, 0, 1];
end
if ~exist('bnorm','var') || isempty(bnorm)
    bnorm = true;
end

% Get size of input data
size_data = size(data);

% Setup size of output mode data
size_data_mode      = size_data;
size_data_mode(dim) = length(mode);

% Set up PCycInc dimensions for broadcasting
idx      = ones(1, length(size_data)); % makes a cell array of '1'
idx(dim) = length(PCycInc);
PCycInc  = reshape(PCycInc, idx);

% Indexing array for selecting mode dimension
idx = repmat({':'}, 1, length(size_data)); % makes a cell array of ':'

% Perform DFT along dim
data_mode = zeros(size_data_mode);
for indMode = 1:length(mode)
    idx{dim}          = indMode; % replace ':' with indices along desired direction
    data_mode(idx{:}) = sum( data .* exp(1i .* mode(indMode) .* PCycInc) , dim );
end

% Normalize by number of phase cycles (mean, rather than sum)
if bnorm
    N = size_data(dim); % Number of phase cycles
    data_mode = (1/N) * data_mode;
end

% if isvector(data)
% 
%     data    = data(:);
%     PCycInc = PCycInc(:);
% 
%     if numel(data) ~= numel(PCycInc)
%         error('MatInput and phi must have the same number of elements!');
%     end
% 
%     data_mode = sum( data(:) .* exp(1i .* mode .* PCycInc(:)) );
% 
%     if bnorm
%         N         = numel(data);
%         data_mode = (1/N) * data_mode;
%     end
% 
% elseif ndims(data) == 4 % 3D image with phase cycles in 4th dimension
% 
%     N  = size(data,4);
% 
%     if numel(PCycInc) ~= N
%         error('size(MatInput,4) and numel(phi) must be the same!');
%     end
% 
%     PCycInc   = reshape(PCycInc,1,1,1,N);
%     data_mode = sum( data .* exp(1i*mode*PCycInc) , 4 );
% 
%     if bnorm
%         data_mode = (1/N) * data_mode;
%     end
% 
% end
