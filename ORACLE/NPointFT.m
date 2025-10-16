function data_mode = NPointFT(data, mode, PCycInc, bnorm)
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

if ~exist('bnorm','var') || isempty(bnorm)
    bnorm = true;
end

if isvector(data)

    data    = data(:);
    PCycInc = PCycInc(:);

    if numel(data) ~= numel(PCycInc)
        error('MatInput and phi must have the same number of elements!');
    end

    N         = numel(data);
    data_mode = (1/N) * sum( data(:) .* exp(1i.*mode.*PCycInc(:)) );

elseif ndims(data) == 4 % 3D image with phase cycles in 4th dimension

    N  = size(data,4);

    if numel(PCycInc) ~= N
        error('size(MatInput,4) and numel(phi) must be the same!');
    end

    PCycInc   = reshape(PCycInc,1,1,1,N);
    data_mode = sum( data .* exp(1i*mode*PCycInc) , 4 );

    if bnorm
        data_mode = (1/N) * data_mode;
    end

end
