% Description: 
% DFT over RF phase increments to obtain Fourier coefficients

function Gp = NPointFT(MatInput, order, phi, bnorm)

if ~exist('bnorm','var') || isempty(bnorm)
    bnorm = true;
end

if isvector(MatInput)

    MatInput = MatInput(:);
    phi      = phi(:);

    if numel(MatInput) ~= numel(phi)
        error('MatInput and phi must have the same number of elements!');
    end

    N  = numel(MatInput);
    Gp = (1/N) * sum( MatInput(:) .* exp(1i.*order.*phi(:)) );

elseif ndims(MatInput) == 4 % 3D image with phase cycles in 4th dimension

    N  = size(MatInput,4);

    if numel(phi) ~= N
        error('size(MatInput,4) and numel(phi) must be the same!');
    end

    phi = reshape(phi,1,1,1,N);
    Gp = sum( MatInput .* exp(1i*order*phi) , 4 );

    if bnorm
        Gp = (1/N) * Gp;
    end

end
