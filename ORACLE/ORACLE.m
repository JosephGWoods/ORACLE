% Description: ORACLE frameworkbSSFPSignal

% This code is for research purposes only.

% Author of function: 
% Nils MJ Plähn, Bern, Switzerland
% E-mail: nils.plaehn@students.unibe.ch
% Department of Diagnostic, Interventional and Pediatric Radiology (DIPR), Inselspital, Bern University Hospital, University of Bern, Switzerland
% Translation Imaging Center (TIC), Swiss Institute for Translational and Entrepreneurial Medicine, Bern, Switzerland


% signal: bSSFP profile
% TR:          Repetition time in seconds
% alpha:       RF excitation angle in rad 

% for more parameter descriptions see: S1_bSSFP_Profile_Generation, S1_Simulation_ORACLE

% optional:

%  phi: RF phase increments can be passed if signal does not start at
%       0° and/or non-sorted phi are used

%  Alaising correction: 
%   The DFT modes {cm1,c0,c1} calculated in the ORACLE_fct
%   function can also be additionally corrected for aliasing correction by
%   simply inserting the function: 
%   [Ab,Bb,zb,xib]=Iterate_DTF2BSSFP_best(c0,cm1,c1,nPC,percent)
%   while c0=Ab, cm1=Bb*conj(zb) and c1 = Ab*zb
%   This takes more time but corrects for eventual aliasing if low nPC are used 

function [T1, T2, theta, M0] = ORACLE(signal, TR, alpha, input, phi, dim)

    if ~exist('dim','var') || isempty(dim)
        dim = length(size(signal)); % dim is assumed to be the last non-zero dimension
    end

    if strcmpi(input,'modes') || strcmpi(input,'mode')

        % Indexing array for selecting mode dimension
        % Get size of input data
        idx = repmat({':'}, 1, length(size(signal))); % makes a cell array of ':'

        idx{dim} = 1;
        cm1 = signal(idx{:});
        
        idx{dim} = 2;
        c0  = signal(idx{:});

        idx{dim} = 3;
        c1  = signal(idx{:});

    else

        %% 1) Pre-Quantification
        % i) uniformly distributed RF phase increments (if not started at 0° can
        %    just pass another "phi" array to this function)
        if isempty(phi)
            NPhaseCycle = size(signal, 4);
            phi         = linspace(0,2*pi,NPhaseCycle+1); 
            phi(end)    = [];
        end

        % ii) VZ depends on the handedness of the coordinate system (Right handed coordinate system has VZ=1 and lefthanded VZ=-1)
        VZ = 1; 
    
        % iii) get central (DFT) modes
        %      Hint: if aliasing correction is desired one can also just pass
        %      the bSSFP modes from the fixed point iteration into this
        %      function if nPC is low
        cm1 = conj(NPointFT(signal, dim, phi, -1.*VZ, true));
        c0  =      NPointFT(signal, dim, phi,  0.*VZ, true);
        c1  =      NPointFT(signal, dim, phi,  1.*VZ, true);

    end

    % iv) Auxilliary definitions 
    z   = c1 ./ c0;
    x   = abs(cm1) ./ abs(c0);
    r   = abs(z);

    %% 2) Quantification - analytical solution function
    % i) T2 quantification
    E2 = (r+x) ./ (1+x.*r);
    T2 = -TR ./ log(E2);

    % ii) T1 quantification
    a = E2 - 2*r + E2.*r.^2;
    b = E2 .* (1-2*E2.*r+r.^2);
    E1 = (a+b*cos(alpha)) ./ (b+a*cos(alpha)); 
    T1 = -TR ./ log(E1);
    
    % iii) PD quantification    
    M0  = abs(c0) + abs(cm1./r);
    M0  = M0 ./ (2*tan(alpha/2)) .* exp(TR./T2/2);

    % iv) B0 inhomogeneity/off-resonance quantification
    %      if df or dB desired by divide by "2*pi*TR" or 
    %      "gamma*TR" respectively -> easy to change 
    theta    = angle(z);

    % gamma = 2*pi*42.577*10^(6); % MHz/T for 1H protons
    % df    = theta/TR/2/pi;% if chemical shift is neglected
    % dB0   = theta/gamma/TR; % if chemical shift is neglected
end

