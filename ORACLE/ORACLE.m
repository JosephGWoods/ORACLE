function [T1, T2, theta, M0] = ORACLE(data, dim, TR, FA, varargin)
% ORACLE - An analytical approach for T1, T2, proton density, and off-resonance mapping with
%          phase-cycled balanced steady-state free precession (bSSFP)
%
% [T1, T2, theta, M0] = ORACLE(data, dim, TR, FA, varargin)
%
% in:
%      data - complex bSSFP data (arbitrary dimensions, but must include a dimension of phase-cycled or mode data along dim)
%      dim  - dimension of the bSSFP data containing the phase cycles or modes
%      TR   - repetition time in seconds
%      FA   - flip angle in radians
%      varargin - optional parameters as name-value pairs:
%           'input':        specifies whether the input data is 'phase-cycles' or 'modes'
%                           default = 'phase-cycles'
%           'phi':          vector of phase cycling angles in radians (required if input = 'phase-cycles')
%                           default = linspace(0, 2*pi, size(data, dim) + 1); phi(end) = [];
%           'mode':         vector of mode indices corresponding to the modes in the input data (required if input = 'modes'
%                           default = [-1, 0, 1]
%           'regularize_r': boolean flag to regularize r based on physically plausible T1 and T2 values
%                           default = true
% out:
%      T1    - longitudinal relaxation time in seconds
%      T2    - transverse relaxation time in seconds
%      M0    - proton density (equilibrium magnetization)
%      theta - off-resonance phase in radians (can be converted to off-resonance frequency by dividing by 2*pi*TR)
%
% Reference:
%   Plähn NMJ, Safarkhanlo Y, Açikgöz BC, et al. ORACLE: An analytical approach for T1, T2, proton density,
%   and off-resonance mapping with phase-cycled balanced steady-state free precession. Magnetic Resonance
%   in Medicine. 2025;93(4). doi:10.1002/mrm.30388
%
% Notes:
%
%   For more parameter descriptions see: S1_bSSFP_Profile_Generation, S1_Simulation_ORACLE.
%
%   Alaising correction: 
%       The bSSFP modes (cm1,c0,c1) can be additionally corrected for aliasing correction by
%       simply inserting the function: 
%           [Ab,Bb,zb,xib] = Iterate_DTF2BSSFP_best(c0,cm1,c1,nPC,percent)
%       where c0=Ab, cm1=Bb*conj(zb) and c1 = Ab*zb
%       This takes more time but corrects for eventual aliasing if low nPC are used
%
% Written by Nils MJ Plähn, Bern, Switzerland
% Department of Diagnostic, Interventional and Pediatric Radiology (DIPR), Inselspital, Bern University Hospital, University of Bern, Switzerland
% Translation Imaging Center (TIC), Swiss Institute for Translational and Entrepreneurial Medicine, Bern, Switzerland

% Parse inputs
p = inputParser;
p.addParameter('input'       , 'phase-cycles', @(x) ischar(x) || isstring(x));
p.addParameter('phi'         , []            , @(x) isnumeric(x) || isempty(x));
p.addParameter('mode'        , [-1, 0, 1]    , @isvector);
p.addParameter('regularize_r', true          , @islogical);
p.parse(varargin{:});
p            = p.Results;
input        = p.input;
phi          = p.phi;
mode         = p.mode;
regularize_r = p.regularize_r;
if isempty(phi)
    phi = linspace(0, 2*pi, size(data, dim) + 1);
    phi(end) = [];
end

% Check input dimensions and values
if strcmpi(input, 'phase-cycles') && size(data, dim) ~= length(phi)
    error('Length of phi must match the size of data along the specified dimension dim.');
end
if strcmpi(input, 'modes') && size(data, dim) ~= length(mode)
    error('Length of mode must match the size of data along the specified dimension dim.');
end
if sort(mode) ~= [-1, 0, 1]
    error('mode must be a vector containing the central mode indices [-1, 0, 1] in any order.');
end

%% Preprocessing

% Extract modes (input = 'modes') or transform phase-cycled data to modes (input = 'phase-cycles')
switch input

    case {'modes', 'mode'}
        % Extract modes from input data
        % Assumes input data is already in mode space, with modes along the the
        % dimension specified by dim, and in the order of the optional mode variable
        % (by default, mode = [-1, 0, 1])

        % Indexing array for selecting mode dimension
        % Get size of input data
        idx = repmat({':'}, 1, length(size(data))); % makes a cell array of ':'

        % Extract -1st mode
        idx{dim} = find(mode==-1);
        cm1 = data(idx{:});
        
        % Extract 0th mode
        idx{dim} = find(mode==0);
        c0  = data(idx{:});

        % Extract 1st mode
        idx{dim} = find(mode==1);
        c1  = data(idx{:});

    case {'phase-cycles','phases','phase'}

        % 1) Pre-Quantification

        % ii) VZ depends on the handedness of the coordinate system
        % (Right handed coordinate system has VZ=1 and lefthanded VZ=-1)
        VZ = 1;
    
        % iii) get central (DFT) modes
        %      Hint: if aliasing correction is desired one can also just pass
        %      the bSSFP modes from the fixed point iteration into this
        %      function if nPC is low
        cm1 = conj(NPointFT(data, dim, 'phi', phi, 'mode', -1.*VZ, 'bsum', true, 'bnorm', true));
        c0  =      NPointFT(data, dim, 'phi', phi, 'mode',  0.*VZ, 'bsum', true, 'bnorm', true);
        c1  =      NPointFT(data, dim, 'phi', phi, 'mode',  1.*VZ, 'bsum', true, 'bnorm', true);

end

% Calculate auxilliary variables
z   = c1 ./ c0;
x   = abs(cm1) ./ abs(c0);
r   = abs(z);

% Regularize r based on physically plausible T1 and T2 values
if regularize_r
    r_min = r_min_boundary(FA, TR);
    r(r<r_min) = r_min;
end

%% Quantification - analytical solution functions

% T2 quantification
E2 = (r + x) ./ (1 + x.*r);
T2 = -TR ./ log(E2);

% T1 quantification
a = E2 - 2*r + E2.*r.^2;
b = E2 .* (1 - 2*E2.*r + r.^2);
E1 = (a + b*cos(FA)) ./ (b + a*cos(FA));
T1 = -TR ./ log(E1);

% PD quantification
M0  = abs(c0) + abs(cm1) ./ r;
M0  = M0 ./ (2*tan(FA/2)) .* exp(TR./T2/2);

% B0 inhomogeneity/off-resonance quantification
%   if df or dB desired by divide by "2*pi*TR" or "gamma*TR" respectively -> easy to change
theta = angle(z);

% gamma = 2*pi*42.577*10^(6); % MHz/T for 1H protons
% df    = theta/TR/2/pi;% if chemical shift is neglected
% dB0   = theta/gamma/TR; % if chemical shift is neglected
