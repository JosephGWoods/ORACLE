function [T1, T2, theta, M0] = run_ORACLE(data, dim, TR, FA, varargin)
% Wrapper function to run ORACLE quantification, check output values are sensible, and optionally plot the resulting T1, T2, theta, and M0 maps.
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
%           'do_plot':      boolean flag to plot the resulting T1, T2, theta, and M0 maps
%                           default = false
%
% out:
%      T1    - longitudinal relaxation time in seconds
%      T2    - transverse relaxation time in seconds
%      M0    - proton density (equilibrium magnetization)
%      theta - off-resonance phase in radians (can be converted to off-resonance frequency by dividing by 2*pi*TR)
%

p = inputParser;
p.addParameter('input'       , 'phase-cycles', @(x) strcmp(x, 'phase-cycles') || strcmp(x, 'modes'));
p.addParameter('phi'         , []            , @(x) isnumeric(x) || isempty(x));
p.addParameter('regularize_r', true          , @islogical);
p.addParameter('do_plot'     , false         , @islogical);
p.parse(varargin{:});
p = p.Results;
input        = p.input;
phi          = p.phi;
regularize_r = p.regularize_r;
do_plot      = p.do_plot;

% Call ORACLE to get T1, T2, theta, and M0 maps from the input signal
[T1, T2, theta, M0] = ORACLE(data, dim, TR, FA, 'input', input, 'phi', phi, 'regularize_r', regularize_r);

% Remove NaN values (e.g. outside of object)
T1(isnan(T1))       = 0;
T2(isnan(T2))       = 0;
theta(isnan(theta)) = 0;
M0(isnan(M0))       = 0;

% Ensure T1 is real
T1(~isreal(T1)) = abs(T1(~isreal(T1)));

% Set negative T1 and T2 values to 0 (they are not physically plausible)
T1(T1<0) = 0;
T2(T2<0) = 0;

% Finally take absolut value of T1
T1 = abs(T1);

if do_plot
    mosaic(T1   , 'color', 'L8', 'window', [0.05,3]              , 'colorbar', true, 'logscale', true );
    mosaic(T2   , 'color', 'L8', 'window', [0.005,1]             , 'colorbar', true, 'logscale', true );
    mosaic(theta, 'color', 'D1', 'window', [-pi, pi]             , 'colorbar', true, 'logscale', false);
    mosaic(M0   , 'color', 'L8', 'window', [0, prctile(M0(:),99)], 'colorbar', true, 'logscale', false);
end

