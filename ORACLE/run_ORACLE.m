function [T1, T2, theta, M0] = run_ORACLE(signal, TR, alpha, bPlot, varargin)

p = inputParser;
p.addParameter('input', 'phasecycles', @(x) ischar(x));
p.addParameter('phi'  , []           , @(x) isnumeric(x) || isempty(x));
p.parse(varargin{:});
p = p.Results;
input = p.input;
phi   = p.phi;

if isempty(phi)
    NPhaseCycle = size(signal,4);
    phi_tmp     = linspace(0,2*pi,NPhaseCycle+1);
    phi         = phi_tmp(1:NPhaseCycle);
end

[T1, T2, theta, M0] = ORACLE_3D(signal, TR, alpha, input, phi);

% Remove NaN values
T1(isnan(T1)) = 0;
T2(isnan(T2)) = 0;
theta(isnan(theta)) = 0;
M0(isnan(M0)) = 0;

T1 = abs(T1);

T1(T1<0) = 0;
T2(T2<0) = 0;

if bPlot
    mosaic(T1, 'color', 'L8', 'window', [0.05,3], 'logscale', true, 'colorbar', true);
    mosaic(T2, 'color', 'L8', 'window', [0.005,1], 'logscale', true, 'colorbar', true);
    mosaic(theta, 'color', 'L8', 'window', [-pi, pi], 'colorbar', true);
    mosaic(M0, 'color', 'L8', 'window', [0, prctile(M0(:),99)], 'colorbar', true);
end

