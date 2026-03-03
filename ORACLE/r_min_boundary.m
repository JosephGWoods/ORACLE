function r_min = r_min_boundary(FA, TR)
% Returns the minimum r boundary value for flip_angle and TR for physically plausibility
% Based on simulations using 20 ms ≤ T1 ≤ 6 s, 5 ms ≤ T2 ≤ 2.5 s, T2 ≤ T1 ≤ 100 T2.
% See eq A.33 in "Simultaneous and robust multi-parameterquantification in magnetic resonance
% imaging", PhD Thesis of Nils Marc Joel Plähn
%
% in:
%      FA - flip angle in radians
%      TR - repetition time in seconds
%
% out:
%      r_min - minimum r boundary value

A = 0.844 * 10^(-176.4 * TR);
B = 318.2 * 10^( 120.7 * TR);

r_min = A ./ (1 + B * tan(FA/2).^2);