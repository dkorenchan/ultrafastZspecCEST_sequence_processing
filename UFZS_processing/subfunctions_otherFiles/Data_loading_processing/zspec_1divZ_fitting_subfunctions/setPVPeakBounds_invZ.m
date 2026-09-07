% setPVPeakBounds_invZ: Sets start points, lower bounds, and upper bounds
%   for Pseudo-Voigt peak fitting in R1*cos^2(theta)*(1/Z-1) space.
%   Bounds are taken from CESTvsW1fit_1divZminus1xR1.m (MT_CEST_fit_QUESP
%   directory) for pools defined there; remaining pools are scaled
%   accordingly for the 1/Z amplitude space.
%
%   NOTE: amplitudes (parameter 1) are in units of R1*(1/Z-1) [s^-1] at
%   the peak center, which are larger than the (1-Z) amplitudes used by
%   the conventional setPVPeakBounds.
%
%   INPUTS:     NONE - Edit values directly in this file to adjust output!
%
%   OUTPUTS:
%       x   -   Struct with .st, .lb, .ub arrays per pool. Parameters:
%               1  Ai        peak amplitude  [R1*(1/Z-1) s^-1]
%               2  alpha     Gaussian proportion
%               3  FWHMl     Lorentzian FWHM  [ppm]
%               4  FWHMrat   Gaussian/Lorentzian FWHM ratio (1 to 2)
%               5  omega_0   offset from water  [ppm]
%               6  phase     zero-order phase  [rad]
%
function x = setPVPeakBounds_invZ()
% Water pool  (from CESTvsW1fit_1divZminus1xR1.m)
x.water.st= [2      .3      1       1   0       0       ];
x.water.lb= [0.1    0       0.01    1   -0.1    0       ];
x.water.ub= [50     1       5       2   0.1     0       ];

% OH pool  (from CESTvsW1fit_1divZminus1xR1.m)
x.OH.st=    [.3     .3      1       1   1.0     0       ];
x.OH.lb=    [0      0       0.2     1   0.1     0       ];
x.OH.ub=    [20     1       5       2   1.4     0       ];

% Amine pool  (from CESTvsW1fit_1divZminus1xR1.m)
x.amine.st= [1      .3      1       1   3.0     0       ];
x.amine.lb= [0      0       1       1   2.8     0       ];
x.amine.ub= [2      1       5       2   3.1     0       ];

% Amide pool  (from CESTvsW1fit_1divZminus1xR1.m)
x.amide.st= [1      .3      1.5     1   3.5     0       ];
x.amide.lb= [0      0       1       1   3.0     0       ];
x.amide.ub= [4      1       5       2   4.0     0       ];

% Trp indole proton pool  (scaled from amide/ppm4pt4 entry)
x.Trp.st=   [1      .3      1.5     1   5.4     0       ];
x.Trp.lb=   [0      0       1       1   5.1     0       ];
x.Trp.ub=   [4      1       5       2   5.7     0       ];

% 4.4 ppm pool  (from CESTvsW1fit_1divZminus1xR1.m)
x.ppm4pt4.st= [1    .3      1.5     1   4.5     0       ];
x.ppm4pt4.lb= [0    0       1       1   4.0     0       ];
x.ppm4pt4.ub= [2    1       5       2   5.5     0       ];

% 7.3 ppm pool  (from CESTvsW1fit_1divZminus1xR1.m)
x.ppm7pt3.st= [1    .3      1.5     1   7       0       ];
x.ppm7pt3.lb= [0    0       1       1   6.5     0       ];
x.ppm7pt3.ub= [2    1       5       2   7.5     0       ];

% 9.8 ppm pool  (from CESTvsW1fit_1divZminus1xR1.m)
x.ppm9pt8.st= [1    .3      1.5     1   10.0    0       ];
x.ppm9pt8.lb= [0    0       1.5     1   9.0     0       ];
x.ppm9pt8.ub= [2    1       5       2   11.0    0       ];

% rNOE pool  (scaled from CEST pool entries; broad negative-ppm peak)
x.NOE.st=   [0.5    .3      2       1   -3.5    0       ];
x.NOE.lb=   [0      0       1       1   -4.5    0       ];
x.NOE.ub=   [3      1       8       2   -2.5    0       ];

% MT pool  (placeholder; MT is typically pre-fitted in Z-space when using
% invZ mode, but bounds are included here for completeness)
x.MT.st=    [0.5    0       30      1   0       0       ];
x.MT.lb=    [0      0       10      1   0       0       ];
x.MT.ub=    [5      1       80      2   1       0       ];
end
