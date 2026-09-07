% setLPeakBounds_invZ: Sets start points, lower bounds, and upper bounds
%   for Lorentzian peak fitting in R1*cos^2(theta)*(1/Z-1) space.
%   Bounds are taken from CESTvsW1fit_1divZminus1xR1.m (MT_CEST_fit_QUESP
%   directory) for pools defined there; remaining pools are scaled
%   accordingly for the 1/Z amplitude space.
%
%   NOTE: amplitudes (parameter 1) are in units of R1*(1/Z-1) [s^-1] at
%   the peak center, which are larger than the (1-Z) amplitudes used by
%   the conventional setLPeakBounds. FWHM ranges are also wider because
%   exchange broadening is more apparent in 1/Z space.
%
%   INPUTS:     NONE - Edit values directly in this file to adjust output!
%
%   OUTPUTS:
%       x   -   Struct with .st, .lb, .ub arrays per pool. Parameters:
%               1  Amplitude     [R1*(1/Z-1) s^-1]
%               2  FWHM          [ppm]
%               3  omega_0       offset from water  [ppm]
%               4  phase         zero-order phase  [rad]
%
function x = setLPeakBounds_invZ()
% Water pool  (from CESTvsW1fit_1divZminus1xR1.m)
x.water.st= [2      3       0       0       ];
x.water.lb= [0.1    0.01    -0.1    0       ];
x.water.ub= [50     10      0.1     0       ];

% OH pool  (from CESTvsW1fit_1divZminus1xR1.m)
x.OH.st=    [.8     3       0.8     0       ];
x.OH.lb=    [0      0.2     0.6     0       ];
x.OH.ub=    [5      10      1.0     0       ];

% Amine pool  (from CESTvsW1fit_1divZminus1xR1.m)
x.amine.st= [1      3       3.0     0       ];
x.amine.lb= [0      0.2     2.5     0       ];
x.amine.ub= [10     10      3.5     0       ];

% Amide pool  (scaled from amine entry above)
x.amide.st= [1      3       3.5     0       ];
x.amide.lb= [0      0.2     3.0     0       ];
x.amide.ub= [10     10      4.0     0       ];

% Trp indole proton pool  (scaled from amine/amide entry)
x.Trp.st=   [0.5    2       5.4     0       ];
x.Trp.lb=   [0      0.2     5.1     0       ];
x.Trp.ub=   [5      10      5.7     0       ];

% 4.4 ppm pool  (scaled from amine entry)
x.ppm4pt4.st= [0.5  2       4.5     0       ];
x.ppm4pt4.lb= [0    0.2     4.0     0       ];
x.ppm4pt4.ub= [5    10      5.0     0       ];

% 7.3 ppm pool
x.ppm7pt3.st= [0.5  2       7.3     0       ];
x.ppm7pt3.lb= [0    0.2     7.0     0       ];
x.ppm7pt3.ub= [5    10      7.5     0       ];

% 9.8 ppm pool
x.ppm9pt8.st= [0.5  2       9.8     0       ];
x.ppm9pt8.lb= [0    0.2     9.0     0       ];
x.ppm9pt8.ub= [5    10      11.0    0       ];

% rNOE pool  (scaled from CEST pool entries; broad negative-ppm peak)
x.NOE.st=   [0.5    3       -3.5    0       ];
x.NOE.lb=   [0      0.2     -4.5    0       ];
x.NOE.ub=   [5      10      -2.5    0       ];

% MT pool  (placeholder; MT is typically pre-fitted in Z-space when using
% invZ mode, but bounds are included here for completeness)
x.MT.st=    [5      30      0       0       ];
x.MT.lb=    [0      10      -1      0       ];
x.MT.ub=    [50     80      1       0       ];
end
