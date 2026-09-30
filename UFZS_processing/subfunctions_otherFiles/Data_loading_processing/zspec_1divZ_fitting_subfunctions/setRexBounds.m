% setRexBounds: Fitting bounds for Rex (Zaiss & Bachert 2013, Eq. 23 terms 1+2).
%   Water uses R2w*sin^2(theta) model:
%       p(1) = R2w      water T2 rate         [s^-1]
%   Each non-water pool has 4 parameters (Gamma approximated with kb>>R2b):
%       p(1) = fb        mole fraction         [dimensionless]
%       p(2) = kb        exchange rate         [s^-1]
%       p(3) = delta_ppm chemical shift        [ppm, relative to water]
%       p(4) = R2b       solute T2 rate        [s^-1]
%
%   INPUTS:     NONE - Edit values directly in this file to adjust output!
%   OUTPUTS:
%       x   -   Struct with .st, .lb, .ub arrays per pool.
%
function x = setRexBounds()
% Water pool  (R2w*sin^2(theta) model: R2w = water T2 relaxation rate [s^-1]; R1w = water T1 realxation rate [s^-1])
x.water.st= [1      0.5];
x.water.lb= [0.0001  0.5];
x.water.ub= [50     0.5  ];

% OH pool
x.OH.st=    [2e-5    500     1.0     50  ];
x.OH.lb=    [0       10      0.1     0   ];
x.OH.ub=    [1e-3    10000   1.4     500 ];

% Amine pool
x.amine.st= [1e-5    4000    3.0      0  ];
x.amine.lb= [0       100      2.8     0   ];
x.amine.ub= [1e-2    20000   3.3      0 ];

% Amide pool
x.amide.st= [1e-5    50      3.5     50  ];
x.amide.lb= [0       10      3.0     0   ];
x.amide.ub= [1e-3    1000    4.0     500 ];

% Trp indole proton pool
x.Trp.st=   [1e-5    500     5.4     50  ];
x.Trp.lb=   [0       10      5.1     0   ];
x.Trp.ub=   [1e-3    10000   5.7     500 ];

% 4.4 ppm pool
x.ppm4pt4.st= [1e-5  500     4.5     50  ];
x.ppm4pt4.lb= [0     10      4.0     0   ];
x.ppm4pt4.ub= [1e-3  10000   5.5     500 ];

% 7.3 ppm pool
x.ppm7pt3.st= [1e-5  500     7.0     50  ];
x.ppm7pt3.lb= [0     10      6.5     0   ];
x.ppm7pt3.ub= [1e-3  10000   7.5     500 ];

% 9.8 ppm pool
x.ppm9pt8.st= [1e-5  500     9.8     50  ];
x.ppm9pt8.lb= [0     10      9.0     0   ];
x.ppm9pt8.ub= [1e-3  10000   11.0    500 ];

% rNOE pool
x.NOE.st=   [1e-5    500     -3.5    50  ];
x.NOE.lb=   [0       10      -4.5    0   ];
x.NOE.ub=   [1e-3    10000   -2.5    500 ];
end
