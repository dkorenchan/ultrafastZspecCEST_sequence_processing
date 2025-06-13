% setSLPeakBounds:  Sets start points, lower bounds, and upper bounds for
%                   super-Lorentzian peak fitting
%
%   INPUTS:     NONE - Edit values directly in this file to adjust output!
%
%   OUTPUTS:
%       x   -   Struct containing start point (.st), lower bound (.lb) and
%               upper bound (.ub) arrays for the four parameters used to
%               fit super-Lorentzian peaks for the pool. The positional index 
%               of each array value is as follows:
%               1     Ai       --  peak amplitude
%               2     "FWHM"   --  peak linewidth (not truly the FWHM, though!!)
%               3     chems    --  displacement
%               4-6   (not currently used)
%
function x = setSLPeakBounds()
% MT pool
x.st=       [0.1    30      0   0   0   0];  % initial value 
x.lb=       [0.0    10      0   0   0   0];  % lower bound
x.ub=       [0.6    50      1   0   0   0];  % upper bound 
end