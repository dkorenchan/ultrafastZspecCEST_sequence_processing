% parseJEOLtextParam: Takes a string containing a value or array of 
% multiple values extracted from a JEOL parameter file and returns it as a
% single value or numerical array
%
%   INPUTS:
%       txtval  -   String containing numeric parameter value(s), either a 
%                   single value followed by a unit (e.g. '2   $$s') or 
%                   multiple values as a comma-separated list within 
%                   brackets along with units (e.g. {0[Hz], 1[Hz], 10[Hz]} )
%
%   OUTPUTS:
%       val     -   Numeric value or array. If single value, it will be 
%                   adjusted to base units based (e.g. if units were $$kHz, 
%                   the value will be returned in Hz)
%       unit    -   String denoting units of value(s) in val. If val is an 
%                   array, this in turn will be a cell array of strings, 
%                   one per value in val 
% 
function [val,unit]=parseJEOLtextParam(txtval)
if contains(txtval,'$$')
    par=strsplit(txtval,' $$');
    val=str2double(par{1});
    % Adjust value based upon prefix (micro, milli, kilo, mega)
    unit=par{2};
    if strcmp(unit(1),'M') %since switch is case-insensitive, use a different 
        %letter if 'M' found
        unit(1)='n';
    end
    switch unit(1) %if a prefix is found, we will take it off of unit and 
        % adjust the value val
        case 'u'
            val=val/1e6;
            unit=unit(2:end);
        case 'm'
            val=val/1e3;
            unit=unit(2:end);
        case 'k'
            val=val*1e3;
            unit=unit(2:end);
        case 'n'
            val=val*1e6;
            unit=unit(2:end);
    end
elseif strcmp(txtval(1),'{')
    vals=extractBetween(txtval,'{','}');
    vals=split(vals,','); %separate out entries
    unit=extractBetween(vals,'[',']'); %extract units
    vals=extractBefore(vals,'['); %isolate numbers
    val=zeros(numel(vals),1);
    for iii=1:numel(vals)
        val(iii)=str2double(vals{iii});
    end
end
end