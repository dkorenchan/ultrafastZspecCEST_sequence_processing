% ReadJEOLjdx: Reads in raw data acquired on a JEOL scanner and outputted
% as a .jdx file, as well as parameter values contained in the file header,
% and returns what is read
%
% INPUTS:
%   aqpath      -   String containing full path to .jdx file (i.e. ends in
%                   .jdx)
%
% OUTPUTS:
%   fid         -   Raw FID data as 1D or 2D array
%   otherdata   -   Structure containing header information
%
function [fid,otherdata]=ReadJEOLjdx(aqpath)
% Read into struct containing data and header information
jstruct=jcampreadJEOL(aqpath);
% Read data in jstruct.Blocks (stored as real, then complex for each 2D
% point) and convert to complex
np1d=numel(jstruct.Blocks(1).YData);
np2d=numel(jstruct.Blocks)/2;
fid=zeros(np2d,np1d);
for iii=1:np2d
    % There appears to be a glitch where the YData in every other block 
    % (imaginary FID component) is multiplied or divided by 10^8! 
    % Correct this below
    if max(abs(jstruct.Blocks(2*iii-1).YData))>1e5
        x=jstruct.Blocks(2*iii-1).YData.*1e-8;
    elseif max(abs(jstruct.Blocks(2*iii-1).YData))<1e-5
        x=jstruct.Blocks(2*iii-1).YData.*1e8;
    else
        x=jstruct.Blocks(2*iii-1).YData;
    end
    if max(abs(jstruct.Blocks(2*iii).YData))>1e5
        y=jstruct.Blocks(2*iii).YData.*1e-8;
    elseif max(abs(jstruct.Blocks(2*iii).YData))<1e-5
        y=jstruct.Blocks(2*iii).YData.*1e8;
    else
        y=jstruct.Blocks(2*iii).YData;
    end    
    fid(iii,:)=x+1i*y;
end
% Read through file again to get header information, store as structure
% otherdata
fin=fopen(aqpath,'r','n','US-ASCII');
while ~feof(fin)
    tline = strtrim(fgetl(fin));  % Read and trim the line
    if startsWith(tline, '##') && contains(tline, '=')
        % Extract header information
        [key, value] = strtok(tline, '=');
        key = strtrim(key(3:end));  % Remove '##' and trim spaces
        value = strtrim(value(2:end));  % Remove '=' and trim spaces
        
        % Convert key to a valid structure field name
        key = matlab.lang.makeValidName(key);
        
        % Convert numeric values if possible
        if ~isempty(value)
            if contains(value,'$$') || strcmp(value(1),'{') %parse out entries 
                %with units specified, or arrays
                if length(value)<60 || strcmp(value(end),'}') %can only 
                    %parse text successfully if it ends in '}'
                    otherdata.(key) = parseJEOLtextParam(value);
                else
                    otherdata.(key) = value;
                end
            else
                numValue = str2double(value);
                if isnan(numValue)
                    otherdata.(key) = value;
                else
                    otherdata.(key) = numValue;
                end
            end
        end

    else %append to previous line as text
        temp=otherdata.(key);
        otherdata.(key)=[temp tline];
        if strcmp(otherdata.(key)(end),'}') %parse out once '}' has been reached
            otherdata.(key) = parseJEOLtextParam(otherdata.(key));
        end
    end
end
end