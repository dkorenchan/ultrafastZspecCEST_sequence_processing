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
    % Correct this below.
    % On the 1st 2D datapoint, ID what the multiplicative factor should be
    % for real and imaginary, then apply the same for all other 2D 
    % datapoints
    if iii==1
        if max(abs(jstruct.Blocks(2*iii-1).YData))>1e5
            Refac=1e-8;
        elseif max(abs(jstruct.Blocks(2*iii-1).YData))<1e-5
            Refac=1e8;
        else
            Refac=1;
        end
        if max(abs(jstruct.Blocks(2*iii).YData))>1e5
            Imfac=1e-8;
        elseif max(abs(jstruct.Blocks(2*iii).YData))<1e-5
            Imfac=1e8;
        else
            Imfac=1;
        end    
    end
    fid(iii,:)=jstruct.Blocks(2*iii-1).YData.*Refac ...
        + 1i*jstruct.Blocks(2*iii).YData.*Imfac;
end

% Unfortunately, the 1E8 correction factor above might be off by a factor 
% of 10.... Check this by seeing if the max abs value is <=100
while max(max(abs(fid)))>100
    if Imfac<1
        fid=real(fid)+1i.*imag(fid)./10;
    elseif Refac<1
        fid=real(fid)./10+1i.*imag(fid);
    else
        fid=fid./10;
    end
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