% ReadTopspinDataPars: Reads in both data, from a 1D or 2D file containing 
% FID or spectral data, and parameter information for a Bruker experiment
%
%   INPUTS:
%       aqpath      -   String containing the full path to the main scan
%                       directory
%       filestr     -   String indicating the type of file to read data
%                       from. Possible options are:
%                       'fid'   -   1D raw FID data
%                       'ser'   -   2D raw FID data
%                       '1r'    -   1D TopSpin-processed spectral data
%                       '2rr'   -   2D TopSpin-processed spectral data
%                       (NOTE: the '1r' and '2rr' data files are assumed to
%                       be in procno 1 (i.e. aqpath/pdata/1/)
%
%   OUTPUTS:
%       data        -   Array containing data read in from Bruker data
%                       file, with the first dimension being each
%                       individual FID/spectrum and the second dimension 
%                       being the time/frequency axis
%       otherdata   -   Struct containing parameter values obtained from
%                       the Bruker parameter files
%   Updated 12/2/25 by DK to be compatible with Topspin versions < 3
%   
function [data,otherdata]=ReadTopspinDataPars(aqpath,filestr)
procpath = fullfile(aqpath,'pdata','1');

% First, let's figure out which version of Topspin is being run, so we know
% how to read in the desired parameters! If below Topspin 3, we will need
% to read in PL rather than PLW
TSvertxt=readParsTopspin(aqpath,'acqus',{'##TITLE'});
TSvertxtsplit=split(TSvertxt,' ');
TSver=str2double(TSvertxtsplit{end}(1));

if TSver<3
    acqusPars=readParsTopspin(aqpath,'acqus',{'##$TD','##$BYTORDA','##$DTYPA','##$D',...
        '##$P','##$CNST','##$SW_h','##$SW','##$SFO1','##$PL'});
else
    acqusPars=readParsTopspin(aqpath,'acqus',{'##$TD','##$BYTORDA','##$DTYPA','##$D',...
        '##$P','##$CNST','##$SW_h','##$SW','##$SFO1','##$PLW'});
end
procsPars=readParsTopspin(procpath,'procs',{'##$BYTORDP','##$DTYPP',...
    '##$SF','##$OFFSET','##$SW_p','##$XDIM','##$SI'});
if sum(strcmp(filestr,{'ser','2rr'}))>0 %also pull in 2D info
    acqu2sPars=readParsTopspin(aqpath,'acqu2s',{'##$TD'});
    proc2sPars=readParsTopspin(procpath,'proc2s',{'##$SF','##$OFFSET','##$SW_p',...
        '##$XDIM','##$SI'});
    td1=str2double(acqu2sPars{1});
    xdim1=str2double(proc2sPars{4});
    si1=str2double(proc2sPars{5});
end

% Obtain values for data size
td2=str2double(acqusPars{1});
xdim2=str2double(procsPars{6});
si2=str2double(procsPars{7});

% DK edits 5/30/24: Added ability to detect big- or little-endian, plus the
% data type
switch filestr
    case {'fid','ser'} %use aqpath to datafile, plus BYTORDA and DTYPA from acqus
        byte_order=str2double(acqusPars{2});
        data_type=str2double(acqusPars{3});
        fname=fullfile(aqpath,filestr);
    case {'1r','2rr'} %use procpath to datafile, plus BYTORDP and DTYPP from procs
        byte_order=str2double(procsPars{1});
        data_type=str2double(procsPars{2});
        fname=fullfile(procpath,filestr);
end

if byte_order==1
    bigendian=true;
elseif byte_order==0
    bigendian=false;
end
if data_type==2
    dtstr='float64';
elseif data_type==0
    dtstr='int32';
end
% end DK edits
% td1pow2=2^ceil(log2(td1));
% fprintf(1,'Reading %i(%i) x %i matrix\n',td1,td1pow2,td2);
%fprintf(1,'xdim1=%i, si1=%i, td1=%i, offset1=%i, sw_p1=%i, sf1=%i\nxdim2=%i,si2=%i, offset2=%i, sw_p2=%i, sf2=%i\n', xdim1, si1, td1, offset1, sw_p1, sf1, xdim2, si2, offset2, sw_p2, sf2);
% Beware: no zero index available, so always index +1
%d=d(2:length(d));
%p=p(2:length(p));
%cnst=cnst(2:length(cnst));
% DK edits 2/19/21: Read in variable parameter lists if present in directory
% DK edits 3/5/21: Added vclist to these
% DK edits 3/30/21: Added spectral width sw to these
% DK edits 5/30/24: Updated to read in for Topspin 4, where lists are in
% folder lists/['va','vd','vp','vc']/
listnames = {'valist','vplist','vdlist','vclist'};
for ii = 1:length(listnames)
    listname = listnames{ii};
    if exist(fullfile(aqpath,listname),'file')
        otherdata.(listname) = ReadVlist(aqpath,listname);
    elseif exist(fullfile(aqpath,'lists',strtok(listname,'list')),'dir')  
        %for newer TopSpin versions, where the file is in its own directory
        fnam = dir(fullfile(aqpath,'lists',strtok(listname,'list')));
        fnam = fnam(3).name;
        otherdata.(listname) = ReadVlist(fullfile(aqpath,'lists',...
            strtok(listname,'list')),fnam);
    else
        otherdata.(listname) = [];
    end
end

% Generate otherdata struct from extracted parameters (NOTE: str2num()
% required for text containing multiple numerical values, to convert 
% properly to an array!)
otherdata.d=str2num(acqusPars{4});
otherdata.p=str2num(acqusPars{5});
otherdata.cnst=str2num(acqusPars{6});
otherdata.sw=str2double(acqusPars{7});
otherdata.sw_p=str2double(acqusPars{8});
otherdata.sfo1=str2double(acqusPars{9});
if TSver<3
    otherdata.pldb=str2num(acqusPars{10});
else
    otherdata.plw=str2num(acqusPars{10});
end

% end DK edits
%flev=fopen(sprintf('%s/%s', procpath, 'level'), 'r');
%[c,count]=fread(flev, inf, 'int32');
%disp('Loading Contour Levels:');
%fclose(flev);
%c=c(3:length(c));
%c
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reading in data file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% will flip later, because of read filling first the columns
% ser=zeros(td2, td1);
%clear(ser);

% disp(['Opening ' fname]);
% DK edits 5/30/24: Updated to use bigendian and dtstr for byte order and
% data type, respectively
if bigendian
  fspec=fopen(fname, 'r','b');
else
  fspec=fopen(fname, 'r','l');
end

switch filestr
    case 'fid'
        data=fread(fspec,inf,dtstr);
        data=data(1:2:(length(data)-1))+1i*data(2:2:(length(data)));
        data=transpose(data);
    case 'ser'
        data=fread(fspec,[td2 td1],dtstr);
        data=transpose(data);
        data=data(:,1:2:(td2-1))+1i*data(:,2:2:td2);
    case '1r'
        data=fread(fspec,inf,dtstr);
        data=transpose(data);
    case '2rr' %there's probably a better way of doing this....
        data=zeros(si1,si2);
        for smxrow=1:xdim1:(si1-xdim1+1)
            for smxcol=1:xdim2:(si2-xdim2+1)
                dummy=fread(fspec,[xdim2 xdim1],dtstr);
                data(smxrow:(smxrow+xdim1-1),smxcol:(smxcol+xdim2-1))=dummy';
            end
        end
        data=data(1:td1,:); %eliminate zeros
end
fclose(fspec);
% % Scale stimmt evtl. um einen Datenpunkt nicht!
% scaleF2=offset2-(sw_p2/si2/sf2)*(1:si2);
% scaleF1=offset1-(sw_p1/si1/sf1)*(1:si1);
% % otherdata=[sf1 sf2];
end