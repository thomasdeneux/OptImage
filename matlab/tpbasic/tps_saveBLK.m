function tps_saveBLK(T,fname)
% function tps_saveBLK(T[,fname])
%---
% Save tps_trial in VDAQ format.

if nargin==0, help tps_saveBLK, return, end

ntrial = length(T);

% check header information (do this first to detect possible errors in the
% headers)
for k=1:ntrial
    if k==1
        header = getheader(T(k));
        header(ntrial) = header;
    else
        header(k) = getheader(T(k));
    end
end

% how to name BLK files?
filemethod = '';
if nargin<2
    answer = questdlg('How do you want BLK file names to be set?','User input', ...
        'Choose a folder and build on original file names','Choose a new base name','Choose a folder and build on original file names');
    if isempty(answer), return, end
    switch answer
        case 'Choose a folder and build on original file names'
            d = fn_getdir('Select folder');
            if ~d, return, end
            filemethod = 'folder';
        case 'Choose a new base name'
            fbase = fn_savefile('*','Select base name');
            filemethod = 'base';
    end
elseif iscell(fname)
    files = fname;
elseif isscalar(T)
    files = {fname};
else
    [d base ext] = fileparts(fname);
    switch ext
        case '.BLK'
            fbase = [d '/' base];
            filemethod = 'base';
        case ''
            prop1 = ['Create a folder ''' base ''' and build on original file names'];
            prop2 = ['Use ''' base ''' as file base'];
            answer = questdlg('How do you want BLK file names to be set?','User input', ...
                prop1,prop2,prop2);
            if isempty(answer), return, end
            if strcmp(answer,prop1)
                d = [d '/' base];
                filemethod = 'folder';
            else
                fbase = [d '/' base];
                filemethod = 'base';
            end
        otherwise
            error 'extension other than ''.BLK'' not handled yet'
    end
end
switch filemethod
    case ''
        % files already defined
    case 'folder'
        s = fn_structedit('prefix','','inbetween','','suffix','','keepextension',false);
        if isempty(s), return, end
        if ~isempty(s.inbetween)
            if s.inbetween(1)=='_', s.inbetween(1)=[]; end
            if s.inbetween(end)~='_', s.inbetween(end+1)='_'; end
        end
        files = {T.file};
        for k=1:ntrial
            [fk extk] = fn_fileparts(files{k},'base','ext','strict');
            if ~s.keepextension, extk = '.BLK'; end
            u = find(fk=='_',1,'last'); if isempty(u), u=0; end
            files{k} = [d '/' s.prefix  fk(1:u) s.inbetween fk(u+1:end) s.suffix extk];
        end
        mkdir(d)
    case 'base'
        files = cell(1,ntrial);
        ndigit = 1+floor(log10(ntrial));
        format = ['_%.' num2str(ndigit) 'i.BLK'];
        for k=1:ntrial
            files{k} = [fbase sprintf(format,k)];
        end
end

% Save
fn_progress('saving file',ntrial)
for k=1:ntrial
    fn_progress(k)
    fid = fopen(files{k},'w');
    writeheader(header(k),fid)
    writedata(T(k),header(k),fid)
end

%---
function s = getheader(T)

s = T.fullinfo;
if ~isfield(s,'checksum_header')
    error 'tps_trial object does not come originally from VDAQ. Some header information will be missing.'
end

% set sizes
s.framewidth = T.nx;
s.frameheight = T.ny;
s.nframesperstim = T.nfr;
s.nstimuli = T.nc;

% set data type to single-precision floating numbers
s.datatype = 14;
s.sizeof = 4;
s.framesize = s.framewidth * s.frameheight * s.sizeof;
s.filesize = s.lenheader + s.framesize * s.nframesperstim * s.nstimuli;

%---
function writedata(T,s,fid)

data = T.dataop;
fseek(fid,s.lenheader,'bof');
fwrite(fid,data,'float32');

%---
function writeheader(s,fid)

%%%%%%%%%%%%%%%%%% BEGIN DEFINITIONS OF VARIABLES %%%%%%%%%%%%%%%%%%

% Note: we save char as uint8, because some special Matlab char type seems
% strange with values actually from 0 to 65535, all those <256 being
% written by fwrite with only one byte, but all those >255 being written
% with 2 bytes!!

% DATA INTEGRITY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fwrite(fid,s.filesize,'int32');
fwrite(fid,s.checksum_header,'int32');
% beginning with the lLen Header field

fwrite(fid,s.checksum_data,'int32');

% COMMON TO ALL DATA FILES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fwrite(fid,s.lenheader,'int32');
fwrite(fid,s.versionid,'float32');
fwrite(fid,s.filetype,'int32');
% RAWBLOCK_FILE          (11)
% DCBLOCK_FILE           (12)
% SUM_FILE               (13)

fwrite(fid,s.filesubtype,'int32');
% FROM_VDAQ              (11)
% FROM_ORA               (12)
% FROM_DYEDAQ            (13)

fwrite(fid,s.datatype,'int32');
% DAT_UCHAR     (11)
% DAT_USHORT    (12)
% DAT_LONG      (13)
% DAT_FLOAT     (14)

fwrite(fid,s.sizeof,'int32');
% e.g. sizeof(long), sizeof(float)

fwrite(fid,s.framewidth,'int32');
fwrite(fid,s.frameheight,'int32');
fwrite(fid,s.nframesperstim,'int32');
fwrite(fid,s.nstimuli,'int32');
fwrite(fid,s.initialxbinfactor,'int32');
% from data acquisition
fwrite(fid,s.initialybinfactor,'int32');
% from data acquisition

fwrite(fid,s.xbinfactor,'int32');
% this file
fwrite(fid,s.ybinfactor,'int32');
% this file

fwrite(fid,s.username,'uint8'); if length(s.username)~=32, error 'bad header', end
fwrite(fid,s.recordingdate,'uint8'); if length(s.recordingdate)~=16, error 'bad header', end
fwrite(fid,s.x1roi,'int32');
fwrite(fid,s.y1roi,'int32');
fwrite(fid,s.x2roi,'int32');
fwrite(fid,s.y2roi,'int32');

% LOCATE DATA AND REF FRAMES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fwrite(fid,s.stimoffs,'int32');
fwrite(fid,s.stimsize,'int32');
fwrite(fid,s.frameoffs,'int32');
fwrite(fid,s.framesize,'int32');
fwrite(fid,s.refoffs,'int32');
fwrite(fid,s.refsize,'int32');
fwrite(fid,s.refwidth,'int32');
fwrite(fid,s.refheight,'int32');

% Common to data files that have undergone some form of
% "compression" or "summing"
% i.e. The data in the current file may be the result of
%      having summed blocks 'a'-'f', frames 1-7
fwrite(fid,s.whichblocks,'uint16'); if length(s.whichblocks)~=16, error 'bad header', end
fwrite(fid,s.whichframes,'uint16'); if length(s.whichframes)~=16, error 'bad header', end

% DATA ANALYSIS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fwrite(fid,s.loclip,'int32');
fwrite(fid,s.hiclip,'int32');
fwrite(fid,s.lopass,'int32');
fwrite(fid,s.hipass,'int32');

fwrite(fid,s.operationsperformed,'uint8'); if length(s.operationsperformed)~=64, error 'bad header', end

% ORA-SPECIFIC%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fwrite(fid,s.magnification,'float32');
fwrite(fid,s.gain,'uint16');
fwrite(fid,s.wavelength,'uint16');
fwrite(fid,s.exposuretime,'int32');
fwrite(fid,s.nrepetitions,'int32');
fwrite(fid,s.acquisitiondelay,'int32');
fwrite(fid,s.interstiminterval,'int32');
fwrite(fid,s.creationdate,'uint8'); if length(s.creationdate)~=16, error 'bad header', end
fwrite(fid,s.datafilename,'uint8'); if length(s.datafilename)~=64, error 'bad header', end
fwrite(fid,s.orareserved,'uint8'); if length(s.orareserved)~=256, error 'bad header', end


if s.filesubtype == 13,   %it's dyedaq file
  
  %  OIHEADER.H
  %  last revised 4.5.97 by Chaipi Wijnbergen for DyeDaq
  %  
  %  DyeDaq-specific
  fwrite(fid,s.includesrefframe , 'int32');     % 0 or 1
  fwrite(fid,s.temp , 'uint8'); if length(s.temp)~=128, error 'bad header', end
  %s.listofstimuli=temp(1:max(find(temp~=0)))';  % up to first non-zero stimulus
  fwrite(fid,s.ntrials , 'int32');
  fwrite(fid,s.scalefactor , 'int32');          % bin * trials
  fwrite(fid,s.cameragain , 'short');         % shcameragain        1,   2,   5,  10
  fwrite(fid,s.ampgain , 'short');            % amp gain            1,   4,  10,  16,
                                             %                    40,  64, 100, 160,
                                             %                    400,1000
  fwrite(fid,s.samplingrate, 'short');       % sampling rate (1/x)
                                             %                     1,   2,   4,   8,
                                             %                     16,  32,  64, 128,
                                             %                     256, 512,1024,2048
  fwrite(fid,s.average ,'short');            % average             1,   2,   4,   8,
                                             %                    16,  32,  64, 128
  fwrite(fid,s.exposuretime , 'short');       % exposure time       1,   2,   4,   8,
                                             %                    16,  32,  64, 128,
                                             %                    256, 512,1024,2048
  fwrite(fid,s.samplingaverage, 'short');    % sampling average    1,   2,   4,   8,
                                             %                    16,  32,  64, 128
  fwrite(fid,s.presentaverage , 'short');
  fwrite(fid,s.framesperstim , 'short');
  fwrite(fid,s.trialsperblock , 'short');
  fwrite(fid,s.sizeofanalogbufferinframes , 'short');
  fwrite(fid,s.cameratrials, 'short');
  fwrite(fid,s.filler,'uint8'); if length(s.filler)~=106, error 'bad header', end;
  
  fwrite(fid,s.dyedaqreserved,'uint8'); if length(s.dyedaqreserved)~=256, error 'bad header', end;
else   % it's not dyedaq specific

  % VDAQ-SPECIFIC%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  fwrite(fid,s.includesrefframe,'int32');
  fwrite(fid,s.listofstimuli,'uint8'); if length(s.listofstimuli)~=256, error 'bad header', end;
  fwrite(fid,s.nvideoframesperdataframe,'int32');
  fwrite(fid,s.ntrials,'int32');
  fwrite(fid,s.scalefactor,'int32');
  fwrite(fid,s.meanampgain,'float32');
  fwrite(fid,s.meanampdc,'float32');
  fwrite(fid,s.vdaqreserved,'uint8'); if length(s.vdaqreserved)~=256, error 'bad header', end;
end    % end of VDAQ specific

% USER-DEFINED%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  fwrite(fid,s.user,'uint8'); if length(s.user)~=256, error 'bad header', end;

% COMMENT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fwrite(fid,s.comment,'uint8'); if length(s.comment)~=256, error 'bad header', end;
fwrite(fid,s.refscalefactor , 'int32');          % bin * trials for reference

%%%%%%%%%%%%%%%%%%  END DEFINITIONS OF VARIABLES  %%%%%%%%%%%%%%%%%%

% Doing ftell(fid) here results in the same as when reading a file with
% oi_headBLK. But there is the same error (see end of oi_headBLK code),
% which is that it is 4 bytes more than s.lenheader!!!


