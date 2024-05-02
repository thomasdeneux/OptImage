function x = tps_readacqsettingsxml(fname)
% function x = tps_readacqsettingsxml(fname)
%---
% first reads xml file and convert its content to a Matlab structure by
% using fn_readxmllabview,
%
% then simplifies a few fields:
% - converts the stimulation paradigm xml representation to a Matlab
%   structure
% - converts the recording channels structure to a cell array

if nargin==0
    fname = fn_getfile('*.xml');
end

% store result in a mat file to save the time to read the XML file
fmat = [fname '.mat'];
okmat = exist(fmat,'file');
version = '1.0'; %#ok<NASGU> % if later we want to change the way to read the file, better save the version which was used
if okmat
    load(fmat,'x')
    return
end

% read xml file
x = fn_readxmllabview(fname);

% stimulation paradigm
if isempty(x.stimulationVI), x.stimulationVI = ''; end
str = x.stimulationparadigmxml;
if ~isempty(str)
    % stupid: needs to write in a file in order to read it directly with
    % fn_readxmllabview
    tmp = 'tmp.xml';
    fid = fopen(tmp,'w');
    fprintf(fid,str);
    fclose(fid);
    x.stimulationparadigmxml = fn_readxmllabview(tmp);
    delete(tmp)
end

% recording channels
s = x.recordingchannels;
x.recordingchannels = {s([s.use]).name};

% save in mat file
save(fmat,'x','version')

