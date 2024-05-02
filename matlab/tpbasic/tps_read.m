function T = tps_read(varargin)
% function trials = tps_read([f][,'header'])
%---
% reads 2-photon data file (CFD, TIF, MPD, 2PLSM, XML or MAT)
% and returns a (multiple if several trials) tps_trial object.
% This function is actually just a link to the appropriate tps_trial
% methods, except that it forces the user to choose a file.

% Input
f = [];
headerflag = false;
for i=1:nargin
    if strcmp(varargin{i},'header')
        headerflag = true;
    else
        f = varargin{i};
    end
end
if isempty(f)
    f = fn_getfile;
end

% Call tps_trial methods
if headerflag
    T = tps_trial.readfileheader(f);
else
    T = tps_trial.readfile(f);
end


