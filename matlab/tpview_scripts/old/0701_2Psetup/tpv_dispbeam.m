function label = tpv_dispbeam(V) 
% function tpv_dispbeam(V) 
% function label = tpv_dispbeam('label') 
%---
% calls tps_dispbeam
%
% See also tps_dispbeam

if nargin==0, V=evalin('base','V'); end

if ischar(V) && strcmp(V,'label')
    label = 'laser beam size';
    return
end

ff = V.file; %#ok<NASGU>
tps_dispbeam
