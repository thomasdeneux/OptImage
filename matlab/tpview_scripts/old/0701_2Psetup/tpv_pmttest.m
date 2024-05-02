function label = tpv_pmttest(V) 
% function tpv_pmttest(V) 
% function label = tpv_pmttest('label') 
%---
% calls tps_pmtcenter
%
% See also tps_pmtcenter

if nargin==0, V=evalin('base','V'); end

if ischar(V) && strcmp(V,'label')
    label = 'PMT test';
    return
end

tps_pmtcenter(strvcat(V.content.trials.file)) %#ok<VCAT>