function label = tpv_inspeck(V) 
% function tpv_inspeck(V) 
% function label = tpv_inspeck('label') 
%---
% calls tps_inspeck
%
% See also tps_inspeck

if nargin==0, V=evalin('base','V'); end

if ischar(V) && strcmp(V,'label')
    label = 'Inspeck beads';
    return
end

assignin('base','ff',strvcat(V.content.trials.file)) %#ok<VCAT>
evalin('base','tps_inspeck') % script, not function!