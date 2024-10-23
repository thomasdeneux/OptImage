function label = script_loadAamirDataAsk(V) 
% function f(V) 
% function label = f('label') 

if nargin==0, V=evalin('base','V'); end
if ischar(V) && strcmp(V,'label')
    label = 'load Aamir data (always ask for folder)';
    return
end

% Write script code below
script_loadAamirData(V, true)