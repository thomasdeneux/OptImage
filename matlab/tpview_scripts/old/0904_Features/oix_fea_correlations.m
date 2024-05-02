function label = oix_fea_correlations(V) 
% function oix_fea_correlations(V) 
% function label = oix_fea_correlations('label') 

if nargin==0, V=evalin('base','V'); end
if ischar(V) && strcmp(V,'label')
    label = 'Correlations';
    return
end

fea_corrgui(V)
