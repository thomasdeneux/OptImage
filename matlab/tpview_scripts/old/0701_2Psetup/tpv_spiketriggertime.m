function label = tpv_spiketriggertime(V)
% function tpv_spiketriggertime(V)
% function label = tpv_spiketriggertime('label')
%---
% do spike-trigger aceraging and display

if nargin==0, V=evalin('base','V'); end

if ischar(V) && strcmp(V,'label')
    label = 'spike-trigger time courses';
    return
end

persistent par
trig = gettriggers(V.content);
F = fieldnames(trig);
spec.trig = F;
spec.nframebefore = 'slider 1 200 1';
spec.nframeafter  = 'slider 1 200 1';
if isempty(par)
    par.trig = F{1};
    par.nframebefore = 3;
    par.nframeafter  = 15;
    par.oddeven = false;
end
if ~fn_ismemberstr(par.trig,F), par.trig = F{1}; end
par = fn_structedit(par,spec);
if isempty(par), return, end
assignin('base','par',par)

% data
a = V.content.signal.data;
a = reshape(a,V.nfr*V.ntrial,V.content.nsel);

% do the triggering
idx = trig.(par.trig);
idx(idx<=par.nframebefore | idx>=V.nfr*V.ntrial-par.nframeafter)=[];
ntrig = length(idx);
multidx = fn_add(idx(:),-par.nframebefore:par.nframeafter);
b = a(multidx,:);
b = reshape(b,ntrig,par.nframeafter+par.nframebefore+1,V.content.nsel);
avg = squeeze(mean(b,1));

% display
figure(1), clf
tidx = (-par.nframebefore:par.nframeafter)*V.dt;
plot(tidx,avg)
