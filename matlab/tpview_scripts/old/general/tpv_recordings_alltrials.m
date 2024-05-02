function label = tpv_recordings_alltrials(V) 
% function f(V) 
% function label = f('label') 

if ischar(V) && strcmp(V,'label')
    label = 'Recordings for all trials';
    return
end

if nargin==0, V=evalin('base','V'); end

T = V.content.trials;
T = T([T.status]=='n');


r = cat(3,T.recording);
recnames = {r(1,:,1).name};
ncond = size(r,1);

s = fn_structedit('recording',{recnames{min(2,length(recnames))} recnames}, ...
    'conditions',{true(1,ncond) ['multcheck' fn_num2str(0:ncond-1,'cell')]});

okrec = strcmp(recnames,s.recording);
r = r(s.conditions,okrec,:);
data = [r.signal];

dt = r(1).dt;
F = V.a4d.F;

hf = figure('numbertitle','off','name',[get(V.hf,'name') ' - ' upper(recnames{okrec})]);
a=fourd(F,data,'plot','in',hf);
a.D.linecol = [0 0 0];
