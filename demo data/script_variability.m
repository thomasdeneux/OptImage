function label = script_variability(V) 
% function f(V) 
% function label = f('label') 

if nargin==0, V=evalin('base','V'); end
if ischar(V) && strcmp(V,'label')
    label = 'variability analysis';
    return
end

% Write script code below

% get data, compute STD
T = V.content.trials;
stimtable = T(1).stimdetails;
okcond = strcmp({stimtable.type},'stim');
okreg = ([T.status]=='n');
x = cat(5,T(okreg).dataop);
x = x(:,:,:,okcond,:);
x = x(:,:,:,:);
x = std(x,0,4);


% display result in optimage
setdata(V,x,'add')