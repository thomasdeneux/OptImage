function label = droso_correlations(V) 
% function f(V) 
% function label = f('label') 

if nargin==0, V=evalin('base','V'); end
if ischar(V) && strcmp(V,'label')
    label = 'correlations (select seed point)';
    return
end

% Write script code below

%% Show timer
watch = fn_watch(V.hf); %#ok<NASGU>

%% Get data
x = V.dataop;
x = fn_imvect(x);

%% Temporal sub-selection
timeselset = V.a4d.SIt.selection.t;
if ~isempty(timeselset) && ~isempty(timeselset.set)
    answer = questdlg('Use temporal selection?');
    switch answer
        case {'' 'Cancel'}
            disp interrupted
            return
        case 'Yes'
            x = x(:,timeselset.set.dataind);
    end
end

%% Transpose 
x = full(x)';

%% Compute correlations

% if isfield(V.trial.user,'droso_correlations')
%     C = V.trial.user.droso_correlations;
% else
disp 'compute correlations'
C = corrcoef(x);
C = reshape(C,[V.nx V.ny 1 1 1 1 V.nx V.ny]);
%     V.trial.user.droso_correlations = C;  % store to avoid re-computing
% end

%% Display

hf = figure('numbertitle','off','integerhandle','off','handlevisibility','off','name','Correlations - with seed point');
G = V.a4d.G;
G.labels(7:8) = G.labels(1:2);
G.units(7:8) = G.units(1:2);
a = fourd(C,V.a4d.G,'in',hf,'2d','proj',[1 2]);
a.D.usercallback = @(D)select_corr_pixel(a.G,hf);
clear watch

%---
function select_corr_pixel(G,hf)

if strcmp(get(hf,'SelectionType'),'open')
    G.ijkl2(7:8) = G.ijkl(1:2);
end

