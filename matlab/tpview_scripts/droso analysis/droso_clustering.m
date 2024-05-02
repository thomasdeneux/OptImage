function label = droso_clustering(V) 
% function f(V) 
% function label = f('label') 

if nargin==0, V=evalin('base','V'); end
if ischar(V) && strcmp(V,'label')
    label = 'clustering';
    return
end

% Write script code below

%% Show timer
watch = fn_watch(V.hf); %#ok<NASGU>

%% Temporal sub-selection
t_idx = [];
timeselset = V.a4d.SIt.selection.t;
if ~isempty(timeselset) && ~isempty(timeselset.set)
    answer = questdlg('Use temporal selection?');
    switch answer
        case {'' 'Cancel'}
            disp interrupted
            return
        case 'Yes'
            t_idx = timeselset.set.dataind;
    end
end
if t_idx
    nfr = length(t_idx);
else
    nfr = V.nfr;
end

%% Get data

disp 'get data'
x = V.dataop;
if t_idx
    x = x(:,:,t_idx);
end
x = fn_float(full(x));

%% Compute correlations

disp 'correlations'
corr = corrcoef(fn_imvect(x)');

%% Clustering

disp 'clustering'
tree = linkage(corr,'weighted');

%% Display

X=fn_control(struct('n__cluster',{100 'stepper 1 1 Inf 1'}),@(s)showclusters(s.n__cluster));
showclusters(X.n__cluster)
clear watch

function showclusters(nclust)

    K = cluster(tree,'maxclust',nclust);
    K = reshape(K,[size(x,1) size(x,2)]);

    fn_figure('Clustering')
    a=fourd(K,'2d',V.a4d.G);
    a.D.cmap = rand(nclust*2,3);

end

end


