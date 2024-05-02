function label = script_functional_connectivity(V) 
% function f(V) 
% function label = f('label') 

if nargin==0, V=evalin('base','V'); end
if ischar(V) && strcmp(V,'label')
    label = 'Functional Connectivity';
    return
end

% Write script code below

%% Get the data
% Current movie or all movies?
answer = questdlg('On which movies do you want to perform the correlation analysis?','Question', ...
    'Current Movie','All Movies','Cancel','All Movies');
if strcmp(answer,'Cancel'), return, end
T = V.content.trials;
if strcmp(answer,'Current Movie')
    T = T(V.ktrial);
else
    T = T([T.status]=='n');
end
data = cat(4,T.data);
%%
data = fn_normalize(data,3,'zscore');
data = fn_reshapepermute(data,{[1 2] [3 4]});

%% Temporal region only?


%% Do one computation per ROI

nsel = V.content.nsel;
if nsel==0
    errordlg('Please select at least one Region of Interest')
    return
end

nrow = ceil(sqrt(nsel));
ncol = ceil(nsel/nrow);

fn_figure('Functional connectivity')
for idxsel = 1:nsel
    ROI = V.content.signal.x(1,idxsel).sel;
    pixelidx = ROI.dataind;
    
    %% Compute the correlation to the reference
    dataref = mean(data(pixelidx,:),1);
    dataref = dataref / sqrt(mean(dataref.^2));
    C = mean(fn_mult(data,dataref),2);
    C = reshape(C,[T(1).nx T(1).ny]);
    
    %% Display
    subplot(nrow,ncol,idxsel)
    imagesc(full(C)')
    fn_drawpoly(ROI)
    axis image
    title(['ROI ' num2str(idxsel)])
end

%% Clustering

[nx ny ntt] = deal(T(1).nx,T(1).ny,size(data,2));
databin = reshape(data,[nx ny ntt]);
databin = fn_bin(databin,4);
[nxbin nybin ntt] = size(databin);
databin = reshape(databin,[nxbin*nybin ntt]);
databin = full(databin);
%%
Z = linkage(databin,'average','correlation');
%
T = cluster(Z,'MaxClust',12);
T = reshape(T,nxbin,nybin);

%% Display
fn_figure('Clustering')
Tlarge = fn_enlarge(T,[nx ny]);
Tlarge = fn_clip(Tlarge,[1 12],fn_colorset('plot12'));
fourd('2dcol',Tlarge,V.a4d.G)
% imagesc(T')
axis image



    
