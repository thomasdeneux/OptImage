function label = droso_correlations_regions(V) 
% function f(V) 
% function label = f('label') 

if nargin==0, V=evalin('base','V'); end
if ischar(V) && strcmp(V,'label')
    label = 'correlations (with regions)';
    return
end

% Write script code below

%% Show timer
watch = fn_watch(V.hf); %#ok<NASGU>

%% Get data
data = V.dataop;
data = fn_imvect(data);

% %% Get the data
% % Current movie or all movies?
% answer = questdlg('On which movies do you want to perform the correlation analysis?','Question', ...
%     'Current Movie','All Movies','Cancel','All Movies');
% if strcmp(answer,'Cancel'), return, end
% T = V.content.trials;
% if strcmp(answer,'Current Movie')
%     T = T(V.ktrial);
% else
%     T = T([T.status]=='n');
% end
% data = cat(4,T.data);
% %%
% data = fn_reshapepermute(data,{[1 2] [3 4]});

%% Temporal sub-selection
timeselset = V.a4d.SIt.selection.t;
if ~isempty(timeselset) && ~isempty(timeselset.set)
    answer = questdlg('Use temporal selection?');
    switch answer
        case {'' 'Cancel'}
            disp interrupted
            return
        case 'Yes'
            data = data(:,timeselset.set.dataind);
    end
end

%% Convert to z-score
data = fn_normalize(data,2,'zscore');

%% Do one computation per ROI

nsel = V.content.nsel;
if nsel==0
    errordlg('Please select at least one Region of Interest')
    return
end

nrow = ceil(sqrt(nsel));
ncol = ceil(nsel/nrow);

fn_figure('Correlations - with regions')
for idxsel = 1:nsel
    ROI = V.content.signal.x(1,idxsel).sel;
    pixelidx = ROI.dataind;
    
    %% Compute the correlation to the reference
    dataref = mean(data(pixelidx,:),1);
    dataref = dataref / sqrt(mean(dataref.^2));
    C = mean(fn_mult(data,dataref),2);
    C = reshape(C,[V.nx V.ny]);
    
    %% Display
    subplot(nrow,ncol,idxsel)
    imagesc(full(C)')
    fn_drawpoly(ROI)
    axis image
    title(['ROI ' num2str(idxsel)])
end




    
