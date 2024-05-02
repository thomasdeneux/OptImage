function label = droso_correlations_inter_regions(V) 
% function f(V) 
% function label = f('label') 

if nargin==0, V=evalin('base','V'); end
if ischar(V) && strcmp(V,'label')
    label = 'correlations (inter regions)';
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
nt = size(data,2);

%% Get average signal for each ROI

nsel = V.content.nsel;
if nsel==0
    errordlg('Please select at least two Region of Interest')
    return
end

dataavg = zeros(nt, nsel);

for idxsel = 1:nsel
    ROI = V.content.signal.x(1,idxsel).sel;
    pixelidx = ROI.dataind;
    
    %% Compute the correlation to the reference
    dataavg(:, idxsel) = full(mean(data(pixelidx,:),1));
end

%% Compute correlations between regions


C = corrcoef(dataavg);

%% Display

fn_figure('Correlations - inter regions')
imagesc(C)
axis image
colormap jet
colorbar 

fn_figure('Correlations - inter regions (graph)')
set(gca,'ydir','reverse')
axis image
axis([0 V.nx 0 V.ny])

for idxsel = 1:nsel
    ROI = V.content.signal.x(1,idxsel).sel.convert('poly2D');
    center = mean(ROI.poly.points,2);
    fn_drawpoly(ROI,'color','k')
    for idx2 = 1:idxsel-1
        ROI2 = V.content.signal.x(1,idx2).sel.convert('poly2D');
        center2 = mean(ROI2.poly.points,2);
        c = C(idxsel,idx2);
        hl = fn_drawpoly([center(:) center2(:)],'linewidth',abs(c)*5);
        if c > 0
            set(hl,'color','b')
        else
            set(hl,'color','r')
        end
    end
end





    
