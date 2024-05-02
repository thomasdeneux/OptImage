function label = automatic_region_level(V) 
% function f(V) 
% function label = f('label') 

if nargin==0, V=evalin('base','V'); end
if ischar(V) && strcmp(V,'label')
    label = 'Create contour region from Picture1';
    return
end

% Write script code below

%% Check whether a spatial filtering was applied

T = V.trial; % current trial
opdef = T.opdef;
operations = {opdef.name};
if ~any(strcmp(operations,'filter') & [opdef.active])
    answer = questdlg('It is recommended to apply some spatial smoothing to the data (data op. > FILTER) before automatic region selection. What do you want to do?','', ...
        'Proceed anyway','Cancel','Cancel');
    if ~strcmp(answer,'Proceed anyway'), return, end
end

%% Get the image displayed in 'Picture1'

a4d = V.a4d;        % structure with information from all the interactive displays
SI = a4d.SI1;       % structure with information of what is displayed in Picture1
image = SI.slice.data;  % the displayed image

%% Parameters

s = fn_structedit( ...
    'mode',     {'maximum'      {'minimum' 'maximum'}   'Look for'}, ...
    'level',    {50             'slider 0 100 1'        'Which % of minimum/maximum?'}, ...
    'baseline', {'baseline = 0' {'baseline = 0', 'baseline = 1', 'minimum pixel', 'maximum pixel', 'pixel average', 'Select a region'} 'Baseline level'}, ...
    'split',    {false          'logical'               'Split multiple components'}); 
if isempty(s), disp interrupted, return, end
level = s.level/100;

%% Get baseline from a region?

switch s.baseline
    case 'baseline = 0'
        baseline = 0;
    case 'baseline = 1'
        baseline = 1;
    case 'minimum pixel'
        baseline = min(image(:));
    case 'maximum pixel'
        baseline = max(image(:));
    case 'pixel average'
        baseline = mean(image(:));
    case 'Select a region'
        hf = figure('name','Please select a region','windowstyle','modal');
        imagesc(image')
        poly = fn_mouse('poly');
        close(hf)
        baselinemask = fn_poly2mask(poly,size(image));
        baseline = mean(image(baselinemask));
end


%% Define the region

% mask of the pixels whose value is below (global maximum * level)
switch s.mode
    case 'minimum'
        m = min(image(:));
        if m>=baseline, waitfor(errordlg('no minimum below baseline found in image: aborting')), return, end
        mask = (image-baseline) <= (m-baseline)*level;
        if all(mask(:)), waitfor(warndlg('all pixels have value below threshold')), end
    case 'maximum'
        M = max(image(:));
        if M<=baseline, waitfor(errordlg('no maximum above baseline found in image: aborting')), return, end
        mask = (image-baseline) >= (M-baseline)*level;
        if all(mask(:)), waitfor(warndlg('all pixels have value above threshold')), end
end

% convert to polygon
poly = fn_mask2poly(mask)';

%% Split polygon if requested

if s.split
    sep = find(isnan(poly(1,:)));
    nroi = length(sep) + 1;
    okroi = true(1,nroi);
    sep = [0 sep size(poly,2)+1];
    sel = cell(1, nroi);
    for k = 1:nroi
        sel{k} = poly(:,sep(k)+1:sep(k+1)-1);
        okroi(k) = ~isempty(sel{k}); % in case there are two NaNs in a row
    end
    nroi = sum(okroi);
    sel = sel(okroi);
else
    nroi = 1;
    sel = {poly};
end

%% Add this region to OptImage

% make a 'selectionND' object out of this polygon
for k = 1:nroi
    sel{k} = selectionND('poly2D',sel{k});
end
sel = [sel{:}];

% add it to OptImage regions
SI = V.a4d.SI1;     % structure with information of what is displayed in Picture1
SI.updateselection('new',[],sel)



