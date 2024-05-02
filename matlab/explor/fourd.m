function a = fourd(varargin)
% function a = fourd('optionflag',optionvalue,...)
%---
% General tool to visualize data
%
% Main concept: the data is usually a n-dimentional array, and only a
% p-dimentional (p<n) section of this data will be displayed. Some
% options allow to define this section, and how to display it. All
% individual displays are linked together, allowing the user to change
% dynamically which part of the data is displayed, select regions, etc.
%
% ex: fn4D('data',rand(10,12,20),'type','2dsim')
%     fn4D(1:5,1:7,1:30,rand(5,7,30))
%---
% Options:
%
% MAIN
% - data    ['data' flag is optional] n-dimentional array
%           the array can be preceded by vectors of coordinates (for
%           example: fn4D(x,y,z,data,...)
%           can also be a VideoReader object, in which case type must be
%           '2d'
% - type    ['type' flag is optional] one of:
%           'plot'    time courses display (section of data along the 1st dimension) [default for 1D data]
%           '2dplot'  time courses display (section of data along the 3rd dimension)
%           'list'    navigation inside 1D data using selection inside a list
%           'slider'  navigation inside 1D data using a scroller
%           'player'  navigation inside 1D data with a looping movie
%           '2d'      image display (2-dimensional section) [default for 2D data]
%           '2dcol'   color image display (3-dimensional section)
%           'snake'
%           'frame'
%           '3d'
%           below are special types that result in multiple displays:
%           '2dsim'   image + time courses displays [default for 3D data]
%           '2dlist'  image + list(s) for dimensions >= 3
%           'plotlist'  time courses + list(s) for dimensions >= 3
%           'mult'    multiple display; 'proj' option must be defined
%                     [default for >=4D data, or when 'proj' option value is a cell array]
%
% GEOMETRIC INFO
% - key         key for linking with other display [default=0], use 'newkey' to use an independent key
% - focus       focus object to synchronize to (can synchronize any data, create with function 'focus')
% - geometry    a geometry object to synchronize to (can synchronize only compatible data, create with 'geometry' or 'rotation(F)' where F is a focus)
% - scale       scale of each dimension [default=1]
% - start       starting point of each dimension [default=1 if unit is '', 'px' or 'frame', 0 otherwise]
% - worlddim    world dimension (in 'focus') to which corresponds each data dimension (in 'geometry') [default=1:nddata]
% - mat         full transformation definition: either {scale start-scale worlddim} or an (1+nd)*(1+nd) array
% - dx          scale of dimensions 1 and 2 (assumed to be x and y)
% - dt          scale of dimension 3 (assumed to be time)
% - tt          full time vector
% - proj        dimensions of the data to base the display on
% - dimsplus    additional displayed dimensions (color channel for image, multiple lines for graphs)
%
% GENERAL GRAPHICS
% - in          handle of graphic container (usually figure, uipanel or axes) [default=current figure/axes]
% - labels      labels of each data dimension [default world dimension labels=x, y, time]
% - units       units of each data dimension [default world dimension units=px, px, frame] 
% - decoration
%
% IMAGE DISPLAY
% - clip
% - clipmode
% - cmap
% - selshow
% - shapemode
%
% PLOT DISPLAY
% - autolinepos
% - ystep
% - movelinegroup
%
% SNAKE
% - snakedx
% 
% LIST DISPLAY
% - averageselection    true [=default] to average over selection           
%
% SLIDER
% - layout      'up', 'down', 'right' or 'left'
%
% OTHER DISPLAY
% - ncol        for 'frame' type
% 
% MISC
% - nmovieframe         number of frames when data is a VideoFrame object

% initialization options
opt = struct( ...
    'key',      0, ...
    'focus',    [], ...
    'geometry', [], ...
    'type',     [], ...
    'mat',      [], ...
    'proj',     [], ...
    'dimsplus', 0, ... default value [] will be set later if needed
    'data',     [], ...
    'in',       [], ...
    'clip',     [], ...
    'clipmode', [], ...
    'cmap',     [], ...
    'ystep',    [], ...
    'dt',       1, ...
    'tt',       [], ...
    'dx',       1, ...
    'scale',    [], ...
    'start',    [], ...
    'worlddim', [], ...
    'ncol',     4, ...
    'selshow',          true, ...
    'seldims',          'xy', ...
    'snakedx',          .8, ...
    'autolinepos',      '', ...
    'movelinegroup',    '', ...
    'linecol',          '', ...
    'navigation',       'pan', ...
    'scrollwheel',      'x', ...
    'scaledisplay',     0, ... default value [] will be set later if needed
    'labels',           [], ...
    'shapemode',        'ellipse', ...
    'units',            [], ...
    'layout',           'auto', ...
    'decoration',       [], ...
    'averageselection', true, ...
    'nmovieframe',      [] ...
    );
opt1 = BuildOptions(opt,varargin{:});
opt2 = ArrangeAndGeometry(opt1);

% no data: we have set up the focus and geometry for next times, this is
% fine we can stop here [doesn't work because no pointer is kept on them
% and they are destroyed!]
if isscalar(opt2) && isempty(opt2.data)
    return
end

% display
a = struct('key',[],'F',[],'G',[],'SI',[],'D',[]);
for i=1:length(opt2)
    a(i) = Display(opt2(i)); 
end

% output?
if nargout==0, clear a, end

%---
function opt = BuildOptions(opt,varargin)
% see help of fn4D

k = 1;
numdata = {};
while k<nargin
    a = varargin{k}; k = k+1;
    if isfigoraxeshandle(a)
        opt.in = a;
    elseif isnumeric(a) || islogical(a)
        numdata{end+1} = a;
    elseif isa(a,'VideoReader')
        opt.data = a;
    elseif isa(a,'focus')
        opt.focus = a;
    elseif isa(a,'geometry')
        opt.geometry = a;
    elseif iscell(a)
        if ischar(a{1})
            opt.type = a;
        else
            opt.proj = a;
        end
    elseif ~ischar(a)
        error argument
    else
        a = lower(a);
        switch a
            case {'list','slider','player','plot','2d','image','2dcol','2dplot','snake', ...
                    'frame','frames','3d','2dsim','2dlist','plotlist','mult'}
                opt.type = a;
            case 'tt'
                tt = varargin{k}; k=k+1;
                if ~isvector(tt) || ...
                        any(abs(diff(tt,2))>3*eps(max(abs(tt))))
                    error('wrong time vector definition')
                end
                opt.tt = tt;
                opt.dt = (tt(end)-tt(1))/(length(tt)-1);
            case 'newkey'
                opt.key = rand;
            case 'colormap'
                opt.cmap = varargin{k}; k = k+1;
            case {'worldim' 'worldims' 'worlddims'}
                opt.worlddim = varargin{k}; k = k+1;
            otherwise
                if isfield(opt,a)
                    opt.(a) = varargin{k}; k=k+1;
                else
                    error('unknown option ''%s''',a)
                end
        end
    end
end

% special: syntax fn4D(x,y,z,...,data,...)
if ~isempty(numdata)
    if ~isempty(opt.data), error 'multiple definitions of ''data''', end
    opt.data = numdata{end};
    if ~isscalar(numdata)
        if ~isempty(opt.mat), error 'multiple definition of ''mat''', end
        nmat = length(numdata)-1;
        mat = {zeros(1,nmat) zeros(1,nmat)}; % scale, translation
        for i=1:nmat
            xi = numdata{i};
            ni = length(xi);
            scale = (xi(ni)-xi(1))/(ni-1);
            mat{1}(i) = scale;
            mat{2}(i) = xi(1)-scale;
        end
        opt.mat = mat;
    end
end

%---
function opt = ArrangeAndGeometry(opt)

% type
% (empty type)
if isempty(opt.type)
    if isempty(opt.proj)
        opt.type = fn_switch(ndims(opt.data),1,'plot',2,'2d',3,'2dsim','2dsim');
    elseif isnumeric(opt.proj)
        opt.type = fn_switch(length(opt.proj),1,'plot',2,'2d',3,'2dsim');
    end
end
% ('mult' type: opt.proj must be defined)
if isempty(opt.type) || (ischar(opt.type) && strcmpi(opt.type,'mult'))
    if isempty(opt.proj) || ~iscell(opt.proj), error '''proj'' option should be a cell array here', end
    ndisplay = length(opt.proj);
    opt.type = cell(1,ndisplay);
    for i=1:ndisplay
        opt.type{i} = fn_switch(length(opt.proj{i}),1,'plot',2,'2d',3,'2dsim');
    end
end
% (other multiple displays; note that lists will automatically be added to multiple displays)
if ischar(opt.type)
    switch lower(opt.type)
        case '2dsim'
            opt.type = {'2d' 'plot'};
        case '2dlist'
            opt.type = {'2d'};
        case 'plotlist'
            opt.type = {'plot'};
        case 'list'
            opt.type = {};
    end
    doaddlist = iscell(opt.type);
else
    doaddlist = false;
end

% attempt to correct row vector data as column vector
if isvector(opt.data) && ~any(opt.proj==2)
    opt.data = opt.data(:);
end
siz = size(opt.data);
siz = siz(1:find(siz~=1,1,'last'));
nddata = length(siz);

% proj and dimsplus
if isempty(opt.proj)
    nd = find(size(opt.data)>1,1,'last');
    if isempty(nd), nd = 1; end
    available = 1:nd;
elseif isnumeric(opt.proj)
    available = opt.proj;
    opt.proj = []; % might be redefined!
    %     nd = sum(opt.proj);
elseif iscell(opt.proj)
    %     nd = sum([opt.proj{:}]);
end
if ~iscell(opt.type)
    opt.type = {opt.type}; 
    if iscell(opt.proj), opt.type = repmat(opt.type,[1 length(opt.proj)]); end
end
opt.type = fn_strrep(opt.type,'image','2d','frames','frame');
ndisplay = length(opt.type);
if ndisplay==1 && strcmp(opt.type{1},'2dplot')
    % first 2D dimensions represent 2D space, do not use them
    available(1:2) = [];
    opt.type{1} = 'plot';
elseif any(strcmp(opt.type,'2dplot'))
    error 'type ''2dplot'' cannot be used in the context of multiple displays'
end
if isempty(opt.proj)    
    % We enter here when 'proj' was either not set, or set as a numerical
    % array, but not when it was set as a cell array; in the later case we
    % assumed that everything is already set correctly, including
    % 'dimsplus'.
    % Here 'proj' will be defined or redefined; note that errors will be
    % generated automatically if 'available' is not long enough.
    opt.proj = cell(1,ndisplay);
    if isequal(opt.dimsplus,0)
        % 'dimsplus' has not been set by user
        opt.dimsplus = repmat({0},1,ndisplay);
    else
        % 'dimsplus' has been set by user
        if ndisplay==1
            if ~iscell(opt.dimsplus), opt.dimsplus = {opt.dimsplus}; end
        elseif isempty(opt.dimsplus)
            opt.dimsplus = repmat({[]},1,ndisplay);
        elseif ~iscell(opt.dimsplus)
            error 'setting of ''dimsplus'' cannot be interpreted correctly in the context of multiple displays'
        end
    end
    for k = 1:ndisplay
        switch lower(opt.type{k})
            case 'list'
                opt.proj{k} = available(1); available(1) = [];
                if k==ndisplay
                    % add as many additional lists as available
                    doaddlist = true;
                end
            case {'slider' 'player'}
                opt.proj{k} = available(1); available(1) = [];
            case 'plot'
                opt.proj{k} = available(1); available(1) = [];
                if isequal(opt.dimsplus{k},0) && ~isempty(available) && ~doaddlist && k==ndisplay
                    opt.dimsplus{k} = available(1);
                    % note that we leave this dimension available... anyway
                    % there are no further display since k==ndisplay
                end
            case '2d'
                if length(available)==1, available(2) = available(1)+1; end
                opt.proj{k} = available(1:2); available(1:2) = [];
            case '2dcol'
                opt.type{k} = '2d';
                if length(available)==1, available(2) = available(1)+1; end
                opt.proj{k} = available(1:2); available(1:2) = [];
                if isequal(opt.dimsplus{k},0)
                    s = size(opt.data);
                    iav = find(s(available)==3,1);
                    if isempty(iav), error 'no dimension has size 3 for ''2dcol'' display', end
                    opt.dimsplus{k} = available(iav);
                    available(iav) = [];
                end
                opt.type{k} = '2d';
            case {'snake' 'frame' 'grid' '3d'}
                if length(available)<3, available(end+1:3) = available(end)+(1:3-length(available)); end
                opt.proj{k} = available(1:3); available(1:3) = [];
            otherwise
                error('unknown type ''%s''',opt.type)
        end
        if isequal(opt.dimsplus{k},0), opt.dimsplus{k} = []; end
    end
    % add lists
    if doaddlist
        navailable = length(available);
        opt.type = [opt.type repmat({'list'},[1 navailable])];
        opt.proj = [opt.proj num2cell(available)];
        opt.dimsplus = [opt.dimsplus repmat({[]},[1 navailable])];
    end
else
    % convert some 'proj' dimensions into 'dimsplus'?
    opt.type = strrep(opt.type,'2dcol','2d');
    if isequal(opt.dimsplus,0)
        opt.dimsplus = cell(1,ndisplay);
    elseif ~iscell(opt.dimsplus) || length(opt.dimsplus)~=ndisplay
        error 'problem with dimsplus format'
    end
    for k=1:ndisplay
        switch(opt.type{k})
            case 'plot'
                if length(opt.proj{k})==2
                    opt.dimsplus{k} = opt.proj{k}(2);
                    opt.proj{k} = opt.proj{k}(1);
                end
            case '2d'
                if length(opt.proj{k})==3
                    opt.dimsplus{k} = opt.proj{k}(3);
                    opt.proj{k} = opt.proj{k}(1:2);
                end
        end
    end
end

% other
if isequal(opt.scaledisplay,0) % && ~isempty(opt.mat)
    opt.scaledisplay = 'xbar'; 
end
if isempty(opt.autolinepos)
    if ~iscell(opt.type) && strcmp(opt.type,'plot') && ~isempty(opt.dimsplus)
        opt.autolinepos = 'cat';
    else
        opt.autolinepos = 'sel';
    end
end
if isempty(opt.movelinegroup), opt.movelinegroup = opt.autolinepos; end
if isempty(opt.linecol), opt.linecol = opt.autolinepos; end
if isempty(opt.clipmode)
    opt.clipmode = 'data';
end

% transform the structure in multiple structure)
if iscell(opt.mat) && ~iscell(opt.mat{1}), opt.mat = {opt.mat}; end
if numel(opt.in)>1 && ~iscell(opt.in), opt.in = num2cell(opt.in); end
if iscell(opt.labels), opt.labels = {opt.labels}; end
if iscell(opt.units), opt.units = {opt.units}; end
F = fieldnames(opt);
C = struct2cell(opt);
FC = [row(F); row(C)];
opt = struct(FC{:});

% case of multiple displays -> divide figure
n = length(opt);
if n>1 && all(fn_isemptyc({opt.in}))
    clf
    % bring the lists to the left
    klist = find(strcmp({opt.type},'list'));
    ord = [klist setdiff(1:n,klist)];
    opt = opt(ord);
    % not all displays have the same width
    w = zeros(1,n);
    for i=1:n, w(i) = fn_switch(opt(i).type,'list',1,'2d',4,'plot',4); end
    wt = sum(w);
    cw = [0 cumsum(w)];
    for i=1:n, opt(i).in = subplot(1,wt,cw(i)+(1:w(i))); end
end

% Focus and Geometry
for k = 1:length(opt)
    optk = opt(k);
    if ~isempty(optk.geometry)
        G = optk.geometry;
        if ~isempty(optk.labels), G.labels = optk.labels; end
        if ~isempty(optk.units), G.units = optk.units; end
    else
        if ~isempty(optk.focus)
            F = optk.focus;
        else
            F = focus.find(optk.key);
            % default labels and units for newborn focus
            if isempty(optk.labels) && isempty(optk.units) && F.nd == 0
                F.labels = {'x' 'y' 'time'};
                F.units = {'px' 'px' 'frame'};
            end 
            optk.focus = F;
        end

        % find a rotation object that matches size and mat
        if ~isempty(optk.mat)
            mat = optk.mat;
        else
            % world dims
            if ~isempty(optk.worlddim)
                worlddim = optk.worlddim;
            else
                worlddim = 1:nddata;
            end
            % scale
            if ~isempty(optk.scale)
                scale = optk.scale;
            else
                if isscalar(optk.dx), optk.dx = optk.dx([1 1]); end
                worldscale = [optk.dx optk.dt];
                worldscale(end+1:max(worlddim)) = 1;
                scale = worldscale(optk.worlddim);
            end
            scale(end+1:nddata) = 1;
            % offset (difficult because default values depends on units!!)
            if ~isempty(optk.start)
                start = optk.start;
            else
                start = zeros(1,0);
            end
            worldunits = F.units;
            worldunits(end:nddata) = {''};
            for d = length(start)+1:nddata
                if ismember(worldunits{worlddim(d)}, {'' 'px' 'frame'})
                    start(d) = 1;
                else
                    start(d) = 0;
                end
            end
            start(end+1:nddata) = 1;
            offset = start-scale;
            % mat
            mat = {scale offset worlddim};
        end
        mat = buildMat(mat,nddata,nddata);

        % does focus already have a child rotation that fits the data sizes and
        % transformation?
        % (if yes, take it)
        c = F.getChildren;
        for i=1:length(c)
            dcomp = intersect(find(c{i}.sizes>1),find(siz>1));
            dmat = [1 1+dcomp];
            if all(c{i}.sizes(dcomp)==siz(dcomp)) ...
                    && isequal(c{i}.mat(1:size(mat,1),dmat), mat(:,dmat))
                optk.geometry = c{i};
                disp 'found matching geometry object'
                break
            end
        end
        % (if not, create a new one)
        if isempty(optk.geometry)
            disp 'create new geometry object'
            optk.geometry = rotation(F,'sizes',siz,'mat',mat);
        end
    end
    opt(k) = optk;
end

%---
function a = Display(opt)

% init structure (focus and geometry are already defined)
a = struct('key',opt.key,'F',opt.focus,'G',opt.geometry,'SI',[],'D',[]);

% Projection and Active Display
switch lower(opt.type)
    case 'plot'
        a.SI = projection(a.G,opt.proj,'dimsplus',opt.dimsplus, ...
            'data',opt.data,'decoration',opt.decoration);
        if ndims(opt.data)>2 && length(opt.proj)>=3
            if ~isempty(opt.in)
                if fn_isfigurehandle(opt.in)
                    hf = opt.in;
                else
                    hf = get(opt.in,'parent');
                end
            else
                hf = gcf;
            end
            clf(hf)
        end
        a.D = activedisplayPlot(a.SI,'in',opt.in, ...
            'ystep',opt.ystep, ...
            'clipmode',opt.clipmode,'clip',opt.clip, ...
            'autolinepos',opt.autolinepos, ...
            'movelinegroup',opt.movelinegroup, ...
            'navigation',opt.navigation,'scrollwheel',opt.scrollwheel);
    case '2d'
        if isa(opt.data,'VideoReader')
            if ~isequal(opt.proj,[1 2])
                error 'projection must be [1 2] for a video object'
            end
            if isempty(opt.nmovieframe)
                a.SI = videoprojection(a.G,opt.data);
            else
                a.SI = videoprojection(a.G,opt.data,opt.nmovieframe);
            end
        else
            a.SI = projection(a.G,opt.proj,'dimsplus',opt.dimsplus,'data',opt.data, ...
                'decoration',opt.decoration);
        end
        a.D = activedisplayImage(a.SI,'in',opt.in, ...
            'clipmode',opt.clipmode,'clip',opt.clip, ...
            'selshow',opt.selshow,'seldims',opt.seldims,'shapemode',opt.shapemode, ...
            'scaledisplay',opt.scaledisplay);
        if ~isempty(opt.cmap), a.D.cmap = opt.cmap; end
    case 'snake'
        a.SI = snake2D(a.G,opt.proj,'dx',opt.snakedx,'data',opt.data, ...
            'dimsplus',opt.dimsplus,'decoration',opt.decoration);
        a.D = activedisplayImage(a.SI,'in',opt.in, ...
            'clipmode',opt.clipmode,'clip',opt.clip, ...
            'selshow',opt.selshow,'seldims',opt.seldims);
    case {'3d' 'frame' 'grid'}
        if strcmpi(opt.type,'3d') && isequal(opt.proj,1:3) ...
                && size(opt.data,4)==3 && isempty(opt.dimsplus)
            opt.dimsplus=4; 
        end
        a.SI = projection(a.G,opt.proj,'data',opt.data, ...
            'dimsplus',opt.dimsplus,'decoration',opt.decoration);
        switch lower(opt.type)
            case '3d'
                a.D = activedisplay3D(a.SI,'in',opt.in, ...
                    'clipmode',opt.clipmode,'clip',opt.clip);
            case 'frame'
                a.D = activedisplayFrames(a.SI,'in',opt.in, ...
                    'clipmode',opt.clipmode,'clip',opt.clip, ...
                    'ncol',opt.ncol);
            case 'grid'
                a.D = activedisplayArray(a.SI,'in',opt.in);
                if ~isempty(opt.cmap), a.D.cmap = opt.cmap; end
        end
    case 'list'
        if ~isscalar(opt.proj), error 'projection should have one dimension for list display', end
        a.SI = projection(a.G,opt.proj); %,'data',opt.data);
        a.D = activedisplayList(a.SI,'in',opt.in,'selmultin',~opt.averageselection);
    case 'slider'
        a.SI = projection(a.G,opt.proj,'decoration',opt.decoration);
        a.D = activedisplaySlider(a.SI,'in',opt.in,'layout',opt.layout);
    case 'player'
        a.SI = projection(a.G,opt.proj,'decoration',opt.decoration);
        a.D = activedisplayPlayer(a.SI,'in',opt.in);
    otherwise
        error('cannot handle type ''%s'' yet',opt.type)
end



%---
function b = isfigoraxeshandle(h)

b = (length(h)==1) && (ishandle(h) || (isnumeric(h) && h>0 && ~mod(h,1)));




