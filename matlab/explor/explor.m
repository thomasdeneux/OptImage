 
classdef explor < interface
    
    properties
        data
        G
    end
    properties (SetObservable = true)
        dodock = false;
    end
    properties (SetAccess = 'private')
        dlist
    end
    
    % Creation
    methods
        function X = explor(a,G)
            opt = struct('defaultList',false,'recentHeaders',struct('hash',{},'value',{},'str',{}));
            hf = figure('integerhandle','off');
            X = X@interface(hf,'DATA EXPLORE',opt);
            
            if nargin==0, a = rand(3,4,5); end
            if ischar(a), a = evalin('base',a); end
            X.data = a;
            
            X.options.defaultList = false; % this option is disabled
            
            init_grob(X)
            interface_end(X)
            init_table(X)
            
            % geometry object 
            if nargin>=2
                X.G = G;
                u = X.grob.table;
                tdata = get(u,'Data');
                cnames = get(u,'ColumnName');
                ilab = find(strcmp(cnames,'Label'),1);
                [tdata{:,ilab}] = deal(G.labels{:});
                iunit = find(strcmp(cnames,'Unit'),1);
                iscale = find(strcmp(cnames,'Scale/Items'),1);
                scale = diag(G.mat(2:end,2:end));
                offset = G.mat(2:end,1);
                start = scale+offset;
                ndim = G.nddata; % (also equal to length(scale) of course)
                for i=1:ndim
                    if iscell(G.units{i})
                        tdata{i,iscale} = fn_strcat(G.units{i},',');
                    else
                        tdata{i,iunit} = G.units{i};
                        tdata{i,iscale} = [num2str(scale(i)) ' [start ' num2str(start(i)) ']'];
                    end
                end
                set(u,'Data',tdata)
            end
        end
        function init_grob(X)
            g = struct;
            % table
            g.table = uitable;
            % headers
            g.hlabel = uicontrol('style','text','string','headers:');
            g.dfl = uicontrol('string','Default','callback',@(u,e)action(X,'defaultlabels'));
            g.closeall = uicontrol('string','close all','callback',@(u,e)action(X,'closeall'));
            g.load = uicontrol('string','Load','callback',@(u,e)loadheaders(X));
            g.save = uicontrol('string','Save','callback',@(u,e)saveheaders(X));
            g.recenth = uicontrol('style','popupmenu','callback',@(u,e)loadrecentheaders(X), ...
                'string',{'(select recent headers)' X.options.recentHeaders.str});
            % display config
            g.dlabel = uicontrol('style','text','string','display:');
            g.dchoice = uicontrol('style','popupmenu','string','(select dimention to display first)');
            g.dock = fn_propcontrol(X,'dodock','checkbox','string','dock figures').hu;
            % go!
            g.ok = uicontrol('string','OK','callback',@(u,e)showdisplay(X,'main'));
            X.grob = g;
        end
        function init_menus(X)
            init_menus@interface(X)
        end
    end
    
    % Table
    methods
        function init_table(X)
            % empty table
            u = X.grob.table;
            set(u,'ColumnName',{'Dim','Size','Label','Unit','Scale/Items','P display','S display' 'List'}, ...
                'ColumnFormat',{'numeric' 'numeric' 'char' 'char' 'char' 'logical' 'logical' 'logical'}, ...
                'ColumnEditable',logical([0 0 1 1 1 1 1 1]))
            [iSI iP iS iL iChar] = deal(5,6,7,8,3:4);
            p = get(u,'pos'); w = p(3);
            widths = {30 55 [] [] [] 55 55 55};
            wavail = w - sum([widths{:}]) - 2;
            idxauto = find(fn_isemptyc(widths));
            [widths{idxauto}] = deal(floor(wavail/length(idxauto)));
            set(u,'ColumnWidth',widths)
            set(u,'RowName',[])
            set(u,'CellEditCallback',@(u,e)action(X,'celledit',e))
            % fill-in the data
            s = size(X.data);
            ndim = length(s);
            tdata = cell(ndim,8);
            for i=1:ndim
                tdata{i,1} = i;
                tdata{i,2} = s(i);
            end
            [tdata{:,iChar}] = deal('');
            [tdata{:,iSI}] = deal('1 [start 1]');
            [tdata{:,[iP iS iL]}] = deal(false);
            set(u,'Data',tdata)
        end
        function action(X,flag,varargin)
            if ~isvalid(X), return, end % can happen when all figures are closed simultaneously
            u = X.grob.table;
            if ~ishandle(u), return, end % can happen when all figures are closed simultaneously
            cnames = get(u,'ColumnName');
            iP = find(strcmp(cnames,'P display'),1);
            iS = find(strcmp(cnames,'S display'),1);
            iL = find(strcmp(cnames,'List'),1);
            tdata = get(u,'Data');
            switch flag
                case 'celledit'
                    e = varargin{1};
                    ename = cnames{e.Indices(2)};
                    switch ename
                        case 'Scale/Items'
                            irow = e.Indices(1);
                            iSI = e.Indices(2);
                            [~, ~, tdata{irow,iSI}] = ReadScalesItems(tdata{irow,iSI},size(X.data,irow),e.PreviousData);
                            set(u,'Data',tdata)
                        case {'P display' 'S display'}
                            % choice between P/S is exclusive
                            iex = setdiff([iP iS],e.Indices(2));
                            [tdata{e.Indices(1),iex}] = deal(false);
                            set(u,'Data',tdata)
                            % update popup menu for display type
                            nshow = sum([tdata{:,iP}]);
                            oldstr = get(X.grob.dchoice,'string');
                            switch nshow
                                case 0
                                    str = '(select dim. for primary display)';
                                case 1
                                    str = {'plot' 'list' 'slider'};
                                case 2
                                    str = {'image'};
                                case 3
                                    str = {'3D' 'frames' 'grid'};
                                otherwise
                                    str = '(too many dim. for primary display)';
                            end
                            if ~isequal(str,oldstr)
                                set(X.grob.dchoice,'string',str,'value',1)
                            end
                        case 'List'
                            if tdata{e.Indices(1),iL}
                                % display list
                                showdisplay(X,'list',e.Indices(1))
                            else
                                % list is supposed to be displayed already,
                                % forbid unchecking the box
                                tdata{e.Indices(1),iL} = true;
                            end
                    end
                    set(u,'Data',tdata)
                case 'defaultlabels'
                    ndim = ndims(X.data);
                    ilab = find(strcmp(cnames,'Label'),1);
                    labels = {'x' 'y' 'time'};
                    [labels{end+1:ndim}] = deal('');
                    for i=1:ndim, tdata{i,ilab} = labels{i}; end
                    iunit = find(strcmp(cnames,'Unit'),1);
                    units = {'px' 'px' 'frame'};
                    [units{end+1:ndim}] = deal('');
                    for i=1:ndim, tdata{i,iunit} = units{i}; end
                    set(u,'Data',tdata)
                case 'closeall'
                    close(findobj('type','figure','tag','DATA DISPLAY'))
                case 'removelists'
                    X.dlist = [];
                    [tdata{:,iL}] = deal(false);
                    set(u,'Data',tdata)
            end
        end
        function loadheaders(X,fname)
            if nargin<2
                fname = fn_getfile('*.csv'); 
                if fname==0, return, end
            end
            % columns of interest
            cnames = get(X.grob.table,'ColumnName');
            columns = {'Dim','Label','Unit','Scale/Items'};
            columns1 = strrep(columns,'/','_');
            ncol = length(columns);
            % read file
            if ischar(fname)
                txt = fn_readtext(fname);
                head = txt{1}; txt(1) = [];
                head = fn_strcut(head,','); nhead = length(head);
                nrow = length(txt);
                c = cell(nrow,nhead);
                for i=1:nrow
                    line = txt{i};
                    line = fn_strcut(line,'"');
                    for j=1:2:length(line)
                        % comma act as separators only outside quote marks
                        line{j} = fn_strcut(line{j},',');
                    end
                    if ~isscalar(line), line = [line{:}]; end
                    c(i,:) = cell2mat(line);
                end
                t = cell2struct(c,head,2);
            else
                % already a structure
                t = fname;
            end
            tnames = fieldnames(t);
            tdata = get(X.grob.table,'Data');
            t = struct2cell(t); t = [t{:}]; % convert structure to cell array
            % check number of dimensions
            if size(t,1)~=size(tdata,1)
                waitfor(warndlg('Header description does not match data number of dimensions'))
            end
            nrow = min(size(t,1),size(tdata,1));
            % fill-in the table
            for i=1:ncol
                id = find(strcmp(columns{i},cnames),1);
                it = find(strcmp(columns1{i},tnames),1);
                if isempty(it), continue, end
                tdata(1:nrow,id) = t(1:nrow,it);
            end
            % check item lists
            iSI = find(strcmp('Scale/Items',cnames),1);
            if isempty(ReadScalesItems(tdata(:,iSI),size(X.data)))
                errordlg 'Some dimension size does not match with number of items in list'
                return
            end
            % finished!
            set(X.grob.table,'Data',tdata);
        end
        function t = saveheaders(X,fname)
            if nargin<2 
                if nargout~=0
                    fname = []; % do not save in file
                else
                    fname = fn_savefile('*.csv');
                    if fname==0, return, end
                end
            end
            cnames = get(X.grob.table,'ColumnName');
            columns = {'Dim','Label','Unit','Scale/Items'};
            columns1 = strrep(columns,'/','_');
            ncol = length(columns);
            keep = zeros(1,ncol);
            for i=1:ncol, keep(i) = find(strcmp(columns{i},cnames),1); end
            tdata = get(X.grob.table,'Data'); % numerics and chars
            c = fn_num2str(tdata(:,keep)); % chars only
            t = cell2struct(num2cell(c,1),columns1,2); % structure, each value is a cell array of nrow values

            % save/output
            if ~isempty(fname)
                % quote strings that contain comma
                for i=1:numel(c), if any(c{i}==','), c{i} = ['"' c{i} '"']; end, end
                % add column headers
                c = [columns1; c];
                % save
                nrow = size(c,1);
                txt = cell(1,nrow);
                for i=1:nrow, txt{i} = fn_strcat(c(i,:),','); end
                fn_savetext(txt,fname)
            end
            if nargout==0
                clear t
            end
        end
        function loadrecentheaders(X)
            k = get(X.grob.recenth,'value')-1;
            if k==0, return, end
            t = X.options.recentHeaders(k).value;
            loadheaders(X,t);
            set(X.grob.recenth,'value',1)
        end
        function saverecentheaders(X)
            t = saveheaders(X); % structure
            nrow = length(t.Label);
            % (for item lists, add the number of items to the description
            % string)
            ScaleItems = t.Scale_Items;
            for i=1:nrow
                nitem = 1+sum(t.Scale_Items{i}==',');
                if nitem>1, ScaleItems{i} = [num2str(nitem) ':' ScaleItems{i}]; end
            end
            c = [t.Label repmat({'('},nrow,1) ScaleItems t.Unit repmat({'), '},nrow,1)]'; c{end} = ')';
            str = [c{:}];
            hash = fn_hash(str);
            s = X.options.recentHeaders;
            kopt = find(strcmp(hash,{s.hash}),1);
            if kopt
                s = s([kopt setdiff(1:end,kopt)]);
            else
                maxitems = 20;
                s = [struct('hash',hash,'value',t,'str',str) s(1:min(end,maxitems-1))];
            end
            X.options.recentHeaders = s;
            saveoptions(X)
            ustr = get(X.grob.recenth,'string');
            set(X.grob.recenth,'string',{ustr{1} s.str},'value',2)
        end
    end
    
    % Display
    methods
        function showdisplay(X,flag,varargin)
            % get data and display information
            u = X.grob.table;
            cnames = get(u,'ColumnName');
            tdata = get(u,'Data');
            nd = size(tdata,1);
            
            iLa = find(strcmp(cnames,'Label'),1);
            iU = find(strcmp(cnames,'Unit'),1);
            iSI = find(strcmp(cnames,'Scale/Items'),1);
            iP = find(strcmp(cnames,'P display'),1);
            iS = find(strcmp(cnames,'S display'),1);
            iL = find(strcmp(cnames,'List'),1);

            % geometry object
            if isempty(X.G) || ~isvalid(X.G)
                %F = focus.find(0);
                %X.G = rotation(F,'mat',mat);
                X.G = geometry('sizes',size(X.data)); % for now, no link between different data
            else
                %X.G.mat = mat; % might have changed...
            end
            
            % update geometry object with header info
            labels = tdata(:,iLa)';
            for i=1:nd, if isempty(labels{i}), labels{i}=sprintf('dim %i',i); end, end
            units = tdata(:,iU)';
            scales_items = tdata(:,iSI);
            [mat listitems] = ReadScalesItems(scales_items,size(X.data));
            for i=1:nd
                if ~isempty(listitems{i}), units{i}=listitems{i}; end
            end
            X.G.grid = [column(mat{1}) column(mat{2})];
            if ~all(fn_isemptyc(labels)), set(X.G,'labels',labels), end
            if ~all(fn_isemptyc(units)), set(X.G,'units',units), end            
            
            
            % display
            switch flag
                case 'main'
                    % new figure
                    hf = figure('tag','DATA DISPLAY','integerhandle','off','numbertitle','off');
                    if X.dodock, set(hf,'windowstyle','docked'), end
                    
                    % main data display
                    proj = find([tdata{:,iP}]);
                    set(hf,'name',fn_strcat(labels(proj),' x '))
                    if ~ismember(length(proj),[1 2 3]), return, end
                    typechoice = get(X.grob.dchoice,'string');
                    type = typechoice{get(X.grob.dchoice,'value')};
                    dimsplus = find([tdata{:,iS}]);
                    a = fourd(X.data,X.G,'type',type,'proj',proj,'dimsplus',dimsplus,'in',hf);
                    switch type
                        case {'2d' '2dcol'}
                    end
                case 'list'
                    % list figure organization
                    if isempty(X.dlist) || ~ishandle(X.dlist.hf)
                        hf = figure('tag','DATA DISPLAY','integerhandle','off','numbertitle','off', ...
                            'name','LISTS', ...
                            'deletefcn',@(u,e)action(X,'removelists'));
                        if X.dodock
                            set(hf,'windowstyle','docked')
                        else
                            s = fn_pixelsize(hf); s(1) = s(2)/2;
                            fn_setfigsize(hf,s)
                        end
                        hu = uicontrol('parent',hf,'style','listbox');
                        X.dlist = struct('hf',hf,'hu',hu,'hsep',[],'xsep',[]);
                    else
                        hf = X.dlist.hf;
                        nlist = length(X.dlist.hu);
                        if nlist==0, error 'no list in list figure', end
                        if ~strcmp(get(hf,'windowstyle'),'docked')
                            s = fn_pixelsize(hf); s(1) = s(1)*(nlist+1)/nlist;
                            fn_setfigsize(hf,s)
                        end
                        col = [1 1 1]*.5;
                        X.dlist.hsep(nlist) = uicontrol('parent',hf, ...
                            'style','frame','enable','off','foregroundcolor',col,'backgroundcolor',col, ...
                            'buttondownfcn',@(u,e)fn_buttonmotion(@()moveLists(X,nlist),'pointer','left'));
                        hu = uicontrol('parent',hf,'style','listbox');
                        X.dlist.hu(nlist+1) = hu;
                        if nlist==1
                            X.dlist.xsep = .5;
                        else
                            X.dlist.xsep = [X.dlist.xsep 1]*nlist/(nlist+1);
                        end
                            
                    end
                    organizeLists(X)
                    
                    % list display
                    dim = varargin{1};
                    SI = projection(X.G,dim);
                    D = activedisplayList(SI,'in',hu);
                    %                 if ~isempty(listitems{lists(i)})
                    %                     D.itemnames = listitems{lists(i)};
                    %                 end
            end
            
            % save headers in the list of recent headers
            saverecentheaders(X)
        end
        function set.data(X,data)
            prevsiz = size(X.data);
%             if ~isequal(prevsiz,[0 0]) && ~isequal(size(data),prevsiz)
%                 error 'new data size does not match previous'
%             end
            X.data = data;
            % update size displays in table
            s = size(data);
            if ~isequal(s,prevsiz) && ~isequal(prevsiz,[0 0]) % prevsize is [0 0] at init
                ndimold = length(prevsiz);
                ndim = length(s);
                u = X.grob.table;
                tdata = get(u,'Data');
                for i=1:min(ndim,ndimold)
                    tdata{i,2} = s(i);
                end
                if ndim>ndimold
                    for i=ndimold+1:ndim
                        tdata{i,1} = i;
                        tdata{i,2} = s(i);
                    end
                    cnames = get(u,'ColumnName');
                    iChar = ismember(cnames,{'Label' 'Unit'});
                    iSI = find(strcmp(cnames,'Scale/Items'),1);
                    iCfg = ismember(cnames,{'P display' 'S display'});
                    iL = find(strcmp(cnames,'List'),1);
                    [tdata{:,iChar}] = deal('');
                    [tdata{:,iSI}] = deal('1 [start 1]');
                    [tdata{:,iCfg}] = deal(false);
                end
                set(u,'Data',tdata)
            end
            % update displays
            if ~isempty(X.G)
                SIs = getChildren(X.G); %#ok<MCSUP>
                for i=1:length(SIs)
                    if ~isempty(SIs{i}.data)
                        SIs{i}.data = data;
                    end
                end
            end
        end
        function organizeLists(X)
            s = X.dlist;
            nlist = length(s.hu);
            xgap = .012; y0 = .01; h = .94;
            wfr = .012; yfr0 = y0; hfr = h;
            edges = [0 s.xsep 1];
            for i=1:nlist
                set(s.hu(i),'units','normalized','pos',[edges(i)+xgap y0 diff(edges(i:i+1))-2*xgap h])
            end
            for i=1:nlist-1
                set(s.hsep(i),'units','normalized','pos',[s.xsep(i)-wfr/2 yfr0 wfr hfr])
            end
        end
        function moveLists(X,i)
            s = X.dlist;
            figsiz = fn_pixelsize(s.hf);
            p = get(s.hf,'currentPoint');
            edges = [0 s.xsep 1];
            xgap = .012; mindist = xgap*2.5;
            x = p(1)/figsiz(1);
            X.dlist.xsep(i) = max(edges(i)+mindist,min(edges(i+2)-mindist,x));
            organizeLists(X)
        end
    end
end

function [mat listitems fulldisplay] = ReadScalesItems(scales_items,siz,previousdata)
% 2 behaviors for 2 types of input:
% - scale: define mat, default listitems (cell array of empty arrays)
% - items list: define listitems (cell array of strings), default mat
% in both cases fulldisplay is defined
% if input cannot interpreted, all outputs are empty

% input
isc = iscell(scales_items);
if ~isc, scales_items = {scales_items}; end
n = length(scales_items);
if nargin>=3 && n>1, error argument, end

% read
mat = {ones(1,n) zeros(1,n)}; % scale, translation
listitems = cell(1,n);
fulldisplay = scales_items;
i = 1;
while i<=n
    a = scales_items{i};
    % scale + start?
    tokens = regexp(a,'^(.*)\[ *start{0,1}(.*)\]$','tokens');
    if ~isempty(tokens)
        tok = tokens{1};
        mat{1}(i) = evalin('base',tok{1});
        start = evalin('base',tok{2});
        mat{2}(i) = start - mat{1}(i);
        fulldisplay{i} = [num2str(mat{1}(i),17) ' [start ' num2str(start,17) ']'];
        i=i+1; continue
    end
    % number?
    try x = evalin('base',a); catch, x = []; end % can also be a string that evaluates to a number
    if isscalar(x)
        mat{1}(i) = x;
        mat{2}(i) = -x;
        fulldisplay{i} = num2str(x,16);
        i=i+1; continue
    elseif isvector(x) && isnumeric(x) && length(x)==siz(i) && max(abs(diff(x,2)))<diff(x(1:2))/1e6
        % vector of values
        if max(abs(diff(x,2)))==0
            dx = diff(x(1:2));
        else
            dx = (x(siz(i))-x(1))/(siz(i)-1);
        end
        mat{1}(i) = dx;
        mat{2}(i) = x(1) - dx;
        fulldisplay{i} = [num2str(mat{1}(i),17) ' [start ' num2str(x(1),17) ']'];
        i=i+1; continue
    end
    % list of items?
    if iscell(x) && all(fn_map(@ischar,x)) && length(x)==siz(i)
        list = x;
    else
        sep = fn_switch(any(a==','),',',' ');
        list = fn_strcut(a,sep);
    end
    if length(list)==siz(i)
        listitems{i} = list;
        fulldisplay{i} = fn_strcat(list,',');
        i=i+1; continue
    end
    % failed to read string
    if nargin>=3
        scales_items{i} = previousdata;
    else
        disp 'could interpret string neither as scale specification or as items list'
        [mat listitems fulldisplay] = deal([]);
        return
    end
end

% output
if ~isc, listitems = listitems{1}; fulldisplay = fulldisplay{1}; end

end


