classdef projection < sliceinfo
    % function P = projection(G,proj,varargin)
    
    % The main point of projection class is its 'slicing' action, i.e. get
    % a section of the data, where:
    % - some dimensions being conserved,
    % - in some specific dimension(s) averaging is performed on specific
    %   regions of interests in specific dimension(s) 
    % - in the remaining dimensions a single index is used
    %
    % Main functions are:
    % - setauxiliary sets general information (which dimensions are
    %   conserved / averaged, etc.)
    % - prepareslicing prepares the data and other variables for slicing
    % - slicing performs the slicing (including loops for multiple regions
    %   of interest)
    
    % communication properties
    properties (SetAccess='private')
        G
    end
    
    % primary properties (i.e. defines how
    properties
        proj = zeros(1,0);
        dimsplus = zeros(1,0);
    end    
    properties (SetAccess='private')
        dataset = struct('name','data','val',[], ...
            'sizesslice',[],'sizesplus',[], ...
            'bgflag','','bgop','','bg',[]);
    end    
    properties (Dependent)
        data
        background
    end   
    properties (SetAccess='private', GetAccess='public')
        % secondary properties (changed in function 'setauxiliary')
        datasizes = [];
        datadims = struct('slice',[],'plus',[],'display',[], ...
            'singl',[], ...
            'orth',[],'avail',[],'nodisplay',[]);
        nodata = true;
        % tertiary properties (changed in function 'prepareslicing')
        slicetool = struct('subsvect',[],'data',[]);
    end
    
    % Constructor, destructor, and events
    methods
        function P = projection(G,proj,varargin)
            fn4D_dbstack
            P = P@sliceinfo(length(proj));
            P.G = G;
            % communication with parent
            addparent(P,G)
            P.datasizes = ones(1,P.G.nddata);
            % set projection - automatic settings and changes happen
            P.proj = proj;
            % set properties
            if nargin==2, return, end
            if ~ischar(varargin{1})
                error('syntax changed for setting dimsplus')
            end
            set(P,varargin{:})
        end        
    end
    
    % Update upon events
    methods
        function updateDown(P,~,evnt)
            fn4D_dbstack(['G2S ' evnt.flag])
            g=P.G;
            switch evnt.flag
                case 'nddata'
                    if length(P.datasizes)<g.nddata, setauxiliary(P), end
                case 'sizes'
                    idx = union(P.datadims.slice,P.datadims.nodisplay);
                    if any(P.datasizes(idx)~=g.sizes(idx)), setauxiliary(P), end
                case 'grid'
                    P.grid = g.grid(P.proj,:);
                case 'labels'
                    P.labels = g.labels(P.proj);
                case 'units'
                    P.units = g.units(P.proj);
                case 'ijkl2'
                    P.ij2 = g.ijkl2(P.proj);
                case 'ijkl'
                    % update slice only if change in datadims.avail
                    % dimensions
                    ijklold = evnt.oldvalue;
                    ijklchg = (g.ijkl~=ijklold);
                    if any(ijklchg(P.datadims.avail))
                        slicing(P,P.datadims.orth,'kl',1)
                    end
                case 'zoom'
                    P.zoom = g.zoom(P.proj,:);
                case 'selection'
                    % evnt.dims is empty when all selections are changed
                    ok = false;
                    if all(ismember(evnt.dims,P.proj)) || isempty(evnt.dims)
                        dims = evnt.dims;
                        if isempty(dims), dims = 1:g.nddata; end
                        ok = true;
                        slicedims = find(ismember(P.proj,dims));
                        setselection(P,slicedims,getselsets(g.selection,P.proj(slicedims)), ...
                            evnt.selflag,evnt.ind,evnt.value)
                    end
                    if all(ismember(evnt.dims,P.datadims.nodisplay)) ...
                             || isempty(evnt.dims)
                         ok = true;
                         slicing(P,evnt.dims,evnt.selflag,evnt.ind,evnt.value)
                    end
                    if ~ok && ~isempty(intersect(evnt.dims,P.datadims.nodisplay))
                        disp('selection not compatible with data or projection')
                    end
            end
        end
        
        function updateUp(P,evnt)
            fn4D_dbstack(['S2G ' evnt.flag])
            g=P.G;
            switch evnt.flag
                case 'labels'
                    g.labels(P.proj) = P.labels;
                case 'units'
                    % do nothing: units can be over-written locally in the
                    % sliceinfo, without changing units in the geometry
                    % parent
                case 'ij2'
                    g.ijkl2(P.proj) = P.ij2;
                case 'zoom'
                    g.zoom(P.proj,:) = P.zoom;
                case 'selection'
                    if strcmp(evnt.selflag,'reset')
                        updateselection(g,[],'reset')
                    else
                        selset = getselset(P.selection,evnt.dims);
                        dims = P.proj(evnt.dims);
                        if isempty(dims), error programming, end
                        setselectiondims(g,dims,selset, ...
                            evnt.selflag,evnt.ind,evnt.value)
                    end
            end
        end
    end
    
    % Slicing
    methods (Access='private')
        function slicing(P,dims,flag,ind,value)
            fn4D_dbstack
            % input
            if P.nodata
                fn4D_dbstack('cancel slicing')
                return
            end
            if nargin<3, dims=[]; flag='all'; end
            if nargin<4, ind=[]; end
            if nargin<5, value=[]; end            
            % shortcuts
            g = P.G; SI = P;
          
            % update datadims.orth (+ auxiliary properties) and/or slicetool
            % if necessary 
            switch flag
                case 'kl'
                    % coordinates in 'available' dimensions have changed
                    prepareslicing(P)
                case 'reset'
                    if ~isequal(dims,[]), error programming, end
                    orths = getdims(g.selection,P.datadims.nodisplay);
                    if ~isequal(orths,{P.datadims.orths}), setauxiliary(P), end
                case {'all' 'new' 'remove'}
                    orths = getdims(g.selection,P.datadims.nodisplay);
                    if ~isequal(orths,{P.datadims.orths}), setauxiliary(P), end
            end
            
            % other info
            orths = P.datadims.orths;
            selexist = ~isempty(orths);
            ndo = length(orths);
            if ndo>1 || (strcmp(flag,'remove') && ~isequal(dims,[orths{:}]))
                % too difficult to do partial update when multiple
                % orthogonal dimension sets!
                if fn_ismemberstr(flag,{'reset'}), error programming, end
                flag = 'all';
            else
                if selexist, orth = orths{1}; else orth = []; end
            end
            
            % slice size
            slicesiz = ones(1,max([2 P.datadims.orth]));
            for i=1:ndo
                orthi = orths{i};
                slicesiz(orthi(1)) = numsel(g.selection,orthi);
            end
            
            % slice
            switch flag
                case 'reset'
                    SI.slice = slicenosel(P);
                case {'all','indices','kl'} 
                    % 'indices' is the case where polygons are the sames
                    % but indices changed because of change in data sizes
                    % 'kl' is the case where some datadims.avail coordinates changed
                    if selexist
                        clear sels slic
                        sels(1,ndo) = selectionND; % pre-allocate
                        for k=[prod(slicesiz) 1:prod(slicesiz)-1] % start with last element for pre-allocation
                            idx = fn_indices(slicesiz,k);
                            for i=1:ndo
                                orthi = orths{i};
                                sels(i) = getsel(g.selection,orthi,idx(orthi(1)));
                            end
                            slic(k) = slicesel(P,sels); %#ok<AGROW>
                        end
                        SI.slice = reshape(slic,slicesiz);
                    elseif strcmp(flag,'kl')
                        updateslice(SI,'change',1,slicenosel(P))
                    else
                        SI.slice = slicenosel(P);
                    end
                case {'new','add','change','affinity'}
                    % in this case, we assume there is one and only one
                    % orthogonal dimension set
                    sel = getsel(g.selection,orth,ind);
                    for i=1:length(ind)
                        slic(i) = slicesel(P,sel(i)); %#ok<AGROW>
                    end
                    % reshape before update to bring the "slicing"
                    % dimension to the first "orthogonal" dimension
                    if ~isscalar(ind)
                        orth1 = orths{1}(1);
                        slic = shiftdim(slic(:),1-orth1);
                    end
                    updateslice(SI,fn_switch(flag,'new','new','change'),ind,slic)
                    if ~isequal(size(SI.slice),slicesiz) % TODO: this is very dirty! updateslice might sometimes fail to set the appropriate size, so...
                        SI.slice = reshape(SI.slice,slicesiz);
                    end
                case {'active','reorder'}
                    % in this case, we assume there is one and only one
                    % orthogonal dimension set
                    updateslice(SI,flag,ind,value)
                case 'remove'
                    % in this case, we assume there is one and only one
                    % orthogonal dimension set
                    if selexist
                        updateslice(SI,flag,ind,value)
                    else
                        SI.slice = slicenosel(P);
                    end
                otherwise
                    error programming
            end         
        end
        
        function prepareslicing(P)
            % shortcut
            d = P.datadims;
            ndo = length(d.orths);
            % collapse dimensions corresponding to the orth selection sets
            for i=1:ndo
                if any(diff(d.orths{i})~=1)
                    error('data dimensions of a selection should be consecutive')
                end
            end
            tool = struct;
            for k=1:length(P.dataset)
                name = P.dataset(k).name;
                data = P.dataset(k).val;
                for i=1:ndo
                    reshapesizes = size(data);
                    reshapesizes(d.orths{i}(1)) = prod(reshapesizes(d.orths{i}));
                    reshapesizes(d.orths{i}(2:end)) = 1;
                    data = reshape(data,reshapesizes);
                end
                tool.(name) = data;
            end
            % prepare subsref structure
            subsvect = cell(1,P.G.nddata);
            for k=d.display,    subsvect{k} = ':'; end
            for k=d.singl,      subsvect{k} = ':'; end
            for k=d.avail,      subsvect{k} = P.G.ijkl(k); end
            for k=[d.orths{:}], subsvect{k} = 1; end % subsvect{d.orth(1)} will be set during slicing
            tool.subsvect = subsvect;
            tool.getall = all(strcmp(subsvect,':'));
            % set property
            P.slicetool = tool;
        end
        
        function slice = slicenosel(P,dobgsub)
            fn4D_dbstack
            % prepare
            tool = P.slicetool; d = P.datadims;
            subsstruct = struct('type','()','subs',{tool.subsvect});
            slice = struct('active',true);
            % data
            for k=1:length(P.dataset)
                s = P.dataset(k);
                subdata = tool.(s.name);
                if ~tool.getall
                    subdata = subsref(subdata,subsstruct);
                end
                subdata = permute(subdata,[d.slice d.plus d.orth d.avail d.singl P.G.nddata+[1 2]]);
                slice.(s.name) = subdata;
            end
            % background subtraction
            if nargin<2 || dobgsub
                slice = slicebgsub(P,slice);
            end
        end
        
        function slice = slicesel(P,sels,dobgsub)
            fn4D_dbstack
            % sel is a selectionND object with indices set; alternatively,
            % it can just be a structure with fields active and dataind

            % prepare
            tool = P.slicetool; d = P.datadims;
            subsvect = tool.subsvect;
            ndo = length(d.orths);
            for i=1:ndo
                subsvect{d.orths{i}(1)} = sels(i).dataind;
            end
            subsstruct = struct('type','()','subs',{subsvect});
            slice = struct('active',prod([sels.active]));
            % data
            for k=1:length(P.dataset)
                s = P.dataset(k);
                subdata = subsref(tool.(s.name),subsstruct);
                for i=1:ndo
                    subdata = nmean(subdata,d.orths{i}(1));
                end
                subdata = permute(subdata,[d.slice d.plus d.orth d.avail d.singl P.G.nddata+1]);
                % check [TODO: remove]
                if size(subdata)~=size(zeros([s.sizesslice s.sizesplus 1])), error programming, end 
                slice.(s.name) = subdata;
            end
            % background subtraction
            if nargin<3 || dobgsub
                slice = slicebgsub(P,slice);
            end
        end
        
        function slice = slicebgsub(P,slice)
            fn4D_dbstack
            for k = 1:length(P.dataset)
                bg = P.dataset(k).bg;
                if isempty(bg), continue, end
                name = P.dataset(k).name;
                switch P.dataset(k).bgop
                    case '/'
                        slice.(name) = slice.(name) ./ bg;
                    case '-'
                        slice.(name) = slice.(name) - bg;
                    otherwise
                        error programming
                end
            end
        end
        
        function eraseslice(P)
            % data became not valid - produce a 'zero' slice of the
            % accurate size
            %---
            % attention! the function does not set P.nodata to true, this
            % should be done in the function which calls eraseslice (this
            % allows eraseslice to be overwritten by superclasses while
            % keeping nodata property private)
            slic = P.slice(1);
            for k=1:length(P.dataset)
                s = P.dataset(k);
                slic.(s.name) = zeros([P.datasizes([P.proj P.dimsplus]) 1],'uint8');
            end
            P.slice = slic;
        end           
    end
    
    % GET/SET data
    methods
        function data = get.data(P)
            fn4D_dbstack
            data = P.dataset(1).val;
        end
        
        function set.data(P,data)
            fn4D_dbstack
            if isempty(data)
                rmdata(P,'data')
                return
            end
            % possible change in G.sizes (in that case, automatic updates
            % happen)
            s = size(data);
            idx = setdiff(find(s>1),P.dimsplus);
            P.G.sizes(idx) = s(idx);
            % set data
            setdata(P,'data',data)
        end
        
        function bg = get.background(P)
            fn4D_dbstack
            bg = P.dataset(1).bg;
        end
        
        function set.background(P,bg)
            fn4D_dbstack
            setbackground(P,'data',bg)
        end
        
        function setdata(P,name,val,bg,op)
            fn4D_dbstack
            % check size
            s = size(val);
            s(end+1:max([P.proj P.dimsplus])) = 1;
            idx = setdiff(find(s>1),P.dimsplus);
            if any(s(idx)~=P.G.sizes(idx)) ...
                    || (any(s(P.proj)~=P.G.sizes(P.proj)))
                error('data size does not match: first set desired data size using setdatasizes(P,s)')
            end
            % which element in data set?
            kdt = find(ismember({P.dataset.name},name));
            if isempty(kdt)
                kdt = length(P.dataset)+1;
                P.dataset(kdt).name = name;
            end
            % set data
            if issparse(val)
                ndval = ndims(val);
                if ~(isequal(P.proj,1:ndval-1) || isequal(P.proj,ndval))
                    disp 'sparse or ndSparse data is not compatible with this projection -> convert to full'
                    val = full(val);
                end
            end
            P.dataset(kdt).val = val;
            % update all auxiliary properties 
            setauxiliary(P)
            % flag for background re-definition
            if nargin<4
                bg = P.dataset(kdt).bgflag;
                if strcmp(bg,'user'), bg = ''; end
            end
            if nargin<5, op = P.dataset(kdt).bgop; end
            % set background (but don't update slice)
            setbackground(P,name,bg,op,false)
            % update slice
            slicing(P,[],'all');
        end
        
        function setbackground(P,name,bg,op,doslicing)
            fn4D_dbstack
            % the 'background' concept is to perform some normalization
            % after slicing the data; typically, a division or a
            % subtraction by the slice obtained from one specific region
            %
            % bg can be one of the flags '', 'avg', 'selfirst', 'sellast',
            % 'selall' 'copydata' or the value of the background itself 
            %---
            % new flag will be one of '', 'avg', 'user', 'copydata'
            if nargin<5, doslicing=true; end
            
            % initializations: check bg, set local flag and new flag,
            % precomputations
            if isempty(bg)
                flag = ''; 
                newflag = flag;
            elseif isnumeric(bg)
                flag = 'user';
                bgval = bg;
                newflag = flag;
            elseif ~ischar(bg)
                error('wrong background argument')
            elseif strcmp(bg,'copydata')
                flag = P.dataset(1).bgflag;
                if strcmp(flag,'user')
                    error('cannot copy the background definition from primary data to other')
                end
                newflag = 'copydata';
            elseif strcmp(bg,'avg')
                flag = 'sel';
                newflag = 'avg';
                % average
                for k=1:length(P.dataset)
                    s = P.dataset(k);
                    subdata = s.val;
                    for i=P.datadims.nodisplay
                        subdata = nmean(subdata,i);
                    end
                    bgval.(s.name) = reshape(subdata,[s.sizesslice s.sizesplus 1 1]);
                end
            elseif fn_ismemberstr(bg,{'selfirst','sellast','selall','avg'})
                flag = 'sel';
                newflag = 'user';
                orth = P.datadims.orth;
                if isempty(orth)
                    error('cannot define a background: no orthogonal selection')
                end
                switch bg
                    case 'selfirst'
                        sel = getsel(P.G.selection,orth,1);
                    case 'sellast'
                        nsel = numsel(P.G.selection,orth);
                        sel = getsel(P.G.selection,orth,nsel);
                    case 'selall'
                        % need to compute indices 'by hand'
                        set = getselset(P.G.selection,orth);
                        ind = cat(1,set.singleset.dataind);
                        ind = unique(ind);
                        sel = struct('active',true,'dataind',ind);
                    case 'avg'
                        updateorthdims(P,P.datadims.nodisplay)
                        ind = 1:prod(P.datasizes(P.datadims.nodisplay));
                        sel = struct('active',true,'dataind',ind);
                        % when data will be changed, it will still be
                        % possible to compute the background 
                        newflag = 'avg'; 
                end
                bgval = slicesel(P,sel,false);
                if strcmp(bg,'avg'), updateorthdims(P,orth), end
            else
                error('wrong background flag')
            end
            
            % which part of data set is concerned?
            if strcmp(name,'data')
                % any additional data whose background calculation is
                % copied from 'data'?
                kbg = 1;
                kplus = find(fn_ismemberstr({P.dataset.bgflag},'copydata'));
                if ~isempty(kplus) && strcmp(flag,'user')
                    error('cannot copy the background definition from primary data to other(s)')
                end
            else
                kbg = find(fn_ismemberstr({P.dataset.name},name));
                if isempty(kbg), error('cannot find name ''%s'' in data set',name), end
                kplus = [];
            end
            
            % set new flag
            P.dataset(kbg).bgflag = newflag;
            
            % compute/set background
            for k=[kbg kplus]
                switch flag
                    case ''
                        P.dataset(k).bgop = '';
                        P.dataset(k).bg = [];
                    case {'user','sel'}
                        if isempty(op)
                            error('background operation (divide or subtract) not defined')
                        elseif ~ismember(op,'/-')
                            error('wrong background operation ''%s''',op)
                        end
                        P.dataset(k).bgop = op;
                        P.dataset(k).bg = bgval.(P.dataset(k).name);
                    otherwise
                        error programming
                end
            end
            
            % update slice?
            if doslicing
                slicing(P,[],'all')
            end
        end
        
        function rmdata(P,name)
            if nargin<2 name = 'all'; end
            if strcmp(name,'data')
                % this 'trick' can be used to avoid slice updating
                P.dataset(1).val = [];
                eraseslice(P);
                P.nodata = true;
            elseif strcmp(name,'all')
                P.dataset(2:end) = [];
                slicing(P,[],'all')
            else
                b = fn_ismemberstr({P.dataset.name},name);
                P.dataset(b) = [];
                slicing(P,[],'all')
            end
        end
    end
    
    % SET proj and dimsplus
    methods
        function set.proj(P,proj)
            fn4D_dbstack
            % checks
            if isequal(P.proj,proj), return, end
            if length(proj)~=P.nd
                error('dimension mismatch between projection and sliceinfo')
            end
            if ~all(diff(proj)>0)
                error('projection dimensions must non-decreasing')
            end
            P.proj = proj;
            % increase in nddata? - automatic updates
            g = P.G; %#ok<MCSUP>
            g.nddata = max(g.nddata,max(proj));
            % update auxiliary properties
            setauxiliary(P)
            
            % many changes; G should not be notified about them
            P.upnotify = false;
            % update grid, labels and units (in P)
            P.grid = g.grid(proj,:);
            P.labels = g.labels(proj);
            P.units = g.units(proj);
            % slice (in P)
            slicing(P,[],'all');
            % projected coordinates (in P)
            P.ij2 = g.ijkl2(P.proj);
            % zoom
            P.zoom = g.zoom(P.proj,:);
            % selection
            selset = getselsets(g.selection,P.proj);
            alldims = {selset.t.dims};
            for i=1:length(alldims)
                setselection(P,alldims{i},getselset(selset,alldims{i}),'all')
            end
            % re-establish communication
            P.upnotify = true;
        end
        
        function set.dimsplus(P,dimsplus)
            fn4D_dbstack
            % checks
            if isequal(P.dimsplus,dimsplus), return, end
            if any(intersect(dimsplus,P.datadims.slice))
                error('dimsplus intersects projection')
            elseif any(intersect(dimsplus,P.datadims.orth))
                error('dimsplus intersects orthogonal selection')
            end
            P.dimsplus = dimsplus;
            % increase in nddata? - automatic updates
            if ~isempty(dimsplus), P.G.nddata = max(P.G.nddata,max(dimsplus)); end
            % update auxiliary properties
            setauxiliary(P)
            
            % update slice
            slicing(P,[],'all');
        end
    end
    
    % Update secondary properties
    methods
        function setauxiliary(P)
            % data sizes
            actualdatasizes = size(P.data);
            actualdatasizes(end+1:P.G.nddata) = 1;
            P.datasizes = actualdatasizes;
            P.datasizes(P.proj) = P.G.sizes(P.proj);
            allsingl = find(P.datasizes==1);
            
            % changes in datadims
            % (display = slice + plus)
            P.datadims.slice = P.proj;
            P.datadims.plus = P.dimsplus;
            if any(ismember(P.proj,P.dimsplus)), error('wrong dimensions'), end
            P.datadims.display = row(union(P.proj,P.dimsplus));
            % (allnodisplay = singl + nodisplay)
            allnodisplay = setdiff(1:P.G.nddata,P.datadims.display);
            P.datadims.singl = row(intersect(allnodisplay,allsingl));
            P.datadims.nodisplay = setdiff(allnodisplay,P.datadims.singl);
            % (nodisplay = orth + avail)
            orths = getdims(P.G.selection,P.datadims.nodisplay);
            orth = sort([orths{:}]);
            if any(diff(orth)==0)
                error 'programming: intersecting selections'
            end
            P.datadims.orths = orths;
            P.datadims.orth = orth;
            P.datadims.avail = setdiff(P.datadims.nodisplay,P.datadims.orth);
            
            % slice sizes
            for k=1:length(P.dataset)
                s = size(P.dataset(k).val); 
                s(end+1:P.G.nddata)=1;
                P.dataset(k).sizesslice = s(P.proj);
                P.dataset(k).sizesplus  = s(P.dimsplus);
            end           
            
            % match between different sizes? 
            % (between primary data and G in dimensions 'proj' and 'nodisplay') 
            idx = union(P.datadims.slice,P.datadims.nodisplay);
            b = any(actualdatasizes(idx)~=P.G.sizes(idx));
            % (between all data in dimensions 'singl' and 'nodisplay') 
            idx = union(P.datadims.singl,P.datadims.nodisplay);
            s0 = P.datasizes(idx);
            for k=1:length(P.dataset)
                s = size(P.dataset(k).val); s(end+1:P.G.nddata) = 1;
                if any(s(idx)~=s0), b=true; end
            end      
            
            % prepare for slicing
            if b
                eraseslice(P)
                P.nodata = true;
            else
                P.nodata = false;
                prepareslicing(P)
            end
        end
    end
    
    % Local menu and its events
    methods
        function initlocalmenu(P,hb)
            initlocalmenu@sliceinfo(P,hb)
            m = P.contextmenu;
            info = get(m,'userdata');
            
            info.hbg(1) = uimenu(m,'label','divide by global average','separator','on', ...
                'callback',@(hu,evnt)setbackground(P,'data','avg','/'));
            info.hbg(2) = uimenu(m,'label','subtract global average', ...
                'callback',@(hu,evnt)setbackground(P,'data','avg','-'));
            info.hbg(3) = uimenu(m,'label','subtract first display', ...
                'callback',@(hu,evnt)setbackground(P,'data','selfirst','-'));
            info.hbg(4) = uimenu(m,'label','subtract last display', ...
                'callback',@(hu,evnt)setbackground(P,'data','sellast','-'));
            info.hbg(5) = uimenu(m,'label','subtract average display', ...
                'callback',@(hu,evnt)setbackground(P,'data','selall','-'));
            info.hbg(6) = uimenu(m,'label','cancel background normalization', ...
                'callback',@(hu,evnt)setbackground(P,'data',''));

            set(m,'userdata',info)
        end
    end
end
