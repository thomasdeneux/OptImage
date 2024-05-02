classdef sliceinfo < fn4Dhandle
    % function SI = sliceinfo(ndim,varargin)
    %---
    % Typing 'sliceinfo' with no argument results in an example,
    % non-functional (since the number of dimensions 'nd' is not set),
    % object. 
    % Some fields of the sliceinfo object can be set directly. The 'slice'
    % field can be set directly, or can be updated by the function
    % 'updateslice(SI,flag,ind,value)'. The 'selection' field must be set
    % by the function 'setselection(SI,selset,flag,ind,value)', or updated
    % by the function 'updateselection(SI,flag,ind,value)'. See the class
    % 'fn4Devent' help [not written yet...] to know more about what flag,
    % ind, value are.
    % The field 'ij' cannot be set directly but is automatically updated
    % by rounding 'ij2' and coercing it within range defined by 'size'.
    % Finally, the fields 'sizes' (size in the first 'nd' dimensions) and
    % 'sizesplus' (size in additional dimensions) cannot be set directly
    % but is automatically updated from the size of the data (in
    % 'slice.data'). 
    
    % properties whose values cannot be changed
    properties (SetAccess='private')
        nd
    end
    % properties whose values can be changed
    properties
        grid = {};
        labels = {};
        units = {};
        
        ij2
        zoom = zeros(0,2);
        decoration
    end
    properties (Dependent)
        slice
    end
    properties (SetAccess='private')
        selection
        sizes
        sizesplus = zeros(1,0);
        ij
    end
    properties (Dependent, SetAccess='private')
        ndplus = 0;
        selectionmarks        
    end  
    properties (SetAccess='protected')
        contextmenu
    end
    properties (Access='private')
        slicecontent = struct('active',1,'data',0);
    end
    
    % technical notes about events:
    % - since ij cannot be set and ij = round(ij2), a 'ij' change
    %   occurs only after a 'ij2' change
    % - a 'sizes' change occurs only after a 'slice' change
    % - a 'sizesplus' change occurs only after a 'slice' change
    
    % INIT
    methods
        function SI = sliceinfo(ndim,varargin)
            fn4D_dbstack
            % number of dimensions is set from the beginning and cannot be
            % changed
            if nargin==0, ndim=1; end
            SI.nd = ndim;
            
            % appropriate initializations
            SI.grid = repmat([1 0],[ndim 1]);
            [SI.labels{1:SI.nd}] = deal('');
            [SI.units{1:SI.nd}] = deal('');
            SI.sizes = ones(1,ndim);
            SI.ij2 = ones(ndim,1); % does not automatically update SI.ij (irrelevant change)
            SI.ij = ones(ndim,1);
            SI.zoom = NaN(ndim,2);
            SI.selection = selectionset(SI.sizes);
            
            % set properties
            if nargin==1, return, end
            set(SI,varargin{:})
       end
    end
    
    % Context menu
    methods
        function initlocalmenu(SI,hb)
            % puts a context menu in the button referenced by hb
            fn4D_dbstack
            delete(get(hb,'uicontextmenu'))
            m = uicontextmenu('parent',get(hb,'parent'));
            SI.contextmenu = m;
            set(hb,'uicontextmenu',m)

            info.si(1) = uimenu(m,'label','sliceinfo in base workspace', ...
                'callback',@(hu,evnt)assignin('base','SI',SI));
            
            set(m,'userdata',info)            
       end
    end
    
    % SET (trigger events only when changes not to identical)
    methods
        function set.ij2(SI,ij2new)
            fn4D_dbstack
            if length(ij2new)~=SI.nd, error('dimension mismatch'), end
            ij2old = SI.ij2;
            SI.ij2 = ij2new(:);
            if irrelevantchange(ij2new,ij2old), return, end % at init, lengths might be different
            % XX: before, the 2 next lines were in opposite order; was
            % there a reason for that?
            
            setij(SI)
            notifycond(SI,fn4Devent('ij2',ij2old))
        end
        
        function setij(SI)
            fn4D_dbstack
            ijold = SI.ij;
            SI.ij = min(max(round(SI.ij2),1),SI.sizes');
            if irrelevantchange(SI.ij,ijold), return, end % at init, lengths might be different
            notifycond(SI,fn4Devent('ij',ijold))
        end
        
        function set.sizes(SI,sizesnew)
            fn4D_dbstack
            % only interest is the 'fn4D_dbstack'!
            SI.sizes = sizesnew;
        end
        
        function set.sizesplus(SI,sizesplusnew)
            fn4D_dbstack
            % only interest is the 'fn4D_dbstack'!
            SI.sizesplus = sizesplusnew;
        end
        
        function slice = get.slice(SI)
            slice = SI.slicecontent;
        end
        
        function set.slice(SI,slicenew)
            fn4D_dbstack
            % shortcut: it is possible to set the 'data' field rather than
            % the complete structure
            if isnumeric(slicenew)
                data = slicenew;
                slicenew = struct('active',true,'data',data);
            end
            SI.slicecontent = slicenew;

            % change in sizes?
            if ~isempty(slicenew)
                siz = datasize(slicenew(1).data,SI.nd); % no trailing zero beyond SI.nd
                sizesnew = siz(1:SI.nd);
                if any(sizesnew~=SI.sizes)
                    sizesold = SI.sizes;
                    SI.sizes = sizesnew;
                    SI.selection = setdatasizes(SI.selection,sizesnew);
                    % ok to trigger 'sizes' event even though ij is not updated
                    % yet
                    notifycond(SI,fn4Devent('sizes',sizesold))
                    % update ij
                    setij(SI)
                end
                
                % in sizesplus?
                sizesplusnew = siz(SI.nd+1:end);
                if ~isequal(sizesplusnew,SI.sizesplus)
                    sizesplusold = SI.sizesplus;
                    SI.sizesplus = sizesplusnew;
                    notifycond(SI,fn4Devent('sizesplus',sizesplusold))
                end
            end
            
            % trigger event only after everything is set correctly
            notifycond(SI,fn4Devent('slice',[],'all'))
        end
        
        function updateslice(SI,flag,ind,value)
            % function updateslice(SI,flag,ind,value)
            %---
            % 'slic' is a structure with at least fields
            % 'active' and 'data' (representing a data slice)
            switch flag
                case 'all'
                    SI.slice = value;
                case {'new','change'}
                    % fill empty index
                    if strcmp(flag,'new') && isempty(ind)
                        ind = length(SI.slicecontent)+(1:length(value)); 
                    end
                    % check no change in size
                    siz = datasize(value(1).data,SI.nd);
                    if ~isequal(siz,[SI.sizes SI.sizesplus]), error('size mismatch'), end
                    % set value: it is a slice, i.e. a structure with
                    % fields 'active' and 'data' at least (it is obviously
                    % not a selectionND object as in 'updateselection'!)
                    nds = ndims(SI.slicecontent);
                    if nds==2, nds = ndims(value); end
                    if nds>2, SI.slicecontent = row(SI.slicecontent); end % otherwise will get error "Attempt to grow array along ambiguous dimension"
                    SI.slicecontent(ind) = value;
                    if nds>2, SI.slicecontent = shiftdim(SI.slicecontent,2-nds); end
                case 'active'
                    [SI.slicecontent(ind).active] = deal(value);
                case 'reorder'
                    perm = value;
                    SI.slicecontent = SI.slicecontent(perm);
                case 'remove'
                    SI.slicecontent(ind) = [];
                    value = []; % last argument might not have been specified
                otherwise
                    % note that flags 'add' and 'affinity' should not
                    % appear here, but be replaced by 'change'
                    error programmin 
            end    
            notifycond(SI,fn4Devent('slice',[],flag,ind,value))
        end
            
        function set.grid(SI,gridnew)
            fn4D_dbstack
            if size(gridnew,1)~=SI.nd, error('dimension mismatch'), end
            gridold = SI.grid;
            SI.grid = gridnew;
            if irrelevantchange(gridnew,gridold), return, end
            notifycond(SI,fn4Devent('grid',gridold))
        end
        
        function set.labels(SI,labelsnew)
            fn4D_dbstack
            if length(labelsnew)~=SI.nd, error('dimension mismatch'), end
            labelsold = SI.labels;
            SI.labels = labelsnew;
            if length(labelsold)~=SI.nd || isequal(labelsnew,labelsold), return, end
            notifycond(SI,fn4Devent('labels',labelsold))
        end
        
        function set.units(SI,unitsnew)
            fn4D_dbstack
            if length(unitsnew)~=SI.nd, error('dimension mismatch'), end
            unitsold = SI.units;
            SI.units = unitsnew;
            if length(unitsold)~=SI.nd || isequal(unitsnew,unitsold), return, end
            notifycond(SI,fn4Devent('units',unitsold))
        end
        
        function set.zoom(SI,zoomnew)
            fn4D_dbstack
            if size(zoomnew,1)~=SI.nd, error('dimension mismatch'), end
            zoomold = SI.zoom;
            SI.zoom = zoomnew;
            if irrelevantchange(zoomnew,zoomold), return, end
            notifycond(SI,fn4Devent('zoom',zoomold))
        end
        
        function updateselection(SI,varargin)
            % function updateselection(SI,[dims,]flag,ind,value) 
            %---
            % set the selection by updating the existing one according to
            % flag and arguments
            fn4D_dbstack
            
            % input
            if isnumeric(varargin{1})
                dims = varargin{1}; 
                varargin(1) = []; 
            elseif strcmp(varargin{1},'reset')
                dims = [];
            else
                dims = 1:SI.nd;
            end
            flag = varargin{1};
            if length(varargin)>=2, ind = varargin{2}; else ind = []; end
            ind = ind(:)';
            if length(varargin)>=3, value = varargin{3}; else value = []; end

            % specific operations
            if fn_ismemberstr(flag,{'new','add','change'}) 
                t = SI.selection.getselset(dims).singleset;
                if isempty(t), flag='new'; end
                if strcmp(flag,'new'), ind = length(t)+(1:length(value)); end
                value = ComputeInd(value,SI.sizes(dims)); % TODO: sure that there are not other existing pointers to value, such that this change in handle object is dangerous?!
                % note that in the case of 'affinity' flag, new indices are
                % computed automatically (see selectionND.selaffinity: this
                % occurs because indices were already computed in the
                % current 'value')
            end
            
            % update selection
            updateselection(SI.selection,dims,flag,ind,value);
            ev = fn4Devent('selection',dims,flag,ind,value);
            notifycond(SI,ev)
        end
        
        function setselection(SI,varargin)
            % function setselection(SI,[dims,]selset,flag,ind,value)
            %---
            % set the selection as 'selset', and triggers event according
            % to flag and arguments
            fn4D_dbstack

            % input
            if isnumeric(varargin{1})
                dims = varargin{1}; 
                varargin(1) = []; 
            else
                selset = varargin{1};
                if isa(selset,'selectionND'), selset = selectionset(selset); end
                if length(selset.t.dims)~=SI.nd, error 'this use is forbidden', end
                dims = 1:SI.nd;
            end
            selset = varargin{1};
            if isa(selset,'selectionND'), selset = selectionset(selset); end
            if isempty(selset.datasizes), setdatasizes(selset,SI.sizes(dims)), end
            if length(varargin)>=2, flag = varargin{2}; else flag = 'all'; end
            if length(varargin)>=3, ind = varargin{3}; else ind = []; end
            ind = ind(:)';
            if length(varargin)>=4, value = varargin{4}; else value = []; end
            
            % set selection
            if selset.singleton
                % set a single selectionset element
                if isconflict(SI.selection,dims)
                    error 'cannot create selections in these dimensions'
                end
                setselection(SI.selection,dims,selset)
            else
                % set a full 'selectionset'
                subdims = getdims(SI.selection,dims);
                for i=1:length(subdims), updateselection(SI.selection,subdims{i},'remove','all'); end
                subidx = getdims(selset);
                for i=1:length(subidx), setselection(SI.selection,dims(subidx{i}),getselset(selset,subidx{i})), end
            end
            notifycond(SI,fn4Devent('selection',dims,flag,ind,value))
        end
        function set.decoration(SI,deco)
            if isa(deco,'selectionND')
                deco = selectionset(deco);
            end
            SI.decoration = deco;
            notifycond(SI,fn4Devent('decoration'))
        end
    end
    
    % GET Dependent
    methods
        function ndplus = get.ndplus(SI)
            ndplus = length(SI.sizesplus);
        end
        function selectionmarks = get.selectionmarks(SI)
            selset = getselset(SI.selection,1:SI.nd);
            selectionmarks = selset.singleset;
        end
    end
        
    % Coordinates conversion
    methods
        function xy = IJ2AX(SI,ij)
            % convert indices (of the data inside the axes) into axes coordinates
            if isnumeric(ij)
                [m n] = size(ij);
                if m~=SI.nd, error('there should be %i axes coordinates',SI.nd), end
                xy = ij .* repmat(SI.grid(:,1),[1 n]) + repmat(SI.grid(:,2),[1 n]);
            else
                mat = [1 0 0; SI.grid(:,2) diag(SI.grid(:,1))];
                switch class(ij)
                    case {'selectionND' 'selectionset'}
                        xy = selaffinity(ij,mat);
                    case 'affinityND'
                        xy = referentialchange(ij,mat);
                    otherwise
                        error programming
                end
            end
        end
        function ij = AX2IJ(SI,xy)
            % convert axes coordinates into indices (of the data inside the axes)
            if isnumeric(xy)
                [m n] = size(xy);
                if m~=SI.nd, error('there should be %i axes coordinates',SI.nd), end
                ij = (xy - repmat(SI.grid(:,2),[1 n])) ./ repmat(SI.grid(:,1),[1 n]);
            else
                mat = [1 0 0; -SI.grid(:,2)./SI.grid(:,1) diag(1./SI.grid(:,1))];
                switch class(xy)
                    case {'selectionND' 'selectionset'}
                        ij = selaffinity(xy,mat);
                    case 'affinityND'
                        ij = referentialchange(xy,mat);
                    otherwise
                        error programming
                end
            end
        end
    end
        
end

%-------------
% BASIC TOOLS
%-------------

function b = irrelevantchange(x,y)

b = (any(size(x)~=size(y)) || all(x(:)==y(:)));

end

function siz = datasize(im,nd)

if iscell(im) % if SI is a 'projection3D' object
    siz = zeros(1,3);
    siz(1:2) = size(im{1});
    siz(3) = size(im{2},1);
else
    siz = size(im);
end

% remove trailing zeros after nd (which could be there if nd<2)
if nd<2 && siz(2)==1, siz=siz(1); end
if nd<1 && siz(1)==1, siz=[]; end

% make siz of length at least nd
siz(end+1:nd) = 1;

end
