classdef selectionset < handle
    % function SET = selectionset([datasizes[,nd|set]])
    % function SET = selectionset(set)
    %---
    % if nd or set is supplied, creates a singleton set, otherwise,
    % creates a non-singleton set
   
    properties (SetAccess='private')
        singleton = false;
        datasizes        
    end
    
    properties (Dependent)
        singleset
    end
    
    properties (SetAccess='private')
        t = struct('dims',{},'set',{});
    end
    
    properties (Dependent, Access='private')
        nd
        nset
    end
    
    % Constructor and utils
    methods
        function SET = selectionset(varargin)
            if nargin==1 && isa(varargin{1},'selectionND')
                % shortcut for input
                varargin = {[] varargin{1}};
            end
            if length(varargin)>=1
                if length(varargin)>2
                    if fn_dodebug
                        error 'syntax for selectionset.selectionset has changed'
                    else
                        varargin{2} = varargin{3}; 
                    end
                end
                SET.datasizes = varargin{1}; 
                % if 'nd' or 'set' is supplied, create a singleton selection set
                if length(varargin)>=2
                    SET.singleton = true;
                    if isnumeric(varargin{2})
                        nd = varargin{2};
                        set = selectionND.empty;
                    else
                        set = varargin{2};
                        if isempty(set)
                            nd = length(SET.datasizes);
                        else
                            nd = set(1).nd;
                        end
                    end
                    if ~isempty(SET.datasizes) && nd~=length(SET.datasizes)
                        error 'mismatch in number of dimensions'
                    end
                    SET.t(1).dims = 1:nd;
                    SET.t(1).set = set;
                end
            end
        end
        
        function access(SET) %#ok<MANU>
            % gives access to private properties for debugging purposes
            keyboard
        end
        
        function SET2 = copy(SET)
            SET2 = selectionset;
            SET2.singleton = SET.singleton;
            SET2.datasizes = SET.datasizes;
            SET2.t = SET.t;
            % no need to copy the selectionND elements of t: the only
            % change which can happen to them through another pointer is a
            % change in data sizes from empty to non-empty, and such change
            % is not problematic (such change can occur only if the SET has
            % no datasizes, and at the moment datasizes are set, either the
            % selectionND elements will be already the correct size, either
            % copies will be created with the new correct size)
        end
        
        function SET2 = singleton2general(SET)
            SET2 = selectionset;
            SET2.singleton = false;
            SET2.datasizes = SET.datasizes;
            SET2.t = SET.t;
            % no need to copy the selectionND elements of t: the only
            % change which can happen to them through another pointer is a
            % change in data sizes from empty to non-empty, and such change
            % is not problematic (such change can occur only if the SET has
            % no datasizes, and at the moment datasizes are set, either the
            % selectionND elements will be already the correct size, either
            % copies will be created with the new correct size)
            if isempty(SET2.t.dims) && fn_debug
                disp 'please help me'
                keyboard
            end
       end
    end
                
    % SET
    methods
        function [SET2 chgdims] = setdatasizes(SET,sizesnew)
            sizesold = SET.datasizes;
            if isequal(sizesold,sizesnew)
                SET2 = SET;
                chgdims = {}; 
                return
            end
            if isempty(sizesold) || all(sizesold==sizesnew(1:length(sizesold)))
                SET2 = SET;
            else
                % more secure to have two distinct copies if sizes do not
                % fit
                SET2 = copy(SET);
            end
            sizesold(end+1:length(sizesnew)) = 0;
            SET2.datasizes(1:length(sizesnew)) = sizesnew;
            sizeschanged = (sizesold~=SET2.datasizes);
            
            % re-compute indices according to which sizes have changed
            chgdims = {};
            for k=1:SET2.nset
                dims = SET2.t(k).dims;
                if any(sizeschanged(dims))
                    fn4D_dbstack(['RECOMPUTE INDICES FOR ' num2str(dims)])
                    chgdims{end+1} = dims; %#ok<AGROW>
                    sizes = SET2.datasizes(dims);
                    for i=1:length(SET2.t(k).set)
                        SET2.t(k).set(i) = ComputeInd(SET2.t(k).set(i),sizes);
                    end
                end
            end
        end
        
        function setselection(SET,dims,SET2)
            % function setselection(SET,dims,SET2)
            fn4D_dbstack
            % checks - note that 'dims' can be empty if t is singleton
            if SET.singleton || ~SET2.singleton
                error('SET must be multiple and SET2 singleton')
            end
            if length(dims)~=SET2.nd, error('dimension mismatch'), end
            if isempty(SET2.datasizes) && ~isempty(SET.datasizes)
                setdatasizes(SET2,SET.datasizes(dims));
            elseif ~(isempty(SET.datasizes) && isempty(SET2.datasizes)) ...
                    && ~isequal(SET.datasizes(dims),SET2.datasizes)
                error('data sizes mismatch')
            end
            % assign
            kset  = DIMS2KSET(SET, dims);
            if ~isempty(kset)
                if isempty(SET2.t.set)
                    SET.t(kset) = [];
                else
                    SET.t(kset).set = SET2.t.set;
                end
            end
        end
        
        function updateselection(SET,dims,flag,ind,value)
            % function updateselection(SET,dims,flag,ind,value)
            fn4D_dbstack
            if strcmp(flag,'reset') % this acts on all dims, contrary to other flags
                if SET.singleton
                    SET.t.set = selectionND.empty;
                else
                    SET.t(:) = [];
                end
                return
            end
            % usually 'dims' is empty if SET is a singleton
            [kset dims] = DIMS2KSET(SET,dims);
            if isempty(kset)
                error 'incompatible selection sets are not handled yet'
            end
            % check ind
            switch flag
                case 'new'
                    if isempty(ind)
                        ind = length(SET.t(kset).set)+(1:length(value));
                    elseif length(ind)~=length(value) || any(ind>length(SET.t(kset).set)+length(value))
                        error('indices are incorrect')
                    end
                case 'all'
                case 'remove'
                    if ischar(ind) 
                        if ~strcmp(ind,'all'), error 'wrong value', end
                        ind = 1:length(SET.t(kset).set);
                    elseif any(ind>length(SET.t(kset).set))
                        disp('some indices to be removed do not exist!')
                        ind(ind>length(SET.t(kset).set)) = [];
                    end
                otherwise
                    if any(ind>length(SET.t(kset).set))
                        error('index out of bound')
                    end
            end
            % precompute
            if fn_ismemberstr(flag,{'add','change','affinity','active'})
                cursel = SET.t(kset).set(ind);
            end
            % compute indices if necessary (note that in the case of
            % 'affinity' action, this will be done automatically in the
            % selaffinity function)
            if fn_ismemberstr(flag,{'new','add','change'})
                if isempty(SET.datasizes)
                    % check on the first value
                    if ~isempty(value(1).datasizes)
                        error('cannot insert a selection with indices in a selection set without indices')
                    end
                else
                    datasizes = SET.datasizes(dims);
                    % check on the first value
                    if ~isempty(value(1).datasizes) && ~isequal(value.datasizes,datasizes)
                        error('programming: data sizes mismatch')
                    end
                    value = ComputeInd(value,datasizes); % TODO: sure that there are not other existing pointers to value, such that this change in handle object is dangerous?!
                end
            end
            % go!
            switch flag
                case 'all'
                    if isempty(value) && ~SET.singleton
                        SET.t(kset) = [];
                    else
                        SET.t(kset).set = value;
                    end
                case 'new'
                    SET.t(kset).set(ind) = value;
                case 'add'
                    SET.t(kset).set(ind) = union(cursel,value);
                case 'change'
                    % note that ind can be a vector: then cursel and value
                    % also are
                    SET.t(kset).set(ind) = substitute(cursel,value);
                case 'affinity'
                    SET.t(kset).set(ind) = selaffinity(cursel,value);
                case 'reorder'
                    SET.t(kset).set = SET.t(kset).set(value);
                case 'active'
                    for k=1:length(ind)
                        cursel(k).active = value(k);
                    end
                    SET.t(kset).set(ind) = cursel;
                case 'remove'
                    SET.t(kset).set(ind) = [];
                    if isempty(SET.t(kset).set) && ~SET.singleton
                        SET.t(kset) = [];
                    end
                case 'reset'
                    if SET.singleton
                        SET.t(kset).set = selectionND.empty;
                    else
                        SET.t(kset) = [];
                    end
            end
        end
    end
    
    % GET dependent, private
    methods
        function nset = get.nset(SET)
            nset = length(SET.t);
        end
        
        function nd = get.nd(SET)
            nd = zeros(1,SET.nset);
            for k=1:SET.nset
                nd(k) = length(SET.t(k).dims);
            end
        end
    end
    
    % GET public
    methods
        function sel = getsel(SET,dims,ind)
            kset = DIMS2KSET(SET,dims,false);
            if isempty(kset)
                sel = selectionND.empty(1,0);
            else
                sel = SET.t(kset).set;
            end
            if nargin>=3, sel = sel(ind); end
        end
        
        function nsel = numsel(SET,dims)
            if nargin==1
                if isempty(SET.t), nsel=0; return, end
                if ~isscalar(SET.t), error 'dims not specified', end
                kset = 1;
            else
                kset = DIMS2KSET(SET,dims,false);
            end
            if isempty(kset)
                nsel = 0;
            else
                nsel = length(SET.t(kset).set);
            end
        end
        
        function b = vide(SET,dims)
            if nargin==1
                b = true;
                for k=1:SET.nset
                    if ~isempty(SET.t(k).set)
                        b = false;
                        return
                    end
                end
            else
                kset = DIMS2KSET(SET,dims,false);
                b = isempty(kset);
                if ~b && isempty(SET.t(kset).set)
                    if fn_dodebug, disp 'found an empty set, should have been removed', end
                    SET.t(kset) = [];
                    b = true;
                end
            end
        end
        
        function b = isconflict(SET,dims)
            % do the dimensions 'dims' intersect with an existing selection
            if SET.singleton, error('''isconflict'' method applies only to non-singleton sets'), end
            for k=1:SET.nset
                dimsk = SET.t(k).dims;
                if ~isempty(intersect(dimsk,dims)) && ~isequal(dimsk,dims)
                    b = true;
                    return
                end
            end
            b = false;
        end
        
        function SET2 = getselset(SET,dims)    
            % function SET2 = getselset(SET,dims)    
            
            % would be absurd that SET is already a singleton set
            if SET.singleton, error('SET must be a multiple selection set'), end
            % SET2 is a singleton set
            if isempty(SET.datasizes)
                sizes=[]; 
            else
                sizes = SET.datasizes;
                sizes(end+1:max(dims)) = 1;
                sizes = sizes(dims); 
            end
            SET2 = selectionset(sizes,length(dims));
            kset = DIMS2KSET(SET,dims,false);
            if ~isempty(kset), SET2.t.set  = SET.t(kset).set; end
        end
        
        function SET2 = getselsets(SET,indims)
            % return sets whose dimensions are included in 'indims'
            % throws error if it encounters one set which dimensions partly
            % intersect 'indims'
            % Attention! the dimensions of each individual set are modified
            % so as to indicate which dimesions WITHIN indims; for example,
            % if SET has subsets with dimensions [1 2], [3 4] and 5 and
            % indims is [1 2 5], then SET2 will have two subsets with
            % dimensions respectively [1 2] and 3!
            
            % would be absurd that SET is a singleton set
            if SET.singleton, error('SET must be a multiple selection set'), end
            
            nsets = length(SET.t);
            takesets = false(1,nsets);
            for k=1:nsets
                kdims = SET.t(k).dims;
                if isempty(SET.t(k).set)
                    if fn_dodebug, disp 'found an empty set, should have been removed', end
                    continue
                end
                if ~isempty(intersect(kdims,indims))
                    if ~isempty(setdiff(kdims,indims))
                        error('found a set which is neither orthogonal to, neither included in ''indims''')
                    end
                    takesets(k) = true;
                end
            end
            
            siz = SET.datasizes;
            if ~isempty(siz), siz = siz(indims); end
            SET2 = selectionset(siz);
            f = find(takesets);
            for k=1:length(f)
                SET2.t(k) = SET.t(f(k));
                SET2.t(k).dims = find(ismember(indims,SET2.t(k).dims));
            end
        end
        
        function nd = numdims(SET)
           if ~SET.singleton, error('use ''getdims'' on non-singleton set'), end 
           nd = length(SET.t.dims);
        end
        
        function dims = getdims(SET,indims)
            if SET.singleton, error('use ''numdims'' on sigleton set'), end
            if nargin == 1
                takesets = true(1,SET.nset);
            else
                takesets = false(1,SET.nset);
                for k=1:SET.nset
                    dimsk = SET.t(k).dims;
                    if ~isempty(intersect(dimsk,indims))
                        if ~isempty(setdiff(dimsk,indims))
                            error('found a set which is neither orthogonal to, neither included in ''indims''')
                        end
                        if isempty(SET.t(k).set)
                            % should not happen since empty sets are
                            % discarded
                            error programming
                        end
                        takesets(k) = true;
                    end
                end
            end
            dims = {SET.t(takesets).dims};
        end
        
        function set = get.singleset(SET)
            if SET.singleton
                set = SET.t.set;
            else
                set = [];
            end
        end
    end
    
    % Operations
    methods
        function SET2 = selsetaffinity(SET,mat,sizes2)
            if nargin<3, sizes2=[]; end
            SET2 = selectionset(sizes2);
            SET2.singleton = SET.singleton;
            % loop on selections
            k2 = 1;
            for k=1:length(SET.t)
                % dims and check mat
                dims = SET.t(k).dims;
                dims2 = find(any(mat(2:end,1+dims),2))';
                if length(dims2)>length(dims)
                    error('wrong affinity matrix')
                elseif length(dims2)<length(dims)
                    % selection in world dimensions not spanned by the data
                    if ~isempty(dims2), disp('OVERLAPPING SELECTION'), end
                    continue
                end
                SET2.t(k2).dims = dims2;
                % affinity and compute indices if necessary
                if isempty(sizes2), siz=[]; else siz = sizes2(dims2); end
                ma = mat([1 1+dims2],[1 1+dims]);
                for i=1:length(SET.t(k).set)
                    SET2.t(k2).set(i) = selaffinity(SET.t(k).set(i),ma,siz);
                end
                k2 = k2+1;
            end
        end
    end
    
    % dims2kset
    methods (Access='private')
        function [kset dims] = DIMS2KSET(SET,dims,createflag)
            % return the index in SET.t corresponding to the set of real
            % world dimensions 'dims'; creates a new index if it does not exist yet,
            % provided this does not conflict with existing indices
            if nargin<3, createflag = true; end
            
            % case of singleton set
            if SET.singleton
                if ~isempty(dims) && ~isequal(dims,SET.t.dims)
                    error('wrong dimensions for access to singleton selection set')
                end
                kset = 1;
                dims = SET.t.dims;
                return
            end
            
            % otherwise dims cannot be empty
            if isempty(dims)
                error('empty dims')
            end
            
            % look for the good set
            kset = find(fn_map(@(d)isequal(d,dims),{SET.t.dims},'array'),1);
            if ~isempty(kset), return, end
            
            % no selection set found corresponding to dims -> create new
            % one if specified -> first check compatibility
            if ~createflag, kset = []; return, end
            for k=1:length(SET.t)
                if ~isempty(intersect(dims,SET.t(k).dims))
                    disp('incompatible selection sets')
                    kset = [];
                    return
                end
            end
            kset = length(SET.t)+1;
            SET.t(kset).dims = dims;
            SET.t(kset).set = selectionND.empty;
        end
    end
    
end

