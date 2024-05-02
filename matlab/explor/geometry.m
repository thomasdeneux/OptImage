classdef geometry < fn4Dhandle
    % defines a geometry (point selected, regions selected) in data
    % coordinates
    
    properties
        nddata = 0;
        sizes = zeros(1,0);
        grid = [];
        labels = {};
        units = {};
        ijkl2 = zeros(0,1);
        zoom = zeros(0,2);
        selection
    end
    
    % dependent properties
    properties (SetAccess='private')
        ijkl = ones(0,1);
    end
    
    % Constructor
    methods
        function G = geometry(varargin)
            fn4D_dbstack
            % do not do the initialization in properties, since this
            % initialization value will change (handle object)!
            G.selection = selectionset([]);    
            % set properties
            if nargin==0, return, end
            set(G,varargin{:})
        end
    end
    
    % SET 
    methods
        function set.nddata(G,nddatanew)
            fn4D_dbstack
            nddataold = G.nddata;
            if nddatanew<=nddataold, return, end
            G.nddata = nddatanew;
            % now update sizes, grid, labels, units, ijkl2, ijkl (but no event triggering)
            [G.upnotify G.downnotify] = deal(false);
            G.sizes(1,nddataold+1:nddatanew) = 1;
            G.grid(nddataold+1:nddatanew,:) = repmat([1 0],[nddatanew-nddataold 1]);
            [G.labels{nddataold+1:nddatanew}] = deal('');
            [G.units{nddataold+1:nddatanew}] = deal('');
            G.ijkl2(nddataold+1:nddatanew) = 0;
            G.zoom(nddataold+1:nddatanew,:) = repmat([-Inf Inf],[nddatanew-nddataold 1]);
            [G.upnotify G.downnotify] = deal(true);
            % notify change only after everything has been set correctly
            notifycond(G,fn4Devent('nddata',nddataold))
        end
        
        function set.sizes(G,sizesnew)
            fn4D_dbstack
            nddatanew = length(sizesnew);
            if nddatanew>G.nddata, G.nddata = nddatanew; end %#ok<*MCSUP>
            sizesold = G.sizes;
            G.sizes(1:nddatanew) = sizesnew;
            if isequal(G.sizes,sizesold), return, end
            % note that the change might be notified before ijkl and
            % selection indices are set correctly
            notifycond(G,fn4Devent('sizes',sizesold))
            % update ijkl (automatic coerce)
            G.ijkl = G.ijkl2;
            % update selection (recompute indices)
            setselectionsizes(G,G.sizes)
        end
        
        function set.grid(G,gridnew)
            fn4D_dbstack
            nddatanew = size(gridnew,1);
            if nddatanew>G.nddata, G.nddata = nddatanew; end
            gridold = G.grid;
            G.grid(1:nddatanew,:) = gridnew;
            if isequal(G.grid,gridold), return, end
            notifycond(G,fn4Devent('grid',gridold))
        end
        
        function set.labels(G,labelsnew)
            fn4D_dbstack
            nddatanew = length(labelsnew);
            if nddatanew>G.nddata, G.nddata = nddatanew; end
            labelsold = G.labels;
            nddataold = length(labelsold);
            for k=1:nddatanew
                if k>nddataold
                    G.labels{k} = labelsnew{k};
                elseif isempty(labelsnew{k})
                    continue
                else
                    if ~isempty(labelsold{k}) && ~strcmp(labelsnew{k},labelsold{k})
                        disp(['Replacing label ''' labelsold{k} ''' by ''' labelsnew{k} ''''])
                    end
                    G.labels{k} = labelsnew{k};
                end
            end
            if isequal(G.labels,labelsold), return, end
            notifycond(G,fn4Devent('labels',labelsold))
        end
        
        function set.units(G,unitsnew)
            fn4D_dbstack
            if ischar(unitsnew), unitsnew = repmat({unitsnew},1,G.nddata); end
            nddatanew = length(unitsnew);
            if nddatanew>G.nddata, G.nddata = nddatanew; end
            unitsold = G.units;
            nddataold = length(unitsold);
            for k=1:nddatanew
                if k>nddataold
                    G.units{k} = unitsnew{k};
                elseif isempty(unitsnew{k})
                    continue
                else
                    if ~isempty(unitsold{k}) && ~isequal(unitsnew{k},unitsold{k})
                        if ischar(unitsold{k}) && ischar(unitsnew{k})
                            disp(['Replacing unit ''' unitsold{k} ''' by ''' unitsnew{k} ''''])
                        else
                            disp 'Replacing units'
                        end
                    end
                    G.units{k} = unitsnew{k};
                end
            end
            if isequal(G.units,unitsold), return, end
            notifycond(G,fn4Devent('units',unitsold))
        end
        
        function set.ijkl2(G,ijkl2new)
            fn4D_dbstack
            nddatanew = length(ijkl2new);
            if nddatanew>G.nddata, G.nddata = nddatanew; end
            ijkl2old = G.ijkl2;
            G.ijkl2(1:nddatanew,1) = ijkl2new(:);
            if isequal(G.ijkl2,ijkl2old), return, end
            % update ijkl (automatic coerce to data range)
            G.ijkl = G.ijkl2;
            % notifycond an event only after everything is set correctly
            notifycond(G,fn4Devent('ijkl2',ijkl2old))
        end 
         
        function set.ijkl(G,ijklnew)
            fn4D_dbstack
            nddatanew = length(ijklnew);
            if nddatanew>G.nddata, error programming; end
            ijklold = G.ijkl;
            G.ijkl = min(max(round(ijklnew),1),G.sizes(1:length(ijklnew))');
            if isequal(G.ijkl,ijklold), return, end
            notifycond(G,fn4Devent('ijkl',ijklold))
        end
         
        function set.zoom(G,zoomnew)
            fn4D_dbstack
            nddatanew = size(zoomnew,1);
            if nddatanew>G.nddata, G.nddata = nddatanew; end
            zoomold = G.zoom;
            G.zoom(1:nddatanew,:) = zoomnew;
            if isequal(G.zoom,zoomold), return, end
            notifycond(G,fn4Devent('zoom',zoomold))
        end
        
        function setselectionall(G,SET)
            fn4D_dbstack
            G.selection = SET;
            notifycond(G,fn4Devent('selection',[],'all'))
        end
        
        function setselectionadd(G,SET)
            fn4D_dbstack
            dims = getdims(SET);
            for k=1:length(dims)
                SUBSET = getselset(SET,dims{k});
                setselectiondims(G,dims{k},SUBSET,'all');
            end
        end
        
        function setselectiondims(G,dims,SET,flag,varargin)
            fn4D_dbstack
            if isa(SET,'selectionND'), SET = selectionset(SET); end
            setselection(G.selection,dims,SET);
            if nargin<4
                flag='all'; 
                varargin = {1:length(SET.singleset) SET.singleset}; 
            end
            if strcmp(flag,'reset'), dims=[]; end % applies to all dims
            notifycond(G,fn4Devent('selection',dims,flag,varargin{:}))
        end
        
        function setselectionsizes(G,datasizes)
            fn4D_dbstack
            [G.selection chgdims] = setdatasizes(G.selection,datasizes);
            for k=1:length(chgdims)
                notifycond(G,fn4Devent('selection',chgdims{k},'indices'))
            end
        end
        
        function updateselection(G,dims,flag,varargin)
            fn4D_dbstack
            updateselection(G.selection,dims,flag,varargin{:});
            ev = fn4Devent('selection',dims,flag,varargin{:});
            notifycond(G,ev)
        end
    end
    
end
