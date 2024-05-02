classdef projection3D < projection
   
    properties (Access='private')
        slicetoolsubsvect
    end
    
    % Constructor, destructor, and events
    methods
        function P = projection3D(G,proj,varargin)
            fn4D_dbstack
            if length(proj)~=3
                error('projection should have 3 elements for a ''projection3D'' object')
            end
            P = P@projection(G,proj);
            % update slices upon change in 'ij'
            hl = connect_listener(P,G,'ij','PostSet',@(m,evnt)updateij(P));
            hl.Recursive = true;
            % set properties
            if nargin==2, return, end
            set(P,varargin{:})
        end
    end
    
    % Slicing
    methods
        function updateij(P)
            updateslice(P,P.datadims.orth,'kl',1)
        end

        function prepareslicing(P)
            fn4D_dbstack
            prepareslicing@projection(P)
            % instead of 1 slice of dimension 3, 3 slices of dimension 2
            % note that subsvect is a cell array of cell arrays
            subsvect = repmat({P.slicetool.subsvect},1,3);
            subsvect{1}{P.proj(3)} = P.G.ijkl(P.proj(3));
            subsvect{2}{P.proj(1)} = P.G.ijkl(P.proj(1));
            subsvect{3}{P.proj(2)} = P.G.ijkl(P.proj(2));
            P.slicetoolsubsvect = subsvect;
        end
        
        function slice = slicenosel(P,dobgsub) %#ok<INUSD>
            fn4D_dbstack
            tool = P.slicetool; d = P.datadims;
            % subsstruct is a 3-element structure
            subsstruct = struct('type','()','subs',P.slicetoolsubsvect);
            slice = struct('active',true);
            idx = {[1 2] [3 2] [1 3]; 3 1 2};
            % data
            for k=1:length(P.dataset)
                s = P.dataset(k);
                for i=1:3
                    subdata = subsref(tool.(s.name),subsstruct(i));
                    subdata = permute(subdata,[d.slice(idx{1,i}) d.plus d.slice(idx{2,i}) d.orth d.avail d.singl]);
                    slice.(s.name){i} = subdata;
                end
            end
        end
        
        function slice = slicesel(P,sel,dobgsub) %#ok<INUSD>
            fn4D_dbstack
            % sel is a selectionND object with indices set; alternatively,
            % it can just be a structure with fields active and dataind
            if isobject(sel) && isempty(sel.datasizes)
                error('indices are not set for selectionND object')
            end
            % prepare
            tool = P.slicetool; d = P.datadims;
            subsvect = tool.subsvect;
            subsvect{d.orth(1)} = sel.dataind;
            subsstruct = struct('type','()','subs',{subsvect});
            slice = struct('active',sel.active);
            % data
            for k=1:length(P.dataset)
                s = P.dataset(k);
                for i=1:3
                    subdata = subsref(tool.(s.name),subsstruct(i));
                    subdata = mean(subdata,d.orth(1));
                    subdata = permute(subdata,[d.slice{idx(i,1)} d.plus d.slice{idx(i,2)} d.orth d.avail d.singl]);
                    slice.(s.name){i} = subdata;
                end
            end
        end
        
        function eraseslice(P)
            if P.nodata, return, end
            % data became not valid - produce a 'zero' slice of the same
            % size as the previous slice
            slic = P.slice(1);
            for k=1:length(P.dataset)
                s = P.dataset(k);
                for i=1:3
                    siz = size(P.slice.(s.name){i});
                    slic.(s.name){i} = zeros(siz);
                end
            end
            P.slice = slic;
        end           
    end
    
    % GET/SET data
    methods
        function setbackground(P,name,bg,op,doupdateslice) %#ok<MANU>
            fn4D_dbstack
            if ~isempty(bg)
                error('no background definition for projection3D object')
            end
        end
    end
    
    % Local menu and its events
    methods
        function initlocalmenu(P,hb)
            initlocalmenu@projection(P,hb)
            m = P.contextmenu;
            info = get(m,'userdata');
            % no background definitions
            delete(info.hbg)
            info = rmfield(info,'hbg');
            set(m,'userdata',info)
        end
    end
    
end


