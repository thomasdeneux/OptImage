classdef agglomerate < sliceinfo
    % function A = agglomerate(SI,agglom,varargin)
    
    % communication properties
    properties (SetAccess='private')
        P       % parent sliceinfo (usually a projection object in fact)
    end
    
    % primary properties
    properties
        agglom = cell(1,0);    % 'agglomerate' the output dimension n to (n-1)
    end    
    
    % Constructor, destructor, and events
    methods
        function A = agglomerate(P,agglom,varargin)
            fn4D_dbstack
            A = A@sliceinfo(length(agglom));
            A.P = P;
            % communication with parent
            addparent(A,A.P)
            % set agglomerate - automatic settings and changes happen
            A.agglom = agglom;
            % set properties
            if nargin==2, return, end
            set(A,varargin{:})
        end
        
    end
    
    % Update upon events and new definition of the agglomeration
    methods
        function set.agglom(A,agglom)
            % check
            if ~isequal(sort([agglom{:}]),1:A.P.nd) %#ok<*MCSUP>
                error('agglomeration definition must be equal to 1:%i once concatenated',A.P.nd)
            end
            % set
            A.agglom = agglom;
            % update child
            p2s_slice(A,'all')
            p2s_labels(A)
            p2s_units(A)
            p2s_ij2(A)
            p2s_zoom(A)
            p2s_selection(A,'all')
        end
        function updateDown(A,~,evnt)
            fn4D_dbstack
            switch evnt.flag
                case 'slice'
                    p2s_slice(A,evnt.selflag,evnt.ind,evnt.value)
                case 'labels'
                    p2s_labels(A)
                case 'units'
                    p2s_units(A)
                case 'ij2'
                    p2s_ij2(A)
                case 'zoom'
                    p2s_zoom(A)
                case 'selection'
                    p2s_selection(A,evnt.selflag,evnt.ind,evnt.value)
            end
        end     
        function updateUp(A,evnt)
            fn4D_dbstack
            switch evnt.flag
                case 'slice'
                    % user has changed slice locally - this is not a good
                    % practice, but any way there is nothing to do about it
                case {'labels' 'units'}
                    % do nothing: labels or units can be over-written
                    % locally in the child sliceinfo, without changing
                    % units in the parent sliceinfo
                case 'ij2'
                    s2p_ij2(A)
                case 'zoom'
                    s2p_zoom(A)
                case 'selection'
                    s2p_selection(A)
            end
        end
    end
    
    % Bi-directional updates
    methods
        function p2s_slice(A,flag,ind,value)
            fn4D_dbstack
          
            if strcmp(flag,'all'), ind=1:length(A.P.slice); end
            switch flag
                case {'all' 'new' 'change'}
                    slic = A.P.slice;
                    perm = [A.agglom{:}];
                    sizp = A.P.sizes;
                    siz = zeros(1,A.nd);
                    for k=1:A.nd
                        siz(k) = prod(sizp(A.agglom{k}));
                    end
                    for i=ind(:)'
                        slic(i).data = reshape(permute(slic(i).data,perm),siz);
                    end
                    updateslice(A,flag,ind,slic)
                case {'active' 'reorder' 'remove'}
                    updateslice(A,flag,ind,value)
                otherwise
                    error programming
            end         
        end
        function p2s_labels(A)
            lab = cell(1,A.nd);
            for k=1:A.nd
                a = A.agglom{k};
                str = A.P.labels{a(1)};
                for j=2:length(a)
                    str = [str '*' A.P.labels{a(j)}]; %#ok<AGROW>
                end
                lab{k} = str;
            end
            A.labels = lab;
        end
        function p2s_units(A)
            unt = cell(1,A.nd);
            for k=1:A.nd
                a = A.agglom{k};
                unt{k} = A.P.units{a(1)};
            end
            A.units = unt;
        end
        function p2s_ij2(A)
            x = zeros(A.nd,1);
            y = A.P.ij2-1;
            for k=1:A.nd
                a = A.agglom{k};
                fact = 1;
                x(k) = 0;
                for j=a
                    x(k) = x(k) + y(j)*fact;
                    fact = fact*A.P.sizes(j);
                end
            end            
            A.ij2 = 1+x;
        end
        function s2p_ij2(A)
            x = A.ij2-1;
            y = zeros(A.P.nd,1);
            for k=1:A.nd
                a = A.agglom{k};
                for j=a
                    sj = A.P.sizes(j);
                    y(j) = mod(x(k),sj);
                    x(k) = (x(k)-y(j))/sj;
                end
            end
            A.P.ij2 = 1+y;
        end
        function p2s_zoom(A)
            z = A.zoom;
            for k=1:A.nd
                a = A.agglom{k};
                if length(a)>1, continue, end % zooming is independent in agglomerated dimensions
                z(k,:) = A.P.zoom(a,:);
            end
            A.zoom = z;
        end
        function s2p_zoom(A)
            z = A.P.zoom;
            for k=1:A.nd
                a = A.agglom{k};
                if length(a)>1, continue, end % zooming is independent in agglomerated dimensions
                z(a,:) = A.zoom(k,:);
            end
            A.P.zoom = z;
        end
        function p2s_selection(A,flag,ind,value) %#ok<MANU,INUSD>
            disp('selections not handled yet')
        end
        function s2p_selection(A,flag,ind,value) %#ok<MANU,INUSD>
            disp('selections not handled yet')
        end
    end
    
    % Local menu and its events
    methods
        function initlocalmenu(A,hb)
            initlocalmenu@sliceinfo(A,hb)
            m = A.contextmenu;
            info = get(m,'userdata');
            
            set(m,'userdata',info)
        end
    end
end

