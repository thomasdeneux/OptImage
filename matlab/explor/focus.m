classdef focus < fn4Dhandle
    % function F = focus(varargin)
    % function F = focus(key)
    
    properties
        nd = 0;
        labels = {};
        units = {};
        
        xyzt = zeros(0,1);
        zoom = zeros(0,2);
        selection 
    end 
    
    % Constructor
    methods
        function F = focus(varargin)
            fn4D_dbstack
            % do not do the initialization in properties, since this
            % initialization value will change (handle object)!
            if mod(nargin,2)
                key = varargin{1};
                varargin(1) = [];
                F = focus.find(key);
            end
            F.selection = selectionset;
            % set properties
            if isempty(varargin), return, end
            set(F,varargin{:})
        end
        function delete(F)
            % unregister
            focus.find(F)
        end
    end
    
    % SET
    methods
        function set.nd(F,ndnew)
            fn4D_dbstack
            ndold = F.nd;
            if ndnew<=ndold, return, end
            F.nd = ndnew;
            % now update xyzt, labels and units (this does not trigger
            % 'xyzt' or 'labels' events)
            F.downnotify = false;
            F.xyzt(ndold+1:ndnew) = 0;
            F.zoom(ndold+1:ndnew,:) = repmat([-Inf Inf],[ndnew-ndold 1]);
            [F.labels{ndold+1:ndnew}] = deal('');
            [F.units{ndold+1:ndnew}] = deal('');
            F.downnotify = true;
            % notify change only after everything has been set correctly
            notifycond(F,fn4Devent('nd',ndold))
        end
        
        function set.labels(F,labelsnew)
            fn4D_dbstack
            ndnew = length(labelsnew);
            if ndnew>F.nd, F.nd = ndnew; end
            labelsold = F.labels;
            ndold = length(labelsold);
            for k=1:ndnew
                if ~isempty(labelsnew{k}) || ndnew>ndold
                    F.labels{k} = labelsnew{k};
                end
            end
            if isequal(F.labels,labelsold), return, end
            notifycond(F,fn4Devent('labels',labelsold))
        end
        
        function set.units(F,unitsnew)
            fn4D_dbstack
            ndnew = length(unitsnew);
            if ndnew>F.nd, F.nd = ndnew; end %#ok<*MCSUP>
            unitsold = F.units;
            ndold = length(unitsold);
            for k=1:ndnew
                if ~isempty(unitsnew{k}) || ndnew>ndold
                    F.units{k} = unitsnew{k};
                end
            end
            if isequal(F.units,unitsold), return, end
            notifycond(F,fn4Devent('units',unitsold))
        end
        
        function set.xyzt(F,xyztnew)
            fn4D_dbstack
            ndnew = length(xyztnew);
            if ndnew>F.nd, F.nd = ndnew; end
            xyztold = F.xyzt;
            F.xyzt(1:ndnew,1) = xyztnew(:);
            if isequal(F.xyzt,xyztold), return, end
            notifycond(F,fn4Devent('xyzt',xyztold))
        end
        
        function set.zoom(F,zoomnew)
            fn4D_dbstack
            ndnew = size(zoomnew,1);
            if ndnew>F.nd, F.nd = ndnew; end
            zoomold = F.zoom;
            F.zoom(1:ndnew,:) = zoomnew;
            if isequal(F.zoom,zoomold), return, end
            notifycond(F,fn4Devent('zoom',zoomold))
        end
        
        function set.selection(F,selectionnew)
            fn4D_dbstack
            F.selection = selectionnew;
            notifycond(F,fn4Devent('selection',[],'all'))
        end
        
        function setselection(F,dims,SET)
            % function setselection(F,dims,SET)
            fn4D_dbstack
            F.downnotify = false; % stupid!?
            setselection(F.selection,dims,SET);
            F.downnotify = true;
            notifycond(F,fn4Devent('selection',dims,'all'))
        end
        
        function updateselection(F,dims,flag,ind,value)
            % function updateselection(F,dims,flag,ind,value)
            fn4D_dbstack
            F.downnotify = false;
            evnt = fn4Devent('selection',dims,flag,ind,value);
            if strcmp(flag,'change')
                % needed for the update of possibly linked matchpoints object(s)
                evnt.oldvalue = getsel(F.selection,dims,ind);
            end
            updateselection(F.selection,dims,flag,ind,value);
            F.downnotify = true;          
            notifycond(F,evnt)
        end
    end
    
    % SET/GET ALL
    methods
        function copyin(F,F2)
            % all the following automatically generate events
            F.nd = F2.nd;
            F.labels = F2.labels;
            F.units = F2.units;
            F.xyzt = F2.xyzt;
            F.zoom = F2.zoom;
            % for selection, we need to generate the event here
            F.selection = F2.selection;
            notifycond(F,'Selection','all')
        end
        
        function F2 = copyout(F)
            F2 = focus;
            F2.copyin(F);
        end
    end
        
    % find object according to key
    methods (Static)
        function obj = find(a)
            persistent reg
            if isempty(reg), reg = struct('key',{},'obj',{}); end
            
            if isobject(a)
                % try to unregister object
                for k=1:length(reg)
                    if reg(k).obj == a
                        reg(k) = [];
                        return
                    end
                end
            else
                % look for existing object
                key = a;
                for k=1:length(reg)
                    if isequal(reg(k).key,key)
                        obj = reg(k).obj;
                        return
                    end
                end
                % create new object if unsuccessful
                obj = focus;
                reg(end+1) = struct('key',key,'obj',obj);
            end
        end
    end
    
end
