classdef fn4Dhandle < hgsetget
   
    % Communication between objects is asymmetric: objects have parents and
    % children. 
    % When an object is changed, its children are notified about the change
    % by a 'ChangeView' event, while the parent(s) is(are) notified by an
    % explicit call to updateUp method.
    % The main interest of this philosophy is that objects do not need to
    % be aware about their children (e.g. watch events notified by their
    % children); note that parent objects are always created before their
    % children.
    % In some exceptional cases it can be usefull however for an object to
    % list its children. This is performed through the 'CheckChildren'
    % event.
    
    properties
        id
    end
       
    properties (SetAccess='protected')
        upnotify = true;
        downnotify = true;
    end
    
    properties (SetAccess='private')
        links = struct('parent',{},'Lview',{},'Lcheck',{}); % parent and 2 listeners to parent 'ChangeView' and 'CheckChildren' events
        childlist = {};
    end
    
    events
        CheckChildren
        ChangeView
    end
    
    % Constructor, destructor, set, save
    methods
        % constructor
        function obj = fn4Dhandle
            obj.id = sprintf('%.2i',floor(rand*100));
            fn4D_dbstack(['CREATE ' obj.id ' ' class(obj)])
        end
        
        % setdisp
        function set(obj,varargin)
            if nargin<3
                desc = setinfo(obj);
                if nargin<2
                    M = metaclass(obj);
                    M = [M.Properties{:}];
                    for k=1:length(M)
                        if ~strcmp(M(k).SetAccess,'public'), continue, end
                        f = M(k).Name;
                        if isfield(desc,f), str = makestr(desc.(f)); else str=[]; end
                        if isempty(str)
                            fprintf('\t%s\n',f)
                        else
                            fprintf('\t%s: %s\n',f,str)
                        end
                    end
                else
                    f = varargin{1};
                    if isfield(desc,f), disp(makestr(desc.(f))), end
                end
            else
                if ~mod(nargin,2), error('Invalid parameter/value pair arguments'), end
                for k=1:2:length(varargin)
                    [obj.(varargin{k})] = deal(varargin{k+1});
                end
            end
        end
        function x = setinfo(obj) %#ok<MANU>
            x = struct;
        end
        
        % destructor: disconnect from parents
        function delete(obj)
            fn4D_dbstack(['DELETE ' obj.id ' ' class(obj)])
            for i=1:length(obj.links)
                A = obj.links(i).parent;
                obj.links(i).parent = [];
                delete(obj.links(i).Lview),  obj.links(i).Lview = [];
                delete(obj.links(i).Lcheck), obj.links(i).Lcheck = [];
                % delete parent if it has no more child: ideally one would
                % like to delete it only if it has no more child AND there
                % is no more Matlab variable pointing to it, but this is
                % impossible because there will always be the pointers
                % inside its own Lview listener
                if isvalid(A) && isempty(getChildren(A))
                    delete(A)
                end
            end
        end
        
        % no save!
        function dum = saveobj(obj) %#ok<MANU,STOUT>
            error('fn4Dhandle objects cannot be saved')
        end
    end
        
    % Virtual methods, must be reimplemented in inheriting classes for
    % proper communications
    methods
        function updateDown(obj,A,ev)
        end
        function updateUp(obj,ev)
        end
    end
    
    % Communication between objects
    methods
        % register a new parent
        function addparent(obj,A)
            obj.links(end+1) = struct('parent',A, ...
                'Lview',    connect_listener(A,obj,'ChangeView',@(u,ev)updateDownNoloop(obj,A,ev)), ...
                'Lcheck',   connect_listener(A,obj,'CheckChildren',@(u,ev)answerCheck(obj,A)));
        end
                
        % notification (both in 'up' and 'down' directions)
        function notifycond(obj,ev)
            fn4D_dbstack(['[event: ' class(obj) num2str(obj.id) ' ' ev.flag ']'])
            % first notify children
            if obj.downnotify
                % send a 'ChangeView' event that will be detected by
                % children; note that the loop preventing mechanism is in
                % function updateDownNoloop below
                notify(obj,'ChangeView',ev)
            else
                fn4D_dbstack('(down notification canceled)')
            end
            % then notify parent(s)
            if obj.upnotify
                % temporarily prevent notifications back down to be
                % treated
                for i=1:length(obj.links), obj.links(i).Lview.Enabled = false; end
                % update parent(s) (directly rather than through an event
                % mechanism)
                try
                    updateUp(obj,ev)
                    for i=1:length(obj.links), obj.links(i).Lview.Enabled = true; end % re-establish communication
                catch ME
                    for i=1:length(obj.links), obj.links(i).Lview.Enabled = true; end % re-establish communication
                    rethrow(ME)
                end
            else
                fn4D_dbstack('(up updating canceled)')
            end
        end
        
        % 'down' updating
        function updateDownNoloop(obj,A,ev)
            % temporarily prevent notifications back up to the parent
            obj.upnotify = false;
            % update object based on event sent by parent
            try 
                updateDown(obj,A,ev);
                obj.upnotify = true; % re-activate up-notifications
            catch ME
                if fn_dodebug
                    disp 'updating object caused an error, try the code below to inspect this object'
                    keyboard
                    typ = class(obj);
                    if strfind(typ,'activedisplay')
                      if ~ishandle(obj.ha)
                          delete(obj)
                      end
                    end
                end
                obj.upnotify = true; % re-activate up-notifications
                rethrow(ME)
            end
        end
        
        function answerCheck(obj,A)
            A.childlist{end+1} = obj;
        end
        
        function c = getChildren(obj)
            obj.childlist = {};
            notify(obj,'CheckChildren')
            c = obj.childlist;
        end
        
        % deprecated
        function disconnect(obj,obj1)
            if fn_dodebug
                disp 'function ''disconnect'' is deprecated, and all the garbage collection of fn4Dhandle objects must be revised'
            end
        end
    end
    
end


function desc = makestr(desc)

if isempty(desc)
    desc = '';
elseif iscell(desc)
    desc = [ '[' sprintf(' %s |',desc{:})];
    desc(end) = ']';
end

end
