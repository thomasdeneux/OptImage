classdef passivedisplayRegister < fn4Dhandle
    % function D = passivedisplay(obj,proplisten,callback)
    %---
    % Registers a passive callback when some properties of object 'obj'
    % change. Callback can be a char array or a function with no argument.
   
    properties (SetAccess = 'private')
        parent
        P2D
        proplisten
    end
    properties
        callback
    end
    
    
    % Constructor
    methods
        function D = passivedisplayRegister(obj,proplisten,callback)
            D.parent = obj;
            % communication with parent
            addparent(D,obj)
            if ~iscell(proplisten), proplisten = {proplisten}; end
            D.proplisten = proplisten;
            D.callback = callback;
        end
    end
    
    % Callback
    methods
        function updateDown(D,~,evnt)
            fn4D_dbstack(['S2D ' evnt.flag])
            if ~fn_ismemberstr(evnt.flag,D.proplisten), return, end
            if ischar(D.callback)
                evalin('base',D.callback)
            else
                feval(D.callback)
            end
        end
    end
end