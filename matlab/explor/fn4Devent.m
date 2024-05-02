classdef fn4Devent < event.EventData
    % function obj = myevent(flag,oldvalue)
    % function obj = myevent(flag,dims,'selflag',ind,value)
    %---
    % returns structure with fields flag, oldvalue, selflag, dims, ind,
    % value
    properties        
        flag
        selflag
        oldvalue
        dims
        ind
        value
    end
    methods
        function obj = fn4Devent(flag,varargin)
            % function obj = myevent(flag,oldvalue)
            % function obj = myevent(flag,dims,'selflag',ind,value)
            obj.flag = flag;
            switch flag
                case {'selection','slice'}
                    obj.dims = varargin{1};
                    obj.selflag = varargin{2};
                    if nargin>3, obj.ind = varargin{3}; end
                    if nargin>4, obj.value = varargin{4}; end
                case 'checkchildren'
                    if nargin>1, obj.value = varargin{1}; else obj.value = struct; end
                otherwise
                    if nargin>1, obj.oldvalue = varargin{1}; end
            end
        end
    end
end