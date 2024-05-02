classdef cliplink < fn4Dhandle
    
    properties
        clip = [0 1];
    end

    events
        ChangeClip
    end
    
    methods
        function set.clip(CL,clip)
            if ~isequal(size(clip),[1 2])
                error('clip must have two values')
            end
            if all(clip==CL.clip), return, end
            CL.clip = clip(:)';
            notify(CL,'ChangeClip')
        end
        
        function delete(obj)
            % unregister object
            cliplink.find(obj)
        end
    end
    
    methods (Static)
        function obj = find(a,clip)
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
                % message if unsuccessful
                disp('could not find object to unregister')
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
                obj = cliplink;
                reg(end+1) = struct('key',key,'obj',obj);
                % initialize the value if specified
                if nargin==2, obj.clip = clip; end
            end
        end
    end
    
end