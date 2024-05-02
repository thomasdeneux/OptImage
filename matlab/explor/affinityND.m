classdef affinityND < handle
    % function mov = affinityND('type',data)
    %---
    % Type can be 'translate2D', 'scale2D'.
   
    properties (SetAccess='private')
        id      % scalar represents the 'identity' of the transformation
        nd
        mat
    end
    
    % Constructor + Load + Display
    methods
        function mov = affinityND(type,data)
            % set id
            mov.id = mod(sum(data(:))+pi+exp(1),1); % pseudo-random number
            
            % number of dimensions
            switch type
                case {'translate1D','scale1D'}
                    mov.nd = 1;
                case {'translate2D','scale2D'}
                    mov.nd = 2;
                otherwise
                    error('unknown affinity type ''%s''',type)
            end
            
            % affinity matrix
            switch type
                case 'translate1D'
                    if ~isscalar(data), error('wrong translation data'), end
                    mov.mat = [1 0 ; data 1];
                case 'scale1D'
                    if ~isscalar(data), error('wrong scaling factor'), end
                    mov.mat = diag([1 data]);
                case 'translate2D'
                    if ~isvector(data) || length(data)~=2, error('wrong translation data'), end
                    mov.mat = [1 0 0; data(:) eye(2)];
                case 'scale2D'
                    if isscalar(data), data = [data data]; end
                    if ~isvector(data) || length(data)~=2, error('wrong scaling data'), end
                    mov.mat = diag([1 data(:)']);
            end
        end
        function disp(impossible_name__) %#ok<MANU>
            warning('off','MATLAB:structOnObject')
            varname = inputname(1);
            eval([varname ' = struct(impossible_name__);'])
            fn_structdisp(varname)
            warning('on','MATLAB:structOnObject')
        end
    end
    
    % Operations 
    methods
        function mov = referentialchange(mov1,mat)
            % function mov = movaffinity(mov1,mat)
            %---
            % scheme prevents headache: 
            %
            %   shape1,ref1   -- mat --> shape1,ref2
            %        |                        |
            %       mov1                     mov
            %        |                        |
            %        V                        V
            %   shape2, ref1  -- mat --> shape2,ref2
            % 
            % that's it!
            mov = mov1; % keep same id and, of course, nd
            mov.mat = mat*mov1.mat*mat^-1;
        end
    end
end

