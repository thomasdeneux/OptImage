classdef rotation < geometry
    % function R = rotation(F,varargin)
    %---
    % automatic link between a 'focus' object and a 'geometry' object
    %
    % The core of a rotation object is its 'mat' property, that defines the
    % affine coordinates transformation to get real world coordinates from
    % data coordinates. It is a (1+m)x(1+n) matrix, where n is the data
    % number of dimensions, and m the real world number of dimensions
    % (n<=m). More precisely, R.mat is equal to [1 0]
    %                                           [t r], where:
    % t is the (mx1) vector defining the translation part, and r is the
    % (mxn) matrix defining the linear part of the transformation. This
    % linear part must defined an orthogonal transformation, i.e. columns
    % of r are orthogonal to each other.
    %
    % When setting R.mat, a number of different formats are accepted
    % - []      empty matrix: R.mat is set to the identity matrix
    % - M       a matrix (the rotation-translation function); if its number
    %           of rows is <1+m or its number of columns is <1+n, it is
    %           appended in the most intuitive way
    % - {rot [trans [perm]]}    a cell array defining separately the
    %                           scaling, translation, and coordinate
    %                           permutation 
    %   rot     a vector of length <=n (scaling)
    %           or a square matrix of size <=n (any linear function)
    %           note that if size or length is <n, coordinate changes for
    %           the missing dimensions are set to the identity function
    %   trans   a vector of length <=n (translation)
    %           if its length is <n, in particular, if it is empty, it is
    %           assumed that there is no translation in the missing
    %           dimensions
    %   perm    a vector of length <=n (permutation)
    %           it says to which dimension of the real world does each data
    %           dimension correspond
    %           if its length is <n, the missing data dimensions are
    %           assumed to correspond to real world dimensions that follow
    %           the higher defined dimension
    
    
    properties (SetAccess='private')
        F
        RotationInMat = false;        
    end
    
    % transformation parameters
    properties
        mat            
        linkselection = true;
    end
    
    % dependent properties
    properties (SetAccess='private')
        mat1
        worlddims
    end
    
    % Constructor, destructor, and events
    methods
        function R = rotation(F,varargin)
            fn4D_dbstack
            R = R@geometry;
            R.F = F;            
            % communication with parent
            addparent(R,F)
            % set properties
            if nargin==1, return, end
            set(R,varargin{:})
        end
    end
    
    % Update upon events
    methods
        function updateDown(R,~,evnt)
            fn4D_dbstack(['F2G ' evnt.flag])
            switch evnt.flag
                case 'nd'
                    % update mat (automatic re-calculations will take place)
                    R.mat = R.mat;
                case 'grid'
                    if fn_dodebug, error 'how can this case happen, focus has no ''grid'' property!', end
                    grid = GridLabels(R);
                    R.grid = grid;
                case 'labels'
                    [grid labels] = GridLabels(R); %#ok<ASGLU>
                    R.labels = labels;
                case 'units'
                    [grid labels units] = GridLabels(R); %#ok<ASGLU>
                    R.units = units;
                case 'xyzt'
                    R.ijkl2 = conv_world2data(R,R.F.xyzt);
                case 'zoom'
                    R.zoom = convzoom_world2data(R,R.F.zoom);
                case 'selection'
                    if ~R.linkselection, return, end
                    if strcmp(evnt.selflag,'all') && isempty(evnt.dims) && fn_dodebug
                        disp 'please check this'
                        keyboard
                    end
                    if strcmp(evnt.selflag,'reset')
                        updateselection(R,[],evnt.selflag)
                    else
                        if ~all(ismember(evnt.dims,R.worlddims))
                            % selection is not in dimensions spanned by data
                            if any(ismember(evnt.dims,R.worlddims))
                                error('overlapping selection')
                            end
                            return
                        end
                        datadims = find(ismember(R.worlddims,evnt.dims));
                        switch evnt.selflag
                            case 'all'
                                if isempty(evnt.dims)
                                    setselectionall(R,selsetaffinity(R.F.selection,R.mat1,R.sizes))
                                else
                                    ok = ismember(evnt.dims,R.worlddims);
                                    if ~any(ok), return, elseif ~all(ok), error 'selection not allowed', end
                                    ndims = length(evnt.dims);
                                    datadims = zeros(1,ndims); for k=1:ndims, datadims(k) = find(R.worlddims==evnt.dims(k),1); end
                                    FSET = getselset(R.F.selection,evnt.dims);
                                    RSET = selsetaffinity(FSET,R.mat1([1 1+evnt.dims],[1 1+datadims]),R.sizes(datadims));
                                    setselectiondims(R,datadims,RSET,'all')
                                end
                            case {'new','add','change'}
                                mat1 = R.mat1([1 1+datadims],[1 1+evnt.dims]); %#ok<PROP>
                                value = selaffinity(evnt.value,mat1); %#ok<PROP>
                                value = ComputeInd(value,R.sizes(datadims));
                                updateselection(R,datadims,evnt.selflag,evnt.ind,value)
                            case 'affinity'
                                mat1 = R.mat1([1 1+datadims],[1 1+evnt.dims]); %#ok<PROP>
                                value = referentialchange(evnt.value,mat1); %#ok<PROP>
                                updateselection(R,datadims,evnt.selflag,evnt.ind,value)
                            otherwise
                                value = evnt.value;
                                updateselection(R,datadims,evnt.selflag,evnt.ind,value)
                        end
                    end
            end
        end
        
        function updateUp(R,evnt)
            fn4D_dbstack(['G2F ' evnt.flag])
            switch evnt.flag
                case 'nddata'
                    % update mat (automatic re-calculations will take place)
                    R.mat = R.mat;
                case 'ijkl2'
                    xyzt = conv_data2world(R,R.ijkl2);
                    R.F.xyzt(R.worlddims) = xyzt(R.worlddims);
                case 'zoom'
                    zoom = convzoom_data2world(R,R.zoom);
                    R.F.zoom(R.worlddims,:) = zoom(R.worlddims,:);
                case 'selection'
                    if ~R.linkselection, return, end
                    switch evnt.selflag
                        case 'indices'
                            % nothing to change in F (no indices for focus
                            % selections)
                        case 'all'
                            if isempty(evnt.dims)
                                error 'not implemented yet' % i am not sure whether this case really happens
                            else
                                SET = getselset(R.selection,evnt.dims);
                                m = R.mat([1 1+R.worlddims(evnt.dims)],[1 1+evnt.dims]);
                                SET = selsetaffinity(SET,m);
                                setselection(R.F,R.worlddims(evnt.dims),SET)
                            end
                        case {'new','add','change'}
                            m = R.mat([1 1+R.worlddims(evnt.dims)],[1 1+evnt.dims]);
                            value = selaffinity(evnt.value,m);
                            updateselection(R.F,R.worlddims(evnt.dims),evnt.selflag,evnt.ind,value)
                        case 'affinity'
                            m = R.mat([1 1+R.worlddims(evnt.dims)],[1 1+evnt.dims]);
                            value = referentialchange(evnt.value,m);
                            updateselection(R.F,R.worlddims(evnt.dims),evnt.selflag,evnt.ind,value)
                        otherwise
                            value = evnt.value;
                            updateselection(R.F,R.worlddims(evnt.dims),evnt.selflag,evnt.ind,value)
                    end
            end
        end
    end
    
    % SET
    methods
        function set.mat(R,mat)
            % function set.mat(R,mat)
            %---
            % Check the help of rotation class for possibles syntaxes for
            % mat
            
            fn4D_dbstack

            % parent focus does not need to know about the changes here
            R.upnotify = false;

            f = R.F; 
            matold = R.mat;
            [R.mat nd nddata] = buildMat(mat,f.nd,R.nddata);
            
            R.mat1 = rotation.invaffinity(R.mat);
            %R.worlddims = find(any(R.mat(2:end,2:end),2))'; % BAD because looses the specific order of dimensions
            if nddata~=size(R.mat,2)-1, error programming, end
            worldd = cell(1,nddata);
            for i=1:nddata, worldd{i}=row(find(R.mat(2:end,1+i))); end
            if ~all(fn_map(@isscalar,worldd)), disp 'beware: rotation handling might still have some bugs', end
            worldd = unique([worldd{:}],'stable');
            if length(worldd)~=nddata, error 'problem with rotation matrix', end
            R.worlddims = worldd;
            R.RotationInMat = any(sum(logical(R.mat(2:end,2:end)))>1);
            
            if nd>f.nd
                % change number of world dimensions
                state = R.downnotify;
                R.downnotify = false; % updateDown must not be executed (would launch a new calculation of R.mat)
                f.nd = nd;
                R.downnotify = state;
            end
            
            if nddata>R.nddata
                % change number of data dimensions
                R.nddata = nddata;
            end
            
            if ~isequal(R.mat,buildMat(matold,f.nd,R.nddata))
                % functional change in mat -> update lot of things
                [R.grid R.labels R.units] = GridLabels(R);
                R.ijkl2 = conv_world2data(R,f.xyzt);
                R.zoom = convzoom_world2data(R,f.zoom);
                if R.linkselection
                    R.upnotify = false;
                    setselectionall(R,selsetaffinity(f.selection,R.mat1,R.sizes)) %#ok<*MCSUP>
                    R.upnotify = true;
                end
            elseif nddata>size(matold,2)-1
                % no functional change but more data dimensions -> update
                % less things
                [R.grid R.labels R.units] = GridLabels(R);
                R.ijkl2 = conv_world2data(R,f.xyzt);
                R.zoom = convzoom_world2data(R,f.zoom);
                if R.linkselection
                    % world dimensions for which we already have selections
                    olddims = find(any(matold(2:end,2:end),2))';
                    olddims = getdims(f.selection,olddims);
                    olddims = [olddims{:}];
                    % additional selections we need to add
                    newdims = setdiff(R.worlddims,olddims);
                    newworldsel = getselsets(f.selection,newdims);
                    setselectionadd(R,selsetaffinity(newworldsel,R.mat1,R.sizes));
                end
            end  
            
            % re-establish notification
            R.upnotify = true;
        end
    end
    
    % Tools (available as public static methods)
    methods (Static)
        function y = affinity(M,x,flag)
            % function y = affinity(M,x[,flag])
            %---
            % Input
            % - x       vertical n-elements vector (or array: operation applied
            %           independently for each column)
            % - M       (m+1,n+1) affinity matrix
            % - flag    checking level:
            %           - 'permissive': the function does not generate an error when x
            %           is not of length n, but either truncate it or pads it with zeros
            %           - 'normal' [default]: the function checks the length of x
            %           - 'strict': the function also checks the first row of M
            %
            % Output
            % - y   vertical m-elements vector (tq [y; 1] = M * [x; 1])
            
            if nargin<3
                flag = 'normal';
            end
            
            [m1 n1] = size(M);
            n = n1-1;
            [mx nx] = size(x);
            if mx~=n
                if strcmp(flag,'permissive')
                    x = x(1:min(mx,n),:);
                    x(mx+1:n,1:nx) = 0;
                else
                    error('x should be one element less than the number of columns of M')
                end
            end
            if strcmp(flag,'strict') && ~(M(1,1)==1 && all(M(1,2:n1)==0))
                error('first row of M should be a 1 followed be 0''s')
            end
            
            x = [ones(1,nx); x];
            y = M * x;
            y = y(2:end,:);
            
        end
        
        %---
        function mat1 = invaffinity(mat)
            % function mat1 = invaffinity(mat)
            %--
            % 'inverts' an affinity matrix
            % we need to have at the end
            % * mat*mat1*[1; 0; ..; 0] = mat*mat1*[1; 0; ..; 0] = [1; 0; ..; 0]
            % * if 'mat' describes an orthogonal projection, so must do 'mat1'
            
            m1 = size(mat,1);
            
            rotation = mat(2:end,2:end);
            translation = mat(2:end,1);
            
            rotation1 = pinv(rotation);
            translation1 = -rotation1*translation;
            
            mat1 = [1 zeros(1,m1-1); translation1 rotation1];
            
        end
    end
end

%-------
% TOOLS
%-------
function [rotgrid rotlabels rotunits] = GridLabels(R)

nddata = R.nddata;
rotgrid = zeros(nddata,2);
rotlabels = cell(1,nddata);
rotunits = cell(1,nddata);

for k=1:nddata
    [dims dum scales] = find(R.mat(2:end,k+1));
    if length(dims)==1
        % moving along dimension k in data exactly corresponds to
        % moving along dimension 'dims' in real world
        rotlabels{k} = R.F.labels{dims};
        rotunits{k} = R.F.units{dims};
        % parameters of the affine correspondance for the coordinates
        rotgrid(k,:) = [scales R.mat(1+dims,1)];
    else
        % dimension k in slice corresponds to a diagonal in real world
        u = unique(R.F.units(dims));
        if isscalar(u), rotunits{k} = u{1}; else rotunits{k} = 'mixed'; end
        % we still can define an affine correspondance to fit the scale
        rotgrid(k,:) = [norm(scales) -norm(scales)];
    end
end

end

%---
function ijkl2 = conv_world2data(R,xyzt)

ijkl2 = rotation.affinity(R.mat1,xyzt);

end

%---
function xyzt = conv_data2world(R,ijkl2)

xyzt = rotation.affinity(R.mat,ijkl2);

end

%---
function worldzoom = convzoom_data2world(R,datazoom)

% conversion too simple if there is any rotation
if R.RotationInMat
    disp('zoom conversion does not work when there is any rotation in ''mat''')
end
% replace Inf by big value to avoid NaNs, then put Inf again...
f = isinf(datazoom);
datazoom(f) = sign(datazoom(f))*1e100;
worldzoom = sort(rotation.affinity(R.mat,datazoom),2);
f = find(any(abs(worldzoom)>1e50,2));
worldzoom(f,:) = repmat([-Inf Inf],[length(f) 1]);

end

%---
function datazoom = convzoom_world2data(R,worldzoom)

% conversion too simple if there is any rotation
if R.RotationInMat && ~all(all(isinf(worldzoom(R.worlddims,:))))
    disp('zoom conversion is erroneous when there is any rotation in ''mat''')
end
% replace Inf by big value to avoid NaNs, then put Inf again...
f = isinf(worldzoom);
worldzoom(f) = sign(worldzoom(f))*1e100;
datazoom = sort(rotation.affinity(R.mat1,worldzoom),2);
f = find(any(abs(datazoom)>1e50,2));
datazoom(f,:) = repmat([-Inf Inf],[length(f) 1]);

end

