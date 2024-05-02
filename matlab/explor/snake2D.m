classdef snake2D < sliceinfo
    % a snake2D object acts as a n-dimensional projection, but later
    % interpolates its data from a curve inside the first 2 dimensions
    
    % communication properties
    properties (SetAccess='private')
        G
        P   % n dimensions - contains the data
    end
    
    % primary properties
    properties
        dx = 1; % pixel unit!
    end    
    properties (SetAccess='private')
        type    % 'line' or 'region'
        snake   % curve
        A       % parameters for sparse interpolation matrix
        tube    % 2 parallel lines in the case of 'region' type
    end
    properties (Dependent, GetAccess='private')
        data
        background
    end
    
    % Constructor, destructor, and events
    methods
        function S = snake2D(obj,proj,varargin)
            fn4D_dbstack
            if isa(obj,'geometry')
                GG = obj;
                % P is an intermediary projection, before snake interpolation
                % will occur
                PP = projection(obj,proj);
            elseif isa(obj,'projection')
                PP = obj;
                GG = obj.G;
                proj = PP.proj;
            end
            if length(proj)<2, error('snake projection should have at least 2 dimensions'), end
            S = S@sliceinfo(length(proj)-1);
            [S.G S.P] = deal(GG,PP);
            % local properties (except data)
            setgrid(S)
            S.labels = [{'snake'} S.P.labels(3:end)];
            setunits(S)
            S.ij2(2:end) = S.P.ij2(3:end);
            S.zoom(2:end,:) = S.P.zoom(3:end,:);
            % local data
            buildsnake(S) % automatic call to interp(S)
            % communication with parent
            addparent(S, S.P)
            % set properties
            if nargin==2, return, end
            set(S,varargin{:})
        end
    end
    
    % Update upon events
    methods
        function updateDown(S,~,evnt)
            fn4D_dbstack
            pa=S.P;
            switch evnt.flag
                case 'grid'
                    setgrid(S)
                case 'labels'
                    S.labels(2:end) = pa.labels(3:end);
                case 'units'
                    setunits(S)
                case 'ij2'
                    if ~isempty(S.snake) && any(pa.ij2(1:2)~=evnt.oldvalue(1:2))
                        S.ij2(1) = project(S,pa.ij2(1:2));
                    end
                    S.ij2(2:end) = pa.ij2(3:end);
                case 'zoom'
                    S.zoom(2:end,:) = pa.zoom(3:end,:);
                case 'slice'
                    interp(S)
                case 'selection'
                    if isequal(evnt.dims,[1 2])
                        switch evnt.selflag
                            case {'new' 'change' 'add'}
                                buildsnake(S,evnt.ind)
                            case {'all' 'remove' 'reset'}
                                buildsnake(S)
                        end
                    end
            end
        end
        
        function updateUp(S,evnt)
            fn4D_dbstack
            pa=S.P;
            switch evnt.flag
                case 'labels'
                    pa.labels(3:end) = S.labels(2:end);
                case 'units'
                    % do nothing: units can be over-written locally in the
                    % sliceinfo, without changing units in the geometry
                    % parent
                case 'ij2'
                    if ~isempty(S.snake) && S.ij2(1)~=evnt.oldvalue(1)
                        pa.ij2(1:2) = merge(S,S.ij2(1));
                    end
                    pa.ij2(3:end) = S.ij2(2:end);
                case 'zoom'
                    pa.zoom(3:end,:) = S.zoom(2:end,:);
            end
        end
    end
    
    % GET/SET data
    methods
        function data = get.data(S)
            fn4D_dbstack
            data = S.P.data;
        end
        
        function set.data(S,data)
            fn4D_dbstack
            S.P.data = data;
        end
        
        function bg = get.background(S)
            fn4D_dbstack
            bg = S.P.bg;
        end
        
        function set.background(S,bg)
            fn4D_dbstack
            S.P.bg = bg;
        end
        
        function setdata(S,varargin)
            fn4D_dbstack
            setdata(S.P,varargin{:})
        end
        
        function setbackground(S,varargin)
            fn4D_dbstack
            setbackground(S.P,varargin{:})
        end
        
        function rmdata(S,name)
            rmdata(S.P,name)
        end
    end
    
    % Snake properties (except snake-specific)
    methods
        function setunits(S)
            setgrid(S)
        end
        
        function setgrid(S)
            units = S.P.units(1:2);
            grid = S.P.grid(1:2,:);
            if strcmp(units{1},units{2}) && grid(1,1)==grid(2,1)
                units = units(1);
                grid = [grid(1,1) -grid(1,1)];
            else
                units = {'pixel'};
                grid = [1 -1];
            end
            S.units = [units S.P.units(3:end)];
            S.grid = [grid; S.P.grid(3:end,:)];
        end
    end
    
    % Snake specific methods (in order of dependence)
    methods
        function set.dx(S,dx)
            if dx==S.dx, return, end
            S.dx = dx;
            buildsnake(S)
        end
        
        function buildsnake(S,ind)
            % default: nothing
            S.snake = [];
            
            % get sel
            sels = S.P.selection.getsel([1 2]);
            if nargin>=2, sels = sels(ind); end
            sel = [];
            for i=1:length(sels)
                if isscalar(sels(i).poly) && strcmp(sels(i).poly.type,'line2D')
                    sel = sels(i);
                    break
                end
            end
            if isempty(sel), eraseslice(S), return, end
            
            % snake definition
            if length(sel.poly)~=1
                error('selection should be a singleton')
            end
            switch sel.poly.type
                %                 case 'point2D'  % collection of points
                %                     error('not implemented yet')
                case 'line2D'   % interpolation along a line
                    S.type = 'line';
                    line = sel.poly.points;
                    % double resampling for accuracy
                    line = resample(line,S.dx/5);
                    S.snake = resample(line,S.dx);
                %                 case 'poly2D'   % interpolation along a tube
                %                     S.type = 'region';
                %                     % look for tube extremities
                %                     a = sel.poly.points;
                %                     a(all(diff(a,1,2)==0)) = [];
                %                     if all(a(:,end)==a(:,1)), a=a(:,1:end-1); end
                %                     b = a(:,[2:end 1]);
                %                     ab = b-a;
                %                     u = ab ./ repmat(sqrt(sum(ab.^2)),[2 1]);
                %                     turn = sum(u(:,[end 1:end-1]).*u(:,[2:end 1]));
                %                     [dum ord] = sort(turn); %#ok<ASGLU>
                %                     ext = ord([1 2]);
                %                     k1 = min(ext); k2 = max(ext);
                %                     tubea = a(:,k1+1:k2);
                %                     tubeb = a(:,[k1:-1:1 end:-1:k2+1]);
                %                     % special resampling
                %                     tubea = resample(tubea,S.dx/5);
                %                     tubeb = resample(tubeb,S.dx/5);
                %                     S.tube = resample({tubea tubeb},S.dx);
                %                     S.snake = mean(S.tube,3); % automatic updates
            end
            
            % parameters for interpolation
            np = size(S.snake,2);
            % (filter)
            ii = floor(S.snake(1,:));
            jj = floor(S.snake(2,:));
            switch S.type
                case 'line'
                    u = mod(S.snake(1,:),1);
                    v = mod(S.snake(2,:),1);
                    AA.p = [1:np 1:np 1:np 1:np];
                    AA.i = [ii ii ii+1 ii+1];
                    AA.j = [jj jj+1 jj jj+1];
                    AA.z = [(1-u).*(1-v) (1-u).*v u.*(1-v) u.*v];
                case 'region'
                    tubea = S.tube(:,:,1);
                    tubeb = S.tube(:,:,2);
                    section = tubeb-tubea;
                    sectlength = sqrt(sum(section.^2));
                    % maximum size for the filter: idx of pixel in the
                    % center
                    icenter = ceil((max(sectlength)+1)/2);
                    nf = 2*icenter;
                    % coordinates of points in the grid - dimensions are
                    % (gridi,gridj,kpoint,xy)
                    [xx yy] = ndgrid((1:nf)-icenter,(1:nf)-icenter,ones(1,np));
                    % coordinates of points in image
                    c = repmat(shiftdim(floor(S.snake'),-2),nf,nf);
                    x = cat(4,xx,yy) + c;
                    % distance to point A, point B and line
                    a = repmat(shiftdim(tubea',-2),nf,nf);
                    b = repmat(shiftdim(tubeb',-2),nf,nf);
                    ax = x-a;
                    bx = x-b;
                    da = sqrt(sum(ax.^2,4));
                    db = sqrt(sum(ax.^2,4));
                    ab = b-a;
                    dab = sqrt(sum(ab.^2,4));
                    u = ab ./ repmat(dab,[1 1 1 2]);
                    cross = u(:,:,:,1).*ax(:,:,:,2) - u(:,:,:,2).*ax(:,:,:,1);
                    dd = abs(cross);
                    % distance to segment
                    dot = sum(u.*ax,4);
                    left = (dot<0);
                    right = (dot>dab);
                    dist = left.*da + (~left & ~right).*dd + right.*db;
                    % filter
                    H = max(1-dist,0);
                    H = H ./ repmat(sum(sum(H)),nf,nf);
                    
                    % parameters for interpolation
                    p = repmat(shiftdim(1:np,-1),nf,nf);
                    f = logical(H);
                    AA.p = p(f);
                    AA.i = c(:,:,:,1)+xx; AA.i = AA.i(f);
                    AA.j = c(:,:,:,2)+yy; AA.j = AA.j(f);
                    AA.z = H(f);
            end
            
            % discard points outside of image
            bad = (AA.i<1 | AA.i>S.P.sizes(1) | AA.j<1 | AA.j>S.P.sizes(2));
            AA.p(bad) = []; AA.i(bad) = []; AA.j(bad) = []; AA.z(bad) = [];
            
            % update data 
            S.A = AA;
            interp(S)
        end
        
        function set.snake(S,snake)
            S.snake = snake;
            if ~isempty(S.snake)
                S.ij2(1) = project(S,S.P.ij2(1:2));
            end
        end
                
        function interp(S)
           if isempty(S.snake), return, end
           slice = S.P.slice;
           if isempty(slice), return, end
           
           % interpolation matrix
           np = size(S.snake,2);
           ni = S.P.sizes(1);
           nj = S.P.sizes(2);
           a = sparse(S.A.p,S.A.i+ni*(S.A.j-1),S.A.z,np,ni*nj);
            
           F = fieldnames(slice);
           for k=2:length(F) % first field, 'active' does not contain data
               f = F{k};
               img = slice.(f);
               s = size(img); 
               nplus = prod(s(3:end));
               img = reshape(img,s(1)*s(2),nplus);
               
               data2 = a*double(img);
               
               slice.(f) = reshape(data2,[np s(3:end)]);
           end
           S.slice = slice;
        end
        
        function i2 = project(S,ij2)
            np = size(S.snake,2);
            
            x = repmat(ij2(:),1,np-1);
            a = S.snake(:,1:np-1);
            b = S.snake(:,2:np);
            ab = b-a;
            ax = x-a;
            bx = x-b;
            dab = sqrt(sum(ab.^2));
            u = ab ./ repmat(dab,2,1);
            dot = ax(1,:).*u(1,:)+ax(2,:).*u(2,:);
            cross = ax(1,:).*u(2,:)-ax(2,:).*u(1,:);
            
            % distance to a segment
            da = sqrt(sum(ax.^2));
            db = sqrt(sum(bx.^2));
            dd = abs(cross);
            left = (dot<0);
            right = (dot>dab);
            dist = left.*da + (~left & ~right).*dd + right.*db;
            
            [dum k] = min(dist);
            i2 = k + dot(k);
        end
        
        function ij2 = merge(S,i2)
            np = size(S.snake,2);
            ij2 = interp1(1:np,S.snake',i2)';
        end

        function eraseslice(S)
            % no snake, or data became not valid - produce a 'zero' slice
            % of a small size
            S.slice = struct('active',true,'data',zeros(1,S.P.sizes(3:end),'uint8'));
        end                   
    end
    
    % Local menu and its events
    methods
        function initlocalmenu(S,hb)
            initlocalmenu@sliceinfo(S,hb)
            m = S.contextmenu;
            info = get(m,'userdata');

            %             info.hbg(1) = uimenu(m,'label','divide by global average','separator','on', ...
            %                 'callback',@(hu,evnt)setbackground(S,'data','avg','/'));
            %             info.hbg(2) = uimenu(m,'label','subtract global average', ...
            %                 'callback',@(hu,evnt)setbackground(S,'data','avg','-'));
            %             info.hbg(3) = uimenu(m,'label','subtract first display', ...
            %                 'callback',@(hu,evnt)setbackground(S,'data','selfirst','-'));
            %             info.hbg(4) = uimenu(m,'label','subtract last display', ...
            %                 'callback',@(hu,evnt)setbackground(S,'data','sellast','-'));
            %             info.hbg(5) = uimenu(m,'label','subtract average display', ...
            %                 'callback',@(hu,evnt)setbackground(S,'data','selall','-'));
            %             info.hbg(6) = uimenu(m,'label','cancel background normalization', ...
            %                 'callback',@(hu,evnt)setbackground(S,'data',''));

            set(m,'userdata',info)
        end
    end
end

%-------
% TOOLS
%-------

function x = resample(x,ds)

if iscell(x)
    longueur1 = [0 cumsum(sqrt(sum(diff(x{1},1,2).^2)))];
    longueur2 = [0 cumsum(sqrt(sum(diff(x{2},1,2).^2)))];
    L = (longueur1(end)+longueur2(end))/2;
else
    dx = diff(x,1,2);
    f = all(dx==0,1);
    x(:,f) = []; dx(:,f) = [];
    longueur = [0 cumsum(sqrt(sum(dx.^2)))];
    L = longueur(end);
end

s = [(0:ds:max(0,L-ds/2)) L];

if iscell(x)
    x = cat(3,interp1(longueur1,x{1}',s,'spline')',interp1(longueur2,x{2}',s,'spline')');
else
    x = interp1(longueur,x',s,'spline')';
end


end




