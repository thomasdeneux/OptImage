function [mat nd nddata] = buildMat(rotation,nd,nddata)
% function [mat nd nddata] = buildMat(rotation,nd,nddata)
%---
% [indices -> real world] transformation matrix
% IT IS LIMITED TO ORTHOGONAL + TRANSLATION
% the first row and column of rotation and rotation1 are dedicated to the
% translation
% the input 'rotation' can have several syntaxes
% note that it can have more rows than columns (injection)
% TODO: add more checks on rotation according to nd

if ~iscell(rotation) && isvector(rotation), rotation = {rotation}; end

if iscell(rotation)
    % scaling/rotation [+ translation [+ permutation]]
    % 1) scaling or rotation -> (1+nddata) square matrix
    rot = rotation{1};
    if isvector(rot)
        rot = diag(rot);
    end
    [m n] = size(rot);
    if m~=n, error('first cell element (rotation) must be a vector or square matrix'), end
    if n>nddata
        nddata = n;
    else
        rot(n+1:nddata,n+1:nddata) = eye(nddata-n);        
    end
    % 2) translation -> update (1+nddata) square matrix
    if length(rotation)<2 || isempty(rotation{2})
        trans = -diag(rot); % coordinates start at zero
    else
        trans = rotation{2};
        if ~isvector(trans), error('wrong rotation format'), end
        trans = rotation{2}(:);
        n = length(trans);
        if n>nddata
            rot(nddata+1:n,nddata+1:n) = eye(n-nddata);
            nddata = n;
        else
            trans(n+1:nddata,1) = -diag(rot(n+1:nddata,n+1:nddata));
        end
    end
    mat0 = [1 zeros(1,nddata); trans rot];
    % 3) permutation: [world dim for first data dim, ...]
    if length(rotation)<3 || isempty(rotation{3})
        perm = 1:nddata;
    else
        perm = rotation{3};
        if ~isvector(perm), error('wrong rotation format'), end
        n = length(perm);
        if n>nddata
            mat0(1+nddata+1:1+n,1+nddata+1:1+n) = eye(n-nddata);
            nddata = n;
        else
            perm(n+1:nddata) = max(perm)+(1:nddata-n);
        end
    end
    nd = max(nd,max(perm));
    mat = zeros(1+nd,1+nddata);
    mat([1 1+perm],:) = mat0;
elseif isempty(rotation)
    nd = max(nd,nddata);
    mat = eye(1+nd,1+nddata);
else
    % rotation + translation; check first row
    mat = rotation;
    [m n] = size(mat); m=m-1; n=n-1;
    % check first row
    if m<n || ~(mat(1,1)==1 && all(mat(1,2:end)==0))
        error('wrong rotation-translation matrix')
    end
    % remove last empty rows
    lastnonempty = find(any(mat,2),1,'last');
    mat = mat(1:lastnonempty,:);    
    m = lastnonempty-1;
    % update nd and nddata
    if n>nddata
        nddata = n;
    else
        mat(1+m+1:1+m+(nddata-n),1+n+1:1+nddata) = eye(nddata-n);
        m = size(mat,1)-1;
    end
    if m>nd
        nd = m;
    else
        mat(1+m+1:1+nd,:) = 0;
    end
end

% Final check: orthogonality
[m n] = size(mat); m=m-1; n=n-1;
if ~(m==nd && n==nddata && mat(1,1)==1 && all(mat(1,2:end)==0))
    error('programming: [indices ->real world] transformation')
else
    transf = mat(2:1+nd,2:1+nddata);
    prods = transf'*transf;
    % orthogonality - any non-zero dot product in the upper triangle of
    % 'prods'?
    if any(triu(prods,1))
        error('[indices ->real world] transformation is not orthogonal')
    end
end


