function xcor = tps_correctshift(x,n,varargin)
% function xcor = tps_correctshift(x,n[,'alternate'][,'full'])

[doalternate dofull] = fn_flags('alternate','full',varargin);
dovalid = ~dofull;

[nx ny nt] = size(x); 

if n>0
    lines1 = 1:2:ny;
    lines2 = 2:2:ny;
else
    lines1 = 2:2:ny;
    lines2 = 1:2:ny;
    n = -n;
end
xcor = zeros(nx-n,ny,nt,'like',x);

if doalternate
    xcor(:,lines1,1:2:nt) = x(n+1:end,lines1,1:2:nt);
    xcor(:,lines2,1:2:nt) = x(1:end-n,lines2,1:2:nt);
    xcor(:,lines1,2:2:nt) = x(1:end-n,lines1,2:2:nt);
    xcor(:,lines2,2:2:nt) = x(n+1:end,lines2,2:2:nt);
else
    xcor(:,lines1,:) = x(n+1:end,lines1,:);
    xcor(:,lines2,:) = x(1:end-n,lines2,:);
end

if dofull
    xcor = [zeros(floor(n/2),ny,nt); xcor; zeros(ceil(n/2),ny,nt)];
end