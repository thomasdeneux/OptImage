function [b delay] = tpv_bidircorr(a,delay,shapeflag)
% function [b delay] = tpv_bidircorr(a)
% function b = tpv_bidirectcorr(a,delay[,'same|valid'])
%---
% correct the misalignment in images due to bidirectional scanning
% default shapeflag is 'valid'

if nargin>=2
    if nargin<3, shapeflag = 'valid'; end
    b = realign(a,delay,shapeflag);
end

%---
function b = realign(a,delay,shapeflag)

[nx ny nt] = size(a);
xx = 1:nx;
xxb = 1+ceil(abs(delay)):nx-ceil(abs(delay));
nxb = length(xxb);

d = floor(delay);
u = delay-d;
v = 1-u;

switch shapeflag
    case 'same'
        b = nan(nx,ny,nt,class(a));
        idx = xxb;
    case 'valid'
        b = zeros(nxb,ny,nt,class(a));
        idx = 1:nxb;
end

if u==0
    b(idx,1:2:end,1:2:end) = a(xxb-delay,1:2:end,1:2:end);
    b(idx,2:2:end,1:2:end) = a(xxb+delay,2:2:end,1:2:end);
    b(idx,1:2:end,2:2:end) = a(xxb+delay,1:2:end,2:2:end);
    b(idx,2:2:end,2:2:end) = a(xxb-delay,2:2:end,2:2:end);
else
    b(idx,1:2:end,1:2:end) = u*a(xxb-d-1,1:2:end,1:2:end) + v*a(xxb-d,1:2:end,1:2:end);
    b(idx,2:2:end,1:2:end) = u*a(xxb+d+1,2:2:end,1:2:end) + v*a(xxb+d,2:2:end,1:2:end);    
    b(idx,1:2:end,2:2:end) = u*a(xxb+d+1,1:2:end,2:2:end) + v*a(xxb+d,1:2:end,2:2:end);
    b(idx,2:2:end,2:2:end) = u*a(xxb-d-1,2:2:end,2:2:end) + v*a(xxb-d,2:2:end,2:2:end);    
end


