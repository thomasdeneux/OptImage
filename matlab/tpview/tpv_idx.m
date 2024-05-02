function x = tpv_idx(data,V,ks,dim)
% function x = tpv_idx(data,V[,k[,dim]])
%---
% average spatially data according to the k^th spatial selection in tpview
%
% Input:
% - data        ND array (at least 3D), first 2 dimensions are spatial
% - V           tpview object
% - k           index of selection in tpview (can be a vector of indices)
%               if it is empty or not specified, all indices are used
% - dim         dimension where to make appear the several selections (in
%               the case of a vector of indices) in the result
%
% Output:
% - x           (N-2)D array if k is a scalar, (N-1)D array if it is a
%               vector

% check
s = size(data);
if V.nx~=s(1) || V.ny~=s(2), error('dimenstion mismatch'), end
selset = V.a4d.SI.selection.singleset;
if nargin<3 || isempty(ks), ks = 1:length(selset); end
if max(ks)>length(selset), error('requested selection does not exist'), end

% about size
if nargin<4, dim = ndims(data)-1; end
newsize = s(3:end);
newsize(dim:end+1) = [1 newsize(dim:end)];

% compute
data = reshape(data,[s(1)*s(2) s(3:end)]);
x = cell(1,length(ks));
for k=ks
    ind = selset(k).dataind;
    x{k} = reshape(mean(data(ind,:),1),newsize);
end
x = cat(dim,x{:});