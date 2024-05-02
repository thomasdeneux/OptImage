function label = tpv_pixelcorreletions(V)
% function tpv_pixelcorreletions(V)
% function label = tpv_pixelcorreletions('label')
%---
% correlations between pixels

if nargin==0, V=evalin('base','V'); end

if ischar(V) && strcmp(V,'label')
    label = 'correlations btw. pixels';
    return
end

idx = getspontaneousindices(V.content);
a = V.dataop(:,:,idx); 
b = reshape(a,V.nx*V.ny,length(idx));
X = corrcoef(b'); 
Y = reshape(X,V.nx,V.ny,V.nx,V.ny); 
figure(1), clf, imagesc(mean(V.data,3)'), colormap gray, fn_imvalue image
assignin('base','Y',Y)
fn_imvalue('register',1,'imagesc(Y(:,1:end-1,j,i)'',[0 .8])')
