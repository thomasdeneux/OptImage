function varargout = tps_motioncorrection(varargin)
% function data = tps_motioncorrection(data[,template][,par])
% function par = tps_motioncorrection('par')
%---

%-------------------------------------------------------------------------%
% ATTENTION!                                                              %
% - i denotes line number, and j column, but the acquisition is           %
%   line per line, which means the pixel in line i, column j is x(j,i)!!! %
% - any point coordinates will be stored, first line then column, i.e.    %
%   p = [i; j]; in other words, image value at point p is x(p(2),p(1))    %
% as a summary, we are more or less using Matlab convention for images    %
% here, which means for example gradient is given as [Ti Tj] = gradient(T)%
%-------------------------------------------------------------------------%

if ischar(varargin{1}) && strcmp(varargin{1},'par')
    par = defaultparameters;
    varargout = {par};
else
    data = varargin{1};
    par = defaultparameters;
    template = mean(data,3);
    for i=2:length(varargin)
        a = varargin{i};
        if isstruct(a)
            par = fn_structmerge(par,a,'strict');
        else
            template = a;
        end
    end
    data = register(data,template,par);
    varargout = {data};
end

%---
function par = defaultparameters

par.interpmethod = 'cubic';   % same as METHOD for interp1 function
par.nlineperpt = 10;     % one point every ~5ms
par.maxmotion = 5;      % pixels
par.scandur = .8;       % perc. of time scanning

%---
function [ii jj] = interpmotion(p,prm)
% interpolate the motion for each pixel (ni*nj values) based on the motion
% of reference points (ipt values)

% p is a npt*2 matrix
% ii and jj are nx vectors
ii = prm.ii(:)+prm.dpixdp*p(:,1);
jj = prm.jj(:)+prm.dpixdp*p(:,2);

%---
function [x dxdp] = resampletemplate(p,prm)

% pixel positions modified by the motion
[ii jj] = interpmotion(p,prm);

% interpolate template at new position (note the extrapolation value 0)
x = interp2(1:prm.ni,1:prm.nj,prm.T,ii,jj,'linear',0);
x = reshape(x,prm.nj,prm.ni);

% derivatives
if nargout>1
    dxdi = interp2(1:prm.ni,1:prm.nj,prm.Ti,ii,jj,prm.interpmethod,0);
    dxdj = interp2(1:prm.ni,1:prm.nj,prm.Tj,ii,jj,prm.interpmethod,0);
    % - dpixdp is a nx*npt matrices (for each pixel, how much its motion
    %   depends from the refence points)
    % - dxdi and dxdj ar ni*nj matrices, i.e. nx vectors -> repmat to
    %   nx*npt matrices 
    % - dxdp is a nx*(npt*2) matrix (p contains first the npt reference
    %   motions along the i-axis, then the npt motions along the j-axis)
    dxdi = repmat(dxdi(:),1,prm.npt);
    dxdj = repmat(dxdj(:),1,prm.npt);
    dxdp = [dxdi.*prm.dpixdp dxdj.*prm.dpixdp];
end

%---
function x = resample(x0,p,prm)

ii = prm.ii(:)-prm.dpixdp*p(:,1);
jj = prm.jj(:)-prm.dpixdp*p(:,2);
x = interp2(1:prm.ni,1:prm.nj,x0,ii,jj,'linear',0);
x = reshape(x,prm.nj,prm.ni);


%---
function [E dEdp d2Edp2] = fittotemplate(x0,p,prm)

[x dxdp] = resampletemplate(p,prm);
diff = x - x0;
badidx = (x==0); % pixels with no value because out of template
diff(badidx)=0;
dxdp(badidx,:)=0;
E = sum(diff(:).^2);
if nargout>1
    % dxdp is a nx*(npt*2) matrix
    % diff is a ni*nj matrix (i.e. a nx vector)
    % dEdp is a (npt*2) vector
    dEdp = 2*dxdp'*diff(:);
    if nargout>2
        % "Gauss-Newton approximation of the Hessian matrix"
        % i.e. ignore terms due to dx2dp2
        d2Edp2 = 2*dxdp'*dxdp;
    end
else 
    disp('fittotemplate called with only one output argument')
end

%---
function [x p] = registerframe(x0,prm,opt)

% motion estimation
p0 = zeros(prm.npt,2);
p = fmincon(@(p)fittotemplate(x0,p,prm),p0,[],[],[],[], ...
    -prm.maxmotion,prm.maxmotion,[],opt);

% resampling
x = resample(x0,p,prm);

%---
function data = register(data,T,par)

% parameters
% (interpolation method, maximum motion)
prm = par;
% (sizes)
[nj ni nfr] = size(data);
prm.ni = ni; prm.nj = nj;
prm.npt = ceil(ni/par.nlineperpt)+1;
prm.ipt = linspace(0,ni,prm.npt);
% (non-integer 'line numbers' for each pixel in image)
ioneline = (.5:nj)'/nj*par.scandur;
prm.ipix = fn_add(ioneline,0:ni-1);
% (precompute indices for function 'interpmotion')
[prm.ii prm.jj] = meshgrid(1:prm.ni,1:prm.nj);
prm.dpixdp = interp1(prm.ipt,eye(prm.npt),prm.ipix(:),prm.interpmethod);
% (template and gradient)
prm.T = T;
[prm.Ti prm.Tj] = gradient(T);

% optimization parameters
opt = optimset('GradObj','on','Hessian','on','Display','none');

% loop on frames
fn_progress('correcting frame',nfr)
for k=1:nfr
    fn_progress(k)
    data(:,:,k) = registerframe(data(:,:,k),prm,opt);
end







        

