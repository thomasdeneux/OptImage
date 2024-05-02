function data = tps_heartbeat(data,varargin)
% function data = tps_heartbeat(data,heartcycle[,par])
% function data = tps_heartbeat(data,dt,freqrange[,par])
% function data = tps_heartbeat(data,heartsignal,dt,freqrange[,par])
% function par = tps_heartbeat('par')
%---
% remove heart beat from signals, or other oscillations that are periodic
% but whose frequency is not completely stable
%
% the phase of the heart can be provided by the user, or inferred from
% heart recordings, or from the data itself; once the phase is known, the
% estimation of heart artefact is reduced to a general linear model
%
% additionally to heart signals, slow drifts can be removed from the data
% 
% Input:
% - data        an nx*ny*nt movie, or an nt*ny array
% - heartcycle  a vector of length nt - the phase of the heart, it must be
%               a signal that increases continuously from 0 to 1 during one
%               heart beat, then jumps back to 0 and repeats like this
% - heartsignal a vector of length nt - heart recording; if neither
%               heartcycle, nor heartsignal is provided, the data itself,
%               averaged over all space, will be used as a heart signal
% - dt          sampling time of both the data and the heart signal
% - freqrange   a 2-elements vector - minimal and maximal frequency of the
%               heart
% - par         a parameters structure with fields
%               .K  number of harmonics for heart regressors (default 5)
%               .N  number of harmonics for additional slow drifts (default 0)
%               .prec precision (the higher, the more flexible the
%                   heartcycle estimation, default .6) 
%               use par = tps_heartbeat('par') to get the default
%               parameters structure
%
% Output:
% - data        data with heart artefact (and possibly slow drifts) removed

% persistent mem_X mem_Y

% default parameters
if nargin==1
    par = defaultpar;
    data = par;
    return
end

% Input
% (data)
if ndims(data)==3
    [nx ny nt] = size(data);
    x = permute(data,[1 3 2]);
else
    [nt nx] = size(data);
    %     ny = nt2/nt;
    %     if mod(ny,1), error programming, end
    x = data;
    if isvector(x), x = x(:); end
end
% (parameters)
if ~isempty(varargin) && isstruct(varargin{end})
    par = fn_structmerge(defaultpar,varargin{end});
    varargin(end) = [];
else
    par = defaultpar;
end
% (heart cycle)
switch length(varargin)
    case 1
        heartcycle = varargin{1};
    case 2
        % compute heart cycle based on data
        [dt freqrange] = deal(varargin{:});
        heartcycle = tps_heartcycle(mean(x,2),dt,freqrange,par.prec);
    case 3
        % compute heart cycle based on provided heart signal
        [heartsignal dt freqrange] = deal(varargin{:});
        heartcycle = tps_heartcycle(heartsignal,dt,freqrange,par.prec);
        if length(heartcycle)~=nt
            error 'not implemented any more, please help me'
        end
    otherwise
        error argument
end

% Regressors based on heart cycle
K = par.K;
A = zeros(nt,2*K); 
for k=1:K
    phase = 2*pi*k*heartcycle;
    A(:,k) = sin(phase); 
    A(:,K+k) = cos(phase); 
end

% Slow regressors
N = par.N;
B = zeros(nt,2*N+1);
B(:,1) = 1; 
for i=1:N
    phase = pi*i*(1:nt)/nt;
    B(:,1+i) = sin(phase);
    B(:,1+N+i) = cos(phase);
end

% All regressors: slow modulation of the heart-based regressors
X = fn_mult(A,permute(B,[1 3 2]));
X = reshape(X,nt,2*K*(2*N+1));
X = fn_normalize(X,1,'-'); % subtract mean

% Regression (last result is kept in memory in order to avoid too many
% useless calculations)

% if isequal(X,mem_X)
%     Y = mem_Y;
% else
%     Y = X*pinv(X);
%     mem_X = X;
%     mem_Y = Y;
% end
% x = x - Y*x;

x = x - X*(X\x);

if ndims(data)==3
    data = reshape(x',nx,ny,nt);
else
    data = x;
end

%---
function par = defaultpar

par.K = 5;
par.N = 0;
par.prec = .6;
