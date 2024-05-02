function out = tps_expfit(data,u)
% function s = tps_expfit(data,n)
% function corrected = tps_expfit(data,s)
%---
% fit data with one constant + n decaying exponentials and remove the
% exponentials
%
% the logic of the function is that the time constants of the exponentials,
% or even the full decay pattern, can be first estimated on an average
% data (first function form; this solves a minimization problem), and then
% used for correction of noisier individual data (second function form;
% this solves only a general linear model)
%
% Input:
% - data    a 2D (1st dimension is time) or 3D or 4D array (3rd dimension
%           is time)
% - n       the number of exponentials for the model
% 
% Input/Output:
% - s       a structure that contains model parameters and optionally also
%           contains fit results
%           when s is an output its fields are:
%           .tau    vector of length n - the time constants of the
%                   exponentials
%           .beta   vector of (1+n) parameters - the weights for each
%                   regressor of the linear model
%           .is3d   logical - true if time is the 3rd dimension, false if
%                   time is the 1st dimension
%           when s is an input, 'tau' and 'is3d' are mandatory fields; if
%           'beta' is present, no general linear model fit is performed on
%           the data, but this set of parameters is used
%
% Output:
% - corrected   the decaying exponentials (but not the constant) have been
%           removed from the data

if nargin<2, u = 2; end
data = double(data);
if isstruct(u) || isa(u,'fn_pointer')
    out = correct(data,u);
else
    out = estimate(data,u);
end

%---
function s = estimate(data,n)

% size -> make a nt x np array
switch ndims(data)
    case 2
        is3d = false;
        if size(data,1)==1, data = data'; end
        [np nt] = size(data);  %#ok<ASGLU>
    case 3
        is3d = true;
        s = size(data);
        nt = s(3);
        np = s(1)*s(2);
        data = reshape(data,np,nt)';
    otherwise
        error 'exponential fitting cannot be performed on 4D data'
end

% initialization of exponents
tau0 = logspace(0,log10(2*nt),n+2);
tau0 = tau0(2:end-1);

% optimization of exponents
FACT = 1e4;
opt = optimset('display','iter','tolfun',1e-3,'algorithm','active-set');
tau = FACT*fmincon(@(x)energy(x*FACT,double(data)),tau0/FACT,[],[],[],[],ones(1,n)/FACT,ones(1,n)*5*nt/FACT,[],opt);

% final fit
[datafit b] = linearfit(data,tau); 
b = fn_div(b,mean(data,1)); % normalize the parameters by the average frame

% display
figure;
plot([mean(data,2) mean(datafit,2)])

% output the parameters for future artefact corrections
s = struct( ...
    'tau',      tau, ...
    'beta',     b, ...
    'is3d',     is3d ...
    );


%---
function [datafit b] = linearfit(data,tau)

[nt np] = size(data); %#ok<NASGU>
tt = (1:nt)';
A = [ones(nt,1) exp(-fn_div(tt,tau))];
b = pinv(A)*data;
b = max(b,0);
datafit = A*b;

%---
function e = energy(tau,data)

datafit = linearfit(data,tau);
d = data-datafit;
e = norm(d(:));
%disp(['tau: ' num2str(tau,'%.4f ') '-> e = ' num2str(e,'%.6f')])

%---
function corrected = correct(data,s)
% note that this removes the bleaching but leaves the average frame
% unchanged

if ~s.is3d, error 'correction not implemented yet for data organized as time x space', end
[nx ny nt nc] = size(data); 

% build the denominator
if isfield(s,'correction')
    correction = s.correction;
else
    data1 = data;
    if s.is3d, data1 = fn_reshapepermute(data1,{3 [1 2 4]}); end
    tt = (1:nt)';
    A = [ones(nt,1) exp(-fn_div(tt,s.tau))];
    if isfield(s,'beta')
        beta = s.beta;
        beta = fn_mult(beta,mean(data1,1)); % rescale by average frame
        if nc>1, beta = repmat(beta,[1 nc]); end
    else
        beta = pinv(A)*data1;
    end
    correction = A(:,2:end) * beta(2:end,:);
    if s.is3d, correction = fn_reshapepermute(correction,[nt nx*ny nc],[2 1 3],[nx ny nt nc]); end
end

% that's it!
corrected = fn_subtract(data,correction); 

