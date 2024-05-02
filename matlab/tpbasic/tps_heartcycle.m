function heartcycle = tps_heartcycle(heartbeat,dt,freqrange,prec)
% function heartcycle = tps_heartcycle(heartbeat,dt[,freqrange|maxfreq[,prec]])
%---
% infers the phase of heart (or other oscillations that are periodic but
% whose frequency is not completely stable) in a signal
%
% Input:
% - heartbeat   vector - signal
% - dt          time interval between successive points
% - freqrange   2-element vector - frequency interval where to locate the
%               peak frequency corresponding to heart in the signal
% - prec        scalar > 0 - the heart signal will be band-passed filtered
%               at cutoff frequencies 'peak frequency'*[2^-prec 2^prec]
%               (default value: 0.6) before getting the phase
%
% Ouput:
% - heartcycle  a vector of same length as the signal, that increases
%               continuously from 0 to 1 during one heart beat, then jumps
%               back to 0 and repeats like this

if nargin==0, help tps_heartcycle, return, end

% Input
if ~isvector(heartbeat), error argument, end
if nargin<3
    freqrange = [0 9]; % default for rat!
elseif isscalar(freqrange)
    freqrange = [0 freqrange];
end
if nargin<4, prec = .6; end
relprec = 2^prec;


heartbeat = heartbeat(:);
nt = length(heartbeat);
fs = 1/dt;

% heart beat cycle
fbeat = fft(heartbeat);
freqstep = fs/nt;
minidx = max(2,1+floor(freqrange(1)/freqstep));
maxidx = min(1+ceil((nt-1)/2),1+ceil(freqrange(2)/freqstep));
[m idx] = max(abs(fbeat(minidx:maxidx))); %#ok<ASGLU> % consider only frequencies < max to avoid harmonics
idx = minidx+idx-1;
rate = (idx-1)*freqstep;
disp(['heart rate ' num2str(rate,'%.1f') ' Hz'])
fbeat([1:round(idx/relprec) round(idx*relprec):end]) = 0; % restrict frequencies in a window around the main heart beat frequency
heartcycle = angle(ifft(fbeat));
heartcycle = mod(heartcycle,2*pi)/(2*pi);


