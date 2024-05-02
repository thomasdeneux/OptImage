function label = tpv_spectrogram(V) 
% function f(V) 
% function label = f('label') 

if nargin==0, V=evalin('base','V'); end
if ischar(V) && strcmp(V,'label')
    label = 'Spectrogram';
    return
end

% Get signals
x = V.getsignals;
deltat = diff(x(1).tidx(1:2));
if strcmp(V.disppar.timesignal,'signal')
    x = [x.data];
else
    x = [x.dataop];
end

% Options
s = fn_structedit('select__signals',{1:size(x,2) 'double'},'frequency__scaling',{'f' {'none' 'f'}});

% Sub-selection of signals
x = x(:,s.select__signals);

% Compute spectrogram
[y freqs yticklog yticklabel] = fn_spectrogram(x,deltat);
switch s.frequency__scaling
    case 'none'
    case 'f'
        y = fn_mult(y,freqs);
    otherwise
        error programming
end
y = permute(y,[3 4 1 5 6 7 2]);

% Scales
G = V.a4d.G;
G.labels{7} = 'frequency';
G.units{7} = '10^x Hz';
deltaf = diff(log10(freqs(1:2)));
G.mat(8,[1 8]) = [log10(freqs(1))-deltaf deltaf];

% Display
figure
a = fourd(y,G,'2d','proj',[3 7], ...
    'cmap','mapgeog','scaledisplay','tick');
set(a.D.ha,'ydir','normal','ytick',yticklog,'yticklabel',yticklabel)
ylabel(a.D.ha,'frequency (Hz)')

% User 
set(a.D,'clipmode','link1')
