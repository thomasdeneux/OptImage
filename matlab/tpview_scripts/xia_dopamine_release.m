function label = xia_dopamine_release(V) 
% function f(V) 
% function label = f('label') 

if nargin==0, V=evalin('base','V'); end
if ischar(V) && strcmp(V,'label')
    label = '(Xia) dopamine release';
    return
end

% Write script code below

%% get data

trials = V.content.trials;
data = {trials.dataop};
ntrial = length(data);

trial_names = cell(1, ntrial);
for i = 1:ntrial
    f = brick.fileparts(trials(i).file, 'name');
    trial_names{i} = strrep(f,'.tiff.bin2.mat','');
end

%% further processing

% remove 1 minute at begin and end because of innacurate temporal high-pass
% in these zones
dt = trials(1).dt;  % should be 0.5 s
nrm = round(60 / dt);

for i = 1:ntrial
    data{i} = data{i}(:,:,1+nrm:end-nrm);
end

%% thresholding

y = zeros(1, ntrial);

thr = 0.01;
for i = 1:ntrial
    x = data{i};
    
    % Divide by median
    m = median(x,3);
    x = x./m;

    % Keep only above threshoold
    x = max(x - thr - 1, 0);

    % Average!
    y(i) = mean(x(:));
end

%% display
brick.figure('Dopamine release')

bar(y)

set(gca,'xtick',1:ntrial,'xticklabel',trial_names,'XTickLabelRotation',30)
ylabel 'average above baseline (DF/F)'
title 'Dopamine release'


