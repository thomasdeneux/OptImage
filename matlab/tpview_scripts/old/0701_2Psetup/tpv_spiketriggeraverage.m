function label = tpv_spiketriggeraverage(V)
% function tpv_spiketriggeraverage(V)
% function label = tpv_spiketriggeraverage('label')
%---
% do spike-trigger aceraging and display

if nargin==0, V=evalin('base','V'); end

if ischar(V) && strcmp(V,'label')
    label = 'spike-trigger averaging';
    return
end

persistent par
trig = gettriggers(V.content);
F = fieldnames(trig);
F1 = [{'all spikes'}; F];
spec.trig = F1;
spec.nframebefore = 'slider 1 50 1';
spec.nframeafter  = 'slider 1 100 1';
spec.descan = 'xslider 0 6 1';
spec.ncol = 'slider 1 10 1';
spec.savedir = 'dir';
if isempty(par)
    par.trig = F{1};
    par.nframebefore = 3;
    par.nframeafter  = 15;
    par.oddeven = false;
    par.descan = [];
    par.ncol = 6;
    par.clip = [.985 1.05];
    par.colormap = 'jet';
    par.applymask = false;
    par.saveresults = false;
    par.savedir = fn_getfile('REP');
    par.prefix = '';
    par.suffix = '';
end
par = fn_structedit(par,spec);
if isempty(par), return, end
assignin('base','par',par)

a = cat(3,V.content.trials.dataop);
if strcmp(par.trig,'all spikes')
    okspike = false(1,length(F));
    for k=1:length(F), okspike(k) = ~isempty(strfind(F{k},'spike')); end
    F1 = F(okspike);
else 
    F1 = {par.trig};    
end
ntrig = length(F1);

% save the results in a tps_trial object and in files
if par.saveresults, T = tps_trial; kfile = 0; swd=pwd; cd(par.savedir), end
for ktrig = 1:ntrig
    f = F1{ktrig};
    idx = trig.(f);    
    
    % % tmp
    % x = fn_triggeravg(a,idx,[par.nframebefore par.nframeafter]);
    % if ~isempty(par.descan)
    %     x = fn_interprows(x,par.descan);
    % end
    % y = squeeze(x(:,1,:))'-1;
    % figure(1), plot(y)
    % assignin('base','x',x)
    % assignin('base','y',y)
    %
    % return
    
    if par.oddeven
        idx = {union(idx(2:4:end),idx(3:4:end)) ...
            union(idx(1:4:end),idx(4:4:end))};
    else
        idx = {idx};
    end
    if isempty(idx{1}), continue, end
    
    % mask the neurons?
    if par.applymask
        mask = logical(getselectionmap(V.content));
    end
    for kidx = 1:length(idx)
        b = fn_triggeravg(a,idx{kidx},[par.nframebefore par.nframeafter]);
        if ~isempty(par.descan)
            b = fn_interprows(b,par.descan);
        end
        figure(kidx), clf
        if par.applymask
            % clip
            if isempty(par.clip)
                m = min(b(:))-1e-6; M = max(b(:))+1e-6;
            else
                m = par.clip(1); M = par.clip(2);
                b = max(m+1e-6,min(M-1e-6,b));
            end
            clip = [m-(M-m) M];
            cmap = [gray(256); jet(256)];
            % mask: neurons in gray levels, neuropil in colors
            [nx ny nt] = size(b);
            b = reshape(b,nx*ny,nt);
            b(mask,:) = b(mask,:)-(M-m);
            b = reshape(b,nx,ny,nt);
        else
            cmap = par.colormap;
            clip = par.clip;
        end
        a4d = fn4D('in',kidx,b,'frames','ncol',par.ncol,'cmap',cmap, ...
            'clipmode','data');
        D = a4d.D;
        if ~isempty(clip)
            D.clip = clip;
        end
        set(D.lines,'color','w')
        %set(D.cross,'color','g','marker','p','markersize',10)
        if par.saveresults
            kfile = kfile+1;
            fname = [par.prefix f '_' num2str(kidx) par.suffix];
            % save in tps_trial object
            T(kfile) = tps_trial(b,V.trial);
            savedata(T(kfile),fname)
            % save in file
            fn_savefig(kidx,[fname '.png'])
        end
    end
    
    if par.saveresults
        fname = [par.prefix 'triggers' par.suffix '.tptrial'];
        save(fname,'T','-MAT')
    end
end

% save the tps_trial object
if par.saveresults, cd(swd), end