function label = tpv_count(V,s) 
% function {par stats} = tpv_count(V[,s]) 
% function label = tpv_count('label') 
%---
% prompt for parameters and call tps_count
% except if a second argument is provided to specify some parameters
%
% Input:
% - V       tpview object
% - s       structure with fields 'LS', 'HS', 'timethr' and 'save'
%
% See also tps_count

persistent par0 ntrial0 file0

if nargin==0, V=evalin('base','V'); end

if ischar(V) && strcmp(V,'label')
    label = 'count spikes/bumps';
    return
end

if isempty(par0) || nargin==2
    par = tps_count('par');
    par.display = true;
    par.save = false;
    par.trials = ['1:' num2str(V.ntrial)];
    ntrial0 = V.ntrial;
    calcperiod = true;
else
    par = par0;
    if ~strcmp(V.content.trials(1).file,file0)
        % probably a new experiment
        par.trials = ['1:' num2str(V.ntrial)];
        calcperiod = true;
        file0 = V.content.trials(1).file;
        ntrial0 = V.ntrial;
    elseif ntrial0<V.ntrial
        % probably new trials have been added
        [start tokens] = regexp(par.trials,'(\d+):(\d+)$','start','tokens');
        if isempty(tokens)
            par.trials = ['1:' num2str(V.ntrial)];
        else
            tokens = tokens{1};
            if str2double(tokens{2})==ntrial0
                par.trials = [par.trials(1:start-1) tokens{1} ':' num2str(V.ntrial)];
            else
                par.trials = [par.trials ' ' num2str(ntrial0+1) ':' num2str(V.ntrial)];
            end
        end
        calcperiod = false;
        ntrial0 = V.ntrial;
    else
        calcperiod = false;
    end
end
par.LS = V.disppar.LS;
par.HS = V.disppar.HS;
par.thr = V.disppar.timethr; if isempty(par.thr), par.thr=0; end
if calcperiod     
    if par.thr==0, decaytime = 0; else decaytime = 1.15; end
    if isfield(V.trial.addinfo,'stim')
        stim = [V.trial.addinfo.stim(1,1) sum(V.trial.addinfo.stim(:,end))];
    else
        % assume a 2s stimulation after a 10s delay
        stim = [10 12];
    end
    nbef = stim(1)-1;
    naft = V.trial.nfr*V.dt - stim(2) - decaytime;
    nsec = floor(min(nbef,naft));
    par.periods = [stim(1)-nsec      nsec/2; ...
        stim(1)-nsec/2               nsec/2; ...
        stim(2)+decaytime            nsec/2; ...
        stim(2)+decaytime+nsec/2     nsec/2];
end   

if nargin<2
    par = fn_structedit(par);
    if isempty(par), return, end
    par0 = par;
    if par.save
        if par.spike
            par.save = [fn_cd('2ps') '/save/tpv_count/'];
        else
            donpil = V.content.signal.opdef.norm_npil_norm(2);
            strmode = fn_switch(donpil,'npilsub','raw');
            fname = V.content.trials(1).file;
            tokens = regexp(fname,'(\d{6}).*/(.*)_','tokens');
            tokens = tokens{1};
            par.save = [fn_cd('2ps') '/save/tpv_count/' ...
                tokens{1} '_' tokens{2} '_' strmode '_'];
        end
    end
    assignin('base','par',par)
else
    par.LS = s.LS;
    par.HS = s.HS;
    par.thr = s.timethr;
    par.save = s.save;
end


stats = tps_count(V.content,par);
if nargout==1, label = struct('par',par,'stats',stats); end


% if all(par.save) && ~par.spike
%     % subtle rename...
%     nameold = [fn_cd('2ps') '/save/tpv_count/' strmode '_result.png'];
%     namenew = [fn_cd('2ps') '/save/tpv_count/result_' strmode '.png'];
%     movefile(nameold,namenew)
% end

