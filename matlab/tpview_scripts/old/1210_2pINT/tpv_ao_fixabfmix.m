function label = tpv_ao_fixabfmix(V) 
% function f(V) 
% function label = f('label') 

if ischar(V) && strcmp(V,'label')
    label = 'Fix recording permutations';
    return
end

perms = struct('trial',{},'permnames',{});
ntrialrec = 0;
for ktrial=1:V.ntrial
   
    Tk = V.content.trials(ktrial);
    r = Tk.recording;
    if isempty(r), continue, end
    ntrialrec = ntrialrec + 1;
    nrec = length(r);
    
    kheart = find(strcmp({r.name},'heart'));
    if isempty(kheart)
        errordlg 'no heart recording in trial, cannot proceed'
        return
    end
    
    % get recording signals
    y = [r.signal];
    % fft
    y1 = abs(fn_normalize(y,1,'-'));
    yf = abs(fft(y1));
    df = 1/(r(1).dt*r(1).n);
    % restrict to >0 and <=10Hz
    freqs = df*(2:ceil(10/df));
    yf = yf(2:ceil(10/df),:);
    % for which signal do we have a strong peak around 7Hz?
    maxfactor = zeros(1,nrec);
    [M ipeak] = max(yf);
    for i = 1:nrec
        if ipeak(i)>5/df
            % peak is between 5 and 10Hz
            fpeak = ipeak(i)*df;
            maxfactor(i) = M(i) / max(yf(freqs<fpeak-.5|freqs>fpeak+.5,i));
        end
    end
    if ~any(maxfactor)
        if fn_dodebug
            nex = ceil(3/r(1).dt);
            figure(1), subplot(211), plot((0:nex-1)*r(1).dt,y(1:nex,:)); subplot(212), plot(freqs,yf);
            disp 'please define kheartreal'
            keyboard
        else
            errordlg(['trial ' num2str(ktrial) ': could not locate a peak between 5 and 10Hz to pick up heart recording, cannot proceed'])
            continue
        end
    else
        [dum kheartreal] = max(maxfactor);
    end
    
    % define permutation
    if kheart~=kheartreal
        nex = ceil(3/r(1).dt);
        figure(1), subplot(211), hl1=plot((0:nex-1)*r(1).dt,y(1:nex,:)); subplot(212), hl2=plot(freqs,yf);
        set([hl1(kheartreal) hl2(kheartreal)],'linewidth',2)
        pause
        perm = fn_mod((1:3)+kheart-kheartreal,3);
        perms(end+1) = struct('trial',ktrial,'permnames',{{r(perm).name}}); %#ok<AGROW>
    end
    
end

% Apply permutation
nperm = length(perms);
if nperm==0
    msgbox('Found no permutation in recordings')
else
    answer = questdlg(['Found permutations in ' num2str(nperm) ' of ' num2str(ntrialrec) ' recordings. Apply them?'], ...
        'question','Yes','No','Yes');
    if strcmp(answer,'Yes')
        for i=1:nperm
            setrecordingname(V.content.trials(perms(i).trial),[],perms(i).permnames)
        end
        V.content.setelectrophys(tps_electrophys(V.content.trials))
    end
end
