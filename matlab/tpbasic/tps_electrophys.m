classdef tps_electrophys < hgsetget
   
    properties (Access='private')
        version = 2.0;
        linedur % old property
        stim    % old property
    end
    properties
        t0 = 0;
        dt = 0;
        n = 0;
        file
        recidx
    end
    properties (Transient) 
        raw
        signal
    end
    properties
        corrpar = struct('name',{},'value',{});
        spikepar
        spikes
    end
    properties (Dependent)
        tidx
    end
    
    % Constructor
    methods
        function E = tps_electrophys(data,fs,stim,linedur)
            % function E = tps_electrophys
            % function E = tps_electrophys(T[,channelname])
            % function E = tps_electrophys(fname|data,fs[,stim[,linedur]])
            %---
            % get trial recording, correct stimulation artefact, guess spikes
            
            % special cases
            if nargin==0, return, end % singleton object with no signal
            if isa(data,'tps_trial')
                T = data;
                for k=1:length(T)
                    E(k) = tps_electrophys; %#ok<AGROW>
                end
                if nargin==2
                    channelname = fs;
                    readfromtrial(E,T,channelname)
                else
                    readfromtrial(E,T)
                end
                return
            end
            
            % general case
            if ischar(data)
                data = cellstr(data); 
            elseif isnumeric(data)
                if isvector(data) data = {data}; else data = num2cell(data,1); end
            end
            ntrial = length(data);
            channel = [];
            for k=1:ntrial
                E(k) = tps_electrophys; %#ok<AGROW>
                E(k).dt = 1/fs; %#ok<AGROW>
                x = data{k};
                if ischar(x)
                    E(k).file = x; %#ok<AGROW>
                    x = fn_readdatlabview(x);
                end
                E(k).n = length(x); %#ok<AGROW>
                if ~isvector(x) && ~isempty(x)
                    if isempty(channel)
                        nchannel = size(x,2);
                        hf = figure;
                        tt = (0:E(k).n-1)*E(k).dt;
                        for i=1:nchannel
                            subplot(nchannel,1,i)
                            plot(tt,x(:,i))
                        end
                        channel = str2double(fn_input('channel__number',0,fn_num2str(1:nchannel,'cell')));
                        close(hf), drawnow
                    end
                    x = x(:,channel);
                end
                E(k).signal = single(x(:)); %#ok<AGROW>
            end
            %okind = ([E.n]>0); WHAT'S THIS!!??
            %E = E(okind);
            if nargin>=3
                if iscell(stim)
                    [E.stim] = deal(stim{:});
                else
                    [E.stim] = deal(stim);
                end
            end
            if nargin>=4
                [E.linedur] = deal(linedur);
            else
                [E.linedur] = deal(0); 
            end
        end
        function readfromtrial(E,T,channelname)
            if nargin<3, channelname = 'electrophysiology'; end
            % read
            for k=1:length(T)
                r = T(k).recording;
                E(k).recidx = find(strcmp({r.name},channelname));
                if isempty(E(k).recidx), continue, end
                if ~isscalar(E(k).recidx), error('cannot handle multiple electrophysiology recording!'), end
                r = r(E(k).recidx);
                E(k).t0 = r.t0;
                E(k).dt = r.dt; 
                E(k).n = r.n; 
                E(k).file = T(k).analogfile;
                E(k).linedur = T(k).linedur;
                if isfield(T(k).addinfo,'stim')
                    E(k).stim = T(k).addinfo.stim;
                end
                E(k).signal = r.signal;
            end
        end
    end
    
    % Get/Set
    methods
        function x = get.raw(E)
            if isequal(E.raw,[]) && ~isempty(E.file)
                [sig E.file] = tps_readrecording(E.file,E.recidx);
                E.raw = sig{1};
            end
            x = E.raw;
        end
        function x = get.signal(E)
            if isequal(E.signal,[])
                removeartefact(E)
            end
            x = E.signal;
        end
        function set.corrpar(E,par)
            E.corrpar = par;
            E.signal = []; %#ok<MCSUP>
        end
        function tidx = get.tidx(E)
            tidx = E.t0 + (0:E.n-1)*E.dt;
        end
    end
    
    % Computations
    methods
        function removeartefact(E)
            % checks
            indtrial = find([E.n]>1);
            if isempty(indtrial), return, end
            DT = unique([E(indtrial).dt]);
            lindur = unique([E(indtrial).linedur]);
            if ~isscalar(DT) || ~isscalar(lindur)
                disp('cannot remove electrophy artefact with different sampling rates or line durations for different trials')
                return
            end
            p = E(indtrial(1)).corrpar;
            for k=indtrial(2:end)
                if ~isequal(E(k).corrpar,p)
                    disp 'cannot remove electrophy artefact with different correction parameters for different trials'
                    return
                end
            end
            
            % init
            [E.signal] = deal(E.raw);                
            
            % apply corrections
            for kcorr = 1:length(E.corrpar)
                name = E(1).corrpar(kcorr).name;
                value = E(1).corrpar(kcorr).value;
                switch name
                    case 'user'
                        for k=indtrial
                            E(k).signal = feval(value,E(k).raw);
                        end
                    case 'xscanning'
                        % X-SCANNING ARTEFACT
                        if lindur
                            lindur = 2*lindur; % group lines 2 by 2 because of double-way scanning
                            nxscan = lindur/DT;
                            if mod(nxscan,1)
                                disp('cannot handle line duration which is not a multiple of recording sampling rate')
                                return
                            end
                            for k=indtrial
                                x = E(k).signal;
                                nblocks = floor(E(k).n/nxscan);
                                x = x(1:nxscan*nblocks);
                                x = reshape(x,nxscan,nblocks); % nxscan x nblocks array, time difference between successive columns is lindur
                                x = fn_filt(x,.5/lindur,'lm',2); % filter at 1/2 seconds rather than just average over blocks (this adjusts for small time shifts between scanning and recording)
                                artefact = x(:);
                                artefact(end+1:E(k).n) = x(end,1:E(k).n-length(artefact));
                                E(k).signal = single(E(k).signal - artefact);
                            end
                        end
                    case 'stim'
                        % STIMULATION ARTEFACT
                        for k=indtrial
                            stim = E(k).stim;
                            nstim = size(stim,2);
                            
                            % checks
                            if isempty(stim), continue, end
                            if nstim<5, error('not enough repetitions of the stimulation to get a good artefact estimation'), end
                            if any(diff(stim(2,:))), error('cannot estimate stimulation artefact with multiple stimulation patterns'), end
                            
                            % window to average over (less than distance btw.
                            % repetitions and no greater than 200ms)
                            window = min(.200,min(diff(E(k).stim(1,:))));
                            
                            % cut blocks
                            nt = floor(window/DT);
                            nbef = floor(.020/DT); % # time instants before stimulation
                            idx = fn_add((-nbef:nt-nbef-1)', 1+round(E(k).stim(1,:)/DT));
                            block = E(k).signal(idx); % nt x nstim array
                            artefact = mean(block,2);
                            artefact = artefact - mean(artefact([1:nbef end-nbef+1:end])); % artefact should be zero just before stim. and away from it
                            
                            % correct signal
                            E(k).signal(idx) = E(k).signal(idx) - repmat(artefact,1,nstim);
                        end
                end
            end
        end
        
        function guessspikes(E,ktrialstart)
            if nargin<2, ktrialstart = 1; end
            % graphics
            hf = figure(250); clf
            set(hf,'name','tps_electrophys.guessspikes','numbertitle','off')
            set(hf,'windowButtonMotionFcn',' ')
            g.hf = hf;
            ha = [axes axes axes];
            g.ha = ha;
            g.hb = axes;
            ntrial = length(E);
            g.u(1) = uicontrol('string','prev spike (<)', ...
                'position',[10 70 120 15],'callback',@(u,e)nextspike(-1));
            g.u(2) = uicontrol('string','next spike (>)', ...
                'position',[10 50 120 15],'callback',@(u,e)nextspike(1));
            g.u(3) = uicontrol('string','remove spike (r)', ...
                'position',[10 30 120 15],'callback',@rmspike);
            g.u(4) = uicontrol('string','next trial ( )', ...
                'position',[10 10 120 15],'callback',@(u,e)nexttrial(1));
            g.list = uicontrol('style','listbox','string',fn_num2str(1:length(E),'cell'), ...
                'callback',@(u,e)settrial(get(u,'val')));
            g.par = uipanel;
            g.slider = uipanel;
            fn_framedesign(g,[fn_fileparts(which('tps_electrophys'),'noext') '_guessspikes'])
            slider = fn_slider(g.slider,'mode','point','width',.1,'callback',@(u,e)movetime(u.value));
            set(hf,'windowKeyPressFcn',@keypress)
            
            s = struct('flag',{'+' {'button' '+' '?' '-' '0'}}, ...
                'LS__ms',{0.1 'xlogslider -2 0'}, ...
                'HS__ms',{1 'xlogslider -.5 2'}, ...
                'thr',{[] 'double'});
            X = fn_control(s,g.par,@(s)computesignal(true));
            
            % variables modified by nested functions
            x = []; tt = []; T = [];
            ktrial = ktrialstart-1;
            spidx = cell(1,ntrial);
            cursp = 0;
            allspikes = [];
            hlthr = []; hl = []; hp = []; hs = [];
            
            % go!
            nexttrial(1)
            
            function nexttrial(inc)
                % change trial
                k = ktrial+inc;
                while k>=1 && k<ntrial && E(k).n<=1
                    k = k+inc;
                end
                if k==0, k=1; elseif k==ntrial+1, k=ntrial; end
                settrial(k)
            end
            
            function settrial(k)
                % set ktrial
                ktrial = k;
                set(g.list,'val',ktrial)
                set(hf,'name',sprintf('guess spikes (trial %i/%i)',ktrial,ntrial))
                % clear
                for i=1:3, cla(ha(i)), end
                cla(g.hb)
                hl=[];
                hp=[];
                % no signal?
                if E(ktrial).n<=1, return, end
                % parameters
                if ~isempty(E(ktrial).spikepar)
                    X.s = E(ktrial).spikepar;
                end
                % spikes already computed
                tt = (0:E(ktrial).n-1)*E(ktrial).dt;
                T = tt(end);
                if size(E(ktrial).spikes,1)==1
                    spidx{ktrial} = find(fn_timevector(E(ktrial).spikes,tt));
                end
                % update slider
                set(slider,'max',T,'width',(2*.6)/T)
                % perform computations and display
                cursp = 1;
                computesignal()
            end
            
            function computesignal(fromcontrol)
                if nargin<1, fromcontrol=false; end
                % signal
                x = E(ktrial).signal(:);
                % band-pass filter
                dtms = E(ktrial).dt*1e3;
                x = fn_filt(x,{X.LS__ms/dtms X.HS__ms/dtms},'b');
                % display signals
                [hl hp hlthr hs] = deal([]);
                nskip = ceil(50/dtms);
                xcut = x(1+nskip:end-nskip);
                ylim = [min(xcut) max(xcut)];
                for i=1:3
                    xlim = get(ha(i),'xlim');
                    h=plot(tt,x,'parent',ha(i),'buttondownfcn',@(u,e)movetime(get(ha(i),'currentpoint')));
                    set(h,'hittest','on')
                    set(ha(i),'xlim',xlim)
                end
                set(ha,'ylim',ylim)
                % initialize threshold?
                if isempty(X.thr)
                    X.thr = max(std(x)*3,max(x)/2);
                end
                % compute spikes
                docomputespikes = (size(E(ktrial).spikes,1)==0); % || docomputespikes;
                computespikes(docomputespikes,~fromcontrol)
            end
            
            function computespikes(docompute,doxlim)
                if nargin<2, doxlim = true; end
                % compute
                if docompute
                    % memorize position of current spike before re-compute
                    if length(E(ktrial).spikes)>=cursp
                        memspike = E(ktrial).spikes(cursp);
                    else
                        memspike = 0;
                    end
                    % compute
                    if X.thr>0, tmp = (x>X.thr); else tmp = (x<X.thr); end
                    tmp = tmp & ~tmp([1 1:end-1]);
                    spidx{ktrial} = find(tmp);
                    % recess time of 1ms
                    d = diff(tt(spidx{ktrial}));
                    f = find(d<.001);
                    spidx{ktrial}(f+1)=[];
                    % new position
                    if ~isempty(spidx{ktrial})
                        d = abs(tt(spidx{ktrial})-memspike);
                        [dum cursp] = min(d); %#ok<ASGLU>
                    end
                end
                % display and save
                displayspikes(doxlim)
            end

            function displayspikes(doxlim)
                if nargin<1, doxlim=true; end
                
                delete([hl hp hlthr hs])
                if isempty(spidx{ktrial})
                    % display threshold
                    for i=1:4
                        if i<=3, hai=ha(i); else hai=g.hb; end
                        hlthr(i) = line([-T 2*T],X.thr*[1 1],'color','r', ...
                            'parent',hai, ...
                            'hittest','on','buttondownfcn',@moveline);
                    end
                    set(ha(1),'xlim',[0 T])
                    % everything else does not exist
                    [hl hp hs] = deal([]);
                else
                    % display all spikes superimposed
                    nframes = round([1 2]*1e-3/E(ktrial).dt);
                    allspikes = fn_triggeravg(x,spidx{ktrial},1,nframes,'noavg');
                    hs = plot((-nframes(1):nframes(2))*E(ktrial).dt*1e3,allspikes,'color',[1 1 1]*.5, ...
                        'parent',g.hb,'buttondownfcn',@showspikeandsave)';
                    set(hs,'hittest','on')
                    % display threshold
                    for i=1:4
                        if i<=3
                            hai = ha(i);
                            tint = [-T 2*T];
                        else
                            hai=g.hb;
                            tint = [-nframes(1) nframes(2)]*E(ktrial).dt*1e3;
                        end
                        hlthr(i) = line(tint,X.thr*[1 1],'color','r', ...
                            'parent',hai, ...
                            'hittest','on','buttondownfcn',@moveline);
                    end
                    set(ha(1),'xlim',[0 T])
                    % display spike markers
                    for i=1:3
                        hl(i) = line(tt(spidx{ktrial}),ones(1,length(spidx{ktrial}))*X.thr, ...
                            'parent',ha(i), ...
                            'linestyle','none','marker','.','color','r', ...
                            'hittest','on','buttondownfcn',@showspikeandsave);
                        hp(i) = line(tt(spidx{ktrial}(1)),X.thr, ...
                            'parent',ha(i), ...
                            'marker','.','color','y', ...
                            'hittest','on','buttondownfcn',@showspikeandsave);
                    end
                end
                % highlight the current spike (if any) and save
                showspikeandsave(doxlim)
            end
            
            function moveline(hl,evnt) %#ok<INUSD>
                fn_moveobject(hl)
                ydata = get(hl,'ydata');
                X.thr = ydata(1);
                set(hlthr,'ydata',ydata)
                computespikes(true)
            end
            
            function showspikeandsave(hl,evnt) %#ok<INUSD>
                % input
                switch nargin
                    case 0
                        doxlim = true;
                    case 1
                        doxlim = hl;
                    case 2
                        doxlim = true;
                        if any(hl==hs)
                            % click on spike time courses -> find which one
                            cursp = find(hl==hs);
                        else
                            % click on spike mark -> find which one
                            hax = get(hl,'parent');
                            p = get(hax,'currentpoint');
                            d = abs(tt(spidx{ktrial})-p(1));
                            [dum cursp] = min(d); %#ok<ASGLU>
                        end
                end
                % show
                if ~isempty(spidx{ktrial})
                    % change view to center on spike
                    idx = spidx{ktrial}(cursp);
                    if doxlim
                        set(ha(2),'xlim',tt(idx)+[-.5 1.5]*.6);
                        set(ha(3),'xlim',tt(idx)+[-.5 1.5]*.02)
                    end
                    set(slider,'value',(tt(idx)-tt(1))/(tt(end)-tt(1)))
                    % show the spike in yellow
                    set(hp,'visible','on','xdata',tt(spidx{ktrial}(cursp)))
                    % highlight this spike in the graph with all spikes
                    % superimposed
                    set(hs,'color',[1 1 1]*.5,'linewidth',1)
                    set(hs(cursp),'color','k','linewidth',2)
                    uistack(hs(cursp),'top')
                elseif doxlim
                    set(hp,'visible','off')
                    set(ha(2),'xlim',tt(1)+[0 2]*.6);
                    set(ha(3),'xlim',tt(1)+[0 2]*.02)
                end
                % save
                E(ktrial).spikes = E(ktrial).t0 + tt(spidx{ktrial});
                E(ktrial).spikepar = X.s;
            end
            
            function movetime(t)
                t = t(1);
                largejump = (t-mean(get(ha(2),'xlim'))>1);
                
                % center ha(2) and ha(3) on t
                set(ha(2),'xlim',t+[-1 1]*.6)
                set(ha(3),'xlim',t+[-1 1]*.02)
                
                if largejump
                    % change current spike to the first visible in ha(2),
                    % or the last before
                    idx = find(E(ktrial).spikes>t-.6 & E(ktrial).spikes<t+.6,1,'first');                    
                    if ~isempty(idx), 
                        % display this new spike and center ha(3) on it
                        cursp = idx;
                        showspikeandsave(false)
                        set(ha(3),'xlim',E(ktrial).spikes(cursp)+[-.5 1.5]*.02)
                    else
                        idx = find(E(ktrial).spikes<=t-.6,1,'last'); 
                        if ~isempty(idx), cursp = idx; showspikeandsave(false), end % this will cause next press on 'next spike' to really go forward rather than backward
                    end
                else
                    idx = find(E(ktrial).spikes>t-.02 & E(ktrial).spikes<t+.02,1,'first');
                    if isempty(idx), idx = find(E(ktrial).spikes<=t-.02,1,'last'); end
                    if ~isempty(idx), cursp = idx; showspikeandsave(false), end % this will cause next press on 'next spike' to really go forward rather than backward
                end
            end
            
            function rmspike(varargin)
                if isempty(cursp), return, end
                spidx{ktrial}(cursp) = [];
                delete(hs(cursp))
                hs(cursp) = [];
                if isempty(spidx{ktrial})
                    cursp = [];
                else
                    cursp = min(cursp,length(spidx{ktrial}));
                end
                set(hl,'xdata',tt(spidx{ktrial}), ...
                    'ydata',ones(1,length(spidx{ktrial}))*X.thr)
                showspikeandsave
            end
            
            function addspike(t)
                [dum idx] = min(abs(tt-t));
                [spidx{ktrial} ord] = sort([idx; spidx{ktrial}]);
                cursp = find(ord==1);
                % show the new spike, but do not change the xlim
                displayspikes(false)
            end
            
            function nextspike(inc)
                cursp = cursp+inc;
                if cursp==0, cursp=1; end 
                if cursp>length(spidx{ktrial}), cursp=length(spidx{ktrial}); end 
                showspikeandsave
            end
            
            function keypress(hfig,evnt) %#ok<INUSL>
                switch evnt.Key
                    case {'rightarrow' 'v'}
                        nextspike(1)
                    case {'leftarrow' 'c'}
                        nextspike(-1)
                    case {'w' 'z'}
                        X.thr = X.thr*(2^(1/10*sign(X.thr)));
                        set(hlthr,'ydata',X.thr*[1 1])
                        computespikes(true)
                    case {'s'}
                        X.thr = X.thr*(2^(-1/10*sign(X.thr)));
                        set(hlthr,'ydata',X.thr*[1 1])
                        computespikes(true)
                    case 'd'
                        nexttrial(1)
                    case 'e'
                        nexttrial(-1)
                    case 'r'
                        rmspike()
                    case {'f' 'g'}
                        if isempty(evnt.Modifier)
                            t = mean(get(ha(2),'xlim'));
                            t1 = t + fn_switch(evnt.Key,'f',-1,1)*1.9*.6;
                        else
                            t = mean(get(ha(3),'xlim'));
                            t1 = t + fn_switch(evnt.Key,'f',-1,1)*1.9*.02;
                        end
                        movetime(t1)
                    case {'a' 'q'}
                        % add a spike
                        if isempty(get(hf,'windowButtonMotionFcn'))
                            set(hf,'windowButtonMotionFcn',' ')
                            return
                        end
                        p = get(ha(3),'currentpoint');
                        addspike(p(1))
                    case 't'
                        % test
                        p = get(ha(3),'currentpoint');
                        disp(p(1))
                        set(hf,'windowButtonMotionFcn',' ')
                end
            end
        end
    end
     
    % Load
    methods (Static)
        function E = loadobj(E)
            if isempty(E), error programming, end
            if E.n==0, return, end
            %             if ~isempty(E.file)
            %                 E.file = tps_locatefile(E.file);
            %                 x = fn_readdatlabview(E.file);
            %                 if isempty(x)
            %                     if E.n>1, error('cannot read data'), end
            %                     return
            %                 end
            if isempty(E.recidx)
                %                     % old version
                %                     nr = size(x,2);
                %                     if nr==1
                %                         E.recidx = 1;
                %                     elseif nr<=3
                %                         hf = figure;
                %                         plot(x)
                prompt = 'Which recording corresponds to electrophysiology?';
                title = 'Question to user';
                choice = num2cell(num2str((1:nr)'));
                answer = questdlg(prompt,title,choice{:},'1');
                E.recidx = str2double(answer);
                close(hf)
                %                     else
                %                         error('cannot decide which recording is electrophysiology')
                %                     end
            end
            %                 E.signal = x(:,E.recidx);
            %                 if length(E.signal)~=E.n, error programming, end
            %                 removeartefact(E)
            %             end
        end
    end
    
    % Misc
    methods
        function access(E)
            keyboard
        end
    end
end

