function label = tpv_crosscorrelation(V) 
% function f(V) 
% function label = f('label') 

if nargin==0, V=evalin('base','V'); end
if ischar(V) && strcmp(V,'label')
    label = 'Cross-correlations';
    return
end

% checks
if ~V.a4d.G.linkselection
    error 'set V.a4d.G.linkselection to true to enable the selection of time intervals'
end

% variables shared with sub-functions
hf = figure('integerhandle','off','numbertitle','off','name','Cross-correlations');
hasignals = [axes axes];
hacorr = axes;
XROI = [];
XTIME = [];


% graphic objects
uicontrol('style','text','string','SOURCES','units','normalized');
panelROI = uipanel;
uicontrol('style','text','string','TIME WINDOWS','units','normalized');
panelTIME = uipanel;
uicontrol('style','pushbutton','string','NEW WINDOW','units','normalized','callback',@(u,e)newtimesegment())
uicontrol('style','pushbutton','string','SAVE','units','normalized','callback',@(u,e)savewindows())
uicontrol('style','pushbutton','string','LOAD','units','normalized','callback',@(u,e)loadwindows())

% graphic objects position
% (evaluate fn_framedesign(gcf,'tpv_crosscorrelation',true); to re-design)
fn_framedesign(hf,'tpv_crosscorrelation');

% ROI selection
sources = ['first optic' 'all optics' setdiff({V.disppar.timerecording.name},{'data' 'sfr'})];
k1 = find(strcmpi(sources,'EEG')); if isempty(k1), k1=2; elseif ~isscalar(k1), k1 = k1(1); end
s = struct( ...
    'signal1',  {sources{k1} sources}, ...
    'signal2',  {'all optics' sources});
XROI = fn_control(s,panelROI,@(x)updateROI());
% TIME selection
controls = struct( ...
    'label',    {'trial' 'condition' 'interval'}, ...
    'style',    'edit', ...
    'type',     'double', ...
    'default',  {1 0 [0 2]}, ...
    'length',   {3 3 8}, ...
    'labellength', 2);
spec = struct('name','WINDOW','controls',controls);
XTIME = fn_supercontrol(panelTIME,spec,@(x)updateTIME());


% SUB-FUNCTIONS
    function newtimesegment()
        timesel = V.a4d.F.selection.getselset(3).t.set;
        if ~isscalar(timesel) || ~isscalar(timesel.poly) || ~strcmp(timesel.poly.type,'line1D')
            errordlg 'please select a single time interval in tpview before pressing NEW WINDOW'
            return
        end
        condarray = tps_readconditionname(V.disppar.timecond,V.trial.nc);
        cond = [condarray{:}];
        if isempty(cond), errordlg 'please display a set of non-calculated conditions', return, end
        time = timesel.poly.points;
        time = time * fn_switch(V.a4d.F.units{3},'s',1,'ms',1e-3);
        if diff(time)<.5, errordlg 'time window should be at least 500ms long', return, end
        XTIME.x(end+1) = struct('name','WINDOW','active',true,'value',{{V.ktrial cond-1 time}});
        [XTIME.specs.controls.default] = deal(V.ktrial,cond-1,time); % values of this new line become default for next ones
        if XTIME.immediateupdate, updateTIME(), end
    end

    function savewindows()
        default = getappdata(hf,'filename'); 
        if isempty(default), default = {}; else default = {default}; end
        fname = fn_savefile('*_crosscorr.mat','Select file where to save cross-correlation windows',default{:});
        if fname==0, return, end
        tokens = regexp(fname,'(_crosscorr|\.mat)*$');
        if ~isempty(tokens), fname = fname(1:tokens-1); end
        fname = [fname '_crosscorr.mat'];
        setappdata(hf,'filename',fname)
        windows = XTIME.x;
        timevisu = V.panels.timevisucontrol.s;
        fn_savevar(fname,windows,timevisu)
    end

    function loadwindows()
        fname = fn_getfile('*_crosscorr.mat');
        if fname==0, return, end
        s = fn_loadvar(fname);
        XTIME.x = s.windows;
    end

    function updateROI()
        crosscorrelation()
    end

    function updateTIME()
        crosscorrelation()
    end

    function crosscorrelation()
        signalorigins = {XROI.s.signal1 XROI.s.signal2};
        windows = XTIME.x([XTIME.x.active]);
        if isempty(windows), for ha=[hasignals hacorr], cla(ha), end, return, end
        
        % Temporal resolution, number of signals
        dtt = zeros(1,2);
        nsig = zeros(1,2);
        for k=1:2
            switch signalorigins{k}
                case 'all optics' 
                    dtt(k) = V.trial.dt_sec;
                    nsig(k) = V.content.nx;
                case 'first optic'
                    dtt(k) = V.trial.dt_sec;
                    nsig(k) = 1;
                otherwise
                    tidx = V.content.getrecording(signalorigins{k},windows(1).value{1}(1)).tidx;
                    dtt(k) = diff(tidx(1:2));
                    nsig(k) = 1;
            end
        end
        [dt klowres] = max(dtt);
        maxlag = ceil(.25/dt); % max time lag of 250ms
        
        % Get the signals
        oldtimecond = V.disppar.timecond;
        setpar(V,'timecond','all conditions') % do this in order to easily get the signals from any condition
        nwin = length(windows);
        signals = cell(2,0);
        for i=1:nwin
            trials = windows(i).value{1};
            for ktrial = trials
                kcond = windows(i).value{2}+1;
                time = windows(i).value{3};
                sigi = cell(2,1);
                for k = [klowres setdiff(1:2,klowres)]
                    % get signal from specified trial
                    switch signalorigins{k}
                        case {'all optics' 'first optic'}
                            % signal from trial
                            switch signalorigins{k};
                                case 'first optic'
                                    x = V.content.getslice(ktrial,1);
                                case 'all optics'
                                    x = V.content.getslice(ktrial);
                            end
                        otherwise
                            x = V.content.getrecording(signalorigins{k},ktrial);
                            if isempty(x) || any(strcmp({x.tag},'empty'))
                                errordlg(['no ' signalorigins{k} ' recording for trial ' num2str(ktrial)])
                                for ha=[hasignals hacorr], cla(ha), end
                                return
                            end
                    end
                    tidx = x(1).tidx;
                    % select condition
                    x = cat(3,x.dataop);
                    x = permute(x(:,kcond,:),[1 3 2]);
                    % filter - it is not nice to do it that way!!!
                    timevisu = V.panels.timevisucontrol.s;
                    LHS = fn_switch(lower(signalorigins{k}), ...
                        {'all optics' 'first optic'},     {timevisu.LS timevisu.HS}, ...
                        'eeg',  {timevisu.LS_eeg timevisu.HS_eeg}, ...
                        'electrophysiology',    {timevisu.LS_elphy timevisu.HS_elphy}, ...
                        {[] []});
                    if ~isempty([LHS{:}])
                        ls = LHS{1}/dtt(k); if isempty(ls), ls=0; end
                        hs = LHS{2}/dtt(k); if isempty(hs), hs=0; end
                        x = fn_filt(x,[ls hs],'m');
                    end
                    % select time window
                    if diff(tidx(1:2))~=dtt(k), error 'time resolution', end
                    if k==klowres
                        % we use this sampling rate
                        % we will also take a wider window to prevent edge
                        % effects
                        % This signal is called the 'lowres' one, and noted y.
                        okfr = find(tidx>=time(1) & tidx<=time(2));
                        ntwin = length(okfr);
                        %                         okfr = okfr(1)-maxlag:okfr(end)+maxlag;
                        if okfr(1)<1 || okfr(end)>size(x,1)
                            errordlg 'window is too close to an edge'
                            for ha=[hasignals hacorr], cla(ha), end
                            return
                        end
                        tidxwin = tidx(okfr)';
                        x = x(okfr,:,:);
                    else
                        % This signal is called the 'interpolated' one, and noted x. 
                        x = interp1(tidx,x,tidxwin);
                    end
                    x = fn_normalize(x,1,'std');
                    sigi{k} = shiftdim(num2cell(x,[1 2]),1);
                end
                % add signals
                signals = [signals [sigi{1}; sigi{2}]]; %#ok<AGROW>
            end
        end
        setpar(V,'timecond',oldtimecond)
        
        % Display the signals
        nwin = size(signals,2);
        for k=1:2
            for i=1:nwin
                sigi = signals{k,i};
                hl=plot((0:size(sigi,1)-1)*dt*1e3,sigi,'parent',hasignals(1));
                if isempty(strfind(signalorigins{k},'optic')), set(hl,'color',[1 1 1]*.4), end
                hold(hasignals(1),'on')
            end
            if k==2, hold(hasignals(1),'off'), end
            fn_axis(hasignals(1),'tight',[1 1.1])
            %ylabel(hasignals(k),['signal' num2str(k)])
        end
        
        % Cross-correlation
        cross = zeros(2*maxlag+1,nsig(1),nsig(2),nwin);
        totallength = 0;
        for k = 1:nwin
            lengthk = size(signals{1,k},1);
            totallength = totallength + lengthk;
            for i=1:nsig(1)
                for j=1:nsig(2)
                    cross(:,i,j,k) = xcov(signals{1,k}(:,j),signals{2,k}(:,i),maxlag,'coef') * lengthk;
                end
            end
        end
        cross = squeeze(sum(cross,4)/totallength);
        
        % Display of correlation
        plot((-maxlag:maxlag)*dt*1e3,cross,'parent',hacorr)
        grid(hacorr,'on')
        fn_axis(hacorr,'tight',[1 1.1])
        ylabel(hacorr,'time lag (ms)')
        ylabel(hacorr,'cross-correlation')

    end

end


