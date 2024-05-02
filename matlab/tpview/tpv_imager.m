classdef tpv_imager < hgsetget
    
    properties
        V
        grob
        params
        paramspec
        curexp
    end
    
    % Init
    methods
        function I = tpv_imager(V)
            I.V = V;
            init_params(I)
            init_grob(I) % position graphic objects and set their properties
            assignin('base','I',I)
        end
        function init_params(I)
            I.params = struct('duration',0,'interval',0,'ncond',1,'ntrial',10,'binning',1);
            I.paramspec = struct('binning','stepper 1 1 Inf','duration','double','interval','double', ...
                'ncond','stepper 1 1 Inf','ntrial','stepper 1 1 Inf');
        end
        function init_grob(I)
            hp = I.V.panels.acquisition;
            delete(get(hp,'children'))
            g = struct;
            g.controls = uipanel('parent',hp,'pos',[.4 .1 .5 .9],'bordertype','none');
            fn_control(I.params,I.paramspec,g.controls,@(s)UpdateParameters(I,s),'nobutton')
            A = .02; W = .3; B = .95; H = .1;
            uicontrol('parent',hp,'units','normalized','pos',[A B-2*H W 2*H], ...
                'string','CAMERA','callback',@(u,e)CameraTool(I))
            uicontrol('parent',hp,'units','normalized','pos',[A B-3*H W H], ...
                'string','4bits colormap','callback',@(u,e)FourBitsColormap(I))
            uicontrol('parent',hp,'units','normalized','pos',[A B-5*H W 2*H], ...
                'string','INIT','callback',@(u,e)PrepareExperiment(I));
            g.startstop = uicontrol('parent',hp,'units','normalized','pos',[A B-7*H W 2*H], ...
                'string','START','callback',@(u,e)StartStopExperiment(I));
            uicontrol('parent',hp,'units','normalized','pos',[A B-8*H W H], ...
                'string','check frame times','callback',@(u,e)CheckFrameTimes(I));
            g.status = uicontrol('parent',hp,'units','normalized','pos',[.02 0 .98 .1], ...
                'style','text','fontname','fixedwidth','foregroundColor','b','horizontalalignment','left');
            I.grob = g;
        end
    end
    
    % Misc
    methods
        function UpdateParameters(I,s)
            I.params = s;
        end
        function status(I,str)
            set([I.grob.status I.V.grob.statusbar],'string',str)
            drawnow
        end
    end
    
    % Camera
    methods
        function CameraTool(I)
            imaqtool
        end
        function FourBitsColormap(I)
            hf = findall(0,'type','figure','tag','IMAQPreviewFigure');
            if isempty(hf), return, end
            cm = gray(16);
            cm = [cm(1:15,:); [1 .8 .8]; ones(256-16,1)*[1 0 0]];
            colormap(hf,cm)
        end
    end
    
    % Experiment
    methods
        function [ok vid nx ny nframe ncond nbit fps] = PrepareExperiment(I)
            ok = false;
            
            % Fix START button if needed
            set(I.grob.startstop,'string','START','enable','on')
            
            % Video input
            vid = iatbrowser.Browser().currentVideoinputObject; % video input object
            if isempty(vid)
                errordlg 'Press ''Camera'' first to adjust camera settings'
                return
            end
            stoppreview(vid)
            src = getselectedsource(vid); % video source object
            
            % Acquisition parameters
            fps = str2double(src.FrameRate);
            nframe = ceil(I.params.duration * fps);
            if nframe==0
                errordlg('Set acquisition duration first!')
                return
            end
            vid.FramesPerTrigger = nframe;
            ncond = I.params.ncond;
            p = vid.ROIPosition;
            xbin = I.params.binning;
            nx = floor(p(3)/xbin); ny = floor(p(4)/xbin);
            nbit = regexp(vid.VideoFormat,'Y(\d*)_','tokens');
            nbit = nbit{1}{1};
            
            % Continue previous experiment
            doinit = isempty(I.curexp);
            if ~isempty(I.curexp) && I.curexp.nacquired>0 && ~I.curexp.ready
                answer = questdlg('Experiment data already present. What do you want to do?','Imager', ...
                    'Continue current experiment','Start new experiment','Cancel','Continue current experiment');
                switch answer
                    case 'Continue current experiment'
                        % nothing to do
                    case 'Start new experiment'
                        doinit = true;
                    case 'Cancel'
                        return
                end
            elseif ~isempty(I.curexp) && I.curexp.nacquired==0 && ~isequal(I.curexp.size,[nx ny nframe ncond])
                % settings have changed: need to re-init
                doinit = true;
            end
            
            % Prompt for file name and prepare for acquisition
            if doinit
                if isempty(I.curexp), I.curexp = struct; end
                I.curexp.ready = false;
                fselect = fn_savefile('*.mat','Select file for saving data');
                I.curexp.fbase = fn_strrep(fselect,'.mat','','_trial001','');
                I.curexp.nacquired = 0;
                I.curexp.size = [nx ny nframe ncond];
                I.curexp.avgdata = zeros(nx,ny,nframe,ncond,'single');
                % make average data already available in tpview
                Tavg = tps_trial(I.curexp.avgdata,[],1/fps);
                Tavg.status = 's';
                I.V.file_open(Tavg)
            else
                if ~isequal(I.curexp.size,[nx ny nframe ncond])
                    errordlg 'Cannot continue experiment: settings have changed'
                    return
                elseif I.V.ntrial~=I.curexp.nacquired+1
                    errordlg '2P VIEW does not have the expected number of trials'
                    return
                end
            end
            
            % Prepare trigger
            if ~isfield(I.curexp,'trigger')                
                s = daq.createSession('ni');
                s.release() % this seems to unreserve the hardware reserved in a previous session 
                warning('off','daq:Session:onDemandOnlyChannelsAdded')
                s.addDigitalChannel('Dev2','Port0/Line0','OutputOnly');
                warning('on','daq:Session:onDemandOnlyChannelsAdded')
                s.outputSingleScan(0);
                I.curexp.trigger = s;
            end
            
            % Ready
            I.curexp.ready = true;
            fname = [I.curexp.fbase '_trial' num2str(I.curexp.nacquired+1,'%.3i') '.mat'];
            fbase = fn_fileparts(fname,'base');
            str = sprintf('%s cond%i/%i: ready to start',fbase,1,ncond);
            I.status(str)
            if nargout==0, clear ok, else ok=true; end
        end
        function StartStopExperiment(I)
            % Start or stop?
            ustartstop = I.grob.startstop;
            if strcmp(get(ustartstop,'string'),'STOP')
                answer = questdlg('Do you want to finish acquiring current trial?','imager','Finish trial','Interrupt now','Cancel','Finish trial');
                switch answer
                    case 'Finish trial'
                        set(ustartstop,'string','stopping','enable','off')
                    case 'Interrupt now'
                        set(ustartstop,'string','interrupting','enable','off')
                end
                return
            end
              
            % Prepare experiment
            [ok vid nx ny nframe ncond nbit] = PrepareExperiment(I);
            if ~ok, return, end
            day2sec = 3600*24;
            
            % Prepare data
            if I.params.binning==1
                data = zeros(nx,ny,nframe,ncond,['uint' nbit]);
            else
                data = zeros(nx,ny,nframe,ncond,'single');
            end
            frametimes = zeros(nframe,ncond);
            
            % Acquisition loop
            set(ustartstop,'string','STOP')
            I.curexp.ready = false; % not ready for a new start any more (confirmation will be asked to user before a new start)
            for kloop = 1:I.params.ntrial
                % Interrupt
                if ismember(get(ustartstop,'string'),{'stopping' 'interrupting'})
                    break
                end
                % File name
                fname = [I.curexp.fbase '_trial' num2str(I.curexp.nacquired+1,'%.3i') '.mat'];
                fbase = fn_fileparts(fname,'base');
                
                % Acquire
                for kcond = 1:ncond
                    if strcmp(get(ustartstop,'string'),'interrupting'), break, end
                    % start acquisition
                    str = sprintf('%s cond%i/%i: acquire',fbase,kcond,ncond);
                    I.status(str)
                    start(vid)
                    tstarted = now()*day2sec;
                    % send trigger
                    I.curexp.trigger.outputSingleScan(1);
                    pause(.2)
                    I.curexp.trigger.outputSingleScan(0);
                    % wait
                    for i=0:10
                        pause(tstarted+(i/10)*I.params.duration-now()*day2sec)
                        I.status([str ' [' repmat('.',1,i) repmat(' ',1,10-i) ']'])
                        if strcmp(get(ustartstop,'string'),'interrupting'), stop(vid), break, end
                    end
                    wait(vid,2)
                    if strcmp(get(ustartstop,'string'),'interrupting'), break, end
                    tdone = tic;
                    str = sprintf('%s cond%i/%i: get data',fbase,kcond,ncond);
                    I.status(str)
                    % get data
                    [x, ~, metadata] = getdata(vid);
                    frametimes(:,kcond) = datenum(cat(1,metadata.AbsTime))*day2sec-tstarted;
                    x = fn_bin(x,I.params.binning);
                    data(:,:,:,kcond) = permute(x,[2 1 4 3]);
                    % wait for next cond
                    if kcond<ncond, pause(I.params.interval - toc(tdone)), end
                end
                if strcmp(get(ustartstop,'string'),'interrupting'), break, end
                I.curexp.nacquired = I.curexp.nacquired + 1;
                nacq = I.curexp.nacquired;
                
                % Save data
                str = sprintf('%s cond%i/%i: save',fbase,kcond,ncond);
                I.status(str)
                Tk = tps_trial(data,I.V.content.trials(1));
                Tk.status = 'n';
                Tk.user.acquisition = struct('frametimes',frametimes, ...
                    'VideoFormat',vid.VideoFormat,'VideoResolution',vid.VideoResolution, ...
                    'ROIPosition',vid.ROIPosition,'binning_post',I.params.binning);
                Tk.savedata(fname)
                
                % Make data available in tpview
                I.curexp.avgdata = I.curexp.avgdata*(1-1/nacq) + single(data)/nacq;
                Tavg = tps_trial(I.curexp.avgdata,I.V.content.trials(1));
                Tavg.status = 's';
                ktrialcur = I.V.ktrial;
                addtrials(I.V.content,[Tk Tavg])
                rmtrial(I.V.content,nacq) % remove the previous average
                if ktrialcur>=nacq-1
                    I.V.content.ktrial = ktrialcur+1; % no automatic display update
                end
                I.V.display_changeview('chgtrial')
                I.V.file_markchange();
                
                % Wait for next trial
                if kloop<I.params.ntrial
                    pause(I.params.interval - toc(tdone))
                end
            end
            I.status('acquisition done')
            set(ustartstop,'string','START','enable','on')
        end
        function CheckFrameTimes(I)
            T = I.V.content.trials;
            T = T([T.status]=='n');
            nframe = T(1).nfr; ncond = T(1).nc; dt = T(1).dt;
            ntrial = length(T);
            frametimes = zeros(nframe,ncond,ntrial);
            for ktrial=1:ntrial, frametimes(:,:,ktrial) = T(ktrial).user.acquisition.frametimes; end
            tt = (1:nframe)'*dt; % theoretical frame time
            fn_figure('Frame times')
            %             plot(tt,tt,'k','linestyle','--')
            %             hold on
            %             plot(tt,frametimes(:,:))
            %             hold off
            ylabel 'actual frame times'
            plot(tt,fn_subtract(frametimes(:,:),tt))
            ylabel 'actual - theoretical'
            xlabel 'theoretical frame times'
            t0 = mean(frametimes(1,:));
            tend = mean(frametimes(nframe,:));
            dtactual = (tend-t0)/(nframe-1);
            title(sprintf('actual values: t0=%.0fms, dt=%.3fms, fs=%.3fHz',t0*1e3,dtactual*1e3,1/dtactual))
        end
    end
end