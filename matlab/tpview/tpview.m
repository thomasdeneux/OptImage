classdef tpview < interface
    % function V = tpview(['optimage',][fname|data][,'noread'])
    
    % Programming notes:
    %
    % DATA DIMENSIONS
    % data dimensions in G are x-y-t-condition-trial-multichannel
    
    % Properties
    % (graphic objects)
    properties (SetAccess='private')
        skin        % '2P VIEW' or 'OptImage'
        modules     % set of flags indicating optional features, will be set at init depending on skin
        % grob      [inherited]
        % menus     [inherited]
        % options   [inherited]
        panels
        a4d
        timer
        precomp = struct; % store here some variables to avoid recomputing several time
    end
    % (content)
    properties (SetAccess='private')
        content % tpv_content object
        internpar
        disppar
        savingpar
    end
    % (shortcuts)
    properties (Dependent, SetAccess='private')
        % (trial information)
        ntrial
        trial
        data
        dataop
        sfr
        file
        fullinfo
        sizes
        nx
        ny
        nfr
        nc
        linedur
        tidx
        nsel
        % (signals)
        signal
    end
    properties (Dependent, SetObservable)
        ktrial
    end
    properties (Dependent)
        addinfo
        type
        scanning
        dx
        dy
        dz
        dt
        t0
    end
    
    % Events
    events
        EventCloseFile
    end
    
    % Constructor / Destructor
    methods
        function V = tpview(varargin)
            % input
            theskin = '2P VIEW'; fnameordata = []; F = []; noread = false;
            for i=1:nargin
                a = varargin{i};
                if ischar(a) 
                    switch a
                        case {'2P VIEW' 'OptImage'}
                            theskin = a;
                        case 'noread'
                            noread = true;
                        otherwise
                            fnameordata = a;
                    end
                elseif isa(a,'focus')
                    F = a;
                else
                    fnameordata = a;
                end
            end
            
            % default options
            opt = struct( ...
                'preferences',struct, ...
                'condcalc',{{}}, ...
                'condmem',{cell(2,0)}, ...
                'favoritefolders',struct('name',{},'path',{}) ...
                );
            
            % init interface
            hf = figure('integerhandle','off','visible','off');
            V = V@interface(hf,theskin,opt);
            %             if ~ok, return, end
            set(hf,'visible','on')
            c = fn_watch(V.hf); %#ok<NASGU>
            
            % skin
            % optional features: 'sfr', 'timevisu', 'acquisition',
            % 'MLspike', 'developer'
            V.skin = theskin;
            switch theskin
                case '2P VIEW'
                    V.modules = {'developer' 'spikes' 'acquisition' 'sfr'};
                otherwise
                    % no module by default
                    V.modules = {};
            end
            
            
            % default content
            V.content = tpv_content;
            
            % init parameters first
            init_parameters(V)
            
            % create and position main graphic objects and containers
            init_display(V)
            interface_end(V)
            
            % initializations
            if ~isempty(F), V.a4d.F = F; end % focus object might be provided by user to synchronize with other window(s)
            init_4d(V) 
            init_slidertrial(V)
            init_header(V)
            init_statusbar(V)
            init_panels(V)
            init_actionbutton(V)
            tpview_keypress(V,'init')
            fn_progress('in',V.grob.statusbar)
            V.timer = timer; %#ok<CPROP>
            
            % set some parameters
            if noread, setpar(V,'readdata',false), end
                
            % load file?
            if ~isempty(fnameordata)
                data_loadfile(V,fnameordata)
            elseif strcmp(theskin,'2P VIEW')
                % it can be useful for developing purposes to show some
                % random data
                x = rand(30,20,100,2,3);
                x = fn_filt(x,15,'l',3);
                x = fn_filt(x,10,'l',[1 2]);
                x = x + randn(size(x))*.01;
                data_loadfile(V,x)
            end
        end
        function delete(V)
            if isfield(V.grob,'time')
                hl = getappdata(V.grob.time,'tpview_display_stimandevents_update');
                delete(hl)
            end
            notify(V,'EventCloseFile')
            delete(V.hf)
            delete(V.timer)
            try delete(V.a4d.Gkeeper), end %#ok<TRYNC>
        end
    end
    
    % INITIALIZATIONS
    methods
        function init_parameters(V)
            internpar0 = tpv_internpar;
            disppar0 = struct( ...
                'im1',          'dataop', ...
                'im1cond',      'C0', ...
                'im1spikecomb',  false, ...
                'im2',          'data', ...
                'im2cond',      'C0', ...
                'im2spikecomb',  false, ...
                'morelist',      struct('SI',{},'SImodel',{},'control',{}), ...
                'time',         'dataop', ... % if several data displayed = name MUST be '_' separated
                'timecond',     'all conditions', ...
                'timesignal',   'signalop', ...
                'timespike',    'spikeboth', ...
                'timefft',      '', ...
                'timerecording', struct('name',{}, 'val',{}), ... % struct('name',{'heart' 'electrophysiology' 'eeg' 'data' 'sfr'}, 'val',{false false false false false})
                'timealltrials',false, ...
                'timealltrialsfilterstim', false, ...
                'timeyzero',    false, ...
                'timerms',      false, ...
                'timermsdcomp', true, ...
                'currentsel',   1, ...
                'curvdist',     -.2, ...
                'showfilter',   true, ...
                'LS',           [], ...
                'HS',           [], ...
                'timethr',      [], ...
                'ncatspike',    5, ...
                'trialsublist', [], ...
                'showstim',     true, ...
                'showevents',   true, ...
                'editevents',   false, ...
                'filelistgroup',true, ...
                'toggle',       'normal' ...
                );
            savingpar0 = struct( ...
                'chg',            false, ...
                'savename',         '' ...
                );
            V.internpar = internpar0;
            V.content.internpar = V.internpar; % share same internal parameters
            V.disppar = disppar0;
            V.savingpar = savingpar0;
            preferences(V,'init')
        end
        function init_panels(V)
            % remember presently selected panel if any
            selobj = get(V.grob.panelselect,'SelectedObject');
            if ~isempty(selobj)
                selidx = (V.panels.select.allobj==selobj);
                selname = V.panels.select.allname{selidx};
            else
                selname = 'files';
            end
            
            % structure template
            V.panels = struct; 
            
            % panel selection
            ps = V.grob.panelselect;
            set(ps,'bordertype','none', ...
                'SelectionChangeFcn',@(hp,evnt)tpview_paneltoggle(V,evnt))
            cpanel = { ...
                'acquisition'   'acquisition'   ismember('acquisition',V.modules)
                'files'         'files'         1
                'dataop'        'data op.'      1
                'signalop'      'signal op.'    1
                'timevisu'      'time visu'     ismember('timevisu',V.modules)
                'selection'     'selection'     1
                'more'          'more'          1
                };
            okpanel = [cpanel{:,3}]; cpanel(~okpanel,:)=[];
            panelnames = row(cpanel(:,1));
            npanel = length(panelnames);
            delete(get(ps,'children'))
            select = struct('allname',{panelnames},'allobj',gobjects(1,npanel));
            for kpanel = 1:npanel
                u = uicontrol('parent',ps, ...
                    'style','togglebutton','units','normalized', ...
                    'position',[(kpanel-1)/npanel 0 1/npanel 1],'string',cpanel{kpanel,2});
                select.(panelnames{kpanel}) = u;
                select.allobj(kpanel) = u;
            end
            V.panels.select = select;
            
            % panels
            p = V.grob.panels;
            delete(get(p,'children'))
            set(p,'bordertype','none');
            for kpanel = 1:npanel
                V.panels.(panelnames{kpanel}) = uipanel('parent',p,'visible','off');
            end
            
            % focus
            for kpanel = 1:npanel
                tpview_panelhide(V,panelnames{kpanel})
            end
            set(V.grob.panelselect,'SelectedObject',select.(selname))
            tpview_panelshow(V,selname)
            
            % content
            if ismember('acquisition',V.modules)
                % Acquisition is handled by a tpv_imager object
                if ~isfield(V.panels,'imager')
                    V.panels.imager = tpv_imager(V);
                else
                    V.panels.imager.init_grob();
                end
            end
            panels_files(V)
            panels_data(V)
            panels_signal(V)
            if ismember('timevisu',V.modules), panels_timevisu(V), end
            panels_selection(V)
            panels_more(V)
        end
        function init_actionbutton(V)
            % action button
            u = V.grob.actionbutton;
            m = uicontextmenu('parent',V.hf);
            set(u,'enable','off', ...
                'uicontextmenu',m,'callback','tpv_actionbutton')
            % create tpv_actionbutton.m if necessary
            fname = fullfile(interface.usercodepath,'tpv_actionbutton.m');
            if ~exist(fname,'file')
                fn_savetext('',fname);
            end
            % context menu
            m1=uimenu(m,'label','enable action button','checked',get(u,'enable'), ...
                'callback',@enabledisable);
            uimenu(m,'label','edit action code', ...
                'callback',@(hu,e)edit('tpv_actionbutton'))
            function enabledisable(hu,e) %#ok<INUSD>
                addpath(interface.usercodepath) % in case path was re-initialized
                val = get(u,'enable');
                val = fn_switch(val,'toggle');
                set(u,'enable',val)
                set(m1,'checked',val)
            end
        end
        function init_menus(V)
            V.menus.items = struct;
            % to prevent the bug in Marseill with 'unstable' menu, create a
            % useless menu (so that there will be always at least one menu)
            m = uimenu('parent',V.hf,'label',V.skin);
            drawnow
            % first menus created by the 'interface' parent class
            init_menus@interface(V)
            % then each menu in an independ function
            V.menus.sep1 = uimenu('parent',V.hf,'label','|','enable','off');
            menus_2pview(V)
            menus_pictures(V)
            menus_moredisplays(V)
            menus_timecourses(V)
            V.menus.sep2 = uimenu('parent',V.hf,'label','|','enable','off');
            menus_trials(V)
            menus_electrophys(V)
            if ismember('spikes',V.modules), menus_spikes(V), end
            V.menus.sep3 = uimenu('parent',V.hf,'label','|','enable','off');
            menus_scripts(V)
            V.menus.sep4 = uimenu('parent',V.hf,'label','|','enable','off');
            menus_help(V)
            % set condition menus
            V.menus.items.ncond = 1; % current setting valid for a unique condition
            data_conditionmenus(V)
            % now delete the useless menu
            drawnow
            delete(m)
        end
    end
    methods (Access='private')
        function init_display(V)
            % figure
            hf = V.hf;
            set(hf, ...
                'tag','2PVIEW main','handlevisibility','off', ...
                'defaultuicontrolfontsize',V.options.preferences.fontsize, ...
                'defaultaxesfontsize',V.options.preferences.textfontsize, ...
                'defaulttextfontsize',V.options.preferences.textfontsize, ...
                'closerequestfcn',@(h,evnt)delete(V));
            
            % structure with handles of graphic objects
            grob.hf = hf;
            
            % axes
            grob.im1 = axes('parent',hf,'units','pixels','xtick',[],'ytick',[]);
            grob.im2 = axes('parent',hf,'units','pixels','xtick',[],'ytick',[]);
            grob.time = axes('parent',hf,'units','pixels');
            
            % sliders
            grob.slidertime = uipanel('parent',hf,'units','pixel');
            grob.slidertrial = uipanel('parent',hf,'units','pixel');
            
            % information display
            grob.header = uicontrol('parent',hf,'style','text');
            grob.statusbar = uicontrol('parent',hf,'style','text');
            grob.timealltrials = uicontrol('parent',hf,'style','text', ...
                'string',fn_switch(V.disppar.timealltrials,'TRIALS','REGIONS'), ...
                'fontweight','bold');
            
            % panels
            grob.panels = uipanel('parent',hf,'units','pixels');
            grob.panelselect = uibuttongroup('parent',hf,'units','pixels');
            
            % action button
            grob.actionbutton = uicontrol('parent',hf,'style','pushbutton');
            
            % skin
            if strcmp(V.skin,'OptImage')
                set(V.hf,'color',fn_colorset('ivory'))
                set([grob.statusbar grob.timealltrials],'backgroundcolor',fn_colorset('ivory'))
            end
            
            % set in V
            V.grob = grob;
        end
        function init_4d(V)
            grob = V.grob;
            
            % structure template
            a = struct( ...
                'F',            [], ...
                'G',            [], ...
                'Gkeeper',      [], ...
                'SI1',          [], ...
                'SI2',          [], ...
                'SIt',          [], ...
                'im1',          [], ...
                'im2',          [], ...
                'time',         [], ...
                'slidertime',   [], ...
                'listeners',    [], ... % listeners to G, SIt
                'slidertrial',          [], ... % this will be set in V.init_slidertrial
                'imdata',               [], ... % this will be set in V.display_image
                'timedisplayfactor',    [] ...  % this will be set in V.display_updatesize
                );
            
            % focus
            if isstruct(V.a4d) && isfield(V.a4d,'F')
                a.F = V.a4d.F;
            else
                a.F = focus('labels',{'x' 'y' 'time' 'condition' 'trial' 'channel'},'units',{'um' 'um' 's'});
            end
            
            % geometry
            a.G = rotation(a.F,'linkselection',true);
            % 'G keeper' points on G, just to prevent G to autodelete
            % because he thinks nobody is pointing on it
            a.Gkeeper = fn4Dhandle;
            addparent(a.Gkeeper,a.G)
            
            % slice info
            a.SI1 = projection(a.G,[1 2],'dimsplus',6); % dimension 6 corresponds to multi-channel
            a.SI2 = projection(a.G,[1 2],'dimsplus',6); % dimension 6 corresponds to multi-channel
            setdata(a.SI1,'data',0);
            setdata(a.SI2,'data',0);
            a.SIt = sliceinfo(1);
            a.SIt.slice = tps_signalx.empty;
            
            % displays
            pref = V.options.preferences;
            scanline = strcmp(V.trial.type,'linescan') || V.content.timeline;
            seldims = fn_switch(scanline,'x','xy');
            a.im1 = activedisplayImage(a.SI1,'in',grob.im1, ...
                'scaledisplay','xbar','clipmode','data', ...
                'seldims',seldims,'shapemode','ellipse', ...
                'channelcolors',[0 1 0; 1 0 0], ...
                'navigation',pref.imnav,'scrollwheel',pref.imscroll ...
                );
            a.im2 = activedisplayImage(a.SI2,'in',grob.im2, ...
                'scaledisplay','xbar','clipmode','data', ...
                'seldims',seldims,'shapemode','ellipse', ...
                'channelcolors',[0 1 0; 1 0 0], ...
                'navigation',pref.imnav,'scrollwheel',pref.imscroll ...
                );
            a.time =  activedisplayPlot(a.SIt,'in',grob.time, ...
                'scaledisplay','ybar', ...
                'slicedisplayfun',@(D,slice)display_timeplot(V,D,slice), ...
                'navigation',pref.tcnav,'scrollwheel',pref.tcscroll ...
                );
            try delete(grob.slidertime), end
            V.grob.slidertime = uipanel('parent',V.hf,'units','pixel','pos',V.options.positions(1).value.slidertime);
            a.slidertime = activedisplaySlider(a.SIt,'in',V.grob.slidertime,'layout','down');
            
            % listeners
            a.listeners = event.listener(a.G,'ChangeView', ...
                @(hs,evnt)display_changeview(V,evnt,'G'));
            a.listeners(2) = event.listener(a.SIt,'ChangeView', ...
                @(hs,evnt)display_changeview(V,evnt,'SIt'));
            
            % that's it!
            V.a4d = a;
        end
        function init_slidertrial(V)
            try delete(V.grob.slidertrial), end
            V.grob.slidertrial = uipanel('parent',V.hf,'units','pixel','pos',V.options.positions(1).value.slidertrial);
            V.a4d.slidertrial = fn_slider(V.grob.slidertrial, ...
                'scrollwheel','default', ...
                'min',1,'max',V.ntrial,'mode','point','layout','down');
            if V.ntrial>1, display_slidertrial(V), end
        end
        function init_header(V)
            set(V.grob.header, ...
                'fontname','monospace','fontsize',8, ...
                'horizontalalignment','left' ...
                )
        end
        function panels_files(V)
            hp    = V.panels.files;
            items = struct;
            
            % panel position
            [W H] = fn_pixelsize(hp);
            
            % create controls
            items.list  = uicontrol('parent',hp,'style','listbox','max', 2, ...
                'position',[0 20 W H-20], ...
                'callback',@(ho,evnt)file_select(V,'files'));
            N = 5;
            items.open  = uicontrol('parent',hp,'string','Open', ...
                'position',[W*0/N 0 W/N 20], ...
                'callback',@(ho,evnt)file_select(V,'open'));
            items.upd   = uicontrol('parent',hp,'string','Update', ...
                'position',[W*1/N 0 W/N 20], ...
                'callback',@(ho,evnt)file_select(V,'update'));
            items.group = uicontrol('parent',hp, ...
                'string',fn_switch(V.disppar.filelistgroup,'Ungroup','Group'), ...
                'position',[W*2/N 0 W/N 20], ...
                'callback',@(ho,evnt)file_select(V,'group'));
            str = {'Favorite folders' V.options.favoritefolders.name 'organize...'};
            items.favorites = uicontrol('parent',hp,'style','popupmenu','string',str, ...
                'position',[W*3/N 0 2*W/N 20], ...
                'callback',@(ho,evnt)file_select(V,'favoritefolders'));
            V.panels.filesitems = items;
            
            % context menu for files
            m = uicontextmenu('parent',V.grob.hf);
            V.panels.filesmenu = m;
            uimenu(m,'label','open average', ...
                'callback',@(hu,evnt)file_select(V,'avg'))
            uimenu(m,'label','open concatenate', ...
                'callback',@(hu,evnt)file_select(V,'cat'))
            uimenu(m,'label','open add', ...
                'callback',@(hu,evnt)file_select(V,'add'))
            m1 = uimenu(m,'label','open special');
            uimenu(m1,'label','bin and open', ...
                'callback',@(hu,evnt)file_select(V,'bin'))
            uimenu(m1,'label','bin and average', ...
                'callback',@(hu,evnt)file_select(V,'bin&avg'))
            uimenu(m1,'label','bin and concatenate', ...
                'callback',@(hu,evnt)file_select(V,'bin&cat'))
            uimenu(m1,'label','average and add', ...
                'callback',@(hu,evnt)file_select(V,'avg&add'))
            uimenu(m1,'label','bin, average and add', ...
                'callback',@(hu,evnt)file_select(V,'bin&avg&add'))
            uimenu(m,'label','rename file(s)','separator','on', ...
                'callback',@(hu,evnt)file_select(V,'rename'))
            uimenu(m,'label','delete file(s)', ...
                'callback',@(hu,evnt)file_select(V,'delete'))
            uimenu(m,'label','edit acq. settings','separator','on', ...
                'callback',@(hu,evnt)file_select(V,'acqsettings'))
            set(V.panels.filesitems.list,'uiContextMenu',m);
            
            % no scroll action when scrolling over the file list
            fn_scrollwheelregister(V.panels.files,'mask')
            
            % list files in current directory
            file_list(V)
        end
        function panels_signal(V)
            hp = V.panels.signalop;
            op = fn_structmerge(tpv_content.signalopdef_default,V.content.signals(1).opdef,'skip');
            spec = tpv_content.signalopdef_spec();
            V.panels.signalcontrol = fn_control(op,@(s)data_signalop(V,s),spec,hp);
        end
        function panels_data(V)
            hp = V.panels.dataop;
            delete(get(hp,'children'))
            
            % get the specification for data operations
            specs = tps_dataopdef.operation_specifications(V);
            V.panels.datacontrol = fn_supercontrol(hp,specs, ...
                @(x)data_opcallback(V,x,V.panels.datacontrol.activechg));
            
            % does the operation apply to all trials or only to current
            % trial?
            siz = fn_pixelsize(hp);
            V.panels.dataitems.linked = uicontrol('parent',hp,'style','checkbox', ...
                'pos',[siz(1)-80 siz(2)-18 78 16], ...
                'string','linked','backgroundcolor',[1 1 1]*.6, ...
                'callback',@(u,e)chglink(u) ...
                ); 
            function chglink(u)
                setdataopdef(V.content,'link',get(u,'value'))
                % update opdef display
                data_dataopdisplay(V)
                % update data display
                imgidx = ~isempty(regexp(V.disppar.im1,'op$')) + 2*~isempty(regexp(V.disppar.im2,'op$')); %#ok<RGXP1>
                if V.internpar.guessspikes, imgidx = 3; end % i don't remember why guessspikes flag requires updating images...
                dotime = any(strfind(V.disppar.time,'op')) || V.internpar.guessspikes;
                display_changeview(V,'datamode',imgidx,dotime,false)
            end
            V.panels.dataitems.mem = uicontrol('parent',hp,'string','M', ...
                'pos',[siz(1)-18 siz(2)-17 16 14], ...
                'backgroundcolor',[1 1 1]*.6, ...
                'callback',@(u,e)data_opmem(V) ...
                );
            
            % now, sets the positions of controls according to the current
            % opdef!
            data_dataopdisplay(V)
        end
        function panels_timevisu(V)
            hp = V.panels.timevisu;
            s = struct( ...
                'curvdist',     {0  'xloglogslider -2 1'    'distance btw. curves'}, ...
                'showfilter',   {true   'logical'           'show filter'}, ...
                'LS',           {[] 'xlogslider -4 1 .01 %.2e [0]'   'low-pass'}, ...
                'HS',           {[] 'xlogslider -4 1 .01 %.2e [0]'   'high-pass'}, ...
                'threshold',    {[] 'xlogslider -6 0 .01 %.2e [1e-3]' 'threshold'}, ...
                'thrdisplay',   {'curve' {'curve' 'line' 'spikes'} 'show as'}, ...
                'LS_eeg',           {[] 'xdouble [.1]'   'LS_eeg'}, ...
                'HS_eeg',           {[] 'xdouble [1]'   'HS_eeg'}, ...
                'LS_elphy',           {[] 'xdouble [1e-4]'   'LS_elphy'}, ...
                'HS_elphy',           {[] 'xdouble [.1]'   'HS_elphy'} ...
                );
            V.panels.timevisucontrol = fn_control(s,@(x)displaysignals(V.a4d.time),hp);
        end
        function panels_selection(V)
            hp = V.panels.selection;
            delete(get(hp,'children'))
            set(hp,'defaultuicontrolunits','normalized')
            items = struct;
            
            % position parameters
            A = .025; B = .35; C = A+B; D = 1-C-A;
            N = 5; H = .12;
            E = (1-N*H)/(N+1); K = (1-E)/N;
            
            % init controls
            uicontrol('parent',hp,'position',[A E+K*2 B H], ...
                'style','text','horizontalalignment','left', ...
                'string','load selection');
            
            items.list = uicontrol('parent',hp,'position',[C E+K*2 D H], ...
                'style','popupmenu','string','blabla', ...
                'callback',@(hu,evnt)selection_load(V));
            
            items.input = uicontrol('parent',hp,'position',[A E+K*1 D H], ...
                'style','edit','horizontalalignment','left', ...
                'backgroundcolor',[.9 .9 .9], ...
                'string','base selection');
            
            uicontrol('parent',hp,'position',[A+D E+K*1 B/2 H], ...
                'style','pushbutton', ...
                'string','save', ...
                'callback',@(hu,evnt)selection_save(V,'save'));
            
            uicontrol('parent',hp,'position',[A+D+B/2 E+K*1 B/2 H], ...
                'style','pushbutton', ...
                'string','delete', ...
                'callback',@(hu,evnt)selection_save(V,'delete'));
            
            V.panels.selitems = items;
            
            % set selection list
            selection_showlist(V)
        end
        function panels_more(V)
            hp = V.panels.more;
            delete(get(hp,'children'))
            set(hp,'defaultuicontrolunits','normalized')
            items = struct;
            
            % position parameters
            A = .025; B = .2; C = .4;
            N = 5; H = .12;
            E = (1-N*H)/(N+1); K = (1-E)/N;
            
            % controls
            items.movie = uicontrol('parent',hp,'position',[A E+K*4 B H], ...
                'style','togglebutton', ...
                ... 'interruptible','off','busyaction','cancel', ...
                'string','MOVIE', ...
                'callback',@(u,evnt)display_movie(V));
            items.speed = uicontrol('parent',hp,'position',[A+B+A E+K*4 C H], ...
                'style','slider', ...
                'min',-.5,'max',3,'value',1, ... % log scale, between 0.3 and 1000 frames/s
                'tooltip','adjust movie speed, use middle or right click to set to real time', ...
                'callback',@(u,evnt)display_movie(V), ...
                'buttondownfcn',@(u,e)display_movie(V,'realtime'));
            items.p1 = uicontrol('parent',hp,'position',[3*A+B+C E+K*4 B H], ...
                'style','checkbox','string','Picture 1 only');
            uicontrol('parent',hp,'position',[A E+K*3 B H], ...
                'string','toggle display', ...
                'callback',@(u,evnt)tpview_toggledisplay(V))
            
            % init the timer
            set(V.timer,'executionmode','fixedRate', ...
                'busymode','error', ...
                'startfcn',@(u,evnt)set(items.movie,'value',1), ...
                'stopfcn',@(u,evnt)moviestop(items.movie,V.grob.statusbar), ...
                'errorfcn',@(u,evnt)stop(u))
            % local definition of stop function
            function moviestop(hmovie,ht)
                set(hmovie,'value',0)
                set(ht,'string','')
            end
            % save items
            V.panels.moreitems = items;
        end
        function init_statusbar(V)
            set(V.grob.statusbar, ...
                'horizontalalignment','left')
        end
    end
    
    % MENUS
    methods
        function menus_2pview(V)
            items = V.menus.items;
            m = V.menus.interface;
            
            % open/save
            m1=uimenu(m,'label','New window...','separator','on');
            uimenu(m1,'label','independent','accelerator','n', ...
                'callback',@(hu,evnt)tpview(V.skin))
            uimenu(m1,'label','linked', ...
                'callback',@(hu,evnt)tpview(V.skin,V.a4d.F))
            uimenu(m,'label','Open...','accelerator','o', ...
                'callback',@(hu,evnt)file_open(V))
            % (open specials)
            m1 = uimenu(m,'label','Open special');
            uimenu(m1,'label','add to trials', ...
                'callback',@(hu,evnt)file_open(V,[],'add'))
            uimenu(m1,'label','average', ...
                'callback',@(hu,evnt)file_open(V,[],'avg'))
            uimenu(m1,'label','concatenate', ...
                'callback',@(hu,evnt)file_open(V,[],'cat'))
            uimenu(m1,'label','bin', ...
                'callback',@(hu,evnt)file_open(V,[],'bin'))
            uimenu(m1,'label','bin and average', ...
                'callback',@(hu,evnt)file_open(V,[],'bin&avg'))
            uimenu(m1,'label','bin and concatenate', ...
                'callback',@(hu,evnt)file_open(V,[],'bin&cat'))
            uimenu(m1,'label','average and add', ...
                'callback',@(hu,evnt)file_open(V,[],'avg&add'))
            uimenu(m1,'label','bin, average and add', ...
                'callback',@(hu,evnt)file_open(V,[],'bin&avg&add'))
            % (save)
            uimenu(m,'label','Save','accelerator','s', ...
                'callback',@(hu,evnt)file_save(V))
            uimenu(m,'label','Save as...', ...
                'callback',@(hu,evnt)file_save(V,[]))
            uimenu(m,'label','Comments', ...
                'callback',@(hu,evnt)file_comments(V))
            
            % repair, debug and bug report
            % (repairs)
            uimenu(m,'label','Auto-repair (light)','separator','on', ...
                'callback',@(hu,evnt)autorepair(V,false))
            uimenu(m,'label','Auto-repair (heavy)', ...
                'callback',@(hu,evnt)autorepair(V,true))
            m1 = uimenu(m,'label','More repairs...');
            uimenu(m1,'label','Clear memory', ...
                'callback',@(hu,evnt)clearmemory(V.content))
            uimenu(m1,'label','Reinit menus', ...
                'callback',@(hu,evnt)init_menus(V))
            uimenu(m1,'label','Reinit panels', ...
                'callback',@(hu,evnt)init_panels(V))
            uimenu(m1,'label','Reinit 4D', ...
                'callback',@(hu,evnt)tpview_reinit_4d(V))
            uimenu(m1,'label','Remove bad lines', ...
                'callback',@(hu,evnt)delete(findobj(V.grob.time,'tag','tpview-plot')))
            uimenu(m1,'label','Reinit trial slider', ...
                'callback',@(hu,evnt)init_slidertrial(V))
            uimenu(m1,'label','Remove watch', ...
                'callback',@(hu,evnt)set(V.hf,'pointer','arrow'))
            uimenu(m,'label','Report bug', ...
                'callback',@(hu,evnt)bugreport())

            % (edit code)
            if ismember('developer',V.modules)
                m1 = uimenu(m,'label','Edit code');
                uimenu(m1,'label','tpview','callback',@(hu,evnt)edit('tpview'))
                uimenu(m1,'label','tpv_content','callback',@(hu,evnt)edit('tpv_content'))
                uimenu(m1,'label','tps_trial','callback',@(hu,evnt)edit('tps_trial'))
                uimenu(m1,'label','tps_signal','callback',@(hu,evnt)edit('tps_signal'))
                uimenu(m1,'label','tps_signalx','callback',@(hu,evnt)edit('tps_signalx'))
                if ismember('acquisition',V.modules)
                    uimenu(m1,'label','tpv_imager','callback',@(hu,evnt)edit('tpv_imager'))
                end
            end
            % (access object)
            if ismember('developer',V.modules)
                uimenu(m,'label','access tpview object', ...
                    'callback',@(hu,evnt)access(V))
            end
            
            % general parameters
            items.tpview.seldotrial = uimenu(m,'label','Trial-specific selection','separator','on', ...
                'callback',@(hu,evnt)setpar(V,'seldotrial','toggle'), ...
                'checked',onoff(V.content.seldotrial));
            
            % read data
            % (remove first one: it is probably useless)
            %             items.tpview.loaddata = uimenu(m,'label','load data at opening','separator','on', ...
            %                 'callback',@(hu,evnt)setpar(V,'loaddata','toggle'), ...
            %                 'checked',onoff(V.internpar.loaddata));
            items.tpview.readdata = uimenu(m,'label','avoid loading data','separator','on', ...
                'callback',@(hu,evnt)setpar(V,'readdata','toggle'), ...
                'checked',onoff(~V.internpar.readdata));
            uimenu(m,'label','read data', ...
                'callback',@(hu,evnt)readalltrials('data'))
            uimenu(m,'label','compute dataop', ...
                'callback',@(hu,evnt)readalltrials('dataop'))
            uimenu(m,'label','compute signals', ...
                'callback',@(hu,evnt)readalltrials('signals'))
            uimenu(m,'label','load all trials', ...
                'callback',@(hu,evnt)readalltrials('all'))
            function readalltrials(flag)
                c = watch(V); %#ok<NASGU>
                switch flag
                    case 'data'
                        fn_progress('reading trial',V.ntrial)
                        for k=1:V.ntrial, fn_progress(k), V.content.trials(k).data; end
                    case 'dataop'
                        fn_progress('operation',V.ntrial)
                        for k=1:V.ntrial, fn_progress(k), operation(V.content.trials(k),'data',[],false); end % don't enlarge possibly binned corrected data
                    case 'signals'
                        fn_progress('compute signals',V.ntrial)
                        for k=1:V.ntrial, fn_progress(k), getslice(V.content,k); end
                    case 'all'
                        fn_progress('load trial',V.ntrial)
                        for k=1:V.ntrial, fn_progress(k), V.ktrial=k; pause(.1), end
                end
            end
            
            % preferences
            uimenu(m,'label','Preferences...','separator','on', ...
                'callback',@(u,e)preferences(V,'edit'))
            
            
            V.menus.items = items;
        end
        function menus_pictures(V)
            hf = V.grob.hf;
            items = V.menus.items;
            datamodelist = {'data','sfr','data_sfr','dataop','dataopmem','sfrop', ...
                'shotnoise','shotnoiseop','data0'}; %,'spikes','shotnoisespikes'};
            if ~ismember('sfr',V.modules)
                datamodelist(~fn_isemptyc(strfind(datamodelist,'sfr')))=[]; 
            end
            nlist = length(datamodelist);
            
            % loop on 2 menus
            for k=1:2
                % parent menu
                flag = ['im' num2str(k)];
                V.menus.(flag) = uimenu('parent',hf,'label',['Picture ' num2str(k)]);
                m = V.menus.(flag);
                % basic data modes
                for i=1:nlist
                    name = datamodelist{i};
                    label = strrep(name,'_','/');
                    items.(flag).(name) = uimenu(m,'label',label, ...
                        'callback',@(hu,evnt)setpar(V,flag,name));
                end
                set(items.(flag).(V.disppar.(flag)),'checked','on')
                %                 % combine with spikes movie
                %                 items.(flag).spikecomb = uimenu(m,'label','combine spikes','separator','on', ...
                %                     'callback',@(hu,evnt)setpar(V,[flag 'spikecomb'],'toggle'));
                % conditions
                items.(flag).condspec  = [];
                items.(flag).condobjs  = [];
            end
            
            V.menus.items = items;
        end
        function menus_timecourses(V)
            hf = V.grob.hf;
            items = V.menus.items;
            % menu
            V.menus.time = uimenu('parent',hf,'label','Time courses');
            m = V.menus.time;
            
            % datamode
            items.time.data = uimenu(m,'label','data', ...
                'callback',@(hu,evnt)setpar(V,'time','data'));
            items.time.dataop = uimenu(m,'label','data op.', ...
                'callback',@(hu,evnt)setpar(V,'time','dataop'));
            items.time.data_dataop = uimenu(m,'label','data/dataop', ...
                'callback',@(hu,evnt)setpar(V,'time','data_dataop'));
            items.time.dataopmem_dataop = uimenu(m,'label','dataopmem/dataop', ...
                'callback',@(hu,evnt)setpar(V,'time','dataopmem_dataop'));
            m1 = uimenu(m,'label','more');
            items.time.dataopmem = uimenu(m1,'label','data op. mem.', ...
                'callback',@(hu,evnt)setpar(V,'time','dataopmem'));
            if ismember('sfr',V.modules)
                items.time.sfr = uimenu(m1,'label','sfr', ...
                    'callback',@(hu,evnt)setpar(V,'time','sfr'));
                items.time.data_sfr = uimenu(m1,'label','data/sfr', ...
                    'callback',@(hu,evnt)setpar(V,'time','data_sfr'));
                items.time.sfrop = uimenu(m1,'label','sfr op.', ...
                    'callback',@(hu,evnt)setpar(V,'time','sfrop'));
                items.time.dataop_sfrop = uimenu(m1,'label','data/sfr op.', ...
                    'callback',@(hu,evnt)setpar(V,'time','dataop_sfrop'));
                items.time.data_shotnoise = uimenu(m1,'label','data/shot noise', ...
                    'callback',@(hu,evnt)setpar(V,'time','data_shotnoise'));
                items.time.dataop_shotnoiseop = uimenu(m1,'label','data/shot noise op.', ...
                    'callback',@(hu,evnt)setpar(V,'time','dataop_shotnoiseop'));
            end
            items.time.shotnoise = uimenu(m1,'label','shot noise', ...
                'callback',@(hu,evnt)setpar(V,'time','shotnoise'));
            items.time.shotnoiseop = uimenu(m1,'label','shot noise op.', ...
                'callback',@(hu,evnt)setpar(V,'time','shotnoiseop'));
            set(items.time.(V.disppar.time),'checked','on')
            
            % spikes display
            if ismember('spikes',V.modules) || ismember('electrophy',{V.disppar.timerecording.name})
                items.time.spikesignal = uimenu(m,'label','signal only','separator','on', ...
                    'callback',@(hu,evnt)setpar(V,'timespike','spikesignal'));
                items.time.spikeboth = uimenu(m,'label','signal+spike', ...
                    'callback',@(hu,evnt)setpar(V,'timespike','spikeboth'));
                items.time.spikespike = uimenu(m,'label','spikes only', ...
                    'callback',@(hu,evnt)setpar(V,'timespike','spikespike'));
            end
            if ismember('spikes',V.modules)
                items.time.spikefit = uimenu(m,'label','signal+spike+fit', ...
                    'callback',@(hu,evnt)setpar(V,'timespike','spikefit'));
                items.time.spikeneuropil = uimenu(m,'label','spikes / neuropil', ...
                    'callback',@(hu,evnt)setpar(V,'timespike','spikeneuropil'));
                set(items.time.(V.disppar.timespike),'checked','on')
            end
            
            % signal mode
            items.time.signal = uimenu(m,'label','signal','separator','on', ...
                'callback',@(hu,evnt)setpar(V,'timesignal','signal'));
            items.time.signalop = uimenu(m,'label','signal op.', ...
                'callback',@(hu,evnt)setpar(V,'timesignal','signalop'));
            items.time.signal_signalop = uimenu(m,'label','signal/signalop', ...
                'callback',@(hu,evnt)setpar(V,'timesignal','signal_signalop'));
            set(items.time.(V.disppar.timesignal),'checked','on')
            
            % recordings display
            recordingnames = lower({V.disppar.timerecording.name});
            if ismember('heart',recordingnames)
                items.time.recheart = uimenu(m,'label','heart','separator','on', ...
                    'callback',@(hu,evnt)setpar(V,'recheart','toggle'), ...
                    'checked',onoff(getrecordingdisp(V,'heart')));
            end
            if ismember('electrophysiology',recordingnames)
                items.time.recelectrophysiology = uimenu(m,'label','electrophysiology', ...
                    'callback',@(hu,evnt)setpar(V,'recelectrophysiology','toggle'), ...
                    'checked',onoff(getrecordingdisp(V,'electrophysiology')));
            end
            if ismember('eeg',recordingnames)
                items.time.receeg = uimenu(m,'label','EEG', ...
                    'callback',@(hu,evnt)setpar(V,'receeg','toggle'), ...
                    'checked',onoff(getrecordingdisp(V,'eeg')));
            end
            missingnames = setdiff(lower(recordingnames),{'heart','electrophysiology','eeg'});
            if ~isempty(missingnames)
                m1 = uimenu(m,'label','more');
            end
            %             items.time.recdata = uimenu(m1,'label','linear data', ...
            %                 'callback',@(hu,evnt)setpar(V,'recdata','toggle'), ...
            %                 'checked',onoff(getrecordingdisp(V,'data')));
            %             items.time.recsfr = uimenu(m1,'label','linear sfr', ...
            %                 'callback',@(hu,evnt)setpar(V,'recsfr','toggle'), ...
            %                 'checked',onoff(getrecordingdisp(V,'sfr')));
            for i=1:length(missingnames)
                recname = ['rec' missingnames{i}];
                items.time.(recname) = uimenu(m1,'label',missingnames{i}, ...
                    'callback',@(hu,evnt)setpar(V,recname,'toggle'), ...
                    'checked',onoff(getrecordingdisp(V,missingnames{i})));
            end
            % fft, linescan display, cells/trials
            items.time.fft = uimenu(m,'label','fft (compute)','separator','on', ...
                'callback',@(hu,evnt)setpar(V,'timefft','togglecompute'), ...
                'checked',fn_switch(V.disppar.timefft,'compute','on','off'));
            items.time.datafft = uimenu(m,'label','fft (from data)', ...
                'callback',@(hu,evnt)setpar(V,'timefft','togglefromdata'), ...
                'checked',fn_switch(V.disppar.timefft,'fromdata','on','off'));
            items.time.line = uimenu(m,'label','linescan', ...
                'callback',@(hu,evnt)setpar(V,'timeline','toggle'), ...
                'checked',onoff(V.content.timeline));
            
            % additional displays
            items.time.yzero = uimenu(m,'label','show zero/one','separator','on', ...
                'callback',@(hu,evnt)setpar(V,'timeyzero','toggle'), ...
                'checked',onoff(V.disppar.timeyzero));
            items.time.rms = uimenu(m,'label','show RMS', ...
                'callback',@(hu,evnt)setpar(V,'timerms','toggle'), ...
                'checked',onoff(V.disppar.timerms));
            items.time.rmsdcomp = uimenu(m,'label','compute RMS based on derivative only', ...
                'callback',@(hu,evnt)setpar(V,'timermsdcomp','toggle'), ...
                'checked',onoff(V.disppar.timermsdcomp));
            
            % regions/trials switch
            items.time.alltrials = uimenu(m,'label','one sel, all trials','separator','on', ...
                'callback',@(hu,evnt)setpar(V,'timealltrials','toggle'), ...
                'checked',onoff(V.disppar.timealltrials));
            items.time.alltrialsfilterstim = uimenu(m,'label','all trials: only with same stim as current', ...
                'callback',@(hu,evnt)setpar(V,'timealltrialsfilterstim','toggle'), ...
                'checked',onoff(V.disppar.timealltrialsfilterstim));
            
            % conditions
            items.time.condspec  = [];
            items.time.condobjs  = [];
            
            % finish
            V.menus.items = items;
        end
        function menus_moredisplays(V)
            hf = V.grob.hf;
            % additional displays
            V.menus.plus = uimenu('parent',hf,'label','More spatial displays');
            m = V.menus.plus;
            uimenu(m,'label','3D view of Picture 1', ...
                'callback',@(hu,evnt)display_more(V,'3d','im1'));
            uimenu(m,'label','3D view of Picture 2', ...
                'callback',@(hu,evnt)display_more(V,'3d','im2'));
            uimenu(m,'label','Frames view of Picture 1','separator','on', ...
                'callback',@(hu,evnt)display_more(V,'frames','im1'));
            uimenu(m,'label','Frames view of Picture 2', ...
                'callback',@(hu,evnt)display_more(V,'frames','im2'));
            uimenu(m,'label','Array view of Picture 1','separator','on', ...
                'callback',@(hu,evnt)display_more(V,'array','im1'));
            uimenu(m,'label','Array view of Picture 2', ...
                'callback',@(hu,evnt)display_more(V,'array','im2'));
            uimenu(m,'label','Movie of Picture 1','separator','on', ...
                'callback',@(hu,evnt)display_more(V,'movie','im1'));
            uimenu(m,'label','Movie of Picture 2', ...
                'callback',@(hu,evnt)display_more(V,'movie','im2'));
            uimenu(m,'label','Line cut from Picture 1','separator','on', ...
                'callback',@(hu,evnt)display_more(V,'linecut','im1'));
            uimenu(m,'label','Line cut from Picture 2', ...
                'callback',@(hu,evnt)display_more(V,'linecut','im2'));
            uimenu(m,'label','Full options...','separator','on', ...
                'callback',@(hu,evnt)display_more(V,'fulloptions'));
        end
        function menus_trials(V)
            hf = V.grob.hf;
            items = V.menus.items;
            V.menus.trials = uimenu('parent',hf,'label','Trials');
            m = V.menus.trials;
            
            
            % headers
            uimenu(m,'label','copy settings from previous experiment...', ...
                'callback',@(hu,evnt)data_copysettings(V))
            uimenu(m,'label','initial binning (all trials)','separator','on', ...
                'callback',@(hu,evnt)setbinning(V,true))
            uimenu(m,'label','initial binning (this trial only)', ...
                'callback',@(hu,evnt)setbinning(V,false))
            uimenu(m,'label','set pixel size and frame duration (all trials)', ...
                'callback',@(hu,evnt)setscales(V,true))
            uimenu(m,'label','set pixel size and frame duration (this trial only)', ...
                'callback',@(hu,evnt)setscales(V,false))
            % functions
            function setbinning(V,doalltrials)
                % check
                if ~V.content.seldotrial && ~doalltrials
                    errordlg 'selection change in specific trials contradicts identical selections mode'
                    return
                end
                % user input
                s0 = struct('xbin',V.trial.xbin,'tbin',V.trial.tbin);
                spec = struct('xbin','stepper 1 1 Inf','tbin','stepper 1 1 Inf');
                s = fn_structedit(s0,spec);
                if isempty(s), return, end % figure closed, action cancelled
                % set trial properties
                ktrials = fn_switch(doalltrials,1:V.ntrial,V.ktrial);
                T = V.content.trials(ktrials);
                chgx = ([T.xbin]~=s.xbin);
                xfactor = s.xbin ./ [T.xbin];
                [T.xbin] = deal(s.xbin);
                [T.tbin] = deal(s.tbin);
                % update selection
                erasedata(V.content,ktrials(chgx))
                if V.content.seldotrial
                    for i=find(chgx)
                        updateselectionhandy(V.content,ktrials(i),'xbinchange',xfactor(i))
                    end
                else
                    if ~isscalar(T) && any(diff(xfactor)), error programming, end
                    updateselectionhandy(V.content,ktrials,'xbinchange',xfactor(1))
                end
                % update display
                display_changeview(V,'data&signals',false) % this will update the display correclty
            end
            function setscales(V,doalltrials)
                s = struct('dx',V.dx,'dy',V.dy,'dz',V.dz,'spatial__unit',V.trial.xunit, ...
                    'frame__duration',V.dt,'temporal__unit',V.trial.tunit);
                spec = struct('dx','double','dy','double','dz','double','spatial__unit','char', ...
                    'frame__duration','double','temporal__unit','char');
                s = fn_structedit(s,spec);
                if isempty(s), return, end % figure closed, action cancelled
                ktrials = fn_switch(doalltrials,1:V.ntrial,V.ktrial);
                T = V.content.trials(ktrials);
                [T.dx] = deal(s.dx);
                [T.dy] = deal(s.dy);
                [T.dz] = deal(s.dz);
                [T.xunit] = deal(s.spatial__unit);
                [T.dt] = deal(s.frame__duration);
                [T.tunit] = deal(s.temporal__unit);
                erasedata(V.content)
                display_changeview(V,'ktrial') % this will update the display correclty
            end
            
            %  add / reorder / remove / filter
            uimenu(m,'label','add average trial...','separator','on', ...
                'callback',@(hu,evnt)data_addtrial(V,'average'));
            uimenu(m,'label','add trial from base workspace...', ...
                'callback',@(hu,evnt)data_addtrial(V,'command'));
            m1 = uimenu(m,'label','reorder trials');
            uimenu(m1,'label','based on experiment and block numbers', ...
                'callback',@(hu,evnt)reordertrials('blocknumber'))
            uimenu(m1,'label','based on file times', ...
                'callback',@(hu,evnt)reordertrials('filetime'))
            uimenu(m1,'label','manual permutation', ...
                'callback',@(hu,evnt)reordertrials('manualperm'))
            uimenu(m,'label','accept/reject trial','accelerator','R', ...
                'callback',@(hu,evnt)action(V,'rejecttrialtoggle'));
            uimenu(m,'label','remove current trial', ...
                'callback',@(hu,evnt)removetrial(false));
            uimenu(m,'label','remove trial + file', ...
                'callback',@(hu,evnt)removetrial(true));
            uimenu(m,'label','remove trials...', ...
                'callback',@(hu,evnt)removetrial('custom'));
            uimenu(m,'label','filter trials...', ...
                'callback',@(hu,evnt)display_trialsublist(V,'set'));
            uimenu(m,'label','unfilter trials', ...
                'callback',@(hu,evnt)display_trialsublist(V,'reset'));
            % functions
            function removetrial(flag)
                if strcmp(flag,'custom')
                    s = fn_structedit('trial__number',1,'remove__files',false);
                    if isempty(s), return, end
                    ktrials = s.trial__number;
                    dormfile = s.remove__files;
                else
                    ktrials = V.content.ktrial;
                    dormfile = flag;
                end
                if dormfile
                    files = {};
                    for s=ktrials
                        [p base] = fileparts(V.content.trials(s).file);
                        if isempty(base), continue, end % no file for this trial
                        tmp = dir([p '/' base '*']);
                        files = [files cellstr([repmat([p '/'],length(tmp),1) char(tmp.name)])]; %#ok<AGROW>
                    end
                    if ~isempty(files)
                        confirmation = questdlg(['Are you sure to delete files?',files], ...
                            'Confirm file deletion','Yes','No','No');
                        if strcmp(confirmation,'No'), return, end
                        for k=1:length(files), delete(files{k}), end
                    end
                end
                rmtrial(V.content,ktrials)
                V.disppar.trialsublist = []; % list became obsolete because trials have been assigned new numbers
                display_changeview(V,'chgtrial')
            end
            function reordertrials(flag)
                T = V.content.trials;
                names = {T.file};
                for k=1:V.ntrial
                    if isempty(names{k})
                        if isfield(T(k).addinfo,'description')
                            names{k} = T(k).addinfo.description;
                        else
                            names{k} = T(k).origin;
                        end
                    end
                end
                switch flag
                    case 'blocknumber'
                        EB = V.ntrial*ones(V.ntrial,2);
                        for k=1:V.ntrial
                            token = regexp(T(k).file,'E(\d*)B(\d*).BLK','tokens');
                            if ~isempty(token)
                                EB(k,:) = [str2double(token{1}{1}) str2double(token{1}{2})];
                            end
                        end
                        [dum perm] = sortrows(EB); %#ok<ASGLU>
                    case 'filetime'
                        datenum = zeros(1,V.ntrial);
                        for k=1:V.ntrial
                            datenum(k) = getfield(dir(T(k).file),'datenum');
                        end
                        [dum perm] = sort(datenum); %#ok<ASGLU>
                    case 'manualperm'
                        if V.ntrial<10
                            def = num2str(1:V.ntrial,'%i ');
                        else
                            def = [num2str(1:9,'%i ') ' 10:' num2str(V.ntrial)];
                        end
                        answer = inputdlg('Permutation:','user input',1,{def});
                        if isempty(answer), return, end
                        str = answer{1};
                        try perm = evalin('base',['[' str ']']); catch, errordlg('could not evaluate string'), return, end %#ok<CTCH>
                        if ~isequal(row(sort(perm)),1:V.ntrial), errordlg 'this is not a permutation', return, end
                end
                permutetrials(V.content,perm)
                % update display
                V.disppar.trialsublist = [];
                display_changeview(V,'ktrial')
            end
            
            % stims
            uimenu(m,'label','edit stims...','separator','on', ...
                'callback',@(hu,evnt)editstims())
            m1 = uimenu(m,'label','set stims');
            uimenu(m1,'label','set stims...', ...
                'callback',@(hu,evnt)markstim('input'));
            uimenu(m1,'label','unmark stims', ...
                'callback',@(hu,evnt)markstim);
            uimenu(m,'label','save stims to file...', ...
                'callback',@(hu,evnt)savestims)
            uimenu(m,'label','load stims from file...', ...
                'callback',@(hu,evnt)loadstims)
            items.showstim = uimenu(m,'label','show stim','checked',onoff(V.disppar.showstim), ...
                'callback',@(hu,evnt)setpar(V,'showstim','toggle'));
            % functions
            function markstim(delay,len,freq,dur)
                if nargin==0
                    stim = [];
                    doalltrial = false;
                elseif strcmp(delay,'input')
                    prompt = {'do for all trials?' 'delay(s)','plateau(ms)','interval(ms)','number'};
                    default = {'0','11','50','200','10'};
                    answer = inputdlg(prompt,'Enter stimulation parameters',1,default);
                    if isempty(answer), return, end
                    doalltrial  = str2double(answer{1});
                    delay       = str2num(answer{2});
                    plateau     = str2num(answer{3});
                    interval    = str2double(answer{4});
                    number      = str2double(answer{5});
                    if isscalar(delay)
                        % train at fix intervals and durations
                        stim = delay + (0:number-1)*(interval/1000);
                        stim(2,:) = (plateau/1000);
                    else
                        % custom starts and durations
                        if isscalar(plateau)
                            plateau = repmat(plateau,1,length(delay));
                        elseif length(delay)~=length(plateau)
                            error('''delay'' and ''plateau'' should be the same length')
                        end
                        stim = [delay; plateau/1000];
                    end
                else
                    stim = delay + (0:len*freq-1)/freq;
                    stim(2,:) = dur;
                    doalltrial = false;
                end
                if doalltrial, ind = 1:V.ntrial; else ind = V.ktrial; end
                for kk=ind
                    V.content.trials(kk).stim = stim;
                end
                display_stimandevents(V)
            end
            function savestims()
                fname = fn_savefile('*.stim','Select file where to select stim definitions');
                if ~fname, return, end
                fname = [fn_fileparts(fname,'noext') '.stim'];
                table = V.content.stimtable.table; %#ok<NASGU>
                save(fname,'table','-MAT')
            end
            function loadstims()
                fname = fn_getfile('*.stim','Select stim definition file');
                if ~fname, return, end
                load(fname,'table','-MAT')
                updatetable(V.content.stimtable,table);
                % update display
                display_changeview(V,'ktrial')
            end
            function editstims()
                hf = msgbox('Edit variable ''STIMTABLE'', then press ''Ok''. Close this figure to Cancel');
                hu = findobj(hf,'type','uicontrol');
                set(hu,'callback',@(u,e)editstimsub(hf))
                stimtable = V.content.stimtable.table;
                F = fieldnames(stimtable);
                ST = [F';  permute(struct2cell(stimtable),[3 1 2])];
                assignin('base','STIMTABLE',ST)
                openvar('STIMTABLE')
            end
            function editstimsub(hf)
                delete(hf)
                ST = evalin('base','STIMTABLE');
                F = ST(1,:);
                ST = ST(2:end,:); if isempty(ST), return, end
                tmp = [F; num2cell(ST,1)];
                ST = struct(tmp{:});
                updatetable(V.content.stimtable,ST)
                file_markchange(V);
                % update display
                display_changeview(V,'ktrial')
            end
            
            % events
            items.showevents = uimenu(m,'label','display events','separator','on', ...
                'callback',@(hu,evnt)setpar(V,'showevents','toggle'), ...
                'checked',onoff(V.disppar.showevents));
            items.editevents = uimenu(m,'label','edit events', ...
                'callback',@(hu,evnt)setpar(V,'editevents','toggle'), ...
                'checked',onoff(V.disppar.editevents));
            
            % analog recordings
            uimenu(m,'label','attach analog recordings','separator','on', ...
                'callback',@(hu,evnt)data_attachanalog(V));
            
            %             % heart
            %             uimenu(m,'label','read heart (all trials)','separator','on', ...
            %                 'callback',@(hu,evnt)setheart(V.content.trials));
            %             uimenu(m,'label','read heart (this trial)', ...
            %                 'callback',@(hu,evnt)setheart(V.content.trial));
            
            % modify trials
            uimenu(m,'label','Make new set with dataop as data','separator','on', ...
                'callback',@(hu,evnt)data_savenewdata(V,'dataop'));
            uimenu(m,'label','Make new set with binned data', ...
                'callback',@(hu,evnt)data_savenewdata(V,'bin'));
            uimenu(m,'label','Make new set with current trial only', ...
                'callback',@(hu,evnt)savecurrenttrialonly());
            uimenu(m,'label','Coregister trials', ...
                'callback',@(hu,evnt)data_savenewdata(V,'coregister'));
            % code
            function savecurrenttrialonly()
                rmi = setdiff(1:V.ntrial,V.ktrial);
                rmtrial(V.content,rmi)
                display_changeview(V,'chgtrial')
                file_save(V,[])
            end
            
            % store items
            V.menus.items = items;
        end
        function menus_electrophys(V)
            if ~ismember('electrophy',{V.disppar.timerecording.name})
                return
            end
            hf = V.grob.hf;
            % electrophys
            V.menus.electrophys = uimenu('parent',hf,'label','Electrophys');
            m = V.menus.electrophys;
            uimenu(m,'label','set electrophysiology signal', ...
                'callback',@(hu,evnt)data_setelectrophy(V))
            uimenu(m,'label','detect true spikes', ...
                'callback',@(hu,evnt)guessspikes(V.content.electrophys,V.ktrial))
        end
        function menus_spikes(V)
            hf = V.grob.hf;
            items = V.menus.items;
            % spikes reconstruction
            V.menus.signals = uimenu('parent',hf,'label','Spikes');
            m = V.menus.signals;
            uimenu(m,'label','compute spikes for this display','accelerator','X', ...
                'callback',@(u,e)data_spikes(V,'current',[],true));
            uimenu(m,'label','compute all spikes','accelerator','A', ...
                'callback',@(u,e)data_spikes(V,'all',[],true));
            
            m1 = uimenu(m,'label','calcium spikes parameters','separator','on');
            uimenu(m1,'label','\default', ...
                'callback',@(u,e)data_spikes(V,'all','calciumdefault',false));
            uimenu(m1,'label','drift', ...
                'callback',@(u,e)data_spikes(V,'all','calciumdrift',false));
            uimenu(m1,'label','custom...', ...
                'callback',@(u,e)data_spikes(V,'all','calciumcustom',false));
            uimenu(m1,'label','active control...', ...
                'callback',@(u,e)data_spikecontrol(V,'calcium'));
            %             uimenu(m,'label','vsd spikes...', ...
            %                 'callback',@(u,e)data_spikecontrol(V,'vsd'));
            
            items.spikes.guessspikes = uimenu(m,'label','auto-spikes','separator','on', ...
                'callback',@(u,e)setpar(V,'guessspikes','toggle'), ...
                'checked',onoff(V.internpar.guessspikes));
            items.spikes.guessspikesall = uimenu(m,'label','auto-spikes all', ...
                'callback',@(u,e)setpar(V,'guessspikesall','toggle'), ...
                'checked',onoff(V.internpar.guessspikesall));
            uimenu(m,'label','(re)compute spikes for this display - no drift','separator','on', ...
                'callback',@(u,e)data_spikes(V,'current','calciumdefault',true));
            uimenu(m,'label','(re)compute spikes for this display - include a drift', ...
                'callback',@(u,e)data_spikes(V,'current','calciumdrift',true));
            uimenu(m,'label','(re)compute spikes for this display - custom parameters', ...
                'callback',@(u,e)data_spikes(V,'current','calciumcustom',true));
            uimenu(m,'label','no spikes for this display','separator','on', ...
                'callback',@(u,e)data_spikes(V,'current','nospike',true))
            uimenu(m,'label','no spikes for this region', ...
                'callback',@(u,e)data_spikes(V,'currentregion','nospike',true))
            uimenu(m,'label','remove invalid spikes', ...
                'callback',@(u,e)data_spikes(V,'all','rminvalid',true))
            uimenu(m,'label','remove all spikes', ...
                'callback',@(u,e)data_spikes(V,'all','rmall',true))
            
            % items to keep track of
            V.menus.items = items;
        end
        function menus_scripts(V)
            % input: package directory (contains user routines)
            
            % default package: 2p analysis
            scriptfolder = fullfile(fileparts(which('optimage')),'tpview_scripts');
            
            % init menu
            hf = V.grob.hf;
            if ~isfield(V.menus,'analysis') || ~any(ishandle(V.menus.analysis))
                V.menus.analysis = uimenu('parent',hf,'label','Scripts');
            end
            m = V.menus.analysis;
            delete(get(m,'children')) % needed when re-scanning
            
            % scan folder and add menus (both for run and edit code)
            medit = uimenu(m,'label','Edit code','separator','on');
            funcol = [.25 .25 .75];
            scanfolder(scriptfolder,m,medit)
            function scanfolder(folder,m1,medit1)
                swd = pwd;
                % sub-folders
                d = dir(folder); d(1:2)=[];
                subfolders = {d([d.isdir]).name};
                for i=1:length(subfolders)
                    m2 = uimenu(m1,'label',subfolders{i});
                    medit2 = uimenu(medit1,'label',subfolders{i});
                    scanfolder(fullfile(folder,subfolders{i}),m2,medit2)
                end
                % m-files
                cd(folder)
                d = dir('*.m');
                mfiles = {d.name};
                for i=1:length(mfiles)
                    fun = mfiles{i}(1:end-2); % remove the '.m' extension
                    if fun(1)=='.'
                        % On MAC all files have a small companion whose
                        % name starts with a '.'
                        continue
                    end
                    try
                        fun = str2func(fun);
                        label = fun('label');
                    catch
                        continue
                    end
                    uimenu(m1,'label',label,'foregroundcolor',funcol, ...
                        'callback',@(u,e)runfile(folder,fun))
                    uimenu(medit1,'label',label,'foregroundcolor',funcol, ...
                        'callback',@(u,e)edit(fullfile(folder,mfiles{i})))
                end
                % go back to previous directory
                cd(swd)
            end
            function runfile(folder,mfile)
                cd(folder)
                feval(mfile,V)
            end
            
            % edit
            uistack(medit,'top') % bring edit button down
            uimenu(m,'label','New script...','callback',@(u,e)newscript())
            function newscript()
                % prompt for label and file name
                scriptname = inputdlg({'Script label' 'file name'},V.skin,1,{'my script' 'script_myscript'});
                if isempty(scriptname), return, end
                label = scriptname{1};
                fscript = scriptname{2}; 
                % modify template
                ftxt = which('script_template.txt');
                txt = fn_readtext(ftxt);
                txt = strrep(txt,'WRITE_SCRIPT_NAME_HERE',fscript);
                txt = strrep(txt,'WRITE_SCRIPT_LABEL_HERE',label);
                % save
                fscript = fullfile(scriptfolder,[fn_fileparts(fscript,'base') '.m']);
                fn_savetext(txt,fscript)
                % re-scan folder
                menus_scripts(V)
                % edit
                edit(fscript)
            end
            
            % organize scripts
            if ispc
                uimenu(m,'label','Organize scripts (in Explorer)','separator','on', ...
                    'callback',['!explorer ' scriptfolder])
            end
            uimenu(m,'label','Re-scan scripts folder','separator',onoff(~ispc), ...
                'callback',@(u,e)menus_scripts(V))
        end
        function menus_help(V)
            % init menu
            hf = V.grob.hf;
            if ~isfield(V.menus,'help') || ~any(ishandle(V.menus.help))
                V.menus.help = uimenu('parent',hf,'label','Help');
            end
            m = V.menus.help;
            delete(get(m,'children')) % needed when re-scanning
            
            doc_file = fullfile( ...
                fileparts(fileparts(which('optimage'))), ...
                'documentation', 'Optimage User Manual.pdf');
            mouse_file = fullfile( ...
                fileparts(fileparts(which('optimage'))), ...
                'documentation', 'optimage mouse actions.pdf');

            function open_pdf(file)
                system(['start "" "' file '"']);
            end
            
            uimenu(m,'label','OptImage doc', ...
                'callback',@(u,e)open_pdf(doc_file))
            
            uimenu(m,'label','Mouse actions', ...
                'callback',@(u,e)open_pdf(mouse_file))
        end
    end
    
    % GET
    methods
        % (trial current and number)
        function ntrial = get.ntrial(V)
            ntrial = length(V.content.trials);
        end
        function ktrial = get.ktrial(V)
            ktrial = V.content.ktrial;
        end
        function trial = get.trial(V)
            trial = V.content.trials(V.content.ktrial);
        end
        % (trial properties)
        function data = get.data(V)
            data = V.content.trials(V.content.ktrial).data;
        end
        function sfr = get.sfr(V)
            sfr = V.content.trials(V.content.ktrial).sfr;
        end
        function file = get.file(V)
            file = V.content.trials(V.content.ktrial).file;
        end
        function fullinfo = get.fullinfo(V)
            fullinfo = V.content.trials(V.content.ktrial).fullinfo;
        end
        function addinfo = get.addinfo(V)
            addinfo = V.content.trials(V.content.ktrial).addinfo;
        end
        function type = get.type(V)
            type = V.content.trials(V.content.ktrial).type;
        end
        function scanning = get.scanning(V)
            scanning = V.content.trials(V.content.ktrial).scanning;
        end
        function sizes = get.sizes(V)
            sizes = V.content.trials(V.content.ktrial).sizes;
        end
        function nx = get.nx(V)
            nx = V.content.trials(V.content.ktrial).nx;
        end
        function ny = get.ny(V)
            ny = V.content.trials(V.content.ktrial).ny;
        end
        function nfr = get.nfr(V)
            nfr = V.content.trials(V.content.ktrial).nfr;
        end
        function nc = get.nc(V)
            nc = V.content.trials(V.content.ktrial).nc;
        end
        function dx = get.dx(V)
            dx = V.content.trials(V.content.ktrial).dx;
        end
        function dy = get.dy(V)
            dy = V.content.trials(V.content.ktrial).dy;
        end
        function dz = get.dz(V)
            dz = V.content.trials(V.content.ktrial).dz;
        end
        function dt = get.dt(V)
            dt = V.content.trials(V.content.ktrial).dt;
        end
        function linedur = get.linedur(V)
            linedur = V.content.trials(V.content.ktrial).linedur;
        end
        function t0 = get.t0(V)
            t0 = V.content.trials(V.content.ktrial).t0;
        end
        function tidx = get.tidx(V)
            tidx = V.content.trials(V.content.ktrial).tidx;
        end
        function dataop = get.dataop(V)
            dataop = V.content.trials(V.content.ktrial).dataop;
        end
        function nsel = get.nsel(V)
            nsel = V.content.nsel;
        end
        % (signals)
        function signal = get.signal(V)
            % get the current signals: needs to compute the data
            signal = V.content.signal;
        end
        % (recording)
        function x = getrecordingdisp(V,flag)
            recidx = find(strcmpi(flag,{V.disppar.timerecording.name}));
            if isempty(recidx)
                error('unknown recording name ''%s''',flag)
            else
                x = V.disppar.timerecording(recidx).val;
            end
        end
    end
    
    % SET
    methods
        % (content)
        function set.ktrial(V,k)
            if k<1 || abs(k-round(k))>1e-4 || k>V.ntrial, error('ktrial must be an integer between 1 and the number of trials'), end
            k = round(k);
            ktrialold = V.content.ktrial;
            
            % simplified procedure if slider is scrolling
            if V.a4d.slidertrial.sliderscrolling
                display_header(V,k)
                return
            end
            
            if (ktrialold==k), return, end
            
            % change value to be inside subset of trials
            sublist = V.disppar.trialsublist;
            if ~isempty(sublist) && ~ismember(k,sublist)
                disp('attempt to change ktrial to a value outside the sublist')
                idx = find(sublist==ktrialold);
                idxnew = fn_switch(k>ktrialold,idx+1,idx-1);
                idxnew = max(1,min(length(sublist),idxnew));
                k = sublist(idxnew);
            end
            % set value (automatic update of signals)
            V.content.ktrial = k;
            % display (also updates the trial slider)
            display_changeview(V,'ktrial')
        end
        function set.addinfo(V,x)
            V.content.trials(V.content.ktrial).addinfo = x;
        end
        function set.type(V,x)
            V.content.trials(V.content.ktrial).type = x;
        end
        function set.scanning(V,x)
            V.content.trials(V.content.ktrial).scanning = x;
        end
        function set.dx(V,x)
            V.content.trials(V.content.ktrial).dx = x;
        end
        function set.dy(V,x)
            V.content.trials(V.content.ktrial).dy = x;
        end
        function set.dz(V,x)
            V.content.trials(V.content.ktrial).dz = x;
        end
        function set.dt(V,x)
            V.content.trials(V.content.ktrial).dt = x;
        end
        function set.t0(V,x)
            V.content.trials(V.content.ktrial).t0 = x;
        end
        % (display parameters)
        function setpar(V,flag,val)
            % function setpar(V,option,value)
            % function setpar(V)
            %---
            % change display parameter value, update the control
            % controling this value, and execute the appropriate action
            %
            % for logical values, value='toggle' result in switching the
            % parameter value from true to false and vice-versa
            %
            % type setpar(V) to get a list of available options
            if nargin<3, setparinfo(V), end
            
            % set new value
            if regexp(flag,'^rec'), recname = flag(4:end); flag = 'rec'; end
            switch flag
                case {'timeline' 'seldotrial'}
                    % in content
                    oldval = V.content.(flag);
                    if strcmp(val,'toggle'), val=~oldval; else val=logical(val); end
                    if val==oldval, return, end
                    if strcmp(flag,'seldotrial') && ~val && V.ntrial>1
                        % canceling the 'trial-specific selections' option
                        % requires check or at least user confirmation
                        linescans = strcmp({V.content.trials.type},'linescan');
                        if any(linescans) && ~all(linescans)
                            errordlg('Selections cannot be identical in all trials because some trials are images and some are line scans')
                            return
                        else
                            answer = questdlg('Are you sure you want to set all trial selections identic to first trial?', ...
                                'Confirmation required','Yes','No','No');
                            if ~strcmp(answer,'Yes'), return, end
                        end
                    end
                    V.content.(flag) = val;
                case {'loaddata' 'readdata' 'guessspikes' 'guessspikesall'}
                    % in internpar
                    oldval = V.internpar.(flag);
                    if strcmp(val,'toggle'), val=~oldval; else val=logical(val); end
                    if val==oldval, return, end
                    V.internpar.(flag) = val;
                case 'rec'
                    % in disppar
                    recidx = find(strcmpi({V.disppar.timerecording.name},recname));
                    if ~isscalar(recidx), error 'could not find recording', end
                    oldval = V.disppar.timerecording(recidx).val;
                    if strcmp(val,'toggle'), val=~oldval; else val=logical(val); end
                    if val==oldval, return, end
                    V.disppar.timerecording(recidx).val = val;
                case 'timefft'
                    % in disppar
                    % val can be '', 'compute', 'fromdata', 'togglecompute'
                    % or 'togglefromdata'
                    oldval = V.disppar.timefft;
                    if strcmp(oldval,val), return, end
                    if strfind(val,'toggle')
                        newval = strrep(val,'toggle','');
                        if strcmp(oldval,newval)
                            % set newval=oldval 'off'
                            V.disppar.timefft = '';
                        else
                            % set newval 'on'!
                            V.disppar.timefft = newval;
                        end
                        val = V.disppar.timefft;
                    end
                otherwise
                    % in disppar
                    oldval = V.disppar.(flag);
                    if islogical(oldval)
                        if strcmp(val,'toggle'), val=~oldval; else val=logical(val); end
                    end
                    if isequal(val,oldval), return, end
                    V.disppar.(flag) = val;
            end
            
            % update control display and execute action
            switch flag
                % (in 'content')
                case 'timeline'
                    % menu mark
                    set(V.menus.items.time.line,'checked',onoff(val))
                    % action
                    display_changeview(V,'timeline')
                case 'seldotrial'
                    % menu mark
                    set(V.menus.items.tpview.seldotrial,'checked',onoff(val))
                    % (in 'internpar')
                case 'loaddata'
                    % menu mark
                    set(V.menus.items.tpview.(flag),'checked',onoff(val))
                    % ensure loaddata => readdata
                    if V.internpar.loaddata && ~V.internpar.readdata
                        V.internpar.readdata = true;
                        set(V.menus.items.tpview.readdata,'checked','off')
                    end
                case 'readdata'
                    % menu mark
                    set(V.menus.items.tpview.(flag),'checked',onoff(~val))
                    % ensure loaddata => readdata
                    if V.internpar.loaddata && ~V.internpar.readdata
                        V.internpar.loaddata = false;
                        set(V.menus.items.tpview.loaddata,'checked','off')
                    end
                    % update display
                    if V.internpar.readdata
                        display_changeview(V,'ktrial')
                    end
                case 'guessspikes'
                    % add spike display if necessary
                    if val && strcmp(V.disppar.timespike,'spikesignal')
                        V.disppar.timespike = 'spikeboth';
                        set(V.menus.items.time.spikesignal,'checked','off')
                        set(V.menus.items.time.spikeboth,'checked','on')
                    end
                    % menu mark
                    set(V.menus.items.spikes.guessspikes,'checked',onoff(val))
                    % action
                    if val, display_changeview(V,'timesignal'), end
                case 'guessspikesall'
                    % menu mark
                    set(V.menus.items.spikes.guessspikesall,'checked',onoff(val))
                    % action: background compute all spikes
                    if val
                        computespikes(V.content,1:V.ntrial,1:V.nsel)
                    end
                    % (in 'disppar')
                case {'im1','im2','time'}
                    % menu marks
                    set(V.menus.items.(flag).(oldval),'checked','off')
                    set(V.menus.items.(flag).(val)   ,'checked','on')
                    % action
                    data_chgmode(V,flag,oldval,val)
                case {'im1spikecomb','im2spikecomb'}
                    inflag = flag(1:3);
                    % menu mark
                    set(V.menus.items.(inflag).spikecomb,'checked',onoff(val))
                    % action
                    display_changeview(V,'datamode',str2double(inflag(3)),false,false)
                case {'timesignal' 'timespike'}
                    % menu marks
                    set(V.menus.items.time.(oldval),'checked','off')
                    set(V.menus.items.time.(val)   ,'checked','on')
                    % action
                    display_changeview(V,'timesignal')
                case 'rec'
                    % menu mark
                    set(V.menus.items.time.(['rec' recname]),'checked',onoff(val))
                    % action
                    display_changeview(V,'timesignal')
                case {'timeyzero' 'timerms' 'timermsdcomp' ...
                        'timealltrialsfilterstim'}
                    % menu marks
                    set(V.menus.items.time.(flag(5:end)),'checked',onoff(val))
                    % action
                    display_changeview(V,'timesignal')
                case {'im1cond','im2cond','timecond'}
                    % menu marks + action
                    inflag = strrep(flag,'cond','');
                    data_chgcondition(V,inflag,oldval,val)
                case 'timefft'
                    % menu mark
                    set(V.menus.items.time.fft,'checked',fn_switch(val,'compute','on','off'))
                    set(V.menus.items.time.datafft,'checked',fn_switch(val,'fromdata','on','off'))
                    % action
                    display_changeview(V,'timefft')
                case 'timealltrials'
                    % menu mark
                    set(V.menus.items.time.(flag(5:end)),'checked',onoff(val))
                    % text mark (for timealltrials)
                    if val
                        set(V.grob.timealltrials,'string','TRIALS', ...
                            'enable','inactive', ...
                            'buttonDownFcn',@(u,e)setpar(V,'currentsel',1+mod(V.disppar.currentsel,V.content.nsel)))
                    else
                        set(V.grob.timealltrials,'string','REGIONS','buttonDownFcn','');
                    end
                    % action
                    display_changeview(V,'timealltrials')
                case 'showstim'
                    % menu mark
                    set(V.menus.items.showstim,'checked',onoff(val))
                    % action
                    display_stimandevents(V)
                case {'showevents' 'editevents'}
                    % more flag changes
                    if V.disppar.editevents && ~V.disppar.showevents
                        switch flag
                            case 'showevents'
                                V.disppar.editevents = false;
                            case 'editevents'
                                V.disppar.showevents = true;
                        end
                    end
                    % menu mark
                    set(V.menus.items.showevents,'checked',onoff(V.disppar.showevents))
                    set(V.menus.items.editevents,'checked',onoff(V.disppar.editevents))
                    % action
                    if V.disppar.editevents
                        set(V.a4d.time,'usercallback',@(D)display_editevents(V))
                    else
                        set(V.a4d.time,'usercallback',[])
                    end
                    display_stimandevents(V)
                case 'currentsel'
                    if val<=0 || val>V.content.nsel || mod(val,1), error('wrong value'), end
                    % this is equivalent to clicking a region, and will
                    % trigger all the desired display updates, as well as
                    % setting V.disppar.currentsel (in display_changeview)
                    updateselection(V.a4d.G,[1 2],'active',val,true)
                case 'filelistgroup'
                    % control display
                    if val, str = 'Ungroup'; else str = 'Group'; end
                    set(V.panels.filesitems.group,'string',str)
                    % action
                    file_list(V)
                otherwise
                    error('wrong display option ''%s''',flag)
            end
        end
        % (indication to user)
        function x = setinfo(V) %#ok<MANU>
            x.ktrial = 'integer';
            x.addinfo = 'structure';
            x.type = {'movie','zstack'};
            x.scanning = {'0','1'};
            x.dx = [];
            x.dy = [];
            x.dz = [];
            x.dt = [];
            x.t0 = [];
        end
        function x = setparinfo(V)
            x.timeline = {'0' '1'};
            x.seldotrial = {'0','1'};
            x.im1 = {'data' 'dataop' 'sfr' 'sfrop' 'shotnoise' 'shotnoiseop'};
            x.im1cond = 'ex: C1/C0';
            x.im2 = {'data' 'dataop' 'sfr' 'sfrop' 'shotnoise' 'shotnoiseop'};
            x.im2cond = 'ex: C1/C0';
            x.time = {'data' 'dataop' 'sfr' 'sfrop' 'shotnoise' 'shotnoiseop' ...
                'dataop_sfrop' 'etc...'};
            x.timecond = 'ex: C1/C0';
            x.timesignal = {'signal' 'signalop' 'signal_signalop'};
            x.timefft = {'' 'compute' 'fromdata'};
            x.timeheart = {'0' '1'};
            x.timeelectrophys = {'0' '1'};
            x.timelineardata = {'0' '1'};
            x.timelinearsfr = {'0' '1'};
            x.timealltrials = {'0' '1'};
            x.LS = [];
            x.HS = [];
            x.timethr = [];
            x.filelistgroup = {'0','1'};
            if nargout==0
                disp(x)
                clear x
            end
        end
    end
    
    % DATA
    methods (Access='private')
        function data_loadfile(V,fname,flag)
            % function data_loadfile(V,fname[,'add|cat'])
            % function data_loadfile(V,data|content)
            if nargin<3, flag=''; end
            if ~fn_ismemberstr(flag,{'','bin','add','avg','cat',...
                    'bin&add','bin&avg','bin&cat','avg&add','bin&avg&add'})
                error('unknown flag ''%''',flag)
            end
            hf = V.grob.hf;
            
            % wait
            c = fn_watch(hf); %#ok<NASGU>
            
            % throw an event for closing file
            notify(V,'EventCloseFile')
            
            % flag
            dobin = any(strfind(flag,'bin'));
            flag = fn_strrep(flag,'bin&','','bin','');
            doadd = any(strfind(flag,'add'));
            flag = fn_strrep(flag,'&add','','add','');
            
            % load file
            okcontentdef = false; oksaved = false; okspecial = false;
            if ischar(fname)
                [p b ext] = fileparts(deblank(fname(1,:)));
                ext = lower(ext(2:end));
                if strfind(ext,'blk'), ext = 'blk'; end
            else
                ext = 'not a file name';
            end
            if strcmp(ext,'tpv')
                % Full data (trials + signals)
                v = load(fname,'-MAT');
                if isfield(v,'s')
                    % old version -> trials in .tptrial file and additional data in variable 's'
                    ftrial = fullfile(p,[b '.tptrial']);
                    T = tps_trial.readheader(ftrial);
                    if V.internpar.loaddata, readata(T), readrecording(T), end
                elseif isfield(v,'version')
                    % also an old version
                    V.content.version = v.version;
                    V.content.trials = v.trials;
                    V.content.signals = v.signals;
                    V.content.ktrial = v.ktrial;
                    if V.internpar.loaddata
                        readdata(V.content.trials)
                        readrecording(V.content.trials)
                    end
                    okcontentdef = true;
                elseif isscalar(fieldnames(v))
                    F = fieldnames(v); f = F{1};
                    V.content = v.(f);
                    if V.internpar.loaddata
                        readdata(V.content.trials)
                        readrecording(V.content.trials)
                    end
                    okcontentdef = true;
                else
                    error('could not load tpv file')
                end
                V.savingpar.savename = fname; oksaved = true;
            elseif isa(fname,'tpv_content')
                % Full data (trials + signals)
                data = fname;
                V.content = data;
                okcontentdef = true;
            elseif fn_ismemberstr(ext,{'mes' 'mesc'})
                % Read trials data from file, prompt for which trials to
                % take
                T = tps_trial.readheader(fname,'prompt');
                if V.internpar.loaddata, readdata(T), readrecording(T), end
                if isempty(T), set(hf,'pointer','arrow'), return, end
            elseif fn_ismemberstr(ext,tps_trial.knownExtensions())
                % Read trials data from file
                if doadd
                    try
                        % this will fail if V.content.trials has some
                        % trials which are the average of several files, in
                        % which case we just don't do the setdiff
                        fname = setdiff(cellstr(fname),{V.content.trials.file});
                    end
                    if isempty(fname)
                        set(hf,'pointer','arrow'), return
                    end
                    fname = char(fname); % go back to char representation
                end
                if dobin
                    if fn_ismemberstr(flag,{'avg' 'cat'})
                        error 'not implemented'
                    end
                    okspecial = true;
                    T = oi_binblocks(fname);
                elseif strcmp(flag,'')
                    T = tps_trial.readheader(fname);
                    if V.internpar.loaddata, readdata(T), readrecording(T), end
                elseif fn_ismemberstr(flag,{'avg' 'cat'})
                    T = tps_trial.readheader(fname,flag);
                    if V.internpar.loaddata, readdata(T), readrecording(T), end
                    okspecial = true;
                else
                    error programming
                end
            elseif isa(fname,'tps_trial')
                % Trials object
                T = fname;
            elseif isnumeric(fname) || iscell(fname)
                % Trials data
                data = fname;
                T = tps_trial(data);
            elseif strcmp(ext,'data')
                % Recording data (LabView format), display it in
                % separate figure and return
                ttotal = V.dt*V.nfr;
                x = fn_readdatlabview(fname); nx = size(x,1);
                fs = nx/ttotal;
                tens = 10^round(log10(fs));
                disp('assume a sampling frequency of 20,000Hz')
                E = tps_electrophys(fname,fs,V.trial.stim,V.trial.linedur);
                assignin('base','E',E)
                figure(180), set(180,'name','recording data','numbertitle','off')
                plot((0:E.n-1)*E.dt,E.signal)
                set(hf,'pointer','arrow'), return
            elseif strcmp(ext,'txt')
                % Text file, open it in editor
                edit(fname)
                set(hf,'pointer','arrow'), return
            else
                if ischar(fname)
                    errordlg(['cannot read file with extension ''' ext ''''])
                else
                    errordlg(['cannot read data of class ''' class(fname) ''''])
                end
                set(hf,'pointer','arrow'), return                
            end
            if ~oksaved
                % New data was loaded, but is not saved yet
                V.savingpar.savename = ''; 
                file_markchange(V);
            end
            
            % Check that optional action was taken
            if (dobin || ismember(flag,{'avg' 'cat'})) && ~okspecial
                warndlg 'option for binning/averaging/concatenating could not be taken into account'
            end
            
            % share internal parameters with content (which itself will
            % share it with every trial)
            V.content.internpar = V.internpar;
            
            % update signals
            if ~okcontentdef
                if doadd
                    oldntrial = V.content.ntrial;
                    addtrials(V.content,T)
                    V.content.ktrial = oldntrial+1;
                    okcontentdef = true;
                else
                    settrials(V.content,T,V.disppar.time,V.disppar.timecond, ...
                        V.panels.datacontrol.x)
                    V.content.seldotrial = V.options.preferences.seldotrial;
                    % selection needs to be set later only, when indices
                    % have been corrected for the new size of the data
                end
                linescans = strcmp({V.content.trials.type},'linescan');
                if any(linescans) && ~all(linescans), V.content.seldotrial = true; end
            end
            selection_showlist(V)
            
            % update display parameters whose values are also stored in
            % tpv_content object (data mode, conditions, recording names)
            data_dispparfromcontent(V)
            
            % figure title
            if isempty(V.savingpar.savename)
                set(V.hf,'name',V.skin)
            else
                set(V.hf,'name',[V.skin ': ' fn_fileparts(fname,'base')])
            end
            
            % update display
            if okcontentdef
                chgclipflag = true; % why should we keep the same clip?!
                display_changeview(V,'data&signals',chgclipflag)
            else
                display_changeview(V,'data')
            end
            
            % background compute spikes if 'guessspikesall' flag is on
            if V.internpar.guessspikesall
                computespikes(V.content)
            end
        end
        function data_copysettings(V,C0)
            if nargin<2
                C0 = fn_getfile('*.tpv',['Select ' V.skin ' file from which to copy settings']); 
                if isequal(C0,0), return, end
            end
            if ischar(C0), C0 = fn_loadvar(C0); end
            if ~isa(C0,'tpv_content'), error 'object to copy settings from must be of class tpv_content', end
            % current data
            T = V.content.trials;
            % get the first trial that is linked
            idx = find([C0.signal.dataopdef.link],1,'first');
            if isempty(idx), idx = 1; end
            T0 = C0.trials(idx);
            
            % Ask user which settings to copy
            settings = cell(0,2);
            % initial binning
            if T0.xbin~=1, settings(end+1,:) = {'xbin' ['initial spatial binning: ' num2str(T0.xbin)]}; end
            if T0.tbin~=1, settings(end+1,:) = {'tbin' ['initial temporal binning: ' num2str(T0.tbin)]}; end
            % pixel size and frame duration
            if ~strcmp(T0.xunit,'px') || T0.dx~=T0.xbin
                settings(end+1,:) = {'dx' ['pixel size: ' num2str(T0.dx) T0.xunit]};
                if T0.dz~=0, settings(end+1,:) = {'dz' ['zstack: image every ' num2str(T0.dz) T0.xunit]}; end
            end
            if ~strcmp(T0.xunit,'frame') || T0.dt~=T0.tbin, settings(end+1,:) = {'dt' ['frame period: ' num2str(T0.dt) T0.tunit]}; end
            % stimuli
            if ~any(diff([T0.nc T.nc]))
                t = T0.stimdetails;
                str = [];
                if ~isempty(t(1).name)
                    str = {t.name};
                elseif ~isempty(t(1).type)
                    str = {t.type};
                elseif ~all(fn_isemptyc({t.stim}));
                    str = cell(1,T0.nc);
                    for i=1:T0.nc, str{i}=['[' fn_chardisplay(t(i).stim) ']']; end
                end
                if ~isempty(str)
                    str = fn_strcat(str,', ');
                    nmax = 100; if length(str)>nmax, str = [str(1:nmax-3) '...']; end
                    settings(end+1,:) = {'stim' ['stim: ' str]};
                end
            end
            % data operations
            if ~isempty(T0.opdef)
                str = upper({T0.opdef.name});
                for i=find(~[T0.opdef.active]), str{i}=['(' str{i} ')']; end
                settings(end+1,:) = {'dataop' ['data operation(s): ' fn_strcat(str,', ')]}; 
            end
            % ask user
            if isempty(settings)
                errordlg('No settings were found that can be copied')
                return
            end
            spec = struct('title',{[] 'label' 'Select which settings to copy:'});
            for i=1:size(settings,1)
                [spec.(settings{i,1})] = deal(true,'logical',settings{i,2});
            end
            s = fn_structedit(spec,'title',V.skin);
            if isempty(s), return, end % user closed window, cancel
            
            % Apply parameters
            F = {'xbin' 'tbin' 'dx' 'dz' 'dt' 'stim' 'dataop'};
            for i=1:length(F), if ~isfield(s,F{i}), s.(F{i})=false; end, end
            % initial binning
            % (make it easy: remove selections)
            if s.xbin
                V.content.updateselection('reset')
                [T.xbin] = deal(T0.xbin);
            end
            if s.tbin, [T.tbin] = deal(T0.tbin); end
            % pixel size and frame duration
            if s.dx
                [T.dx T.dy] = deal(T0.dx);
                [T.xunit] = deal(T0.xunit);
                if s.dz, [T.dz] = deal(T0.dz); end
            end
            if s.dt
                [T.dt] = deal(T0.dt);
                [T.tunit] = deal(T0.tunit);
            end
            % stimuli
            if s.stim
                [T.stimtable] = deal(T0.stimtable);
                ids = T0.stimid;
                for i=1:length(T), T(i).setstim(ids); end
            end
            % data operations
            if s.dataop
                [T.opdef] = deal(T0.opdef);
                % check whether some indices might need to be updated
                pb = {};
                for i=1:length(T0.opdef)
                    op = T0.opdef(i);
                    if fn_ismemberstr(op.name,{'time' 'mask' 'cameranoise' 'heart' 'motion'})
                        kchoice = fn_switch(op.name,'time',2,1);
                        kind = kchoice+1; str = 'indices';
                        if strcmp(op.value{kchoice},str) && isnumeric(op.value{kind})
                            pb{end+1} = upper(op.name);
                        end
                    end
                end
                if ~isempty(pb)
                    waitfor(warndlg(['Operation(s) ' fn_strcat(pb,', ') ' involve using the data of specific spatial region(s). ' ...
                        'Such region(s) from previous experiment are probably unaccurate in the new experiment and need to be redefined.']))
                end
            end
            % update display
            if s.xbin || s.dataop
                display_changeview(V,'data&signals',false)
            else
                display_changeview(V,'ktrial')
            end
        end
        function data_addtrial(V,flag)
            % add trials, either through computation or from the base
            % workspace
            persistent spec
            if isempty(spec)
                spec = struct( ...
                    'dataflag', {'data'     {'data' 'dataop' 'data+sfr'}    'Based on'}, ...
                    'idx',      {1:V.ntrial 'double'    'Filter indices'}, ...
                    'idxout',   {[]         'double'    'Filter out indices'}, ...
                    'idlist',   {[]         'double'    'Stim id'}, ...
                    'rejected', {false      'logical'   'Keep rejected trials'}, ...
                    'dobadavg', {false      'logical'   'Average different conditions together'}, ...
                    'donorm',   {false      'logical'   'Normalize by the inter-trial variability'}, ...
                    'desc',     {''         'char 15'   'Description of new trial'} ...
                    );
            end
            c = watch(V); %#ok<NASGU>
            switch flag
                case 'average'
                    % prompt user
                    s = fn_structedit(spec);
                    if isempty(s), return, else spec(1) = s; end
                    desc = s.desc;
                    dosfr = strcmp(s.dataflag,'data+sfr') && V.content.trials(1).sfrchannel;
                    if dosfr, s.dataflag = 'data'; end
                    % filtering trial indices
                    ind = 1:V.ntrial;
                    if ~isempty(s.idx), ind = intersect(ind,s.idx); end %#ok<*ST2NM>
                    if ~isempty(s.idxout), ind = setdiff(ind,s.idxout); end
                    if ~isempty(s.idlist)
                        % note that filtering according to stim id will
                        % automatically discard multi-condition trials
                        ok = false(1,length(ind));
                        for i=1:length(ind)
                            id = V.content.trials(ind(i)).stimid;
                            ok(i) = isscalar(id) && ismember(id,s.idlist);
                        end
                        ind = ind(ok);
                    end
                    if ~s.rejected
                        ok = ([V.content.trials(ind).status]~='r');
                        ind = ind(ok);
                    end
                    if isempty(ind), errordlg('No trial fit conditions'), return, end
                    % prepare for normalization
                    if s.donorm
                        data = V.content.trials(ind(1)).(s.dataflag);
                        m = mean(data(:));
                        if m>-.1 && m<.1
                            datamean = 0;
                        elseif m>.9 && m<1.1
                            datamean = 1;
                        else
                            errordlg 'could not determine whether data is near zero or near one'
                            return
                        end
                        if dosfr, error 'normalization not implemented for sfr yet', end
                    end
                    % subset of trials
                    T0 = V.content.trials(ind);
                    % check compatibility
                    siz = cat(1,V.content.trials(ind).sizes);
                    if any(any(diff(siz(:,1:2),1,1))), errordlg 'frame sizes do not match', return, end
                    if any(diff(siz(:,3)))
                        answer = questdlg('number of frames do not match, truncate all trials to the minimum number of frames?', ...
                            '','Yes','Cancel','Cancel');
                        if ~strcmp(answer,'Yes'), return, end
                        if s.donorm, errordlg 'not compatible with ''donorm'' option', return, end
                    end
                    siz = siz(1,:);
                    ncond = siz(4);
                    stimids = cat(1,T0.stimid); % ntrial * ncond
                    ismultstim = any(any(diff(stimids,1,1)));
                    if ncond>1 && ismultstim, errordlg 'multi-condition trials but conditions do not fit', return, end
                    if any(diff([T0.sfrchannel])), errordlg 'only part of trials have an sfr channel', return, end
                    % average
                    ustimids = unique(stimids(:));
                    if ~ismultstim || s.dobadavg
                        data = trialaverage(T0,s.dataflag,V.internpar.floattype,true);
                        if dosfr, sfr = trialaverage(T0,'sfr',V.internpar.floattype,true); else sfr=[]; end
                        if s.donorm
                            datav = std(cat(5,T0.(s.dataflag)),1,5);
                            if datamean==0
                                data = fn_div(data,datav);
                            else
                                data = 1 + fn_div(data-1,datav);
                            end
                        end
                    else
                        nstim = length(ustimids);
                        siz(4) = nstim;
                        data = zeros(siz,V.internpar.floattype);
                        if dosfr, sfr = zeros(siz); else sfr = []; end
                        fn_progress('averaging condition',nstim)
                        for k=1:nstim
                            fn_progress(k)
                            indk = (stimids == ustimids(k));
                            data(:,:,:,k) = trialaverage(T0(indk),s.dataflag,V.internpar.floattype,false);
                            if dosfr, sfr(:,:,:,k) = trialaverage(T0(indk),'sfr',V.internpar.floattype,false); end
                            if s.donorm
                                datav = std(cat(5,T0(indk).(s.dataflag)),1,5);
                                if datamean==0
                                    data(:,:,:,k) = fn_div(data(:,:,:,k),datav);
                                else
                                    data(:,:,:,k) = 1+fn_div(data(:,:,:,k)-1,datav);
                                end
                            end
                        end
                        fn_progress end
                    end
                case 'command'
                    answer = inputdlg( ...
                        {'Add new trial from base workspace. Variable containing data:','Variable containing sfr','Description'}, ...
                        'Add trial', ...
                        1, ...
                        {'' '' 'user-defined'});
                    if isempty(answer), return, end
                    [name1 name2 desc] = deal(answer{:});
                    data = evalin('base',name1);
                    if isempty(name2), sfr = []; else sfr = evalin('base',name2); end
                    ismultstim = false;
                    T0 = V.trial;  % copy trial info from current trial
                otherwise
                    error('unknown flag ''%s''',flag)
            end
            T = tps_trial(data,T0(1)); % T = tps_trial(data,sfr,T0(1));
            T.opdef = T0(end).opdef;
            %T.sfrchannel = ~isempty(sfr);
            if isfield(T.addinfo,'comment'), T.addinfo = rmfield(T.addinfo,'comment'); end
            T.addinfo.description = desc;
            T.status = 's'; % special trial
            if ismultstim
                if s.dobadavg
                    % valid stim only if the non-empty ones are all equal
                    stim = getglobalstim(T.stimtable,ustimids);
                    setstim(T,struct('stim',stim,'name','undefined'))
                else
                    setstim(T,ustimids)
                end
            end
            addtrials(V.content,T)
            if ~isempty(V.disppar.trialsublist), V.disppar.trialsublist(end+1) = V.ntrial; end
            display_slidertrial(V)
            V.ktrial = V.ntrial;
        end
        function data_chgmode(V,inflag,dataflagold,dataflag)
            if strcmp(dataflagold,dataflag), return, end
            dotime = strcmp(inflag,'time');
            doimg = ~dotime;
            % change clipping?
            if doimg
                inidx = str2double(inflag(3));
                chgclipflag = false(1,2);
                chgclipflag(inidx) = ~fn_ismemberstr([dataflagold ' ' dataflag], ...
                    {'data shotnoise','shotnoise data', ...
                    'dataop shotnoiseop','shotnoiseop dataop'});
            else
                inidx = 0; chgclipflag = false;
            end
            % change in signals - automatic erase of old data
            if dotime, setdatamode(V.content,dataflag), end
            % update display
            display_changeview(V,'datamode',inidx,dotime,chgclipflag)
        end
        function data_conditionmenus(V,forceflag)
            % this function sets the conditions menus, but starts by
            % checking that the conditions displayed exist indeed
            forceflag = nargin>1 && strcmp(forceflag,'force');
            items = V.menus.items;
            F = {'im1','im2','time'};
            
            % any change?
            if V.ntrial==0, nc=1; else nc=V.nc; end
            st = unique([V.content.trials.stimtable]); 
            if isscalar(st) && ~any(diff(fn_itemlengths({V.content.trials.stimid}))) && ~any(row(diff(cat(1,V.content.trials.stimid))))
                st = st.table; 
                stids = [st.id];
                ids = V.content.trials(1).stimid;
                condnames = cell(1,V.nc);
                for k=1:V.nc, condnames{k} = st(stids==ids(k)).name; end
                donames = ~any(fn_isemptyc(condnames)) && (length(unique(condnames))==length(condnames));
            else
                st = []; 
                donames = false;
            end
            if nc==items.ncond && (~isfield(items,'stimtable') || isequal(items.stimtable,st)) && ~forceflag, return, end
            items.stimtable = st;
            items.ncond = nc;
            
            % check which conditions are displayed: conditions that are not
            % defined are replaced by 'C0'
            if nc>1
                dispcond = {V.disppar.im1cond V.disppar.im2cond V.disppar.timecond};
                ok = validcondition(dispcond,nc);
                [dispcond{~ok}] = deal('C0');
                [V.disppar.im1cond V.disppar.im2cond V.disppar.timecond] = deal(dispcond{:});
            end
            
            % delete existing menu entries
            for i=1:3
                f = F{i};
                delete(items.(f).condobjs)
                delete(items.(f).condspec)
                items.(f).condobjs  = [];
                items.(f).condspec  = [];
            end
            
            % no need for condition menus?
            if nc==1, V.menus.items=items; return, end
            
            % cell with each single condition
            condsimp = cell(1,nc);
            for k=1:nc, condsimp{k} = ['C' num2str(k-1)]; end
            % cell with pre-defined calculations
            % load from settings file, but only consider those which use
            % existing conditions
            condcalc = V.options.condcalc;
            accept = validcondition(condcalc,nc);
            condcalc = condcalc(accept);
            ncalc = length(condcalc);
            % almost the same for recently used calculations
            condmem = V.options.condmem;
            if size(condmem,1)==1
                % handle previous version
                condmemname = condmem;
            else
                condmemname = condmem(2,:);
                condmem = condmem(1,:);
                for i=1:length(condmem)
                    if isempty(condmemname{i}), condmemname{i}=condmem{i}; end
                end
            end
            accept = validcondition(condmem,nc);
            condmem = condmem(accept);
            condmemname = condmemname(accept);
            nmem = length(condmem);
            
            for i=1:3
                f = F{i};
                m = V.menus.(f);
                set(m,'enable','on')
                condobjs = [];
                % memory conditions
                for k=1:nmem
                    % note that the label of these menu items will
                    % change....
                    condobjs(end+1) = uimenu(m,'label',condmemname{k}, ...
                        'callback',@(hu,evnt)setpar(V,[f 'cond'],condmem{k})); %#ok<AGROW>
                end
                if nmem>0, set(condobjs(1),'separator','on'), end
                % simple conditions
                m1 = uimenu(m,'label','simple conditions');
                condspec(1) = m1;
                for k=1:nc
                    if donames && ~isempty(condnames{k})
                        namek = [condsimp{k} ' (' condnames{k} ')']; 
                    else
                        namek = condsimp{k};
                    end
                    condobjs(end+1) = uimenu(m1,'label',namek, ...
                        'callback',@(hu,evnt)setpar(V,[f 'cond'],condsimp{k})); %#ok<AGROW>
                end
                % calculated conditions
                m1 = uimenu(m,'label','calculated conditions');
                condspec(2) = m1;
                for k=1:ncalc
                    condobjs(end+1) = uimenu(m1,'label',condcalc{k}, ...
                        'callback',@(hu,evnt)setpar(V,[f 'cond'],condcalc{k})); %#ok<AGROW>
                end
                condspec(3) = uimenu(m,'label','new calculated condition...', ...
                    'callback',@(hu,evnt)data_conditionnewcalc(V,f));
                % all conditions
                if i==3
                    condobjs(end+1) = uimenu(m,'label','all conditions', ...
                        'callback',@(hu,evnt)setpar(V,[f 'cond'],'all conditions')); %#ok<AGROW>
                end
                % store names and handles
                items.(f).condspec  = condspec;
                items.(f).condobjs  = condobjs;
                % check the appropriate item
                condnames = [condmem condsimp condcalc fn_switch(i==3,'all conditions',{})];
                set(condobjs(strcmp(condnames,V.disppar.([f 'cond']))),'checked','on')
            end
            V.menus.items = items;
        end
        function data_conditionnewcalc(V,inflag)
            % update calculated conditions list
            ncalcmax = 20;
            [newcalc newname] = tps_definecondcalc(V.trial);
            if isempty(newcalc), return, end % figure close by user
            condcalc = V.options.condcalc;
            condcalc(strcmp(condcalc,newcalc)) = [];
            V.options.condcalc = [newcalc condcalc(1:min(end,ncalcmax-1))];
            % saving options and updating condition menus will be done by
            % data_chgcondition
            V.disppar.([inflag 'cond']) = newcalc;
            data_chgcondition(V,inflag,[],newcalc,newname);
        end
        function data_chgcondition(V,inflag,oldcond,newcond,newname) %#ok<*INUSL>
            % currently, the 'oldcond' argument is unused...
            % update recently-used conditions
            nmemmax = 5;
            condmem = V.options.condmem;
            imem = find(strcmp(condmem(1,:),newcond));
            if isempty(imem)
                condmem = [{newcond; newcond} condmem(:,1:min(end,nmemmax-1))];
            else
                condmem = condmem(:,[imem(1) setdiff(1:end,imem)]);
            end
            if nargin>=5, condmem{2,1} = newname; end
            V.options.condmem = condmem;
            saveoptions(V)
            % update condition menus
            data_conditionmenus(V,'force')
            % update display
            dotime = strcmp(inflag,'time');
            if dotime, setdatacond(V.content,newcond), end
            display_changeview(V,'datamode',fn_switch(dotime,0,str2double(inflag(3))),dotime,false)
        end
        function data_dataopdisplay(V)
            V.panels.datacontrol.x = V.trial.opdef;
            set(V.panels.dataitems.linked,'value',V.content.signals(1).dataopdef(V.ktrial).link);
            drawnow
        end
        function data_signalopdisplay(V)
            V.panels.signalcontrol.s = V.content.signals(1).opdef;
            drawnow
        end
        function data_opcallback(V,op,activechg) % data_opcallback(V,flag)
            % reads the controls, sets the description of data operation
            % in signals and trials, and updates display
            % reads the control
            
            op = tps_dataopdef(op);
            
            % Special operations
            chgop = false;
            for k=1:length(op)
                value = op(k).value;
                switch op(k).name
                    case 'space'
                        kchoice = 2; kind = 3;
                        if strcmp(value{kchoice},'current selection')
                            sel = V.a4d.SIt.selection.getsel(1);
                            op(k).value{kchoice} = 'frames';
                            if isempty(sel)
                                errordlg('no temporal selection')
                                ind = [];
                            else
                                ind = sel(1).dataind;
                            end
                            op(k).value{kind} = ind;
                            chgop = true;
                        end
                    case {'time' 'mask' 'cameranoise' 'heart' 'motion'}
                        kchoice = fn_switch(op(k).name,'time',2,1);
                        kind = kchoice+1;
                        ind = [];
                        switch value{kchoice}
                            case 'first selection'
                                sel = V.content.signals(1).x(V.ktrial,1).sel;
                                ind = sel.dataind;
                            case 'last selection'
                                sel = V.content.signals(1).x(V.ktrial,end).sel;
                                ind = sel.dataind;
                            case {'out selections' 'out selection'}
                                sel = [V.content.signals(1).x(V.ktrial,:).sel];
                                sel = sel([sel.active]);
                                indsel = cat(1,sel.dataind);
                                ind = setdiff(1:V.nx*V.ny,indsel);
                            case 'user'
                                answer = inputdlg('Type command that evaluates in the base workspace to either a list of pixel indices or a logical mask', ...
                                    V.skin,1,{''});
                                str = answer{1};
                                ok = ~isempty(str);
                                if ok
                                    try 
                                        ind = evalin('base',str);
                                        V.trial.data0(ind); % does this expression evaluate without error?
                                    catch
                                        ok = false; 
                                    end
                                end
                                if ~ok
                                    waitfor(errordlg('Command does not evaluate correctly to a list of indices or a logical mask'))
                                    ind = '[]'; % empty mask
                                end
                        end
                        if ~isempty(ind)
                            op(k).value{kchoice} = 'indices';
                            op(k).value{kind} = ind;
                            chgop = true;
                        end
                end
            end
            if chgop, V.panels.datacontrol.x = op; end
            
            % set content: automatic update of signals
            linked = get(V.panels.dataitems.linked,'value');
            setdataopdef(V.content,op,linked)
            
            % update display if necessary
            if activechg
                imgidx = ~isempty(regexp(V.disppar.im1,'op$')) + 2*~isempty(regexp(V.disppar.im2,'op$')); %#ok<RGXP1>
                if V.internpar.guessspikes, imgidx = 3; end % i don't remember why guessspikes flag requires updating images...
                dotime = any(strfind(V.disppar.time,'op')) || V.internpar.guessspikes;
                display_changeview(V,'datamode',imgidx,dotime,false)
            end
        end
        function data_opmem(V)
            % memorize current operation in all the linked trials
            
            % apply to all linked trials if the current trial is linked
            if V.content.signals(1).dataopdef(V.ktrial).link
                ktrials = find([V.content.signals(1).dataopdef.link]);
            else
                ktrials = V.content.ktrial;
            end
            
            % memorize
            memorizeop(V.content.trials(ktrials));
            
            % remove non up-to-date signals - TODO: further checks to
            % remove only what is necessary!
            dotime = any(strfind(V.disppar.time,'opmem'));
            if dotime, erasedata(V.content,ktrials), end
            
            % update display if necessary
            imgidx = ~isempty(strfind(V.disppar.im1,'opmem'))+2*~isempty(strfind(V.disppar.im2,'opmem'));
            if V.internpar.guessspikes, imgidx = 3; end % i don't remember why guessspikes flag requires images to be updated...
            display_changeview(V,'datamode',imgidx,dotime,false)
        end
        function data_signalop(V,op)
            setsignalopdef(V.content,op)
            display_changeview(V,'datamode',0,true,false)
        end
        function data_dispparfromcontent(V)
            signal = V.content.signals(1);
            % data mode for time display
            if ~strcmp(signal.datamode,V.disppar.time)
                V.disppar.time = signal.datamode;
            end
            % condition for time display
            if ~strcmp(signal.datacond,V.disppar.timecond)
                newcond = signal.datacond;
                V.disppar.timecond = newcond;
                % what follows is copied from data_chgcondition
                nmemmax = 5;
                condmem = V.options.condmem;
                imem = find(strcmp(condmem(1,:),newcond));
                if isempty(imem)
                    condmem = [{newcond; newcond} condmem(:,1:min(end,nmemmax-1))];
                else
                    condmem = condmem(:,[imem(1) setdiff(1:end,imem)]);
                end
                V.options.condmem = condmem;
                saveoptions(V)
            end
            % analog recordings
            for k=1:V.ntrial
                r = V.content.trials(k).recording;
                if ~isempty(r), break, end
            end
            rnames = unique({r.name});
            missing = ~ismember(lower(rnames),lower({V.disppar.timerecording.name}));
            if any(missing)
                V.disppar.timerecording = [V.disppar.timerecording struct('name',rnames(missing),'val',repmat({false},1,sum(missing)))];
                init_menus(V) % reinit menus so that new recording names will appear in the recordings list
            end
            % reinit menus, to make the above changes appear correctly
            init_menus(V)
        end
    end
    methods
        function data_attachanalog(V,files)
            if nargin==1
                attachanalogfile(V.content.trials)
            else
                attachanalogfile(V.content.trials,files)
            end
            setelectrophys(V.content,tps_electrophys(V.content.trials));
            data_dispparfromcontent(V)
            display_changeview(V,'timesignal')
        end
        function data_setelectrophy(V,channelname)
            if nargin<2
                for k=1:V.ntrial
                    r = V.content.trials(k).recording;
                    if ~isempty(r), break, end
                end
                rnames = unique({r.name});  
                if isempty(r), errordlg 'This experiment does not seem to have any recording', return, end
                channelname = fn_input('Select electrophy channel',rnames{1},rnames);
                if isempty(channelname), return, end
            end
            E = tps_electrophys(V.content.trials,channelname);
            setelectrophys(V.content,E);
            data_dispparfromcontent(V)
            display_changeview(V,'timesignal')
        end
        function data_spikes(V,whatflag,parflag,docompute)
            % function data_spikes(V,'current|all',parflag,docompute)
            %---
            % compute (or erase) spikes according to 'flag' and update display
            % possible flags are 'current' [default], 'trialdrift', 'all',
            % 'alldrift', 'erase' and 'eraseall'
            
            % input
            if nargin<4, error('not enough input arguments'), end
            
            % start watch
            c = watch(V); %#ok<NASGU>
            
            % which spikes to compute/erase
            switch whatflag
                case 'currentregion'
                    trials  = 1:V.ntrial;
                    regions = V.disppar.currentsel;
                case 'current'
                    if V.disppar.timealltrials
                        if strcmp(parflag,'nospike'), disp('must be in ''trial'' view mode for this action'), return, end
                        trials  = 1:V.ntrial;
                        regions = V.disppar.currentsel;
                    else
                        trials = V.ktrial;
                        regions = 1:V.content.nx;
                    end
                case 'all'
                    trials  = 1:V.ntrial;
                    regions = 1:V.content.nx;
                otherwise
                    error('unknown flag ''%s''',whatflag)
            end
            
            % update spike parameters and compute spikes if necessary
            if ischar(parflag) && any(strfind(parflag,'rm'))
                erasespikes(V.content,trials,regions,strcmp(parflag,'rmall'))
            else
                computespikes(V.content,trials,regions,parflag,docompute)
            end
            
            % update display
            display_changeview(V,'spikes')
        end
        function data_spikecontrol(V,flag)
            % default par
            switch flag
                case 'calcium'
                    spikefun = @tps_mlspikes;
                case 'vsd'
                    spikefun = @tps_vsdspikes;
            end
            s = fn_structmerge(struct('fun',spikefun),feval(spikefun,'par'));
            defaultpar = V.content.signals(1).spikepar;
            if ~isempty(defaultpar) && isequal(defaultpar.fun,spikefun)
                s = fn_structmerge(s,defaultpar,'skip');
            end
            
            % control
            fn_control(s,@(p)setandcomp(p));
            function setandcomp(p)
                computespikes(V.content,1:V.ntrial,1:V.content.nx,p,false)
                data_spikes(V,'current',p,true)
            end
        end
        function data_savenewdata(V,flag)
            % user input for new file names
            fname = fn_fileparts(V.savingpar.savename,'noext');
            switch flag
                case 'dataop'
                    extm = '.corr';
                    if ~isempty(fname), fname = [fname ' - corrected.tpv']; end
                case 'bin'
                    s = fn_structedit('xbin',{2 'stepper 1 1 100' 'spatial binning'},'tbin',{1 'stepper 1 1 100' 'temporal binning'});
                    if isempty(s), return, end
                    xbin = s.xbin;
                    tbin = s.tbin;
                    extm = ['.bin' num2str(xbin)];
                    if ~isempty(fname), fname = [fname ' - bin' num2str(xbin) '.tpv']; end
                case 'coregister'
                    if ~isempty(fname), fname = [fname ' - coregistered.tpv']; end
                otherwise
                    error 'unknown flag'
            end
            fname = fn_savefile('*.tpv','Select name for file and folder',fname);
            if ~fname, return, end
            fbase = fn_fileparts(fname,'noext');
            fname = [fbase '.tpv'];
            savedir = [fbase '_data'];
            if ~exist(savedir,'dir'), mkdir(savedir), end
            T = V.content.trials;
            % some sub-selection for 'dataop' and 'bin'; in the case of
            % 'coregister', there will be another user input inside
            % function tps_register
            if ~strcmp(flag,'coregister')
                s = struct( ...
                    'extension',        {[extm '.mat']   'char'}, ...
                    'save__rejected__trials',   {false  'logical'}, ...
                    'save__special__trials',    {false  'logical'} ...
                    );
                s = fn_structedit(s);
                if isempty(s), return, end
                % (keep only requested trials)
                if ~s.save__rejected__trials
                    % remove rejected trials
                    T([T.status]=='r') = [];
                end
                if ~s.save__special__trials
                    % remove special trials
                    T([T.status]=='s') = [];
                end
                ntr = length(T);
                if ntr==0, return, end
            end
            c = watch(V); %#ok<NASGU>
            
            % new tps_trial object
            if strcmp(flag,'coregister')
                % dedicated m-file for this
                % (first get the current time selection)
                sel = V.a4d.SIt.selection.getsel(1);
                if isempty(sel)
                    ind = [];
                else
                    ind = sel(1).dataind;
                end
                Top = tps_register(T,fbase,{V.ktrial ind});
                if isempty(Top), return, end % empty output signals an error
            else
                Top = T;
                fn_progress('trial',ntr)
                for k=ntr:-1:1
                    fn_progress(1+ntr-k)
                    switch flag
                        case 'dataop'
                            data = single(T(k).dataop);
                            Top(k) = tps_trial(data,T(k));
                        case 'bin'
                            data = fn_float(T(k).data);
                            data = fn_bin(data,[xbin xbin tbin]);
                            data = single(data);
                            Top(k) = tps_trial(data,T(k));
                            Top(k).dx = T(k).dx*xbin;
                            Top(k).dy = T(k).dy*xbin;
                            Top(k).dt = T(k).dt*tbin;
                        case 'coregister'
                            [shift e datareg] = fn_register(T(k).data,par); %#ok<ASGLU>
                            Top(k) = tps_trial(datareg,T(k));
                    end
                    ftrial = [savedir '/' fn_fileparts(T(k).file,'base') s.extension];
                    savedata(Top(k),ftrial)
                end
            end
            
            % Save tpv_content object and load into new window
            V2 = tpview(V.skin);
            settrials(V2.content,Top,'data',V.content.datacond,tps_dataopdef.empty(1,0),copy(V.a4d.SI1.selection))
            file_save(V2,fname)
            display_changeview(V2,'data&signals',false)
        end
    end
    
    % FILE
    methods
        function file_open(V,fname,flag)
            if nargin<2 || isempty(fname)
                [dataext tpstrial_prompt] = tps_trial.knownExtensions;
                dataextstar = fn_map(dataext,@(s)['*.' s]);
                dataextstr = fn_strcat(dataextstar,';');
                fname = fn_getfile( ...
                    [{['*.tpv;' dataextstr],['All ' V.skin ' Files']}; ...
                    {'*.tpv', [V.skin ' Analysis']}; ...
                    tpstrial_prompt; ...
                    {'*.*', 'All Files (*.*)'}], ...
                    ['Open ' V.skin ' File']);
            end
            if isequal(fname,0), return, end
            if nargin<3, flag=''; end
            data_loadfile(V,fname,flag)
        end
        function file_save(V,fname)
            if nargin<2
                fname = V.savingpar.savename;
            end
            if isempty(fname)
                % suggest a name
                files = char(V.content.trials.file);
                if ~isempty(files)
                    f = find(var(files),1,'first');
                    if isempty(f), f = size(files,2)+1; end
                    [d fname] = fileparts(files(1,1:f-1));
                    f = find(fname=='_',1,'last');
                    if f, fname = fname(1:f-1); end
                else
                    fname = '';
                end
                [fname filterindex] = fn_savefile( ...
                    {'*.tpv' [V.skin ' analysis']; '*.BLK*' 'Optical Imaging data'; '*.tiff' 'Tiff images'}, ...
                    'Save data in',fname);
                if ~fname, return, end
                [d base ext] = fileparts(fname);
                switch filterindex
                    case 1
                        saveformat = 'tpv';
                        if ~isempty(ext) && ~strcmpi(ext,'.tpv'), error('bad extension'), end
                        fname = fullfile(d,[base '.tpv']);
                    case 2
                        saveformat = 'BLK';
                    case 3
                        saveformat = 'tiff';
                end
            else
                [d base ext] = fileparts(fname);
                switch ext
                    case {'' '.tpv'}
                        saveformat = 'tpv';
                    case {'.BLK' '.tiff'}
                        saveformat = ext(2:end);
                    otherwise
                        errordlg('unknown file format')
                        return
                end
            end
            
            c = watch(V); %#ok<NASGU>
            
            % save/remove trials whose data is not saved
            idxbad = false(1,V.ntrial);
            for i=1:V.ntrial, idxbad(i) = isempty(V.content.trials(i).file); end
            idxbad = find(idxbad);
            if ~isempty(idxbad)
                if ~strcmp(saveformat,'tpv')
                    errordlg('Some temporary trials exist, save first in ''tpv'' file format')
                end
                answer = questdlg('Save or Remove temporary trials?', ...
                    'confirmation needed','Save','Remove','Cancel','Cancel');
                switch answer
                    case 'Cancel'
                        return
                    case 'Remove'
                        rmtrial(V.content,idxbad)
                    case 'Save'
                        T = V.content.trials(idxbad);
                        nf = length(idxbad);
                        descs = cell(1,nf);
                        for i=1:nf, ai = T(i).addinfo; if isstruct(ai) && isfield(ai,'description'), descs{i} = ai.description; else descs{i}=''; end, end
                        okdesc = (length(unique(descs))==nf);
                        prompt = cell(1,nf);
                        for i=1:nf, prompt{i} = ['trial ' num2str(idxbad(i))]; end
                        default = cell(1,nf);
                        for i=1:nf
                            if okdesc
                                desc = regexprep(descs{i},'[^a-zA-Z0-9_-,;.()\[\]]','-');
                            else
                                desc = num2str(idxbad(i),'trial%.2i');
                            end
                            default{i} = [base '_' desc];
                        end
                        name='Enter file names';
                        numlines=1;
                        answer = inputdlg(prompt,name,numlines,default);
                        if isempty(answer), return, end
                        dd = [d '/' base '_userdata/'];
                        if ~exist(dd,'dir'), mkdir(dd), end
                        drawnow
                        for i=1:nf, savedata(T(i),[dd answer{i}]), end
                end
            end
            
            % save tpview data
            drawnow
            switch saveformat
                case 'tpv'
                    x = V.content;  %#ok<NASGU>
                    crd = V.internpar.readdata;
                    V.internpar.readdata = false; % don't read data, in particular if data0 field of trials are empty, save them empty rather than computing them
                    save(fname,'x','-MAT')
                    V.internpar.readdata = crd;
                    disp saved
                    V.savingpar.savename = fname;
                    V.savingpar.chg  = false;
                case 'BLK'
                    tps_saveBLK(V.content.trials,[d '/' base])
                case 'tiff'
                    fn_saveimg(cat(3,V.data),fname)
            end
            
            % update display
            fn_getfile('REP',fileparts(V.savingpar.savename)) % make current directory the one containing the file
            file_list(V,true)
            if strcmp(saveformat,'tpv')
                display_header(V)
                set(V.hf,'name',[V.skin ': ' fn_fileparts(fname,'base')])
            end
        end
        function file_list(V,keepvalue)
            hlist = V.panels.filesitems.list;
            if nargin<2, keepvalue=false; end
                      
            % list directory
            wd = fn_getfile('REP');
            list = dir(wd);
            
            % order directories according to names
            isdir = cat(1,list.isdir);
            d = find(isdir);
            
            % order files according to types
            f = find(~isdir);
            nf = length(f);
            ext = lower(fn_fileparts({list(f).name},'extnodot'));
            [~, ord] = sort(2*~strcmp(ext,'tpv')+~ismember(ext,tps_trial.knownExtensions()));
            f = f(ord);
            list = list([d; f]);
            
            % group together files which match the same prototype
            % 'patter_number.ext&description'
            if V.disppar.filelistgroup
                oldpattern = [];
                groupstart = 0;
                k = 0;
                % attention! no for loop because k can be changed inside
                % the loop
                while k<=length(list)
                    k = k+1;
                    if k<=length(list)
                        fname = list(k).name;
                        f1 = find(fname=='_' | fname=='-',1,'last');
                        f2 = find(fname=='.',1,'first');
                        x = str2double(fname(f1+1:f2-1));
                        % file follows the prototype pattern_number.ext?
                        if ~isempty(f1) && ~isempty(f2) && ~isnan(x)
                            pattern = [fname(1:f1) '#.' fname(f2+1:end)];
                        else
                            pattern = [];
                        end
                    else
                        % reached end
                        pattern = [];
                    end
                    % close previous group?
                    if groupstart && ~strcmp(pattern,oldpattern)
                        if k>=groupstart+1
                            list(groupstart).name = strvcat(list(groupstart:k-1).name); %#ok<VCAT>
                            list(groupstart).name = [oldpattern ' [' num2str(k-groupstart) ' files]'];
                            list(groupstart+1:k-1) = [];
                            k = groupstart+1;
                        end
                        groupstart = 0;
                    end
                    % start of new group?
                    if ~groupstart && ~isempty(pattern)
                        groupstart = k;
                        oldpattern = pattern;
                    end
                end
            end
            
            % display list
            if ischar(keepvalue)
                value = find(strcmp(keepvalue,{list.name}));
                if isempty(value), value = 1; end
                listboxtop = max(1,value-5);
            elseif keepvalue
                listboxtop = get(hlist,'listboxtop');
                listboxtop = min(length(list),listboxtop);
                value = get(hlist,'value');
                value = min(length(list),value(1));
            else
                listboxtop = 1;
                value = 1;
            end
            
            if isempty(list), str = {}; else str = {list.name}; end
            set(hlist, ...
                'userdata',struct('wd',wd,'list',list), ...
                'string',str, ...
                'listboxtop',listboxtop,'value',value)
        end
        function file_select(V,flag)
            % flag marks which control was used to call file_select
            hf = V.grob.hf;
            hlist = V.panels.filesitems.list;
            userdata = get(hlist,'userdata');
            if userdata.wd(end)~='/'; userdata.wd(end+1) = '/'; end
            klist = get(hlist,'value');
            if klist>length(userdata.list), klist=1; end
            fil = userdata.list(klist);
            fname = char(fil.name);
            fnamecomp = [repmat(userdata.wd,[size(fname,1) 1]) fname];
            
            % which action?
            switch flag
                case 'update'
                    file_list(V,true)
                case 'group'
                    setpar(V,'filelistgroup',~V.disppar.filelistgroup); % automatic update
                case {'files','open','avg','cat','add','bin', ...
                        'bin&add','bin&avg','bin&cat','avg&add','bin&avg&add'}
                    opennow = ~strcmp(flag,'files')...
                        || strcmp(get(hf,'selectiontype'),'open');
                    if ~opennow, return, end
                    if fil(1).isdir
                        keepvalue = false;
                        switch fil.name
                            case '.'
                                fnamecomp = userdata.wd;
                            case '..'
                                fnamecomp = fn_fileparts(userdata.wd,'path');
                                keepvalue = fn_fileparts(userdata.wd,'name'); % trick to get back in the parent directory and place selection and view on the current directory
                        end
                        fn_getfile('REP',fnamecomp)
                        file_list(V,keepvalue)
                    else
                        if fn_ismemberstr(flag,{'files' 'open'})
                            flag=''; 
                        end
                        data_loadfile(V,fnamecomp,flag)
                    end
                case 'acqsettings'
                    tps_trial.EditAcqSettings(fnamecomp(1,:))
                case 'delete'
                    confirmation = questdlg({'Are you sure to delete files?',fname}, ...
                        'Confirm file deletion','Yes','No','No');
                    if strcmp(confirmation,'Yes')
                        swd = pwd;
                        cd(userdata.wd)
                        fnamec = cellstr(fname);
                        delete(fnamec{:});
                        cd(swd)
                        file_list(V,true)
                    end
                case 'rename'
                    fnamec = cellstr(fname);
                    nf = length(fnamec);
                    prompt=cell(1,nf); %num2cell(char('a'+(0:nf-1)));
                    name='Rename files';
                    numlines=1;
                    defaultanswer=fnamec;
                    answer = inputdlg(prompt,name,numlines,defaultanswer);
                    if isempty(answer), return, end
                    swd = pwd;
                    cd(fn_getfile('REP'))
                    for i=1:nf
                        if ~strcmp(fnamec{i},answer{i})
                            movefile(fnamec{i},answer{i})
                        end
                    end
                    cd(swd)
                    file_list(V,true)
                case 'favoritefolders'
                    u = V.panels.filesitems.favorites; % popup menu
                    favorites = V.options.favoritefolders;
                    str = get(u,'string');
                    sel = get(u,'value'); set(u,'value',1)
                    changed = false;
                    if sel==1
                        % 'Favorite folders' -> nothing to do
                        return
                    elseif sel<length(str)
                        % Go to selected folder
                        kfolder = sel-1;
                        if ~exist(favorites(kfolder).path,'dir')
                            answer = questdlg(['Folder ''' favorites(kfolder).path ''' does not exist.' ...
                                'What do you want to do?'], V.skin, ...
                                'Select new location', 'Remove entry', 'Cancel', 'Cancel');
                            switch answer
                                case 'Select new location'
                                    p = fn_getdir(['Select new location for folder ''' favorites(kfolder).name '''.']);
                                    if isequal(p,0), return, end
                                    name = inputdlg('Folder nickname',V.skin,1,{favorites(kfolder).name});
                                    if isempty(name) || isempty(name{1}), return, end                                    
                                    favorites(kfolder).name = name{1};
                                    favorites(kfolder).path = p;                                    
                                    changed = true;
                                    fn_getfile('REP',favorites(kfolder).path)
                                    file_list(V)
                                case 'Remove entry'
                                    favorites(kfolder) = [];
                                    changed = true;
                                case 'Cancel'
                                    return
                            end
                        else
                            fn_getfile('REP',favorites(kfolder).path)
                            file_list(V)
                        end
                    else
                        % 'organize...'
                        [favorites changed] = fn_listorganize(favorites,'name',@newfavoritefolder);
                    end

                    if changed
                        V.options.favoritefolders = favorites;
                        saveoptions(V)
                        str = {'Favorite folders' favorites.name 'organize...'};
                        set(u,'string',str)
                    end
                otherwise
                    error programming
            end
            
            % sub-function for adding new favorite folder
            function s = newfavoritefolder
                s = [];
                p = fn_getdir('Select location of favorite folder.');
                if isequal(p,0), return, end
                name = inputdlg('Folder nickname',V.skin,1,{fn_fileparts(p,'name')});
                if isempty(name) || isempty(name{1}), return, end
                s = struct('name',name{1},'path',p);
            end
        end
        function file_comments(V)
            % function file_comments(V)
            %---
            % open in the editor a file where to write comments on the
            % experiment
            if isempty(V.savingpar.savename), errordlg 'Save tpv file first', return, end
            fname = '';
            if ~isempty(V.content.docfile) && exist(V.content.docfile,'file')
                fname = V.content.docfile;
            end
            if isempty(fname)
                fbase = V.savingpar.savename;
                d = [fn_fileparts(fbase,'noext') '_userdata/'];
                if ~exist(d,'dir'), mkdir(d), end
                fname = [d fn_fileparts(fbase,'base') '_comments.txt'];
                % create empty file
                fid=fopen(fname,'w'); fclose(fid);
                % save
                V.content.docfile = fname;
                file_save(V)
            end
            edit(fname)
        end
        function file_markchange(V)
            if V.savingpar.chg, return, end % already marked as changed
            V.savingpar.chg = true;
            display_header(V)
        end
    end
    
    % SELECTION
    methods
        function selection_showlist(V)
            storesignals = V.content.signals(2:end);
            if isempty(storesignals)
                names = {};
            else
                names = {storesignals.name};
            end
            set(V.panels.selitems.list, ...
                'value',1, ...
                'string',[{'Choose selection to load'} names])
        end
        function selection_save(V,flag)
            % locate entry in collection
            selname = get(V.panels.selitems.input,'string');
            if isempty(selname)
                msgbox('please give a name to selection')
                return
            end
            if isscalar(V.content.signals)
                idx = [];
            else
                idx = find(strcmp(selname,{V.content.signals(2:end).name}))+1;
            end
            % save or delete
            switch flag
                case 'save'
                    % create new entry in V.content.signals?
                    if isempty(idx)
                        idx = length(V.content.signals)+1;
                    else
                        answer = questdlg('overwrite existing selection?');
                        if ~strcmp(answer,'Yes'), return, end
                    end
                    % save
                    storesignals(V.content,idx,selname);
                case 'delete'
                    confirmation = questdlg(['Are you sure to delete selection ''' selname '''?'], ...
                        'Confirm selection deletion','Yes','No','No');
                    if strcmp(confirmation,'No'), return, end
                    V.content.signals(idx) = [];
                otherwise
                    error programming
            end
            %             % save everything in the current *.tpv file
            %             file_save(V)
            % update list display
            selection_showlist(V)
        end
        function selection_load(V)
            % change current selection from value stored in
            % V.signals(2:end)
            idx = get(V.panels.selitems.list,'value');
            if idx==1, return, end
            c = watch(V); %#ok<NASGU>
            % use routine of tpv_content class (update opdef in trials, etc.)
            loadsignals(V.content,idx);
            % update display of data and signal operations
            data_dataopdisplay(V)
            data_signalopdisplay(V)
            % update display parameters whose value are stored in signals
            data_dispparfromcontent(V)
            % update data display
            dochgopdef = true; % it became too difficult to check whether operation definition in trials has changed
            chgimg = dochgopdef && any(strfind([V.disppar.im1 V.disppar.im2],'op'));
            display_changeview(V,'signals',chgimg)
            % update default saving name
            set(V.panels.selitems.input,'string',V.content.signals(idx).name)
        end
    end
    
    % DATA DISPLAY
    methods
        function display_changeview(V,flag,varargin)
            % function display_changeview(V,'data&signals',chgclipflag)
            % function display_changeview(V,'data')
            % function display_changeview(V,'ktrial')
            % function display_changeview(V,'signals',chgimg)
            % function display_changeview(V,'datamode',imgidx,dotime,chgclipflag)
            % function display_changeview(V,evnt,source)
            
            a = V.a4d;
            if isa(flag,'fn4Devent')
                evnt=flag;
                flag=evnt.flag;
            end
            % show watch, deactivate listeners
            c = watch(V); %#ok<NASGU>
            if ~fn_ismemberstr(flag,{'timeline'})
                % upon 'timeline' change, complex updates occur, in
                % particular the selections are reset, and this needs to be
                % detected by the geometry listener
                [a.listeners.Enabled] = deal(false);
            end
            % handle event
            switch flag
                % change in data
                case 'data&signals'
                    % menus and controls
                    data_conditionmenus(V)
                    V.disppar.trialsublist = [];
                    display_slidertrial(V)
                    display_header(V,V.content.ktrial)
                    set(V.menus.items.tpview.seldotrial,'checked',onoff(V.content.seldotrial))
                    set(V.menus.items.time.line,'checked',onoff(V.content.timeline))
                    % operations
                    data_dataopdisplay(V)
                    data_signalopdisplay(V)
                    % selection
                    chgsiz = display_updatesize(V);
                    chgclip = varargin{1};
                    display_image(V,true,chgclip|chgsiz)
                    % selection (after image, to have the correct size in V.a4d.SI1 and SI2)
                    display_updateselection(V) % get SI1.selection from signals
                    updateselection(V.content,'ij',[],V.a4d.G.ijkl(1:2)) % set current pixel in content (and get signal at this location if no selection)
                    V.disppar.currentsel = V.content.nx;
                    % time
                    timeold = V.disppar.time;
                    timenew = V.content.datamodes;
                    timenew = fn_switch(isempty(timenew{2}),timenew{1},[timenew{1} '_' timenew{2}]);
                    if ~strcmp(timeold,timenew)
                        V.disppar.time = timenew;
                        set(V.menus.items.time.(timeold),'checked','off')
                        set(V.menus.items.time.(timenew),'checked','on')
                    end
                    display_time(V)
                    display_stimandevents(V)
                    % more
                    display_more_update(V)
                case 'data'
                    % menus and controls
                    data_conditionmenus(V)
                    V.disppar.trialsublist = [];
                    display_slidertrial(V)
                    display_header(V,V.content.ktrial)
                    set(V.menus.items.tpview.seldotrial,'checked',onoff(V.content.seldotrial))
                    % operations
                    data_dataopdisplay(V)
                    data_signalopdisplay(V)
                    % image
                    chgsiz = display_updatesize(V);
                    display_image(V,true,chgsiz)
                    % selection
                    setselection(V.content,a.SI1.selection)
                    updateselection(V.content,'ij',[],V.a4d.SI1.ij)
                    V.disppar.currentsel = max(1,min(V.disppar.currentsel,V.content.nsel));
                    % time
                    display_time(V)
                    display_stimandevents(V)
                    % more
                    display_more_update(V)
                case {'ktrial' 'chgtrial'}
                    if strcmp(flag,'chgtrial'), display_slidertrial(V), end
                    data_conditionmenus(V)
                    if isempty(V.disppar.trialsublist)
                        kslider = V.ktrial;
                    else
                        kslider = find(V.disppar.trialsublist==V.ktrial);
                        if ~isscalar(kslider), error('selected trial is not inside sub-list'), end
                    end
                    set(V.a4d.slidertrial,'value',kslider)
                    display_header(V,V.content.ktrial)
                    data_dataopdisplay(V)
                    chgsiz = display_updatesize(V);
                    if V.content.seldotrial
                        display_updateselection(V) % get SI1.selection from signals
                    end
                    display_image(V,true,chgsiz)
                    if ~V.disppar.timealltrials || ~strcmp(flag,'ktrial')
                        display_time(V)
                        display_stimandevents(V)
                    end
                    display_more_update(V)
                case 'signals'
                    set(V.menus.items.tpview.seldotrial,'checked',onoff(V.content.seldotrial))
                    set(V.menus.items.time.line,'checked',onoff(V.content.timeline))
                    display_updateselection(V) % get SI1.selection from signals
                    V.disppar.currentsel = V.content.nx;
                    updateselection(V.content,'ij',[],V.a4d.SI1.ij)
                    chgimg = varargin{1};
                    if chgimg
                        display_image(V,true,true)
                        display_more_update(V)
                    end
                    display_time(V)
                case 'datamode'
                    [imgidx dotime chgclipflag] = deal(varargin{:});
                    if imgidx, display_image(V,true,chgclipflag,imgidx), display_more_update(V), end
                    if dotime, display_time(V), display_stimandevents(V), end
                    % change in display options
                case 'timesignal'
                    display_time(V)
                case 'timealltrials'
                    D = V.a4d.time;
                    %                     switch V.disppar.timealltrials
                    %                         case true
                    %                             if strcmp(D.linecol,'sel')
                    %                                 D.linecol = fn_colorset(V.disppar.currentsel);
                    %                             end
                    %                         case false
                    %                             if isnumeric(D.linecol)
                    %                                 D.linecol = 'sel';
                    %                             end
                    %                     end
                    display_time(V,'reset')
                case 'timefft'
                    display_updatesize(V)
                    display_time(V,'redisplay')
                case 'timeline'
                    % special: view time courses as if from line scan
                    if ~strcmp(V.type,'movie'), return, end
                    display_updatesize(V)
                    display_image(V,true)
                    display_updateselection(V) % get SI1.selection from signals
                    display_time(V)
                    % 'more' displays became bad
                    display_more_update(V,'disconnect')
                case 'spikes'
                    if ~strcmp(V.disppar.timespike,'spikesignal')
                        display_time(V)
                    end
                    if any([strcmp({V.disppar.im1 V.disppar.im2},'spikes') ...
                            V.disppar.im1spikecomb V.disppar.im2spikecomb])
                        display_image(V,true)
                        display_more_update(V)
                    end
                case 'stim'
                    display_stimandevents(V)
                    % events from displays
                case 'ijkl2'
                    % source is G
                    chg = (a.G.ijkl2~=evnt.oldvalue);
                    if chg(3)
                        a.SIt.ij2 = a.G.ijkl2(3);
                    end
                case 'ijkl'
                    % source is G
                    chg = (a.G.ijkl2~=evnt.oldvalue);
                    if any(chg(1:2))
                        updateselection(V.content,'ij',[],a.G.ijkl(1:2))
                        if ~V.content.nsel
                            % a fake 'one-pixel' selection has been created
                            display_time(V)
                        end
                    end
                case 'ij2'
                    % source is SIt
                    if ~strcmp(V.disppar.timefft,'compute')
                        a.G.ijkl2(3) = a.SIt.ij2;
                    end
                case 'ij'
                case 'zoom'
                    source = varargin{1};
                    if ~strcmp(V.disppar.timefft,'compute')
                        switch source
                            case 'G'
                                a.SIt.zoom = a.G.zoom(3,:);
                            case 'SIt'
                                a.G.zoom(3,:) = a.SIt.zoom;
                        end
                    end
                case 'selection'
                    source = varargin{1};
                    if strcmp(evnt.selflag,'reset')
                        % update through fn4D objects
                        switch source 
                            case 'G'
                                updateselection(a.SIt,'reset')
                            case 'SIt'
                                updateselection(a.G,[],'reset')
                        end
                        % additional display updates
                        V.disppar.currentsel=1;
                        updateselection(V.content,'reset')
                        display_time(V,'reset')
                        display_more_update(V)
                    elseif strcmp(source,'G') && any(evnt.dims==1) % can be 1 (linescan mode) or [1 2] (normal mode)
                        % update index of current selection
                        switch evnt.selflag
                            case 'reorder'
                                perm = evnt.value;
                                if isscalar(evnt.ind)
                                    V.disppar.currentsel = find(perm==evnt.ind);
                                else
                                    V.disppar.currentsel = V.content.nsel;
                                end
                            case {'remove','reset','all'}
                                seldimsnum = a.im1.seldims-'x'+1;
                                V.disppar.currentsel = numsel(a.G.selection,seldimsnum);
                                if V.disppar.currentsel==0, V.disppar.currentsel=1; end
                            case 'indices'
                                % no change in V.disppar.currentsel
                            otherwise
                                V.disppar.currentsel = evnt.ind(end);
                        end
                        % replace abstract representation of the change by
                        % the result to avoid recomputing things
                        doselshift = strcmp(evnt.selflag,'affinity') && length(evnt.ind)==V.content.nsel && V.content.seldotrial;
                        if doselshift, mov = evnt.value; end
                        if fn_ismemberstr(evnt.selflag,{'add' 'affinity'})
                            evnt.selflag = 'change';
                            evnt.value = a.G.selection.getselset(evnt.dims).singleset(evnt.ind);
                        end
                        % update internal saving of selections
                        if doselshift
                            % (selection shift)
                            % update the memory of general shift
                            % between trials; note that this will replace
                            % the selection of the next trials that have
                            % the same shift value by the selection of the
                            % current trial, even if initially the
                            % selections are not the same!!!
                            updateselectionshift(V.content,mov,evnt.value);
                        else
                            % (selection)
                            if strcmp(evnt.selflag,'all'), evnt.value = V.a4d.G.selection.getselset([1 2]); end
                            updateselection(V.content,evnt.selflag,evnt.ind,evnt.value)
                        end
                        % update display
                        display_time(V,evnt.selflag,evnt.ind,evnt.value) % TODO: update only indices which changed -> from tpv_content.operation
                    elseif (strcmp(source,'SIt') || (strcmp(source,'G') && isequal(evnt.dims,3))) && ~strcmp(V.disppar.timefft,'compute')
                        switch source
                            case 'SIt'
                                % do not transmit to G when in 'MOVIE'
                                % mode!!!
                                if get(V.panels.moreitems.movie,'value')
                                    % update movie time range
                                    display_movie(V)
                                else
                                    updateselection(a.G,3,evnt.selflag,evnt.ind,evnt.value)
                                end
                            case 'G'
                                updateselection(a.SIt,evnt.selflag,evnt.ind,evnt.value)
                        end
                   end
                   file_markchange(V);
                case 'slice'
                    % 'slice' event serves only to refresh the display(s)
                    % attached to the 'sliceinfo' object that emits it, but
                    % does not imply any other update in synchronized objects
                case {'sizes' 'sizesplus' 'units' 'nddata' 'labels'}
                otherwise
                    if fn_dodebug
                        fprintf('change in fn4D object: case ''%s'' not handled by display_changeview\nplease check\n',flag)
                        keyboard
                    end
            end
            [a.listeners.Enabled] = deal(true);
        end
        function display_slidertrial(V)
            if isempty(V.disppar.trialsublist)
                nmax = V.ntrial;
                k = V.content.ktrial;
            else
                nmax = length(V.disppar.trialsublist);
                k = find(V.disppar.trialsublist==V.ktrial);
                if ~isscalar(k), error programming, end
            end
            set(V.a4d.slidertrial,'max',nmax, ...
                'sliderstep',[1 1]./max(1,(nmax-1)), ...
                'value',k, ...
                'callback',@(hu,evnt)trialcallback(V,get(hu,'value')))
        end
        function hl = display_timeplot(V,D,slice)
            % SPECIAL FUNCTION: displays one selection in the time display
            % window this function is called by
            % activedisplayPlot/slicedisplay, D is an activedisplayPlot
            % object
            
            ischannel = ~isempty(slice.tag);
            if ismember('timevisu',V.modules)
                timevisu = V.panels.timevisucontrol.s;
                if isempty(timevisu.curvdist), timevisu.curvdist = 0; end
            else
                timevisu = struct('curvdist',0,'showfilter',false,'LS',[],'HS',[], ...
                    'threshold',[],'thrdisplay','curve','LS_eeg',[],'HS_eeg',[],'LS_elphy',[],'HS_elphy',[]);
            end
            
            % data: can be multiple
            hl = zeros(2,100); nl = 0; % pre-allocate more than necessary
            ndata = 1 + any([V.disppar.time V.disppar.timesignal]=='_');
            ncond = size(slice.dataop,2);
            for kdata=1:ndata
                % show signal and/or spike?
                if isempty(V.precomp.spikes{kdata})
                    spikes = [];
                    yfit = [];
                else
                    spikes = slice.(V.precomp.spikes{kdata});
                    yfit = slice.(V.precomp.spikefit{kdata});
                end
                showsignal = fn_switch(V.disppar.timespike, ...
                    {'spikesignal' 'spikeboth' 'spikefit'},true, ...
                    'spikespike',false, ...
                    'spikeneuropil',isempty(slice.spikepar) || isempty(slice.spikepar.fun));
                showspike = fn_switch(V.disppar.timespike, ...
                    'spikesignal',false, ...
                    {'spikespike' 'spikeboth' 'spikeneuropil' 'spikefit'},~isempty(spikes));
                showfit = strcmp(V.disppar.timespike,'spikefit');
                
                % vertical shift
                ydec = ((kdata-1)*timevisu.curvdist);
                
                if showsignal
                    % data
                    if ischannel
                        y = slice.dataop;
                    else
                        y = slice.(V.precomp.datamode{kdata});
                    end
                    ny = length(y);
                    
                    % tidx: use the one in slice?
                    if V.precomp.useslicetidx(kdata) || length(D.tidx)~=length(y)
                        tixs = slice.tidx(:);
                    else
                        tixs = D.tidx(:) / V.a4d.timedisplayfactor;
                    end
                    if isscalar(tixs)
                        deltat = 0;
                    else
                        deltat = diff(tixs(1:2)); % unit = second
                    end
                    
                    % filtering
                    LHS = fn_switch(lower(slice.tag), ...
                        '',     {timevisu.LS timevisu.HS}, ...
                        'eeg',  {timevisu.LS_eeg timevisu.HS_eeg}, ...
                        'electrophysiology',    {timevisu.LS_elphy timevisu.HS_elphy}, ...
                        {[] []});
                    if ~isempty([LHS{:}])
                        % note that we want the mean to be unchanged in order to
                        % still apply the same threshold -> hence the 'z' flag
                        ls = LHS{1}/deltat; if isempty(ls), ls=0; end
                        hs = LHS{2}/deltat; if isempty(hs), hs=0; end
                        if 0 %~ls
                            % special: Butterworth filter
                            [b a] = butter(1,min(1/hs,.999),'high');
                            yf = filter(b,a,y,b(2)*mean(y)); % assume a flat signal at mean(y) before y
                        elseif hs
                            yf = fn_filt(y,[ls hs],'mz');
                        else
                            yf = fn_filt(y,[ls hs],'m');
                        end
                    else
                        yf = y;
                    end
                    if timevisu.showfilter, y = yf; end
                    
                    % fft
                    if ~isempty(V.disppar.timefft)
                        % note that we have already set V.a4d.SIt to show
                        % accurate frequencies in display_updatesize(V)
                        switch V.disppar.timefft
                            case 'compute'
                                % compute fft
                                y    = abs(fft(y));
                                y(1) = 0;
                            case 'fromdata'
                                y = abs(y);
                        end
                        tixHz = (0:ny-1)/ny * (1/deltat); % unit = Hz
                        tix = tixHz / V.a4d.timedisplayfactor; % display unit for frequency
                    else
                        tix = tixs * V.a4d.timedisplayfactor; % display unit for time
                    end
                    
                    % display line at y=0 or y=1
                    if V.disppar.timeyzero && kdata==1
                        m = min(y); %M = max(y);
                        y0 = fn_switch(m(1)<=0 || m(1)>1, 0, 1);
                        if ~isempty(y0)
                            h = line(tix,y0*ones(length(tix),1), ...
                                'parent',D.ha,'tag','tpview-timeplot', ...
                                'linestyle','--');
                            hl(1,nl+1) = h;
                            hl(2,nl+1) = ydec;
                            nl = nl+1;
                        end
                    end
                    
                    % display data
                    if isempty(y)
                        h = line(tix([1 end]),[0 0],'parent',D.ha,'tag','tpview-timeplot')';
                    else
                        if issparse(y) 
                            if nnz(y)==0
                                tix = tix([1 end]);
                                y = [0; 0];
                            elseif nnz(y)<length(y)/10
                                % make the vectors significantly smaller, but
                                % this requires some work as adjacent zero
                                % values need to be included
                                [idx, ~, v] = find(y);
                                n = length(idx);
                                d = diff(idx);
                                idx = [idx-1 idx idx+1]';
                                values = zeros(3,n); values(2,:) = v';
                                ok = true(3,n);
                                ok(1,1) = (idx(1)>1);
                                ok(1,2:end) = (d>2);
                                ok(3,1:end-1) = (d>1);
                                ok(3,end) = (idx(end)<length(y));
                                tix = tix([1; idx(ok); end]);
                                y = [0; values(ok); 0];
                            end
                        end
                        h = line(tix,full(y),'parent',D.ha,'tag','tpview-timeplot')';
                    end
                    hl(1,nl+(1:ncond)) = h;
                    hl(2,nl+(1:ncond)) = ydec;
                    nl = nl+ncond;
                    
                    % RMS
                    if V.disppar.timerms
                        if V.disppar.timermsdcomp
                            rms = sqrt(var(diff(y))/2);
                        else
                            rms = [sqrt(var(diff(y))/2); sqrt(var(y))];
                        end
                        m = mean(y(:,1)); if m>1e-4, rms = rms/m; end
                        tpos = tix(ceil(end/40));
                        ypos = double(max(y(1:ceil(end/5),:)));
                        for kcond=1:ncond
                            str = ['RMS: ' num2str(rms(1,kcond),'%.2e')];
                            if size(rms,1)>1, str = [str ' / ' num2str(rms(2,kcond),'%.2e')]; end %#ok<AGROW>
                            hl(1,nl+kcond) = text(tpos,ypos(kcond),0,str, ...
                                'parent',D.ha,'fontweight','bold');
                        end
                        hl(2,nl+(1:ncond)) = ydec;
                        nl = nl+ncond;
                    end
                    
                    % thresholding
                    if ~isempty(timevisu.threshold) && ~ischannel
                        y0 = mean(yf(:));
                        thr = timevisu.threshold + fn_switch(y0<.5,0,1);
                        for kcond=1:ncond
                            ok = (yf(:,kcond)>thr);
                            if any(ok)
                                nl = nl+1;
                                switch timevisu.thrdisplay
                                    case 'curve'
                                        ythr = y(:,kcond);
                                        ythr(~ok) = NaN;
                                        hl(1,nl) = line(tix,ythr,'parent',D.ha,'linewidth',2,'tag','tpview-timeplot');
                                    case 'line'
                                        % make isolated spikes double
                                        yk = y(:,kcond);
                                        ok = ok | (ok([1 1:end-1]) & ~ok([1 2 1:end-2]));
                                        y0 = mean(yk); m = prctile(y(:,kcond),3);
                                        d = y0 + 1.5*(m-y0);
                                        ythr = ones(ny,1)*d;
                                        ythr(~ok) = NaN;
                                        hl(1,nl) = line(tix,ythr,'parent',D.ha,'linewidth',7,'tag','tpview-timeplot');;
                                    case 'spikes'
                                        yk = y(:,kcond);
                                        y0 = mean(yk);
                                        d = y0-prctile(yk,3);
                                        ypos = y0+[-1.4 -1.1]*d';
                                        thrspikes = fn_timevector(ok,deltat,'times');
                                        hl(1,nl) = fn_spikedisplay(thrspikes*V.a4d.timedisplayfactor, ...
                                            ypos,'parent',D.ha,'tag','tpview-timeplot');
                                end
                                hl(2,nl) = ydec(kcond);
                            end
                        end
                    end
                end
                
                % spikes
                if showspike
                    y = slice.(V.precomp.datamode{1})(:,1);
                    if strcmp(slice.tag,'electrophysiology')
                        y0 = max(y);
                        nl=nl+1;
                        hl(1,nl) = line(spikes*V.a4d.timedisplayfactor,ones(1,length(spikes))*y0, ...
                            'linestyle','none','marker','*', ...
                            'parent',D.ha,'tag','tpview-timeplot');
                    else
                        y0 = mean(y);
                        d = y0-prctile(y,3);
                        ypos = y0+[-1.4 -1.1]*d';
                        nl = nl+1;
                        hl(1,nl) = fn_spikedisplay(spikes*V.a4d.timedisplayfactor,ypos,'linewidth',2, ...
                            'parent',D.ha,'tag','tpview-timeplot');
                    end
                    hl(2,nl) = ydec;
                end
                if showfit && ~isempty(slice.spikepar) && ~isempty(yfit)
                    nl = nl+1;
                    hl(1,nl) = line(tix,yfit(:,1),'linewidth',2,'parent',D.ha,'tag','tpview-timeplot')';
                    hl(2,nl) = ydec;
                    if size(yfit,2)>=2
                        nl = nl+1;
                        hl(1,nl) = line(tix,yfit(:,2),'linestyle','--','parent',D.ha,'tag','tpview-timeplot')';
                        hl(2,nl) = ydec;
                    end
                end
            end
            
            % remove what was allocated too much
            hl = hl(:,1:nl);
            
            % color 
            % (channel data)
            if ischannel
                set(hl(1,:),'color',[1 1 1]*.5)
                for i=1:size(hl,2)
                    setappdata(hl(1,i),'color',[1 1 1]*.5)
                end
            end
            %             % (condition) TODO: this is really a hack, try to do it better
            %             if strcmp(D.linecol,'cat') && isvector(y)
            %                 col = fn_colorset(slice.kcond);
            %                 set(hl(1,:),'color',col)
            %                 for i=1:size(hl,2)
            %                     setappdata(hl(1,i),'color',col)
            %                 end
            %             end
        end
        function display_more(V,flag,im)
            persistent s0
            
            % full options
            dofull = strcmp(flag,'fulloptions');
            if dofull
                spec = struct( ...
                    'type',     {'frames' {'3d' 'frames' 'array' 'movie'} 'type'}, ...
                    'var1',     {'time' {'time' 'trial' 'condition'} 'third dimension'}, ...
                    'var1sel',  {[] 'double' 'sub-range'}, ...
                    'var2',     {'' {'' 'trial' 'condition'} 'fourth dimension'}, ...
                    'var2sel',  {[] 'double' 'sub-range'}, ...
                    'im',       {'im1' {'default' 'im1' 'im2'} 'other parameters as in'} ...
                    );
                if isempty(s0), s0 = spec(1); end
                spec = spec(2:3);
                s = fn_structedit(s0,spec);
                if isempty(s), return, end
                s0 = s;
            else
                s = struct( ...
                    'type',     flag, ...
                    'var1',     'time', ...
                    'var1sel',  [], ...
                    'var2',     '', ...
                    'var2sel',  [], ...
                    'im',       im ...
                    );
            end
            
            % display from which we will copy clipping and color map
            kim = fn_switch(s.im,'default',0,'im1',1,'im2',2);
            if kim
                D = V.a4d.(s.im);
                s.data = V.disppar.(s.im);
                if fn_ismemberstr(s.data,{'dataspikes' 'shotnoisespikes'})
                    errodlg 'not implemented yet'
                    return
                end
                s.cond = V.disppar.([s.im 'cond']);
            else
                errordlg 'not implemented yet'
                return
            end
            
            % data
            if ~dofull
                SImodel = ['SI' num2str(kim)];
                data = V.a4d.(SImodel).data;
            else
                switch [s.var1 s.var2]
                    case 'time'
                        data = getcondition(V.trial,s.data,s.cond);
                    case {'timecondition' 'condition'}
                        conds = fn_switch(s.var1,'time',s.var2sel,'condition',s.var1sel);
                        if isempty(conds)
                            condname = 'all conditions';
                        else
                            condname = fn_strcat(conds-1,'C',',C','');
                        end
                        data = getcondition(V.trial,s.data,condname);
                    case {'timetrial' 'trial'}
                        trials = fn_switch(s.var1,'time',s.var2sel,'trial',s.var1sel);
                        if isempty(trials)
                            trials = 1:V.ntrial;
                        end
                        ntrial = length(trials);
                        data = getcondition(V.content.trials(trials(1)),s.data,s.cond);
                        data(1,1,1,1,ntrial) = 0;
                        for i = 2:ntrial
                            data(:,:,:,1,i) = getcondition(V.content.trials(trials(i)),s.data,s.cond);
                        end
                    otherwise
                        errordlg('this case was not implemented yet!') % missing trialcondition and conditiontrial
                end
                if isempty(data), errordlg('no data'), return, end
            end
            % display
            switch s.type
                case {'3d' 'frames' 'framescond' 'array' 'linecut'}
                    Flag = fn_switch(s.type,'3d','3D','frames','Frames','framescond','Conditions','array','Array','linecut','Line cut');
                    hf = figure('Tag','2PVIEW display','menubar','none', ...
                        'numbertitle','off','integerhandle','off', ...
                        'color',get(V.hf,'color'), ...
                        'name',[V.skin ' - ' Flag]);
                    proj = fn_switch(s.var1,'time',[1 2 3],'condition',[1 2 4],'trial',[1 2 5]);
                    dimsplus = fn_switch(s.var2,'',[],'condition',4,'trial',5);
                    SI = projection(V.a4d.G,proj,'dimsplus',dimsplus,'data',data);
                    cmap = fn_switch(strcmp(D.cmap,'user'),D.cmapval,D.cmap);
                    switch s.type
                        case '3d'
                            activedisplay3D(SI,'in',hf, ...
                                'clipmode','data','clip',D.clip, ...
                                'cmap',cmap ...
                                );
                        case 'frames'
                            activedisplayFrames(SI,'in',hf, ...
                                'clipmode','data','clip',D.clip, ...
                                'cmap',cmap, ...
                                'ncol',5,'xbin',D.binning ...
                                );
                        case 'array'
                            activedisplayArray(SI,'in',hf, ...
                                'xbin',5 ...
                                );
                        case 'linecut'
                            set(V.a4d.(s.im),'shapemode','segment')
                            SNK = snake2D(SI,[]);
                            activedisplayImage(SNK,'in',hf, ...
                                'clipmode','data','clip',D.clip, ...
                                'cmap',cmap ...
                                );
                            msgbox(['Select a line segment in Picture ' num2str(kim)])
                    end
                    % possibility to read the current data displayed in
                    % tpview
                    if ~dofull
                        m = uimenu('parent',hf,'label','Auto-update: off','userdata',false, ...
                            'callback',@(h,e)display_more_update(V,SI));
                        V.disppar.morelist(end+1) = struct('SI',SI,'SImodel',SImodel,'control',m);
                    end
                case 'movie'
                    fn_movie(squeeze(data),'clip',D.clip,'cmap',D.cmap);
                otherwise
                    error('unknown flag ''%s'' for additional display',s.type)
            end
        end
        function display_more_update(V,SI)
            morelist = V.disppar.morelist;
            if isempty(morelist), return, end
            dodisconnect = (nargin==2 && ischar(SI) && strcmp(SI,'disconnect'));
            % update control
            if nargin==2 && ~dodisconnect
                morelist = morelist(SI==[morelist.SI]);
                autoupdate = ~get(morelist.control,'userdata');
                set(morelist.control,'userdata',autoupdate,'label',['Auto-update: ' onoff(autoupdate)])
            end
            % remove invalid displays (for example, when figure has been
            % closed)
            ok = isvalid([morelist.SI]) & ishandle([morelist.control]);
            badcontrols = [morelist(~ok).control];
            delete(badcontrols(ishandle(badcontrols)))
            morelist = morelist(ok);
            V.disppar.morelist = morelist;
            controls = [morelist.control];
            % 'disconnect' flag: set all auto-update values to off
            if dodisconnect
                set([morelist.control],'userdata',false,'label','Auto-update: off')
                return
            end
            % update display(s)
            autoupdate = cell2mat(fn_get(controls,'userdata','cell'));
            for i = row(find(autoupdate))
                set(morelist(i).SI,'data',V.a4d.(morelist(i).SImodel).data)
            end
        end
        function display_movie(V,flag)
            persistent framelist kframe nframe
            
            dorealtime = (nargin>=2) && strcmp(flag,'realtime');
            isoff = strcmp(V.timer.Running,'off');
            stop(V.timer)
            ht = V.grob.statusbar;
            set(ht,'string','')
            items = V.panels.moreitems;
            G = V.a4d.G; % doing V.a4d.G.ijkl2 = ... would cause V.a4d.G to be [] temporarily, and would cause an error
            p1only = get(items.p1,'value');
            
            % stop
            if ~get(items.movie,'value')
                if isoff, return, end % Movie was not on
                % note that 'Picture1 only' mode can disrupt the synchrony
                % between graphs: we force a display update
                G.ijkl2(3) = framelist(kframe);
                % also, the time selection might not be synchronized
                setselectiondims(V.a4d.G,3,V.a4d.SIt.selection.getselset(1))
                return
            end
            
            % time selection -> restric range
            if ~vide(V.a4d.SIt.selection)
                framelist = V.a4d.SIt.selection.getsel(1,1).dataind;
                % remove temporal selection in G so that individual (rather
                % than average) frames will be displayed
                [V.a4d.listeners.Enabled] = deal(false); % disconnect G from SIt (selection will remain in SIt)
                updateselection(V.a4d.G,3,'remove',1)
                [V.a4d.listeners.Enabled] = deal(true);
            else
                framelist = 1:V.nfr;
            end
            nframe = length(framelist);
            
            % set speed
            if dorealtime
                if isempty(V.trial.dt_sec)
                    errordlg 'set frame duration to play movie at real time speed'
                    return
                end
                speed = 1/V.trial.dt_sec; % use real-time frames/s
                set(V.panels.moreitems.speed,'value',log10(speed))
            else
                speed = 10^get(items.speed,'value');
            end
            str0 = sprintf('movie %.1ffr/s',speed);
            if ~isempty(V.trial.dt_sec)
                fact = num2str(speed*V.trial.dt_sec,'%.1g');
                if strcmp(fact,'1')
                    str0 = [str0 ' (realtime)'];
                else
                    str0 = [str0 ' (' fact ' x realtime)'];
                end
            end
            set(ht,'string',str0)
            
            % current frame
            t0 = []; 
            countprev = 0;
            kframe = find(framelist>=G.ijkl(3),1,'first');
            if isempty(kframe), kframe=1; end
            
            % set and start
            set(V.timer,'Period',round(1e3/speed)/1e3, ...
                'ExecutionMode','FixedRate','BusyMode','drop', ...
                'timerfcn',@moviestep)
            start(timer('TimerFcn',@(u,e)start(V.timer),'StartDelay',0.1)) %#ok<CPROP> % needed to avoid errors (let Matlab breathe)
            
            % local definition of moviestep function
            function moviestep(~,~)
                % first exec: just init counter
                if isempty(t0)
                    t0 = tic;
                    return
                end
                % frame
                count = round(toc(t0)*speed);
                if count-countprev==1
                    set(ht,'string',str0)
                else
                    set(ht,'string',[str0 ' [missing frames!]'])
                end
                if ~p1only
                    % current frame might have been changed by user
                    kframe = find(framelist>=G.ijkl(3),1,'first');
                    if isempty(kframe), kframe=1; end
                end
                kframe = fn_mod(kframe + (count-countprev),nframe);
                % update display
                if p1only
                    % only change image in Picture 1!
                    frame = V.a4d.SI1.data(:,:,framelist(kframe));
                    set(V.a4d.im1.img,'cdata',permute(frame,[2 1 3]))
                else
                    G.ijkl2(3) = framelist(kframe);
                end
                countprev = count;
                drawnow
            end
        end
        function display_trialsublist(V,flag)
            switch flag
                case 'set'
                    types = {V.trial.stimtable.table.type};
                    [types{fn_isemptyc(types)}] = deal('');
                    s = struct( ...
                        'label',    {'' 'label' 'filter trials by:'}, ...
                        'trials',   {1:V.ntrial     'double'    'trials'}, ...
                        'type',     {'' unique(['' types]) 'type'}, ...
                        'id',       {unique([V.content.trials.stimid]) 'double' 'id(s)'}, ...
                        'normal',   {true 'logical' 'keep normal trials'}, ...
                        'rejected', {true 'logical' 'keep rejected trials'}, ...
                        'special',  {true 'logical' 'keep special trials'} ...
                        );
                    s = fn_structedit(s);
                    if isempty(s)
                        return
                    else
                        ok = true(1,V.ntrial);
                        for i=1:V.ntrial
                            Ti = V.content.trials(i);
                            if ~isempty(s.trials)
                                ok(i) = ok(i) && ismember(i,s.trials);
                            end
                            if isscalar(Ti.stimid)
                                % filtering on id and type only for
                                % single-condition trials
                                if ~isempty(s.id)
                                    ok(i) = ok(i) && ismember(Ti.stimid,s.id);
                                end
                                if ~isempty(s.type)
                                    ok(i) = ok(i) && strcmp(Ti.stimdetails.type,s.type);
                                end
                            end
                            ok(i) = ok(i) ...
                                && fn_switch(Ti.status,'n',s.normal,'r',s.rejected,'s',s.special);
                        end
                        if ~any(ok), errordlg('No trial satisfies condition'), return, end
                        V.disppar.trialsublist = find(ok);
                        if ~ok(V.ktrial), V.ktrial = V.disppar.trialsublist(1); end
                    end
                case 'reset'
                    V.disppar.trialsublist = [];
                otherwise
                    error('unknown flag ''%s''',flag)
            end
            display_slidertrial(V)
        end
    end
    methods (Access='private')
        function display_header(V,ktrial)
            if nargin<2
                ktrial = V.content.ktrial;
            end
            trial = V.content.trials(ktrial);
            
            % header content
            header = struct('ktrial',ktrial,'type',trial.type,'size',trial.sizes);
            dt = trial.dt_sec;
            if ~isempty(dt), header.frametime = fn_switch(dt<1,[num2str(trial.dt_sec*1000) 'ms'],[num2str(dt) 's']); end
            stimd = trial.stimdetails;
            if isfield(stimd,'name') && isscalar(stimd) && ~isempty(stimd.name)
                header.stim = stimd.name;
            else
                header.stim = [stimd.id];
            end
            if ~isempty(trial.addinfo)
                F = fieldnames(trial.addinfo);
                for k=1:length(F)
                    f = F{k};
                    header.(f) = trial.addinfo.(f);
                end
            end
            
            % background color indicates trial status and whether data is
            % saved
            switch trial.status
                case 'n' % normal trial
                    col = fn_switch(isempty(trial.file),'w','default');
                case 'r' % rejected trial
                    col = fn_switch(isempty(trial.file),[1 .7 .7],'r');
                case 's' % special trial
                    col = fn_switch(isempty(trial.file),[.9 .9 .4],[.9 .7 .5]);
            end
            set(V.grob.header,'backgroundcolor',col)
            
            % prepare title
            if isempty(trial.file)
                folder = ''; ftrial = ''; ext = '';
            else
                [d ftrial ext] = fileparts(trial.file(1,:));
                [d folder] = fileparts(d);
            end
            if size(trial.file,1)>1, ftrial = [ftrial '+...']; end
            ftrial = [ftrial ext];
            if ~isempty(trial.preprocess.op), ftrial = [ftrial ' - corrected']; end
            changedflag = fn_switch(V.savingpar.chg,'*','');
            if ~isempty(V.savingpar.savename)
                [d fmat] = fileparts(V.savingpar.savename);
                [d folder] = fileparts(d);
                str1 = {folder [fmat '.tpv' changedflag] ftrial};
            else
                str1 = {[' ' changedflag] folder ftrial};
            end
            
            % prepare information
            F = fieldnames(header);
            nf = length(F);
            str2 = {};
            for k=1:nf
                f = F{k};
                % special cases
                switch f
                    case 'PMT1'
                        % show two PMT values together to save space
                        f = 'PMTs';
                        a = sprintf('%.0f (%.0f)',header.PMT1,header.PMT2);
                    case 'PMT2'
                        continue
                    otherwise
                        a = header.(f);
                end
                switch class(a)
                    case 'char'
                        val = a;
                    case 'double'
                        if isvector(a) && length(a)<=5
                            val = num2str(a(:)');
                        else
                            s = size(a);
                            val = sprintf('[%ix%i',s(1),s(2));
                            for i=3:length(s), val = [val 'x' num2str(s(i))]; end %#ok<AGROW>
                            val = [val ' double]']; %#ok<AGROW>
                        end
                    otherwise
                        s = size(a);
                        val = sprintf('[%ix%i',s(1),s(2));
                        for i=3:length(s), val = [val 'x' num2str(s(i))]; end %#ok<AGROW>
                        val = [val ' ' class(a) ']']; %#ok<AGROW>
                end
                str2{end+1} = [f ': ' val]; %#ok<AGROW>
            end
            
            % display
            set(V.grob.header,'string', ...
                [str1 {' '} str2])
            
        end
        function chgsiz = display_updatesize(V)
            a = V.a4d;
            
            % does the new data look similar to the one currently displayed
            chgsiz = any(V.sizes(1:2)~=a.G.sizes(1:2));
            oldscale = a.G.grid(1:2,1);
            
            % geometric information
            type = V.type;
            if V.content.timeline, type = 'linescan'; end
            T = V.trial;
            xunit = T.xunit;
            
            % conversion between frame / second / actual unit displayed
            % (get the frame time in second... or attempt to do it)
            tlabel = 'time';
            tunit = 's';
            fact = 1;                           % conversion from second to second
            deltat = T.dt_sec;                  % conversion from frame to second
            if isempty(deltat)
                if ~strcmp(T.tunit,'frame')
                    errordlg(['unknown time unit ''' T.tunit ''', please update code'])
                end
                deltat = V.internpar.defaultdt;
                tunit = 'frames';
                fact = 1/deltat;                % conversion from second to frame
            end
            t0 = T.t0 * (deltat/T.dt);          % t0 in seconds
            if strcmp(type,'zstack')
                tlabel = 'z';
                tunit = xunit;
                fact = abs(T.dz)/deltat;        % factor for conversion of second to xunit
            end
            % (special: linescan)
            if strcmp(type,'linescan')
                deltat = deltat/T.ny;
            end
            % (display another unit than s?)
            if strcmp(tunit,'s') 
                tstop = T.dt_sec*T.nfr;
                if tstop<1 || (deltat<=.002 && tstop<10)
                    tunit = 'ms';
                    fact = 1e3;     % factor for conversion of second to millisecond
                elseif tstop>180*60
                    tunit = 'h';
                    fact = 1/3600;  % factor for conversion of second to hour
                elseif tstop>600
                    tunit = 'min';
                    fact = 1/60;    % factor for conversion of second to minute
                end
            end
            % (store factor)
            V.a4d.timedisplayfactor = fact;     % conversion from second to display unit
            tscale = deltat*fact;               % conversion from frame to display unit
            t0 = t0*fact;
            
            % set units and scales
            % (starting values)
            switch type
                case {'movie' 'zstack'}
                    SIproj = 1:2;
                    scale = [T.dx T.dy tscale]';
                    translate = [-T.dx/2 -T.dy/2 t0-tscale];
                    labels = {'x' 'y' tlabel};
                    units = {xunit xunit tunit};
                    siz = T.sizes(1:3);
                    seldims = 'xy';
                case 'linescan'
                    SIproj = [1 3];
                    scale = [T.dx 1 tscale]';
                    translate = [-T.dx/2 -1/2 t0-tscale];
                    labels = {'x' '' tlabel};
                    units = {xunit '' tunit};
                    siz = [T.nx 1 T.ny*T.nfr];
                    seldims = 'x';
            end
            % (special: fft)
            if ~isempty(V.disppar.timefft)
                flabel = 'f';
                funit = fn_switch(tunit,'s','Hz','ms','kHz',[tunit '^-1']);
                fs = 1/tscale; % sampling frequency in display unit
                maxfreq = fs;
                fscale = maxfreq/(siz(3)-1); % go from 0Hz to maxfreq in nt steps
                [tlabel tunit tscale] = deal(flabel,funit,fscale);
                t0 = 0;
                if strcmp(V.disppar.timefft,'fromdata')
                    [labels{3} units{3} scale(3)] = deal(flabel,funit,fscale);
                    translate{3} = -tscale;
                end                
            end
            
            % update display
            [a.listeners.Enabled] = deal(false);
            set([a.im1 a.im2],'seldims',seldims) % important to do it first: this will reset selections if needed
            set(a.G,'sizes',siz,'mat',{scale translate}) % 'left side' is on zero
            set([a.SI1 a.SI2],'proj',SIproj)
            a.F.labels(1:3) = labels;
            a.F.units(1:3)  = units;
            set(a.SIt,'grid',[tscale t0-tscale],'labels',{tlabel},'units',{tunit})
            [a.listeners.Enabled] = deal(true);
        end
        function display_updateselection(V)
            % this function updates the selection displayed according to
            % the information stored about signals
            
            
            a = V.a4d;
            [a.listeners.Enabled] = deal(false);
            signal = V.content.signals(1);
            if V.content.nsel==0
                sel = selectionND.empty(1,0);
            else
                sel = [signal.x(V.ktrial,:).sel];
            end
            % replace selection set in G
            scanline = (strcmp(V.type,'linescan') || V.content.timeline);
            if scanline
                siz = V.nx;
                seldimsnum = 1;
            else
                siz = [V.nx V.ny];
                seldimsnum = [1 2];
            end
            selsetnew = selectionset(siz,sel);
            setselectiondims(a.G,selsetnew.t.dims,selsetnew,'all')
            
            [a.listeners.Enabled] = deal(true);
        end
        %         function display_scaleselection(V,zoomratio)
        %             error 'function not valid any more'
        %             % change selection in G and in SI if scale has changed
        %             if all(zoomratio==1), return, end
        %             a = V.a4d;
        %             % TODO: update selectionset.selsetaffinity method so that it
        %             % can accept an affinityND object instead of a matrix to define
        %             % the affinity; then, there is no need to do the tricky things
        %             % below
        %             t = a.SI.selection.getselset([1 2]).singleset;
        %             scal = affinityND('scale2D',zoomratio);
        %             t = selaffinity(t,scal);
        %             selsetnew = selectionset(a.SI.sizes,2,t);
        %             setselectiondims(a.G,1:2,selsetnew,'all')
        %             setselection(a.SI,selsetnew,'all')
        %         end
        function display_image(V,recompute,chgclipflag,imgidx)
            if nargin<2, recompute=false; end
            if nargin<3, chgclipflag=false; end
            if nargin<4, imgidx=[1 2]; end
            if isequal(imgidx,3), imgidx=[1 2]; end
            a = V.a4d;
            % which data is displayed
            D = struct('in',{},'data',{},'cond',{},'spk',{});
            p = V.disppar;
            if any(imgidx==1)
                dataflag = fn_strcut(p.im1,'_'); % e.g. data becomes {data}, data_sfr becomes {data sfr}
                D = [D struct('in',1,'data',dataflag,'cond',p.im1cond,'spk',p.im1spikecomb)];
            end
            if any(imgidx==2)
                dataflag = fn_strcut(p.im2,'_'); % e.g. data becomes {data}, data_sfr becomes {data sfr}
                D = [D struct('in',2,'data',dataflag,'cond',p.im2cond,'spk',p.im2spikecomb)];
            end
            idata = {find([D.in]==1) find([D.in]==2)};    % 1x2 cell array: which part of this data for image 1 and for image 2
            ndata = length(D);
            % set data
            if recompute
                % compute the needed data
                imdata = cell(1,ndata);
                scanline = strcmp(V.trial.type,'linescan') || V.content.timeline;
                for k=1:ndata
                    dk = D(k);
                    % show spikes
                    showspike = fn_ismemberstr(dk.data,{'spikes' 'shotnoisespikes'});
                    combspike = dk.spk;
                    if ~showspike
                        imdata{k} = getcondition(V.trial,dk.data,dk.cond);
                        if isempty(imdata{k})
                            imdata{k} = zeros(V.nx,V.ny,fn_switch(scanline,V.nfr,1),'uint8');
                        elseif scanline && strcmp(dk.data,'data0')
                            imdata{k} = repmat(imdata{k},[1 1 V.nfr]);
                        end
                    end
                    if showspike || combspike
                        kcond = cell2mat(tps_readconditionname(dk.cond));
                        spikedataflag = fn_switch( ...
                            any(strfind(dk.data,'shotnoise')),'shotnoise', ...
                            any(strfind(dk.data,'sfr')),'', ...
                            'data');
                        if showspike
                            if isscalar(kcond) && ~isempty(spikedataflag)
                                imdata{k} = spikemovie(V.content,V.ktrial,kcond,spikedataflag);
                            else
                                imdata{k} = zeros(V.nx,V.ny,'uint8');
                            end
                        elseif isscalar(kcond) && ~isempty(spikedataflag)
                            % combine imdata{k} with spikes
                            imdata{k} = spikemovie(V.content,V.ktrial,kcond,spikedataflag,imdata{k});
                        end
                    end
                    % scan line? -> reshape
                    if V.content.timeline || strcmp(V.trial.type,'linescan')
                        [nx ny nt] = size(imdata{k});
                        imdata{k} = reshape(imdata{k},[nx 1 ny*nt]);
                    end
                end
                % update display
                for kdisplay=1:2
                    if ~ismember(kdisplay,imgidx), continue, end
                    if isscalar(idata{kdisplay})
                        % we need to separate this case, because doing 
                        % 'datak = cat(6,imdata{idata{kdata}});' would
                        % stupidly duplicate a potentially large array!
                        datak = imdata{idata{kdisplay}};
                    else
                        datak = cat(6,imdata{idata{kdisplay}}); % dimension 6 corresponds to multi-channel
                    end
                    SI = fn_switch(kdisplay,1,V.a4d.SI1,2,V.a4d.SI2);
                    setdata(SI,'data',datak)
                end
            end
            if isscalar(chgclipflag), chgclipflag(2)=chgclipflag; end
            if any(imgidx==1) && chgclipflag(1), autoclip(a.im1), end
            if any(imgidx==2) && chgclipflag(2), autoclip(a.im2), end
        end
        function display_time(V,flag,ind,value)
            % this functions updates both:
            % - the signals information storage (V.content.signals(1))
            %   [note that this part can be erroneous if trials are of
            %   different sizes, since the selection data sizes are not
            %   updated!]
            % - the display (by changing V.a4d.SIt.slice)
            
            if nargin<2, flag=''; end % no change in selections - display only
            SIt = V.a4d.SIt;
            
            % precompute some variables (see display_timeplot function)
            V.precomp.scanline = strcmp(V.type,'linescan');
            V.precomp.datamode{1} = fn_switch(V.disppar.timesignal, ...
                'signalop','dataop', ...
                'data');
            V.precomp.datamode{2} = fn_switch(V.disppar.timesignal, ...
                'signal_signalop','dataop', ...
                'signalop','data2op', ...
                'data2');
            V.precomp.spikes = fn_switch(V.disppar.timesignal, ...
                'signal_signalop',{'' 'spikes'}, ...
                {'spikes' 'spikes2'} ...
                );
            V.precomp.spikefit = fn_switch(V.disppar.timesignal, ...
                'signal_signalop',{'' 'spikefit'}, ...
                {'spikefit' 'spikefit2'} ...
                );
            V.precomp.useslicetidx = fn_switch( ...
                strcmp(V.type,'zstack') || strcmp(V.disppar.timesignal,'signal'), [0 0], ...
                strcmp(V.disppar.timesignal,'signal_signalop'), [0 1], ...
                [1 1]);
            
            % special case: flag 'redisplay' indicates that the content of
            % the slice has not changed (what might have changed however is
            % the parameters specific to function display_timeplot)
            if strcmp(flag,'redisplay')
                displaysignals(V.a4d.time)
                return
            end
            
            % go!
            if V.disppar.timealltrials
                % current selection, all trials: probably need to recompute
                % everything
                cursel = V.disppar.currentsel;
                if isempty(cursel)
                    % no selection -> no display
                    SIt.slice(:)=[];
                    title(V.grob.time,'')
                else
                    ktrials = V.disppar.trialsublist;
                    if isempty(ktrials), ktrials=0; end % use 0 to say 'all trials'
                    slice = getslice(V.content,ktrials,cursel);
                    
                    % show only trials with same stimulation as current
                    % trial
                    if V.disppar.timealltrialsfilterstim
                        if ktrials
                            ntr = length(ktrials);
                            T = V.content.trials(ktrials);
                        else
                            ntr = V.ntrial;
                            T = V.content.trials;
                        end
                        ok = false(1,ntr);
                        id = V.trial.stimid;
                        if ~isscalar(id)
                            if fn_dodebug, disp 'improve this', end
                            id = id(1); 
                        end
                        for k=1:ntr
                            idk = T(k).stimid;
                            ok(k) = isscalar(idk) && idk==id;
                        end
                        slice = slice(ok);
                    end
                    
                    % note that there is no sense in displaying
                    % channel/recording data of one specific trial
                    if length(SIt.slice)==length(slice) && isequal(size(slice(1).data,2),SIt.sizesplus)
                        updateslice(SIt,'change',1:length(slice),slice);
                    else
                        SIt.slice = slice;
                    end
                    
                    % change the color on time courses to that of the current selection
                    D = V.a4d.time;
                    %if isnumeric(D.linecol), D.linecol = fn_colorset(cursel); end
                    if strcmp(flag,'reset'), V.a4d.time.autolineposition('reset'), end
                    
                    % title
                    ttl = fn_fileparts(V.savingpar.savename,'base');
                    if ~isempty(ttl), ttl = [ttl ', ']; end
                    if V.disppar.timealltrialsfilterstim
                        if isempty(V.disppar.trialsublist)
                            ttl = [ttl 'All trials for condition ID ' num2str(id) ', '];
                        else
                            ttl = ['Subselection of trials for condition ID ' num2str(id) ', '];
                        end
                    else
                        if isempty(V.disppar.trialsublist)
                            ttl = [ttl ' All trials, '];
                        else
                            ttl = [ttl 'Subselection of trials, '];
                        end
                    end
                    ttl = [ttl 'region ' num2str(cursel)];
                    title(V.grob.time,ttl,'interpreter','none')
                end
            else
                % all selections, current trial
                switch flag % TODO: case of no selection!!
                    case {'' 'remove' 'reset' 'all'}
                        slice = getslice(V.content,V.ktrial);
                        r = V.disppar.timerecording;
                        if any([r.val])
                            names = {r([r.val]).name};
                            slice = [slice getrecording(V.content,names)];
                        end
                        if length(SIt.slice)==length(slice) ...
                                && fn_sizecompare(size(slice(1).data),[SIt.sizes SIt.sizesplus])
                            updateslice(SIt,'change',1:length(slice),slice);
                        else
                            SIt.slice = slice;
                            if strcmp(flag,'reset'), V.a4d.time.autolineposition('reset'), end
                        end
                    case {'change' 'add' 'affinity' 'new'}
                        slice = getslice(V.content,V.ktrial,ind);
                        if strcmp(flag,'new')
                            nnew = length(ind);
                            nadd = sum([V.disppar.timerecording.val]);
                            ntot = V.content.nx + nadd;
                            % remove the time courses of a single point
                            if ntot<length(SIt.slice)+nnew
                                updateslice(SIt,'remove',1)
                            end
                            % add new time courses
                            updateslice(SIt,'new',[],slice)
                            % reorder to maintain recordings below
                            if nadd>0
                                npre = V.content.nx - nnew;
                                updateslice(SIt,'reorder',[],[1:npre npre+nadd+(1:nnew) npre+(1:nadd)])
                            end
                        else
                            updateslice(SIt,'change',ind,slice)
                        end
                    case 'active'
                        updateslice(SIt,flag,ind,value)
                    case 'reorder'
                        % the permutation might need to take into account
                        % the additional recording time courses
                        nval = length(value);
                        nslice = length(SIt.slice);
                        if nslice>nval, value(nval+1:nslice) = nval+1:nslice; end
                        updateslice(SIt,flag,ind,value)
                    case 'indices'
                        % can happen when creating external objects
                    otherwise
                        error programming
                end
                
                % title
                ttl = fn_fileparts(V.savingpar.savename,'base');
                if ~isempty(ttl), ttl = [ttl ', ']; end
                ttl = [ttl V.disppar.timecond];
                title(V.grob.time,ttl,'interpreter','none')
            end
        end
        function display_stimandevents(V)
            %             if ~isvalid(V), return, end % happens at program exit, because Matlab is resetting all objects before deleting the, which causes the ylim listener to react
            ha = V.grob.time;
            % initialization: re-display each time ylim will change
            if ~isappdata(ha,'tpview_display_stimandevents_update')
                hl = addlistener(ha,'YLim','PostSet',@(h,e)display_stimandevents(V));
                setappdata(ha,'tpview_display_stimandevents_update',hl)
            end
            % delete previous display
            delete(findobj(ha,'tag','tpview_display_stim'))
            delete(findobj(ha,'tag','tpview_display_event'))
            % display stims
            if V.disppar.showstim
                stim = getstimbycond(V.trial,V.disppar.timecond);
                stim = stim * V.a4d.timedisplayfactor; % convert seconds to time display units
            else
                stim = [];
            end
            if ~isempty(stim)
                tps_displaystim(stim,ha, ...
                    'color',[1 1 1]*.8,'tag','tpview_display_stim','hittest','off')
            end
            % display events
            cond = singlecond(V.disppar.timecond);
            if V.disppar.showevents && isscalar(cond)
                eventlist = V.trial.eventlist;
                nevent = length(eventlist);
                for k=1:nevent
                    okcond = (eventlist(k).cond==cond);
                    eventk = eventlist(k).time(okcond) * V.a4d.timedisplayfactor;
                    if isempty(eventk), continue, end
                    tps_displaystim(eventk,ha, ...
                        'color',fn_colorset(k),'tag','tpview_display_event', ...
                        'buttondownfcn',@(h,e)display_editevents(V,k,h),'hittest',onoff(V.disppar.editevents))
                end
            end
        end
        function display_editevents(V,k,hobj)
            eventlist = V.trial.eventlist;
            switch get(V.hf,'selectiontype')
                case 'open'
                    % create new
                    cond = singlecond(V.disppar.timecond);
                    if ~isscalar(cond), return, end
                    k = 1;
                    if isempty(eventlist)
                        eventlist(1).name = '';
                        eventlist(1).cond = [];
                        eventlist(1).time = [];
                    end
                    p = get(V.grob.time,'currentpoint');
                    t = p(1) / V.a4d.timedisplayfactor;
                    eventlist(k).cond(:,end+1) = cond;
                    eventlist(k).time(:,end+1) = t;
                case 'normal'
                    % move
                    dp = fn_moveobject(hobj);
                    dt = dp(1) / V.a4d.timedisplayfactor;
                    eventlist(k).time = eventlist(k).time + dt;
                case 'alt'
                    % delete
                    eventlist(k) = [];
            end
            V.trial.eventlist = eventlist;
            display_stimandevents(V)
        end
    end
    
    % GENERAL DISPLAY
    methods
        function chgframepositions(V,varargin)
            % after positions have been changed, update the panel positions
            chgframepositions@interface(V,varargin{:})
            if ~isempty(V.panels), init_panels(V), end % V.panels is empty when chgframepositions is called the first time by interface_end()
        end
        function tpview_paneltoggle(V,evnt)
            % evnt can be a structure with fields OldValue and NewValue
            % and values control handles, or the name of the panel to show
            select = V.panels.select;
            % names of previous and new panels
            if isstruct(evnt) || isobject(evnt)
                kold = (evnt.OldValue==select.allobj);
                oldname = select.allname{kold};
                knew = (evnt.NewValue==select.allobj);
                newname = select.allname{knew};
            else
                oldobj = get(V.grob.panelselect,'SelectedObject');
                kold = (oldobj==select.allobj);
                oldname = select.allname{kold};
                newname = evnt;
                set(select.(newname),'value',1)
            end
            % change focus
            if strcmp(oldname,newname), return, end
            tpview_panelhide(V,oldname)
            tpview_panelshow(V,newname)
        end
        function tpview_reinit_4d(V)
            init_4d(V)
            init_slidertrial(V)
            display_changeview(V,'data&signals',true)
        end
        function tpview_redisplay(V)
            chgclip = false;
            display_changeview(V,'data&signals',chgclip)
        end
        function tpview_toggledisplay(V)
            im1 = V.grob.im1;
            im2 = V.grob.im2;
            time = V.grob.time;
            pos1 = get(im1,'pos');
            pos2 = get(im2,'pos');
            pos3 = get(time,'pos');
            switch V.disppar.toggle
                case 'normal'
                    set(im1,'pos',pos3)
                    set(im2,'pos',pos1)
                    set(time,'pos',pos2)
                    V.disppar.toggle = 'image';
                case 'image'
                    set(im1,'pos',pos2)
                    set(im2,'pos',pos3)
                    set(time,'pos',pos1)
                    V.disppar.toggle = 'normal';
                otherwise
                    error programming
            end
        end
        function tpview_keypress(V,e)
            % specials
            if ischar(e)
                switch e
                    case 'init'
                        objs = [V.hf findobj(V.hf,'type','uicontrol','-not','style','edit')'];
                        objs = setdiff(objs,V.panels.filesitems.list);
                        set(objs,'keyPressFcn',@(h,evnt)tpview_keypress(V,evnt))
                        return
                    case 'shortcutslist'
                        return
                    otherwise
                        error('unknown flag ''%s''',e)
                end
            end
            
            % action
            switch(e.Key)
                case {'pageup' 'period' 'pagedown' 'slash' 'backquote'}
                    % change trial
                    step = fn_switch(e.Key,{'pageup' 'period'},-1,1);
                    list = V.disppar.trialsublist;
                    if isempty(list), list = 1:V.ntrial; end
                    kcur = find(V.ktrial==list);
                    kcur = min(length(list),max(1,kcur+step));
                    V.ktrial = list(kcur);
                case {'rightarrow' 'leftarrow' 'uparrow' 'downarrow'}
                    % move selections
                    if V.content.nsel==0, return, end
                    if isempty(e.Modifier)
                        dx = 1;
                    elseif strcmp(e.Modifier{1},'control')
                        dx = .25;
                    elseif strcmp(e.Modifier{1},'shift')
                        dx = 4;
                    elseif strcmp(e.Modifier{1},'command')
                        % Mac 'apple' key
                        return
                    else
                        error('unknown key modifier ''%s''',e.Modifier{1})
                    end
                    dp = dx * fn_switch(e.Key, ...
                        'rightarrow',[1 0], ...
                        'leftarrow', [-1 0], ...
                        'uparrow',   [0 -1], ...
                        'downarrow', [0 1] ...
                        );
                    % use low-level function directly to the sliceinfo
                    % object to define translation
                    SI1 = V.a4d.SI1;
                    mov = affinityND('translate2D',dp);
                    ijmov = AX2IJ(SI1,mov);
                    updateselection(SI1,'affinity',1:V.content.nsel,ijmov)
                case {'0' 'control' 's'}
                    % ignore
                case 'qwfe' % 'x' already used to compute spikes!!
                    % action button
                    if strcmp(get(V.grob.actionbutton,'enable'),'on')
                        tpv_actionbutton
                    end
                case 'r'
                    % reject trial
                    action(V,'rejecttrialtoggle')
                case {'c' 'd'}
                    % move to next/previous cell - note that it is also
                    % possible to click on the 'TRIALS' button
                    step = fn_switch(e.Key,'c',1,'d',-1);
                    setpar(V,'currentsel',fn_mod(V.disppar.currentsel+step,V.content.nsel))
                otherwise
                    if fn_dodebug, disp(e), end
            end
        end
    end
    methods (Access='private')
        function tpview_panelshow(V,name)
            set(V.panels.select.(name),'foregroundcolor',zeros(1,3))
            set(V.panels.(name),'visible','on')
        end
        function tpview_panelhide(V,name)
            set(V.panels.select.(name),'foregroundcolor',ones(1,3)*.4)
            set(V.panels.(name),'visible','off')
        end
    end
    
    % MISC
    methods
        function access(V) %#ok<MANU>
            keyboard
        end
        function autorepair(V,doall)
            if nargin<2, doall = false; end
            % clear memory
            if doall
                clearmemory(V.content)
            end
            % make sure some content properties are set correctly
            autorepair(V.content) % clean all signals (erasedata(V.content))
            setdatamode(V.content,V.disppar.time)
            setdatacond(V.content,V.disppar.timecond)
            setdataopdef(V.content,V.panels.datacontrol.x, ...
                get(V.panels.dataitems.linked,'value'))
            % repair 4D
            if doall
                init_4d(V)
            else
                a = V.a4d;
                % repair listeners
                delete(a.listeners)
                a.listeners = event.listener(a.G,'ChangeView', ...
                    @(hs,evnt)display_changeview(V,evnt,'G'));
                a.listeners(2) = event.listener(a.SIt,'ChangeView', ...
                    @(hs,evnt)display_changeview(V,evnt,'SIt'));
                V.a4d = a;
                % solve problem with selections
                nsel = length(a.SI1.selectionmarks);
                if V.content.nsel~=nsel || V.content.nx~=fn_switch(nsel,nsel,1)
                    % problem in the number of selections -> cancel all
                    % selections
                    updateselection(a.G,[],'reset')
                    updateselection(a.SI1,'reset')
                    updateselection(a.SI2,'reset')
                    updateselection(a.SIt,'reset')
                    updateselection(V.content,'reset')
                end
            end
            init_slidertrial(V)
            % update display
            delete(findobj(V.grob.time,'tag','tpview-timeplot')) % these are lines created just before an error occured (without error, tag should have been replaced by fn4D_lines)
            if doall
                display_changeview(V,'data&signals',true)
            else
                display_changeview(V,'signals',false)
                displaysignals(a.time)
            end
            % remove watch
            set(V.hf,'pointer','arrow')
        end
        function preferences(V,flag,varargin)
            % function preferences(V,'init') loads preferences at init
            % function preferences(V,'edit') open GUI to edit preferences
            % function preferences(V,'reload') set new preferences
            
            % update preferences
            [s spec] = tpv_internpar.default;
            switch flag
                case {'init' 'reload'}
                    s = fn_structmerge(s,V.options.preferences,'skip');
                case 'edit'
                    s = fn_structedit(V.options.preferences,spec);
                otherwise
                    error programming
            end
            if isempty(s), return, end
            if ~isscalar(s), error programming, end
            V.options.preferences = s;
            if ~strcmp(flag,'reload'), saveoptions(V), end
            
            % update parameters according to preferences
            s = V.options.preferences;
            doinit = strcmp(flag,'init'); % at init, not all actions are needed
            c = watch(V); %#ok<NASGU>
            % (display)
            if ~doinit
                set(V.hf,'defaultuicontrolfontsize',s.fontsize, ...
                    'defaultaxesfontsize',V.options.preferences.textfontsize, ...
                    'defaulttextfontsize',V.options.preferences.textfontsize)
                set(findobj(V.hf,'type','uicontrol'),'fontsize',s.fontsize)
                set(findobj(V.hf,'type','axes'),'fontsize',s.textfontsize)
                set(findobj(V.hf,'type','text'),'fontsize',s.textfontsize)
            end
            % (general behavior)
            V.internpar.memorysize = s.memorysize;
            V.internpar.floattype = s.floattype;
            memorypool.setmaxmem(num2str(s.memorysize,'%.1fGB'))
            if doinit, V.internpar.loaddata = s.loaddata; else setpar(V,'loaddata',s.loaddata), end
            if doinit, V.content.seldotrial = s.seldotrial; end
            % (data handling)
            V.internpar.defaultdt = s.defaultdt;
            % (navigation)
            if ~doinit
                a = V.a4d;
                set([a.im1 a.im2],'navigation',s.imnav,'scrollwheel',s.imscroll)
                set(a.time,'navigation',s.tcnav,'scrollwheel',s.tcscroll)
            end
            if ~doinit, end
        end
        function loaddefaultoptions(V)
            % load default options
            loaddefaultoptions@interface(V)
            % update what needs to be updated
            init_menus(V)
            init_panels(V)
            preferences(V,'reload')
        end
    end
    
    % USER
    methods
        function action(V,flag,varargin)
            % write here new methods with no need to re-compile the class
            switch flag
                case 'changeview'
                    display_changeview(V,varargin{:});
                case 'rejecttrialtoggle'
                    switch V.trial.status
                        case 'n'
                            V.trial.status = 'r';
                        case 'r'
                            V.trial.status = 'n';
                        case 's'
                            errordlg 'cannot reject a special trial'
                            return
                    end
                    display_header(V)
                case 'addtriggertrial'
                    kframetrig = cell(1,V.ntrial);
                    T = V.content.trials;
                    for k=1:V.ntrial
                        if isempty(T(k).eventlist), continue, end
                        ttrig = T(k).eventlist(1).time;
                        kframetrig{k} = round(ttrig/T(k).dt);
                    end
                    alltrig = [kframetrig{:}];
                    ntrig = length(alltrig);
                    nbef = min(alltrig);
                    naft = T(1).nfr-max(alltrig);
                    data = 0;
                    for k = 1:V.ntrial
                        for i=1:length(kframetrig{k})
                            cond = T(k).eventlist.cond(i);
                            ktrig = kframetrig{k}(i);
                            data = data + T(k).data(:,:,ktrig-nbef+1:ktrig+naft,cond)/ntrig;
                        end
                    end
                    setdata(V,data,'add')
                    V.trial.eventlist = struct('name','','cond',1,'time',[nbef*V.trial.dt_sec; 0]);
                    display_stimandevents(V)
                otherwise
                    error('unknown action flag ''%s''',flag)
            end
        end
        function updatedisplay(V)
            display_changeview(V,'data&signals',false)
        end
        function updatedataop(V,trials)
            % function updatedataop(V,trials|'all'|'current')
            %---
            % force dataop to be recomputed for specified trials, or all of
            % them, or the current trial [default]
            if nargin<2
                trials = V.ktrial;
            elseif ischar(trials)
                trials = fn_switch(trials,'all',1:V.ntrial,'current',V.ktrial);
            end
            erasedata(V.content.trials(trials))
            erasedata(V.content,trials)
            
            % update display
            imgidx = ~isempty(strfind(V.disppar.im1,'op')) + ~isempty(strfind(V.disppar.im2,'op'));
            dotime = any(strfind(V.disppar.time,'op'));
            display_changeview(V,'datamode',imgidx,dotime,false)
        end
        function setdata(V,varargin)
            % setdata(V,data[,sfr][,'set|add|rep']) [default is 'set']
            % setdata(V,data[,sfr],'rep',ind)
            % setdata(V,sfr,'sfr'[,ind])
            % setdata(V,'perm',perm)
            % setdata(V,'rm',ind)
            %---
            % data can be either a numerical array, or a tps_trial object
            % sfr must be a numerical array
            
            if nargin==1 || isequal(varargin{1},'help'), help tpview.setdata, return, end
            
            % input
            T = tps_trial.empty; data = []; sfr = []; sfrflag = false; flag = 'set'; ind = [];
            k = 1;
            while k<nargin
                a = varargin{k};
                if isa(a,'tps_trial')
                    T = a;
                elseif isa(a,'tpv_content')
                    data_loadfile(V,a)
                    return
                elseif isnumeric(a) && isvector(a)
                    ind = a;
                elseif isnumeric(a)
                    if isempty(data)
                        data = a;
                    else
                        sfr = a;
                    end
                elseif iscell(a)
                    if length(a)==2
                        data = a{1};
                        if numel(a)==2, sfr = a{2}; end
                    elseif isempty(data)
                        data = a;
                    else
                        sfr = a;
                    end
                elseif ischar(a)
                    switch a
                        case 'sfr'
                            sfrflag = true;
                            flag = 'rep';
                        case {'set' 'add' 'rep' 'perm' 'rm'}
                            flag = a;
                        otherwise
                            if exist(a,'file')
                                data_loadfile(V,a)
                                return
                            else
                                error('unknown flag ''%s''',a)
                            end
                    end
                else
                    error argument
                end
                k = k+1;
            end
            
            switch flag
                case {'set' 'add' 'rep'}
                    if isempty(T)
                        % data/sfr
                        if size(data,4)>1, data = num2cell(data,1:3); end
                        if ~iscell(data), data = {data}; end
                        s = size(data{1});
                        if any(s(1:2)~=V.trial.sizes(1:2))
                            disp('warning: data size does not fit current data')
                        end
                        if ~isempty(sfr)
                            if size(sfr,4)>1, sfr = num2cell(sfr,1:3); end
                            if ~iscell(sfr), sfr = {sfr}; end
                            if numel(sfr)~=numel(data) || any(size(sfr{1})~=s), error('data and sfr differ in size'), end
                        else
                            sfr = cell(1,numel(data));
                        end
                        % sfr only?
                        if sfrflag
                            sfr = data;
                            data = {V.content.trials(ind).data};
                        end
                        % create tps_trial object
                        for k=1:numel(data)
                            T(k) = tps_trial(data{k},sfr{k},V.trial);
                            T(k).addinfo.description = 'user-defined';
                        end
                    end
                    % indices
                    if strcmp(flag,'rep')
                        if ~isempty(ind)
                            if length(ind)~=numel(sfr), error('number of indices does not fit number of data'), end
                        elseif numel(sfr)==1
                            ind = V.ktrial;
                        elseif numel(sfr)==V.ntrial;
                            ind = 1:V.ntrial;
                        else
                            error('indices to replace are not defined');
                        end
                    end
                    % add trials
                    ntrialold = V.ntrial;
                    if isempty(V.content.signal.x)
                        % program was just started!
                        settrials(V.content,T,'data','all conditions')
                        display_changeview(V,'data')
                    else
                        addtrials(V.content,T)
                        V.disppar.trialsublist = []; % hey fais pas chier avec le filtrage des trials!!
                        % rearrange trials
                        switch flag
                            case 'add'
                                V.content.ktrial = ntrialold+1;
                            case 'set'
                                rmtrial(V.content,1:ntrialold)
                                V.content.ktrial = 1;
                            case 'rep'
                                indfinal = ntrialold+(1:length(ind));
                                perm = 1:ntrialold;
                                perm(ind) = indfinal;
                                perm = [perm ind];
                                permutetrials(V.content,perm)
                                rmtrial(V.content,indfinal)
                        end
                        % update display
                        display_changeview(V,'chgtrial')
                    end
                case 'permute'
                    perm = ind(:)';
                    if ~isequal(sort(perm),1:V.ntrial), error('argument is not a permutation'), end
                    permutetrials(V.content,perm)
                    V.disppar.trialsublist = [];
                    display_changeview(V,'chgtrial')
                case 'rm'
                    rmtrial(V.content,ind)
                    V.disppar.trialsublist = []; % list became obsolete because trials have been assigned new numbers
                    display_changeview(V,'chgtrial')
                otherwise
                    error('unknown action flag ''%s''',flag)
            end
        end
        function x = getsignals(V,varargin)
            % function x = getsignals(V)            -> all signals
            % function x = getsignals(V,'current')  -> currently displayed signals
            % function x = getsignals(V,k)          -> kth currently displayed signal 
            % function x = getsignals(V,ktrial,ksel)-> specify trial and sel number
            %                                          (use [] for all trials or all cells) 
            
            if nargin>1 && strcmp(varargin{1},'help'), help tpview.getsignals, return, end
            
            switch length(varargin)
                case 0
                    ktrial = []; ksel = [];
                case 1
                    k = varargin{1};
                    if fn_ismemberstr(k,{'cur' 'current'}), k = []; end
                    if V.disppar.timealltrials
                        ktrial = k; ksel = V.disppar.currentsel;
                    else
                        ktrial = V.ktrial; ksel = k;
                    end
                case 2
                    [ktrial ksel] = deal(varargin{:});
                otherwise
                    error argument
            end
            if isempty(ktrial), ktrial = 1:V.ntrial; end
            if isempty(ksel), ksel = 1:V.content.nx; end
            x = getslice(V.content,ktrial,ksel);
        end
    end
end

% SUB-FUNCTIONS

function ok = validcondition(conds,nc)

if ischar(conds), conds = {conds}; end
ncond = length(conds);
ok = true(1,ncond);
for k=1:ncond
    if isempty(conds{k})
        warning 'empty condition, please check this' %#ok<WNTAG>
        ok(k) = false;
        continue
    end
    if strcmp(conds{k},'all conditions'), continue, end % this special case is ok
    match = regexp(conds{k},'\d*','match');
    if isempty(match) && fn_dodebug, error programming, end
    ok(k) = all(str2num(char(match))<nc);
end

end

function markers = SpikeStyles(n)

markers = '*d^po';
nmark = length(markers);
markers = markers(mod(0:n-1,nmark)+1);

end

function mousescroll(V,n)
% if several trials, scroll trials, otherwise scroll time
if V.ntrial>1
    list = V.disppar.trialsublist;
    if isempty(list), list = 1:V.ntrial; end
    kslider = find(V.ktrial==list);
    kslider = min(length(list),max(1,kslider+n));
    V.ktrial = list(kslider);
else
    SIt = V.a4d.SIt;
    SIt.ij2 = min(SIt.sizes,max(1,SIt.ij2 + n));
end
end

function trialcallback(V,k)
if isempty(V.disppar.trialsublist)
    V.ktrial = k;
else
    V.ktrial = V.disppar.trialsublist(k);
end
end

function tpview_reinit(V)

% keep content
okcontent = (V.ntrial>1) || ~isequal(V.data,0);
if okcontent
    content = V.content;
end
% restart the program
V = tpview(V.skin);
% reload content
if okcontent
    data_loadfile(V,content)
end

end

function cond = singlecond(condname)

cond = regexp(condname,'^C(\d+)$','tokens');
if isempty(cond)
    cond = [];
else
    cond = str2double(cond{1});
end

end

function c = watch(V)

% start watch
c = fn_watch(V.hf);
% progress in status bar
if isfield(V.grob,'statusbar'), fn_progress('in',V.grob.statusbar), end

end

function x = trialaverage(T,dataflag,floattype,doprogress)

nfr = min([T.nfr]);
ntrial = length(T);
x = zeros(T(1).nx,T(1).ny,nfr,T(1).nc,floattype);
if doprogress, fn_progress('trial',ntrial), end
for k=1:ntrial
    if doprogress, fn_progress(k), end
    xk = fn_float(T(k).(dataflag));
    if T(k).nfr>nfr, xk = xk(:,:,1:nfr,:); end
    x = x+xk; clear xk
end
x = x/ntrial;

end

% BUG REPORT
function bugreport()

email = 'thomas.deneux@unic.cnrs-gif.fr';
subject = '[Optimage bug report]';
version = fn_readtext('tpview_version');
body = {'' 'Dear Thomas,' '' ...
    ['I have encountered the following bug with OptImage version ' version '.'] '' ...
    'Matlab error message:' '[please copy-paste red lines from the Matlab command]:' '' ...
    'This occured in the following circumstances:' '[please explain briefly]' '' ...
    'Thanks'};
for i=1:length(body), if isempty(body{i}), body{i}=' '; end, end
web(['mailto:' email '?subject=' subject fn_strcat(body,'&body=','&body=')])

end

