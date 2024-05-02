classdef activedisplayPlot < fn4Dhandle  %#ok<*MCSUP>
    
    properties 
        axis = [0 0];
        clip = [0 0];
        clipmode = 'data';
        movelinegroup = 'sel';  % 'sel', 'cat', or 'none'
        moveyonly = true;
        linecol = 2; % a dimension, 0 for 'category'
        colorset = fn_colorset('plot12'); % a n*3 array, can be set however with a string that either is accepted as a color set name by fn_colorset, or is a color map
        scaledisplay = ''; % 'ybar' or ''
        movescale = true;
        baseline = 'auto'; % 'mean', '0', '1' or 'auto'
        slicedisplayfun = @slicedisplaydefault; % hl = slicedisplayfun(D,slice) - for 'buttondownfcn' property, use app-def data 'ksel' and 'icat' - don't use 'color' property but create 'color' app-def data - don't use 'tag' property
        menustayvisible = false;
        usercallback
        lineusercallback
    end
    properties (SetObservable=true)
        scrollwheel = 'x'; % '' ['off' tolerated], 'x', 'y' or 'xy' ['on' tolerated]
        navigation = 'zoom'; % 'zoom' or 'pan'
    end
    
    properties (SetAccess='private')
        ha
        hf
        scalebar % scale bar and text
        tidx = 0;
    end
    
    properties (Access='private')
        hu
        ha2
        slider
        buttons
        hplot = {}; % nsel x (nlinpersel) cell array of vector arrays
        yxdec = {};  % nsel x (3*nlinpersel) cell array of arrays - 1 = y, set in slicedisplay; 2,3 = x,y, set in moveline and reset in autolineposition
        activeidx = []; % nsel vector array
        hline
        htsel
        hdeco = struct('lines',[]);
        menu        
        menuitems
        ystep
        oldclip = [0 1];
        oldselectiontype = 'normal';
        autopos = 1;  % vector of slice dimensions that are dispatched; use 0 for categories (i.e. several lines that appear in a single slice item)
    end
    
    properties (Dependent, SetAccess='private')
        signals
    end
    
    properties (Dependent)
        autolinepos  % reflects D.autopos
    end
    
    properties (SetAccess='private')
        CL
        SI
        C2D
        listenaxpos
        listenparentcolor
    end
    
    % Constructor and destructor
    methods
        function D = activedisplayPlot(varargin)
            fn4D_dbstack

            % options for initialization
            opt = struct( ...
                'clip',         [], ...
                'in',           [], ...
                'ystep',        []);
            if nargin==0 || ~isobject(varargin{1})
                D.SI = sliceinfo(1);
                for i=1:3, D.SI.slice(i) = struct('active',true,'data',rand(20,2)); end
                [opt optadd] = fn4D_parseInput(opt,varargin{:});
            else
                D.SI = varargin{1};
                [opt optadd] = fn4D_parseInput(opt,varargin{2:end});
            end
            
            % type check
            if D.SI.nd~=1 || D.SI.ndplus>2
                error('activedisplayPlot class can display only one-, two-, or three-dimensional data slices')
            end
            
            % figure and axes
            if ~isempty(opt.in)
                if ~ishandle(opt.in) && mod(opt.in,1)==0 && opt.in>0
                    figure(opt.in)
                end
                switch get(opt.in,'type')
                    case 'figure'
                        D.hf = opt.in;
                        figure(opt.in)
                        D.ha = gca;
                    case 'axes'
                        D.ha = opt.in;
                        D.hf = get(opt.in,'parent');
                    otherwise
                        error('bad handle')
                end
            else
                D.ha = gca;
                D.hf = get(D.ha,'parent');
            end
            cla(D.ha,'reset')
            if isempty(get(D.hf,'Tag')), set(D.hf,'Tag','used by fn4D'), end
            
            % axes options
            set(D.ha,'NextPlot','ReplaceChildren','YLimMode','manual')
            
            % time marker line
            D.hline = line('parent',D.ha,'color','black', ...
                'HitTest','on','ButtonDownFcn',@(hl,evnt)movebar(D,hl));
            
            % scale bar (and listener for re-positioning upon axes resize)
            D.scalebar(1) = line('Parent',D.ha,'Color','k','visible','off', ...
                'linewidth',3);
            D.scalebar(2) = text('Parent',D.ha,'Color','k','visible','off', ...
                'horizontalalignment','right','verticalalignment','middle');
            D.listenaxpos = fn_pixelsizelistener(D.ha,D,@(h,evnt)displayscalebar(D));
            fn4D_enable('off',D.listenaxpos)
            set(D.scalebar,'DeleteFcn',@(h,evnt)fn4D_enable('on',D.listenaxpos))
            
            % buttons (bottom-up or local) - previous to displaygrid which
            % will change axis and hence try to move the controls
            D.hu = uicontrol('parent',D.hf,'CallBack',@(u,evnt)chgtime(D,'-1'),'String','<');
            fn_controlpositions(D.hu(1),D.ha,[0 1], [0 -17 18 18]);
            D.hu(2) = uicontrol('parent',D.hf,'CallBack',@(u,evnt)chgtime(D,'+1'),'String','>');
            fn_controlpositions(D.hu(2),D.ha,[0 1], [18 -17 18 18]);
            D.hu(3) = uicontrol('parent',D.hf,'CallBack',@(u,evnt)autolineposition(D,'toggle'),'String','A');
            fn_controlpositions(D.hu(3),D.ha,[1 1], [-35 -17 18 18]);
            D.hu(4) = uicontrol('parent',D.hf,'CallBack',@(u,evnt)set(D,'linecol','toggle'),'String','C');
            fn_controlpositions(D.hu(4),D.ha,[1 1], [-17 -17 18 18]);
            D.hu(5) = uicontrol('parent',D.hf,'CallBack',@(u,evnt)autolineposition(D,'+'),'String','+');
            fn_controlpositions(D.hu(5),D.ha,[1 1], [-17 -35 18 18]);
            D.hu(6) = uicontrol('parent',D.hf,'CallBack',@(u,evnt)autolineposition(D,'-'),'String','-');
            fn_controlpositions(D.hu(6),D.ha,[1 1], [-17 -53 18 18]);
            set(D.hu,'fontsize',8)
            
            % axes on the side serve to generate "outside" mouse events
            D.ha2 = axes('parent',D.hf,'buttondownfcn',@(u,e)Mouse(D,'outsidedown'));
            fn_controlpositions(D.ha2(1),D.ha,[0 0 1 0], [0 -15 0 15]);
            D.ha2(2) = axes('parent',D.hf,'buttondownfcn',@(u,e)Mouse(D,'outsideleft'));
            fn_controlpositions(D.ha2(2),D.ha,[0 0 0 1], [-15 0 15 0]);
            D.ha2(3) = axes('parent',D.hf,'buttondownfcn',@(u,e)Mouse(D,'outsideboth'));
            fn_controlpositions(D.ha2(3),D.ha,[0 0 0 0], [-15 -15 15 15]);
            set(D.ha2,'handlevisibility','off')
            D.listenparentcolor = connect_listener(D.hf,D,'Color','PostSet',@(u,e)updateSecondaryAxesColor());
            updateSecondaryAxesColor()
            function updateSecondaryAxesColor()
                try %#ok<TRYNC>
                    col = get(D.hf,'color');
                    set(D.ha2,'color',col, ...
                        'xtick',[],'xcolor',col,'ytick',[],'ycolor',col)
                end
            end
            uistack(D.ha2,'bottom')
            
            % sliders
            D.slider = fn_slider('parent',D.hf,'mode','area', ...
                'visible','off','scrollwheel','on', ... % set 'scrollwheel' after 'visible' so that scroll map won't need to be recalculated
                'callback',@(u,evnt)chgzoom(D,'slider',u));
            fn_controlpositions(D.slider(1),D.ha,[0 1 1 0], [36 -8 -70 9]);
            D.slider(2) = fn_slider('parent',D.hf,'mode','area', ...
                'visible','off','scrollwheel','on', ...
                'callback',@(u,evnt)chgclip(D,'slider',u));
            fn_controlpositions(D.slider(2),D.ha,[1 0 0 1], [-8 0 9 -62]);
            
            % more buttons
            D.buttons = uicontrol('Parent',D.hf, ...
                'backgroundcolor',[.5 0 0],'foregroundcolor',[.5 0 0]);
            fn_controlpositions(D.buttons(1),D.ha,[1 1], [-13 -65 12 12]);
            %             fn_controlpositions(D.buttons(1),D.ha,[1 1], [-17 -62 9 9]);
            %             D.buttons(2) = uicontrol('Parent',D.hf,'hittest','off', ...
            %                 'backgroundcolor',[.8 .8 0],'foregroundcolor',[.8 .8 0]);
            %             fn_controlpositions(D.buttons(2),D.ha,[1 1], [-8 -62 9 9]);
            %             D.buttons(3) = uicontrol('Parent',D.hf, ...
            %                 'backgroundcolor',[.8 .8 0],'foregroundcolor',[.8 .8 0]);
            %             fn_controlpositions(D.buttons(3),D.ha,[1 1], [-1 -62 2 9]);
            set(D.buttons,'style','frame','enable','off')
            initlocalmenu(D)
            %             initlocalmenu(D.SI,D.buttons(3))
            
            % clipping
            if ~isempty(opt.clip), D.clip = opt.clip; end 
            D.oldclip = D.clip;
            D.clipmode = 'data';
            
            % ystep
            if ~isempty(opt.ystep), D.ystep = opt.ystep; end
            
            % display
            displaygrid(D)
            displaydata(D,'all')
            displaylabels(D)
            displayzoom(D)
            displayline(D)
            displayselection(D)
            displaydecoration(D)

            % axes event (bottom-up)
            set(D.ha,'ButtonDownFcn',@(ha,evnt)Mouse(D))
            fn_scrollwheelregister(D.ha,@(n)Scroll(D,n),~isempty(D.scrollwheel))

            % communication with parent
            addparent(D,D.SI)
            
            % trick to make reset of axes trigger object deletion
            line('parent',D.ha,'visible','off','deletefcn',@(x,y)delete(D))            
            
            % set more properties
            if ~isempty(optadd)
                set(D,optadd{:})
            end
        end
        
        function initlocalmenu(D)
            fn4D_dbstack
            hb = D.buttons(1);
            delete(get(hb,'uicontextmenu'))
            m = uicontextmenu('parent',D.hf);
            D.menu = m;
            set(hb,'uicontextmenu',m)
            
            % keep menu visible
            info.pin = uimenu(m,'label','keep menu visible','checked',onoff(D.menustayvisible), ...
                    'callback',@(u,e)set(D,'menustayvisible',~D.menustayvisible));
            function xx(fun,varargin)
                if D.menustayvisible, set(D.menu,'visible','on'), end
                feval(fun,varargin{:})
            end
            
            % clipping and y offset
            m1 = uimenu(m,'label','clipping mode','separator','on');
            info.clip(1) = uimenu(m1,'label','slice', ...
                'callback',@(u,e)xx(@()set(D,'clipmode','slice')));
            info.clip(2) = uimenu(m1,'label','data', ...
                'callback',@(u,e)xx(@()set(D,'clipmode','data')));
            info.clip(3) = uimenu(m1,'label','link1', ...
                'callback',@(u,e)xx(@()set(D,'clipmode','link1')));
            info.clip(4) = uimenu(m1,'label','link2', ...
                'callback',@(u,e)xx(@()set(D,'clipmode','link2')));            
            b = fn_ismemberstr({'slice','data','link1','link2'},D.clipmode);
            set(info.clip(b),'checked','on') 
            info.usrclip = uimenu(m,'label','user clip...', ...
                'callback',@(u,e)xx(@()set(D,'clip',fn_input('clip',D.clip))));
            m1 = uimenu(m,'label','baseline');
            info.baseline(1) = uimenu(m1,'label','auto', ...
                'callback',@(u,e)xx(@()set(D,'baseline','auto')));
            info.baseline(2) = uimenu(m1,'label','signal mean', ...
                'callback',@(u,e)xx(@()set(D,'baseline','mean')));
            info.baseline(3) = uimenu(m1,'label','0', ...
                'callback',@(u,e)xx(@()set(D,'baseline','0')));
            info.baseline(4) = uimenu(m1,'label','1', ...
                'callback',@(u,e)xx(@()set(D,'baseline','1')));
            info.baseline(5) = uimenu(m1,'label','specific time...', ...
                'callback',@(u,e)xx(@()set(D,'baseline','meantime')));
            baseflag = D.baseline; if iscell(baseflag), baseflag = baseflag{1}; end
            b = fn_ismemberstr({'auto','mean','0','1','meantime'},baseflag);
            set(info.baseline(b),'checked','on') 
            
            % resets
            info.sel(1) = uimenu(m,'label','reset signals display','separator','on', ...
                'callback',@(u,e)xx(@()displaydata(D,'all')));
            info.sel(2) = uimenu(m,'label','reset selection display', ...
                'callback',@(u,e)xx(@()displayselection(D)));
            info.sel(3) = uimenu(m,'label','reset selection', ...
                'callback',@(u,e)xx(@()updateselection(D.SI,'reset')));
            
            % line moving
            info.moveline(1) = uimenu(m,'label','move lines separately','separator','on', ...
                'callback',@(u,e)xx(@()set(D,'movelinegroup','')));      
            info.moveline(2) = uimenu(m,'label','move lines per selection', ...
                'callback',@(u,e)xx(@()set(D,'movelinegroup','sel')));      
            info.moveline(3) = uimenu(m,'label','move lines per category', ...
                'callback',@(u,e)xx(@()set(D,'movelinegroup','cat')));  
            info.moveyonly = uimenu(m,'label','move only in y','checked',onoff(D.moveyonly), ...
                'callback',@(u,e)xx(@()set(D,'moveyonly',~D.moveyonly)));
            b = fn_ismemberstr({'','sel','cat'},D.movelinegroup);
            set(info.moveline(b),'checked','on') 
            
            % navigation
            fn_propcontrol(D,'scrollwheel', ...
                {'menu' {'' 'xy' 'x' 'y'} {'off' 'on' 'x-axis only' 'y-axis only'}}, ...
                {'parent',m,'label','scroll wheel','separator','on'});
            fn_propcontrol(D,'navigation', ...
                {'menu' {'zoom' 'pan'}}, ...
                {'parent',m});
            
            % decorations
            % (features sub-menu)
            m1 = uimenu(m,'label','features','separator','on');
            uimenu(m1,'label','vertical line', ...
                'callback',@(u,e)xx(@()fn_setpropertyandmark(D.hline,'visible',u,'toggle')))
            uimenu(m1,'label','xgrid', ...
                'callback',@(u,e)xx(@()fn_setpropertyandmark(D.ha,'xgrid',u,'toggle')))
            uimenu(m1,'label','ygrid', ...
                'callback',@(u,e)xx(@()fn_setpropertyandmark(D.ha,'ygrid',u,'toggle')))
            uimenu(m1,'label','ticks...', ...
                'callback',@(u,e)xx(@()gridspacing()))
            uimenu(m1,'label','ticks and grid color...', ...
                'callback',@(u,e)xx(@()gridcolor()))           
            function gridspacing
                s = fn_structedit(struct( ...
                    'x',    {[] 'double' ['x (' D.SI.units{1} ')']}, ...
                    'y',    {[] 'double' 'y (y units)'}, ...
                    'ynum', {false 'logical' 'y: line numbers'} ...
                    ));
                if isempty(s), return, end
                if s.x
                    set(D.ha,'xtick',s.x*(ceil(D.tidx(1)/s.x):floor(D.tidx(end)/s.x)))
                elseif s.x==0 % i. e., s.x is not empty but equal to 0
                    set(D.ha,'xtickmode','auto')
                end
                if s.ynum && ~isempty(D.autopos)
                    if D.autopos(1)==0
                        nl = max(fn_map(@length,D.hplot)); % organize first by categories
                    else
                        nl = size(D.SI.slice,D.autopos(1));
                    end
                    ytick = -nl:-1;
                    yticklabel = fn_num2str(-ytick,'cell');
                    set(D.ha,'ytick',ytick,'yticklabel',yticklabel)
                elseif s.y
                    if ~isempty(D.autopos)
                        ytick = 0:s.y/D.ystep:D.oldclip(2);
                        yticklabel = cell(1,length(ytick));
                        oky = ismember(ytick,1:D.oldclip(2)-1);
                        yticklabel(oky) = fn_num2str(ytick(oky),'cell');
                        set(D.ha,'ytick',ytick,'yticklabel',yticklabel) %,'ytickmode','auto')
                    else
                        ytickstart = ceil(D.oldclip(1)/s.y)*s.y;
                        set(D.ha,'ytick',ytickstart:s.y:D.oldclip(2),'yticklabelmode','auto');
                    end
                elseif s.y==0 % the last alternative would be s.y empty -> nothing to do
                    set(D.ha,'ytickmode','auto')
                end
            end
            function gridcolor
                col = uisetcolor(get(D.ha,'xcolor'));
                set(D.ha,'xcolor',col,'ycolor',col)
            end
            % (scale bar)
            info.scalebar = uimenu(m,'label','scale bar','checked',onoff(strcmp(D.scaledisplay,'ybar')), ...
                'callback',@(u,e)xx(@()set(D,'scaledisplay',onoff(D.scaledisplay,'','ybar',''))));

            uimenu(m,'label','duplicate in new figure','separator','on', ...
                'callback',@(u,e)xx(@()duplicate(D)));
            uimenu(m,'label','display object ''D'' in base workspace', ...
                'callback',@(u,e)xx(@()assignin('base','D',D)));
            uimenu(m,'label','export signals to Matlab...', ...
                'callback',@(u,e)xx(@()exportsignals(D)))
            
            % programming
            uimenu(m,'label','reinit menu','separator','on', ...
                'callback',@(u,e)xx(@()initlocalmenu(D)))

            D.menuitems = info;
        end
        function delete(D)
            cla(D.ha,'reset')
            obj = [D.hu D.ha2 D.buttons D.menu];
            delete(obj(ishandle(obj)))
            delete(D.slider(isvalid(D.slider)))
            delete(D.listenaxpos(ishandle(D.listenaxpos)))
            delete(D.listenparentcolor(ishandle(D.listenparentcolor)))
        end
    end
    
    % Display
    methods (Access='private')
        function displayscalebar(D)
            fn4D_dbstack
            if ~strcmp(D.scaledisplay,'ybar')
                set(D.scalebar,'visible','off')
                fn4D_enable('off',D.listenaxpos)
                return
            else
                set(D.scalebar,'visible','on')
                fn4D_enable('on',D.listenaxpos)
            end
            % find a nice size for bar: in specialized function
            barsize = BarSize(D);
            % label - no unit yet
            label = num2str(barsize);
            % size in axes
            if ~isempty(D.autopos) && ~isempty(D.ystep)
                barsize = barsize/D.ystep; 
            end                
            % positions
            barorigin = fn_coordinates(D.ha,'b2a',[55 12]','position');
            barpos = [barorigin barorigin+[0 barsize]'];
            textpos = mean(barpos,2) + ...
                fn_coordinates(D.ha,'b2a',[-5 0]','vector');
            % set properties
            set(D.scalebar(1),'xdata',barpos(1,:),'ydata',barpos(2,:))
            set(D.scalebar(2),'position',textpos,'string',label)
            if D.movescale
                set(D.scalebar,'hittest','on','buttondownfcn', ...
                    @(hobj,evnt)fn_moveobject(D.scalebar)) %,'latch'))
            else
                set(D.scalebar,'hittest','off')
            end
        end
        
        function x = BarSize(D)
            % minimal/maximal size (25 pix) in axes coordinates
            xmin = fn_coordinates(D.ha,'b2a',[0 25],'vector');
            xmin = xmin(2);
            
            % desired size in axes coordinates
            ysiz = diff(get(D.ha,'ylim'));
            if ~isempty(D.autopos)
                x = D.ystep*max(xmin,min(2/3,ysiz/5));
                if isempty(x)
                    % no signal displayed, scale does not mean anything
                    x = .1;
                end
            else
                x = max(xmin,ysiz/5);
            end
            
            % round to a nice value
            x10 = 10^floor(log10(x));
            x = x / x10;
            vals = [1 2 5];
            f = find(x+1e-3>=vals,1,'last'); % +1e-3 is to avoid error due to numerical approximation
            x = vals(f) * x10;            
        end
        
        function displaygrid(D)
            fn4D_dbstack
            D.tidx = IJ2AX(D.SI,1:D.SI.sizes(1));
            D.axis = [D.tidx(1) D.tidx(end)];
        end
           
        function displaylabels(D)
            fn4D_dbstack
            str = D.SI.labels{1};
            if ~isempty(D.SI.units{1})
                str = [str ' (' D.SI.units{1} ')'];
            end
            xlabel(D.ha,str)
        end
        
        function displaydata(D,flag,ind,value)
            fn4D_dbstack
            %XX minimal update!
            
            slice = D.SI.slice;
            
            % no data -> no display
            if isempty(slice)
                delete(findobj(D.ha,'Tag','fn4D_line'))
                D.hplot = {}; D.yxdec = {}; D.activeidx = []; 
                return
            end

            nsel = numel(slice);
            if nsel==0, error programming, end
            
            % specific actions
            if isempty(flag), flag='all';
            elseif strcmp(flag,'remove') 
                % remove all?
                if isempty(ind)
                    % nothing to do
                    return
                elseif length(ind)==length(D.hplot)
                    flag='all'; ind=1; 
                end
            elseif strcmp(flag,'changeall')
                flag = 'change';
            elseif strcmp(flag,'new') && ~isscalar(ind)
                flag = 'all';
            end
            
            % display time courses
            switch flag
                case {'all','reset'}
                    delete(findobj(D.ha,'Tag','fn4D_line'))
                    D.hplot = cell(1,nsel); D.yxdec = cell(1,nsel);
                    for k=1:nsel
                        [D.hplot{k} dec] = slicedisplay(D,slice(k));
                        D.yxdec{k} = [dec; zeros(2,length(dec))];
                        slicenum(D,k)
                    end
                    D.activeidx = cumsum(logical([slice.active]));
                    autolineposition(D,fn_switch(flag,'all','set','reset','reset'))
                case 'new'
                    if ind==1 && ~isempty(D.hplot), delete(D.hplot{1}), end % start of selection mode
                    [D.hplot{ind} dec] = slicedisplay(D,slice(ind));
                    D.yxdec{ind} = [dec; zeros(2,length(dec))];
                    slicenum(D,ind)
                    D.activeidx(ind) = sum(logical([slice(1:ind).active]));
                    autolineposition(D,'set',ind)
                case {'change','add','kl','user'}
                    if nargin<3 || isempty(ind), ind = 1:nsel; end
                    for k=ind
                        hold = D.hplot{k}; delete(hold(ishandle(hold)))
                        [D.hplot{k} dec] = slicedisplay(D,slice(k));
                        % smart change of D.yxdec: keep information about
                        % curves moved by user
                        nl = length(dec);
                        if nl<size(D.yxdec{k},2)
                            D.yxdec{k} = D.yxdec{k}(:,1:nl);
                        elseif nl>size(D.yxdec{k},2)
                            D.yxdec{k}(:,end+1:nl)=0;
                        end
                        D.yxdec{k}(1,:) = dec;
                        slicenum(D,k) % set the color
                    end
                    autolineposition(D,'set',ind) % set the position
                case 'active'
                    % update activeidx
                    prevactiveidx = D.activeidx;
                    D.activeidx = cumsum(logical([slice.active]));
                    % line visibility
                    for k=ind
                        set(D.hplot{k},'visible',onoff(slice(k).active))
                    end
                    % re-dispatch iff no line has been moved by user
                    allyxdec = [D.yxdec{:}];
                    if ~any(row(allyxdec([2 3],:)))
                        % re-dispatch
                        autolineposition(D)
                    else
                        % keep same positions,
                        % update yxdec according to change into activeidx
                        ddy = D.activeidx - prevactiveidx;
                        for k=find(ddy)
                            D.yxdec{k}(3,:) = D.yxdec{k}(3,:) - ddy(k);
                        end
                    end
                case 'reorder'
                    perm = value;
                    D.hplot = D.hplot(perm);
                    D.yxdec  = D.yxdec(perm);
                    for k=1:nsel
                        slicenum(D,k)
                    end
                    prevactiveidx = D.activeidx;
                    D.activeidx = cumsum(logical([slice.active]));
                    % re-dispatch iff no line has been moved by user
                    allyxdec = [D.yxdec{:}];
                    if ~any(row(allyxdec([2 3],:)))
                        % re-dispatch
                        autolineposition(D)
                    else
                        % keep same positions,
                        % update yxdec according to change into activeidx
                        ddy = D.activeidx - prevactiveidx;
                        for k=find(ddy)
                            D.yxdec{k}(3,:) = D.yxdec{k}(3,:) - ddy(k);
                        end
                    end
                case {'remove'}
                    delete(D.hplot{ind})
                    D.hplot(ind) = [];
                    D.yxdec(ind) = [];
                    % no re-dispatch - TODO
                    autolineposition(D,'set',[]) % this only changes D.oldclip
                    for k=1:nsel
                        slicenum(D,k)
                    end
                    prevactiveidx = D.activeidx;
                    prevactiveidx(ind) = [];
                    D.activeidx = cumsum(logical([slice.active]));
                    ddy = D.activeidx - prevactiveidx;
                    for k=find(ddy)
                        D.yxdec{k}(3,:) = D.yxdec{k}(3,:) - ddy(k);
                    end
            end
            
        end
            
        function slicenum(D,ksel)
            hl = D.hplot{ksel};
            if isempty(hl), return, end
            if ~isappdata(hl(1),'color')
                %                 if ~isscalar(D.linecol)
                %                     set(hl,'color',D.linecol)
                %                 else
                ncolor = size(D.colorset,1);
                if D.linecol == 0
                    for i=1:length(hl)
                        set(hl(i),'color',D.colorset(1+mod(i-1,ncolor),:))
                    end
                else
                    siz = size(D.SI.slice); siz(end+1:D.linecol) = 1;
                    idx = fn_indices(siz,ksel); i = idx(D.linecol);
                    set(hl,'color',D.colorset(1+mod(i-1,ncolor),:))
                end
            end
            for i=1:length(hl)
                setappdata(hl(i),'ksel',ksel)
            end
        end
        
        function [hl ydec] = slicedisplay(D,slice)
            % line handles and vertical shift
            hl = feval(D.slicedisplayfun,D,slice);
            switch size(hl,1)
                case 1
                    ydec = zeros(1,length(hl));
                case 2
                    ydec = hl(2,:);
                    hl = hl(1,:);
                otherwise
                    error('wrong output for user-defined slice display function')
            end
            % additional properties
            if slice.active
                visible = 'on';
            else
                visible = 'off';
            end
            % handle existing buttondown function
            for i=1:length(hl)
                oldbtdwnfcn = get(hl(i),'buttondownfcn');
                setappdata(hl(i),'actdispPlot_oldbtdwnfcn',oldbtdwnfcn);
            end
            set(hl,'visible',visible,'Tag','fn4D_line', ...
                'buttondownfcn',@(hlin,evnt)moveline(D,hlin))
            % attached data: note that ksel is set (and changed) in 'displaydata'
            for i=1:length(hl)
                try
                    xdatai = get(hl(i),'xdata'); 
                    ydatai = get(hl(i),'ydata'); 
                catch
                    p = get(hl(i),'pos');
                    xdatai = p(1);
                    ydatai = p(2);
                end
                setappdata(hl(i),'xdata',xdatai)
                setappdata(hl(i),'ydata',ydatai)
                setappdata(hl(i),'icat',i)
            end
        end
        
        function hl = slicedisplaydefault(D,slice)
            if isempty(slice.data), hl=[]; return, end
            [dum nplus] = size(slice.data); %#ok<ASGLU>
            if isfield(slice,'tidx'), tt = slice.tidx; else tt = D.tidx; end
            %markers = PlotStyles(nplus);
            hl = zeros(1,nplus);
            for j=1:nplus
                y = slice.data(:,j);
                hl(j) = line(tt,y,'parent',D.ha); % 'marker',markers(j));
            end
        end
        
        function displayline(D)
            fn4D_dbstack
            t = IJ2AX(D.SI,D.SI.ij2);
            set(D.hline,'xdata',[t t]);
            % make sure that the line is in the foreground
            uistack(D.hline,'top')
        end
        
        function displayzoom(D)
            fn4D_dbstack
            zoom = IJ2AX(D.SI,D.SI.zoom);
            % set axis value - automatic range check
            D.axis = zoom;
        end
            
        function displayselection(D)
            fn4D_dbstack
            delete(findobj(D.ha,'tag','ActDispIm_Sel'))
            nsel = length(D.SI.selectionmarks);
            if nsel>1
                disp('only one selection along plot x axis can be shown')
                nsel = 1;
            end
            if nsel==1 && D.SI.selectionmarks(1).active
                poly = D.SI.selectionmarks(1).poly;
                ax = axis(D.ha); %#ok<CPROP>
                for k=1:length(poly)
                    tsel = IJ2AX(D.SI,poly(k).points);
                    switch poly(k).type
                        case 'point1D'
                            line(tsel,ones(1,length(tsel))*ax(3), ...
                                'marker','.','linestyle','none', ...
                                'color','r','parent',D.ha,'HitTest','off', ...
                                'tag','ActDispIm_Sel');
                        case 'line1D'
                            for i=1:2:length(tsel)
                                line(tsel([i i+1]),[ax(3) ax(3)], ...
                                    'linestyle','-','linewidth',5, ...
                                    'color','r','parent',D.ha,'HitTest','off', ...
                                    'tag','ActDispIm_Sel');
                            end
                        otherwise
                            error programming
                    end
                end
            end
            D.htsel = [];
        end
        
        function displaydecoration(D)
            delete(findobj(D.ha,'Tag','fn4D_deco'))
            deco = D.SI.decoration;
            if isempty(deco), return, end
            
            D.hdeco = struct('lines',[]);
            for k=1:length(deco.t)
                if ~isequal(deco.t(k).dims,1)
                    disp 'only decoration if the first dimension can be displayed'
                    continue
                end
                sel = deco.t(k).set;
                for i=1:length(sel)
                    poly = sel(i).poly;
                    for j=1:length(poly)
                        switch poly(j).type
                            case 'point1D'
                                points = poly(j).points;
                                hl = fn_spikedisplay(points,D.clip, ...
                                    'color',[1 1 1]*.5,'hittest','off','tag','fn4D_deco');
                                for idx=1:length(hl), uistack(hl(idx),'bottom'), end
                                D.hdeco.lines = [D.hdeco.lines hl];
                            otherwise                                
                                disp 'only decoration of type ''point1D'' can be displayed'
                        end
                    end
                end
            end
        end
        
        function updatedecoration(D)
            hl = D.hdeco.lines;
            if ~isempty(hl)
                ycur = get(hl(1),'ydata');
                ylim = D.clip;
                if ycur(1)>ylim(1) || ycur(2)<ylim(2)
                    for idx=1:length(hl)
                        ydata = get(hl(idx),'ydata');
                        ydata(1:3:end) = min(ydata(1),ylim(1));
                        ydata(2:3:end) = max(ydata(2),ylim(2));
                        set(hl(idx),'ydata',ydata)
                    end
                end
            end
        end
    end
    methods
        function autolineposition(D,flag,ind)
            % function autolineposition(D,'set|reset|toggle|+|-',ind)
            fn4D_dbstack
            
            % Input
            if nargin<2, flag = 'set'; end
            if nargin<3, ind = 1:length(D.hplot); end % all selections
            resetstep = (strcmp(D.clipmode,'slice') && fn_ismemberstr(flag,{'set','reset'}));
            clipignoreprev = strcmp(D.clipmode,'slice') || fn_ismemberstr(flag,{'toggle' 'reset'});

            % Some variables
            updatestep = 0;
            ncat = max(fn_map(@length,D.hplot));
            % davailable gathers dimensions that can be used
            % for dispatching, including possibly 0 for
            % 'categories'
            davailable = find(size(D.SI.slice)>1);
            if ncat>1, davailable(end+1) = 0; end
            
            % Button press actions
            switch flag
                case {'set' 'reset'}
                    % nothing to do
                case 'toggle'
                    allyxdec = [D.yxdec{:}];
                    okxydec = ~any(row(allyxdec([2 3],:)));
                    if isempty(D.autopos)
                        nl = 1;
                    elseif D.autopos(1)==0
                        nl = ncat;
                    else
                        nl = size(D.SI.slice,D.autopos(1));
                    end
                    okclip = all(D.axis==D.tidx([1 end])) && ...
                        all(D.clip==D.oldclip) && (isempty(D.autopos) || all(D.oldclip==[-(nl+.6) -.4]));
                    if okxydec && okclip
                        % change mode only if lines are at their default
                        % position and clipping has a correct value
                        
                        % we need to implement a logical way to enumerate
                        % all the possible ordered subsets of davailable;
                        % this is tricky! the best is to use an example: 
                        % if davailable = [1 3 0],
                        % the possible subsets are, [], 1, 1 3, 1 3 0, 1 0,
                        % 1 0 3, 3, 3 1, 3 1 0, 3 0, 3 0 1, 0, 0 1, 0 1 3,
                        % 0 3, 0 3 1
                        D.autopos = intersect(D.autopos,davailable,'stable'); % remove dimensions that are no longer available
                        if isempty(davailable)
                            D.autopos = [];
                        elseif length(D.autopos)<length(davailable)
                            dm = setdiff(davailable,D.autopos,'stable');
                            D.autopos = [D.autopos dm(1)];
                        else
                            % find the rightmost (if any) element that can
                            % be 'increased'
                            D.autopos(end) = [];
                            while ~isempty(D.autopos)
                                elem = D.autopos(end);
                                D.autopos(end) = [];
                                dm = setdiff(davailable,D.autopos,'stable');
                                idx = find(dm==elem);
                                if idx<length(dm)
                                    D.autopos(end+1) = dm(idx+1);
                                    break
                                end
                            end
                        end
                    end
                    for k=1:length(D.yxdec), D.yxdec{k}(2:3,:) = 0; end
                    % also reset the x-axis
                    D.axis = D.tidx([1 end]);
                case '+'
                    updatestep=-1;
                case '-'
                    updatestep=+1;
            end
            
            % dispatch lines or not
            if ~isempty(D.autopos)
                
                % 'DISPATCHED' MODES
                
                % guess an optimal step
                autostep = 0;
                for k=1:length(D.hplot)
                    hl = D.hplot{k};
                    for i=1:length(hl)
                        ydata = getappdata(hl(i),'ydata');
                        autostep = max(autostep,max(ydata)-min(ydata));
                    end
                end
                autostep = double(autostep);
                
                % use this automatic step?
                % if specified, or if the current step is more than 10
                % times larger or smaller than the autostep (this probably
                % means that the data has changed importantly in magnitude,
                % and that it is wise to replace the current step by the
                % automatic one)
                if isempty(D.ystep) || resetstep || (~updatestep && ~isnan(autostep/autostep) && abs(log10(D.ystep/autostep))>2)
                    D.ystep = autostep;
                    newstep = true;
                else
                    newstep = false;
                end
                if isnan(D.ystep/D.ystep), D.ystep=1; end % this includes D.ystep being zero, Inf or NaN
                
                % round ystep to a nice value
                if newstep || updatestep
                    y10 = 10^floor(log10(D.ystep));
                    y = D.ystep / y10;
                    vals = [.75 1 1.5 2 3 4 5 7.5 10];
                    f = find(y*1.1>vals,1,'last');
                    
                    % update if specified
                    f = f+updatestep;
                    D.ystep = vals(f) * y10;
                end
                
                %                 % title
                %                 title(D.ha,['lines separated by ' num2str(D.ystep)])
                
                % set line positions according to ystep
                nd = max([ndims(D.SI.slice) max(D.autopos)]); % e.g. D.autopos could still have a 3, whereas D.SI.slice is no longer of dimension 3; works also when D.autopos is empty
                siz = [ncat size(D.SI.slice)]; siz(end+1:nd+1) = 1;
                decs = zeros(1,1+nd);
                offsetdec = 1; % for the case where D.autopos is empty
                for i=1:length(D.autopos)
                    if i==1
                        curdec = 1;
                        offsetdec = 0;
                    else
                        si = siz(1+D.autopos(i));
                        curdec = curdec / (si+.5);
                        offsetdec = offsetdec - curdec*(si+1)/2;
                    end
                    decs(1+D.autopos(i)) = curdec;
                end
                
                for k=ind
                    hl = D.hplot{k};
                    if isempty(hl), continue, end
                    % removing of a constant: use the data of the first
                    % line
                    idx = row(fn_indices(siz(2:end),k));
                    for i=1:length(hl)
                        ydata = fn_float(getappdata(hl(i),'ydata'));
                        baseflag = D.baseline; if iscell(baseflag), baseflag = baseflag{1}; end
                        switch baseflag
                            case 'auto'
                                ymin = min(ydata); ymax = max(ydata); closethr = .05;
                                if (ymin<=0 && ymax>=0) || (min(abs([ymin ymax]))/max(abs([ymin ymax])) < closethr)
                                    yoffset = 0;
                                elseif ymin<=1 && ymax>=1 || ((ymin-1)/(ymax-1) < closethr)
                                    yoffset = 1;
                                else
                                    yoffset = mean(ydata(~isnan(ydata) & ~isinf(ydata)));
                                end
                            case 'mean'
                                yoffset = mean(ydata(~isnan(ydata) & ~isinf(ydata)));
                            case '0'
                                yoffset = 0;
                            case '1'
                                yoffset = 1;
                            case 'meantime'
                                meantime = D.baseline{2};
                                idxedge = (meantime-D.SI.grid(2))/D.SI.grid(1);
                                idxt = max(1,round(idxedge(1))):min(length(ydata),round(idxedge(2)));
                                yoffset = mean(ydata(idxt));
                        end
                        kdec = offsetdec + sum([i idx].*decs);
                        ydata = ydata/D.ystep + (-yoffset/D.ystep - kdec + sum(D.yxdec{k}([1 3],i)));
                        fn_set(hl(i),'ydata',ydata) % fn_set can deal not only with lines but also rectangles, etc.
                        if ~D.moveyonly
                            xdata = getappdata(hl(i),'xdata');
                            fn_set(hl(i),'xdata',xdata + D.yxdec{k}(2,i))
                        end
                    end
                end
                
                % clip
                if isempty(D.autopos)
                    nl = 1;
                elseif D.autopos(1)==0
                    nl = ncat;
                else
                    nl = size(D.SI.slice,D.autopos(1));
                end
                mM = [-(nl+.6) -.4];
                prevclip = D.clip;
                if clipignoreprev
                    D.oldclip = mM;
                    D.clip = mM; % automatic update of scale bar
                elseif mM(1)<D.oldclip(1) || mM(2)>D.oldclip(2)
                    % need to extend the current clipping, this also
                    % updates automatically D.oldclip
                    D.clip = [min(mM(1),D.oldclip(1)) max(mM(2),D.oldclip(2))];
                elseif newstep || updatestep
                    % update scale bar
                    displayscalebar(D)
                end
                
                % ticks: keep only integer numbers
                if any(D.clip~=prevclip) 
                    set(D.ha,'ytickmode','auto')
                    ytick = get(D.ha,'ytick');
                    ytick = ytick(ismember(ytick,-nl:-1));
                    yticklabel = fn_num2str(-ytick,'cell');
                    set(D.ha,'ytick',ytick,'yticklabel',yticklabel)
                end
                
            else
                
                % 'ALIGNED' MODE
                
                doset = ~updatestep; % if '+' or '-' has been selected, action only on the clipping
                %                 if doset, title(D.ha,''), end
                % compute min and max while replacing lines if necessary
                if length(ind)==length(D.hplot)
                    % apply on all signals -> estimate limits from scratch
                    m = Inf; M = -Inf;
                else
                    % apply on some signals only -> only extend limits if needed
                    if updatestep, error programming, end
                    m = D.clip(1); M = D.clip(2);
                end
                for k=ind
                    hl = D.hplot{k};
                    for i=1:length(hl)
                        ydata = getappdata(hl(i),'ydata') + sum(D.yxdec{k}([1 3],i));
                        yy = ydata(~isinf(ydata) & ~isnan(ydata)); % note that yy can be empty (e.g. if all ydata are NaN)
                        m = min([m min(yy)]);
                        M = max([M max(yy)]);
                        if doset, fn_set(hl(i),'ydata',ydata), end
                        if ~D.moveyonly
                            xdata = getappdata(hl(i),'xdata');
                            fn_set(hl(i),'xdata',xdata+D.yxdec{k}(2,i))
                        end
                    end
                end
                mM = [m M];
                if any(isinf(mM)) || any(isnan(mM))
                    if ~isinf(m) && ~isnan(m), mM = [m m+1]; elseif ~isinf(M) && ~isnan(M), mM = [M-1 M]; else mM = [0 1]; end
                end
                % update clipping
                prevclip = D.clip;
                if updatestep
                    c = mean(D.clip); d = D.clip(2)-c;
                    d = d * 2^(updatestep*1/2);
                    D.clip = c + [-1 1]*d;
                    D.oldclip = fn_minmax('minmax',mM,D.clip);
                else
                    if clipignoreprev
                        D.oldclip = mM;
                        D.clip = mM; % automatic update of scale bar
                    elseif mM(1)<D.oldclip(1) || mM(2)>D.oldclip(2)
                        % need to extend the current clipping, this also
                        % updates automatically D.oldclip
                        D.clip = fn_minmax('minmax',mM,D.oldclip);
                    end
                end
                % no weird stuff with yticks!
                if any(D.clip~=prevclip)
                    set(D.ha,'ytickmode','auto','yticklabelmode','auto')
                end
            end
            
        end        
    end
    
    % Update routines
    methods (Access='private')
    end
    
    % GET/SET clip
    methods
        function set.clip(D,clip)
            fn4D_dbstack
            clip = double(clip);
            clipold = D.clip;
            if isempty(clip), return, end % occurs when pressing 'user clip...' menu but later canceling
            if all(clip==clipold), return, end
            if diff(clip)<=0
                if diff(clip)==0
                    clip = clip + [-1 1];
                else
                    clip = [0 1];
                end
            end
            D.clip = clip;
            set([D.ha D.ha2(2)],'YLim',clip);
            % if the new clip extends beyond oldclip, update oldclip value
            D.oldclip = [min(clip(1),D.oldclip(1)) max(clip(2),D.oldclip(2))];
            % propagate change in the case of clip link
            if fn_ismemberstr(D.clipmode,{'link1','link2'})
                D.CL.clip = clip;
            end
            % update the selection display (red marks should stay in the bottom)
            if clip(1)~=clipold(1) && ~isempty(D.SI.selectionmarks)
                displayselection(D)
            end
            % change slider parameters
            set(D.slider(2),'value',D.clip)
            slidercliphide(D)
            % update scale bar (TODO: only one call when set.axis +
            % set.clip)
            displayscalebar(D)
            % update decorations
            updatedecoration(D)
        end
        
        function set.oldclip(D,clip)
            if all(clip==D.oldclip), return, end
            D.oldclip = clip;
            % coerce the current clip to fit inside oldclip
            if D.clip(1)<=clip(2) && D.clip(2)>=clip(1)
                D.clip = [max(clip(1),D.clip(1)) min(clip(2),D.clip(2))];
            end
            % change slider parameters
            set(D.slider(2),'minmax',D.oldclip)
            slidercliphide(D)
            % change size of time-marking line
            set(D.hline,'ydata',clip)
        end
        
        function slidercliphide(D)
            if all(D.clip==D.oldclip)
                set(D.slider(2),'visible','off')
            else
                set(D.slider(2),'visible','on')
            end
        end
        
        function set.clipmode(D,clipmode)
            fn4D_dbstack
            if strcmp(clipmode,D.clipmode), return, end
            if fn_ismemberstr(D.clipmode,{'link1','link2'})
                % cancel previous cliplink and listener
                disconnect(D,D.CL), delete(D.C2D)
            end
            D.clipmode = clipmode;
            switch clipmode
                case 'slice'
                    kitem = 1;
                    autolineposition(D)
                case 'data'
                    kitem = 2;
                    % no change in current clipping
                case {'link1','link2'}
                    kitem = 2+str2double(clipmode(5));
                    D.CL = cliplink.find(clipmode,D.clip);
                    D.clip = D.CL.clip;
                    D.C2D = connect_listener(D.CL,D,'ChangeClip', ...
                        @(cl,evnt)clipfromlink(D,D.CL));
                otherwise
                    D.clipmode = slice;
                    error('wrong clip mode')
            end
            % check mark in uicontextmenu
            set(D.menuitems.clip,'checked','off')
            set(D.menuitems.clip(kitem),'checked','on')
        end
        
        function clipfromlink(D,CL)
            D.clip = CL.clip;
        end
    end
       
    % GET/SET 
    methods
        function signals = get.signals(D)
            fn4D_dbstack
            signals = cat(2,D.SI.slice.data);
        end        
        function set.axis(D,axis)
            fn4D_dbstack
            if all(axis==D.axis), return, end
            % range check
            if length(D.tidx)==1
                D.axis = D.tidx([1 1]);
                set([D.ha D.ha2(1)],'xLim',D.tidx+[-1 1])
            else
                D.axis = [max(D.tidx(1),axis(1)) min(D.tidx(end),axis(2))];
                if diff(D.axis)<=0
                    disp('new axis is outside of range - reset axis')
                    D.axis = D.tidx([1 end]);
                end
                set([D.ha D.ha2(1)],'xLim',D.axis);
            end
            % change slider parameters
            set(D.slider(1),'value',D.axis)
            slideraxishide(D)
            % update scale bar (TODO: only one call when set.axis +
            % set.clip)
            displayscalebar(D)
        end        
        function slideraxishide(D)
            if all(D.axis==D.tidx([1 end]))
                set(D.slider(1),'visible','off')
            else
                set(D.slider(1),'visible','on')
            end
        end       
        function set.tidx(D,tidx)
            if isequal(tidx,D.tidx), return, end
            D.tidx = tidx;
            set(D.slider(1),'minmax',tidx([1 end]))
            % re-set axis (automatic coerce)
            D.axis = D.axis;
        end
        function set.scaledisplay(D,flag)
            if ~fn_ismemberstr(flag,{'ybar',''})
                error('wrong value for ''scaledisplay'' property')
            end
            D.scaledisplay = flag;
            displayscalebar(D)
            set(D.menuitems.scalebar,'checked',onoff(strcmp(flag,'ybar')))
        end
        function set.movescale(D,b)
            D.movescale = b;
            displayscalebar(D)
        end
        function set.slicedisplayfun(D,fun)
            if isempty(fun)
                D.slicedisplayfun = @slicedisplaydefault;
            elseif ischar(fun) && strcmp(fun,'update')
                % function has not changed, but its output might be
                % different due to parameters stored elsewhere: update the
                % display 
            elseif isa(fun,'function_handle')
                D.slicedisplayfun = fun;
            else
                error('''fun'' is not a function handle')
            end
            % re-display data: special flag 'changeall' -> complete reset
            % only if the number of curves has changed
            displaydata(D,'changeall')
        end        
        function autolineppos = get.autolinepos(D)
            autolineppos = D.autopos;
        end        
        function set.autolinepos(D,autopos)
            if islogical(autopos) % old old version
                autopos = fn_switch(autopos,1,[]);
            elseif ischar(autopos) % old version
                autopos = fn_switch(autopos,'',[],'sel',1,'cat',0,'selcat',[1 0],'catsel',[0 1]);
            end
            if strcmp(autopos,D.autopos), return, end
            D.autopos = autopos;
            autolineposition(D,'set')
        end        
        function set.movelinegroup(D,flag)
            if strcmp(D.movelinegroup,flag), return, end
            % update check marks
            switch flag
                case ''
                    kitem = 1;
                case 'sel'
                    kitem = 2;
                case 'cat'
                    kitem = 3;
                otherwise
                    error('flag from line-move grouping must be one of '''', ''sel'', ''cat''')
            end
            set(D.menuitems.moveline,'checked','off')
            set(D.menuitems.moveline(kitem),'checked','on')
            % set property
            D.movelinegroup = flag;
        end
        function set.moveyonly(D,b)
            % set property
            D.moveyonly = b;
            % update mark
            set(D.menuitems.moveyonly,'checked',onoff(b))
        end
        function set.linecol(D,val)
            oldval = D.linecol;
            if strcmp(val,'toggle')
                % davailable gathers dimensions that can be used
                % for coloring, including possibly 0 for
                % 'categories'
                davailable = find(size(D.SI.slice)>1);
                ncat = max(fn_map(@length,D.hplot));
                if ncat>1, davailable(end+1) = 0; end
                if isempty(davailable)
                    val = 2;
                else
                    idx = find(davailable==oldval);
                    if isempty(idx) || idx==length(davailable), idx=0; end
                    val = davailable(idx+1);
                end
            end
            % old version
            if ischar(val)
                val = fn_switch(val,'sel',1,'cat',0);
            end
            if isempty(val) || isequal(oldval,val), return, end
            %             % update check marks
            %             oldflag = fn_switch(isnumeric(oldval),'custom',oldval);
            %             curflag = fn_switch(isnumeric(val),'custom',val);
            %             set(D.menuitems.linecol.(oldflag),'checked','off')
            %             set(D.menuitems.linecol.(curflag),'checked','on')
            % set property
            D.linecol = val;
            % update display
            displaysignals(D)
        end
        function set.colorset(D,val)
            % convert string to set of colors
            if ischar(val)
                try
                    try 
                        val = fn_colorset(val);
                    catch
                        val = feval(val);
                    end
                catch
                    error 'string cannot be evaluated to a valid color set'
                end
            end
            % checks
            if ~isnumeric(val) || size(val,2)~=3 || ~all(row(val>=0 & val<=1))
                error 'color set must be a numerical array with 3 rows, with values between 0 and 1'
            end
            % set property
            if isequal(val,D.colorset), return, end
            D.colorset = val;
            % update display
            displaysignals(D)
        end
        function x = setinfo(D) %#ok<MANU>
            x.axis = '2-elements vector';
            x.autolinepos = 'vector';
            x.linecol = 'scalar (dimension or 0 for ''category'')';
            x.movelinegroup = {'sel' 'cat' ''};
            x.scaledisplay = {'ybar' ''};
            x.movescale = {'0' '1'};
            x.slicedisplayfun = 'function hl = slicedisplayfun(D,slice) - for ''buttondownfcn'' property, use app-def data ''ksel'' and ''icat'' - don''t use ''color'' property but create ''color'' app-def data - don''t use ''tag'' property';
            x.clip = '2-elements vector';
        end
        function set.scrollwheel(D,val)
            val = fn_switch(val,'on','xy','off','',val);
            if ~fn_ismemberstr(val,{'' 'x' 'xy' 'y'})
                error '''scrollwheel'' property must be '''', ''x'', ''y'' or ''xy'''
            end
            oldval = D.scrollwheel;
            if strcmp(val,oldval), return, end
            % update value
            D.scrollwheel = val;
            % update callback
            if isempty(val)
                fn_scrollwheelregister(D.ha,'off')
            else
                fn_scrollwheelregister(D.ha,'on')
            end 
            % (marks are automatically updated by the fn_propcontrol object)
        end
        function set.navigation(D,val)
            fn_ismemberstr(val,{'zoom' 'pan'},'doerror')
            D.navigation = val;
            % (marks are automatically updated by the fn_propcontrol object)
        end
        function set.baseline(D,val)
            meantime = [];
            if isnumeric(val)
                if isscalar(val)
                    val = num2str(val);
                elseif length(val)==2
                    meantime = val;
                    val = 'meantime';
                else
                    error 'wrong baseline value'
                end
            end
            kitem = find(strcmp(val,{'auto' 'mean' '0' '1' 'meantime'}));
            if isempty(kitem), error 'wrong baseline value', end
            if strcmp(val,'meantime') 
                if isempty(meantime)
                    meantime = fn_mouse(D.ha,'xsegment','Select baseline time segment');
                end
                val = {'meantime' meantime};
            end
            if isequal(D.baseline,val), return, end
            D.baseline = val;
            set(D.menuitems.baseline,'checked','off')
            set(D.menuitems.baseline(kitem),'checked','on')
            displaysignals(D)
        end
        function set.menustayvisible(D,val)
            D.menustayvisible = val;
            set(D.menuitems.pin,'checked',onoff(val))
            if val
                set(D.menu,'pos',get(D.hf,'currentpoint'),'visible','on')
            end
        end
    end
    
    % Events (bottom-up: mouse)
    methods (Access='private')
        function Mouse(D,flag)
            fn4D_dbstack
            % different selection types are:
            % - point with left button          -> change time
            % - area with left button           -> zoom and clip to region
            % - double-click with left button   -> zoom and clip reset
            %   (or point with middle button outside of axis)
            % - point/area with middle button   -> add point/area to current selection
            % - double-click with middle button -> cancel current selection
            %   (or point with middle button outside of axis)
            % - point/area with right button    -> add new selection
            % - double-click with right button  -> cancel all selections
            %   (or point with right button outside of axis)
            
            % shortcuts
            hb = D.ha;
            si = D.SI;
            
            % flag and selection type
            if nargin>=2
                outtype = strrep(flag,'outside','');
                flag = 'outside';
                selectiontype = get(D.hf,'selectiontype');
                outsideleft = fn_switch(outtype,{'left' 'both'},true,'down',false);
                outsidedown = fn_switch(outtype,{'down' 'both'},true,'left',false);
            else
                flag = get(D.hf,'selectiontype');
                if strcmp(flag,'open')
                    selectiontype = D.oldselectiontype;
                else
                    selectiontype = flag;
                end
                % special case - click outside of axis
                ax = axis(hb); %#ok<CPROP>
                point =  get(hb,'CurrentPoint'); point = point(1,[1 2]);
                if (point(1)<ax(1) || point(1)>ax(2) || point(2)<ax(3) || point(2)>ax(4))
                    flag = 'outside';
                    outsideleft = (point(1)<ax(1) || point(1)>ax(2));
                    outsidedown = (point(2)<ax(3) || point(2)>ax(4));
                end
            end
            D.oldselectiontype = selectiontype;
            
            switch flag
                case 'normal'                           % CHANGE VIEW AND/OR MOVE CURSOR
                    if strcmp(D.navigation,'zoom') || all(D.axis==D.tidx([1 end]) & D.clip==D.oldclip)
                        rect = fn_mouse(hb,'rect-^');
                    else
                        rect = [];
                    end
                    if isempty(rect)
                        moved = pan(D);                 % pan
                        if moved
                            return
                        else
                            p = get(D.ha,'currentpoint'); 
                            rect = [p(1,1:2) 0 0]; 
                        end
                    end
                    if any(rect(3:4))                   % zoom in
                        newaxis = [rect(1)+[0 rect(3)]; rect(2)+[0 rect(4)]];
                        % change zoom (automatic display update)
                        if rect(3)>diff(D.axis)/100
                            si.zoom = AX2IJ(si,newaxis(1,:));
                        end
                        % and change clip (local action if no link clipmode)
                        if rect(4)>diff(D.clip)/100
                            D.clip = newaxis(2,:);
                        end
                    else                                % change t
                        si.ij2 = AX2IJ(si,rect(1));
                    end
                case 'extend'                           % EDIT CURRENT SELECTION
                    rect = fn_mouse(hb,'rect-');
                    flag = fn_switch(vide(si.selection),'new','add');
                    if rect(3)==0                       % select a time point
                        updateselection(si,flag,1, ...
                            selectionND('point1D',AX2IJ(si,rect(1))));
                    else                                % select a temporal segment
                        updateselection(si,flag,1, ...
                            selectionND('line1D',AX2IJ(si,[rect(1) rect(1)+rect(3)])));
                    end
                case 'alt'                              % CHANGE SELECTION
                    rect = fn_mouse(hb,'rect-');
                    flag = fn_switch(vide(si.selection),'new','change');
                    if rect(3)==0                       % select a time point
                        updateselection(si,flag,1, ...
                            selectionND('point1D',AX2IJ(si,rect(1))));
                    else                                % select a temporal segment
                        updateselection(si,flag,1, ...
                            selectionND('line1D',AX2IJ(si,[rect(1) rect(1)+rect(3)])));
                    end
                case 'open'                             % MISC
                    switch selectiontype
                        case 'normal'    
                            if ~isempty(D.usercallback) % user callback
                                feval(D.usercallback,D)
                            else
                                % reset zoom
                                si.zoom = [-Inf Inf];   % zoom reset
                                % and clip
                                D.clip = D.oldclip;
                            end
                        otherwise
                            return
                    end
                case 'outside'                          % RESET SELECTION
                    switch selectiontype
                        case 'normal'                   % zoom reset
                            rect = fn_mouse(hb,'rect-');
                            if all(rect(3:4))           % zoom in
                                newaxis = [rect(1)+[0 rect(3)]; rect(2)+[0 rect(4)]];
                                if outsidedown
                                    % change zoom 
                                    si.zoom = AX2IJ(si,newaxis(1,:));
                                end
                                if outsideleft
                                    % change clip
                                    D.clip = newaxis(2,:);
                                end
                            else                        % zoom reset
                                % here we look "where exactly outside"
                                if outsidedown
                                    % reset zoom
                                    si.zoom = [-Inf Inf];
                                end
                                if outsideleft
                                    % reset clip
                                    D.clip = D.oldclip;
                                end
                            end
                        case {'extend','alt'}           % unselect all regions
                            if ~vide(si.selection), updateselection(si,'remove',1); end
                        otherwise
                            return
                    end
                otherwise
                    error programming
            end
        end
        function Scroll(D,nscroll)
            si = D.SI;
            zoomfactor = 1.5^nscroll;
            p = get(D.ha,'currentpoint'); 
            p = p(1,1:2);
            % x-axis
            if any(D.scrollwheel=='x')
                xaxis = [max(.5,si.zoom(1)) min(si.sizes+.5,si.zoom(2))];
                origin = AX2IJ(D.SI,p(1));
                newzoom = origin + zoomfactor*(xaxis-origin);
                newzoom = fn_minmax('maxmin',newzoom,[.5 si.sizes+.5]);
                si.zoom = newzoom;
            end
            % y-axis
            if any(D.scrollwheel=='y')
                D.clip = p(2) + zoomfactor*(D.clip-p(2));
            end
        end
        function movebar(D,hl)
            fn4D_dbstack
            if ~strcmp(get(D.hf,'selectiontype'),'normal')
                % execute callback for axes
                Mouse(D)
                return
            end           
            fn_buttonmotion({@movebarsub,D,hl},get(D.ha,'parent'),'doup');
        end
        function moved = pan(D)
            set(D.hf,'pointer','hand')
            p0 = get(D.ha,'currentpoint'); p0 = p0(1,1:2)';
            xlim0 = get(D.ha,'xlim'); ylim0 = get(D.ha,'ylim');
            moved = false;
            fn_buttonmotion(@chgaxis,D.hf)
            function chgaxis
                moved = true;
                p = get(D.ha,'currentpoint'); p = p(1,1:2)';
                mov = p0-p;
                sidedist = D.tidx([1 end]) - D.axis;
                %mov(1) = max(sidedist(1),min(sidedist(2),mov(1)));
                set(D.ha,'xlim',get(D.ha,'xlim')+mov(1),'ylim',get(D.ha,'ylim')+mov(2))
            end
            D.axis = get(D.ha,'xlim')+eps;
            D.clip = get(D.ha,'ylim')+eps;
            set(D.hf,'pointer','arrow')
            if moved, D.SI.zoom = AX2IJ(D.SI,D.axis); end
        end
    end   
    % (but this guy must be public to be accessible by fn_buttonmotion!)
    methods
        function p = movebarsub(D,hl) %#ok<INUSD>
            p = get(D.ha,'currentpoint'); p = p(1);
            D.SI.ij2 = AX2IJ(D.SI,p);
            %set(hl,'xdata',[p p]);
        end
    end
    
    % Events (bottom-up/local: buttons)
    methods (Access='private')
        function chgtime(D,flag)
            fn4D_dbstack
            D.SI.ij2 = D.SI.ij2 + str2double(flag);
        end
        
        function chgzoom(D,flag,hu)
            fn4D_dbstack
            ax = D.axis;
            mov = ax(2)-ax(1);
            switch flag
                case '-1'
                    ax = ax-min(mov,ax(1)-D.tidx(1));
                case '+1'
                    ax = ax+min(mov,D.tidx(end)-ax(2));
                case 'slider'
                    ax = get(hu,'value');
                otherwise
                    error programming
            end
            if hu.sliderscrolling
                set([D.ha D.ha2(1)],'xlim',ax)
            else
                D.SI.zoom = AX2IJ(D.SI,ax);
            end
        end
        
        function chgclip(D,flag,hu)
            fn4D_dbstack
            if ~strcmp(flag,'slider'), error programming, end
            newclip = get(hu,'value');
            if hu.sliderscrolling
                set([D.ha D.ha2(2)],'ylim',newclip)
            else
                D.clip = newclip;
            end
        end
        
        function duplicate(D,hobj)
            if nargin<2, hobj=figure; end
            D2 = activedisplayPlot(D.SI,'in',hobj);
            D2.clip = D.clip;
            D2.autopos = D.autopos;
            D2.movelinegroup = D.movelinegroup;
            D2.linecol = D.linecol;
            D2.slicedisplayfun = D.slicedisplayfun;
        end
    end
  
    % Events (local: move lines)
    methods (Access='private')
        function moveline(D,hl)
            fn4D_dbstack
            ksel = getappdata(hl,'ksel');
            icat = getappdata(hl,'icat');
            % not left-click -> execute callback for line or axes, or user
            % callback
            if ~strcmp(get(D.hf,'selectiontype'),'normal')
                if strcmp(get(D.hf,'selectiontype'),'open') && ~isempty(D.lineusercallback)
                    feval(D.lineusercallback,D,row(fn_indices(size(D.SI.slice),ksel)),icat)
                else
                    fun = getappdata(hl,'actdispPlot_oldbtdwnfcn');
                    if ~isempty(fun), fn_evalcallback(fun,hl), else Mouse(D), end
                    return
                end
            end           
            % which lines to move?
            switch D.movelinegroup
                case ''
                    hls = D.hplot{ksel}(icat);
                case 'sel'
                    hls = D.hplot{ksel}';
                case 'cat'
                    hls = [];
                    for k=1:length(D.hplot)
                        if length(D.hplot{k})>=icat
                            hls = [hls; D.hplot{k}(icat)]; %#ok<AGROW>
                        end
                    end
            end
            % point and data from which we start moving
            p0 = get(D.ha,'currentpoint'); p0 = p0(1,1:2);
            xydata0 = fn_get(hls,'xdata','ydata');
            [dx dy] = fn_buttonmotion({@movelinesub,D,hls,p0,xydata0,D.ha},get(D.ha,'parent')); 
            % end: did the line move at all?
            if any([dx dy])
                % (save time course in base workspace)
                if length(xydata0)==1, tmp=xydata0(1).ydata; else tmp={xydata0.ydata}; end
                assignin('base','currentlineydata',tmp);
                % (update 'yxdec')
                switch D.movelinegroup
                    case ''
                        if ~D.moveyonly, D.yxdec{ksel}(2,icat) = D.yxdec{ksel}(2,icat)+dx; end
                        D.yxdec{ksel}(3,icat) = D.yxdec{ksel}(3,icat)+dy;
                    case 'sel'
                        if ~D.moveyonly, D.yxdec{ksel}(2,:) = D.yxdec{ksel}(2,:)+dx; end
                        D.yxdec{ksel}(3,:) = D.yxdec{ksel}(3,:)+dy;
                    case 'cat'
                        for k=1:length(D.yxdec)
                            if size(D.yxdec{k},2)>=icat
                                if ~D.moveyonly, D.yxdec{k}(2,icat) = D.yxdec{k}(2,icat)+dx; end
                                D.yxdec{k}(3,icat) = D.yxdec{k}(3,icat)+dy;
                            end
                        end
                end
            else
                % (no move, then try an event in the line or parent axes)
                fun = getappdata(hl,'actdispPlot_oldbtdwnfcn');
                if ~isempty(fun)
                    fn_evalcallback(fun,hl)
                else
                    % wait just a little bit, so a double-click on the line
                    % can be detected before the time cursor (black
                    % vertical line) is covering it
                    if ~isempty(D.lineusercallback), pause(.1), end
                    % clicked with left button and no move -> change t
                    D.SI.ij2 = AX2IJ(D.SI,p0(1));
                end
            end
        end
    end      
    % (but this guy must be public to be accessible by fn_buttonmotion!)
    methods
        function [dx dy] = movelinesub(D,hl,p0,xydata0,ha)
            p = get(ha,'currentpoint'); p = p(1,1:2);
            dx = p(1) - p0(1);
            dy = p(2) - p0(2);
            for i=1:length(hl)
                if D.moveyonly
                    fn_set(hl(i),'ydata',xydata0(i).ydata+dy)
                else
                    fn_set(hl(i),'xdata',xydata0(i).xdata+dx,'ydata',xydata0(i).ydata+dy)
                end
            end
            ax = D.clip;
            p2 = p(2);
            if p2<ax(1), D.clip = [p2 ax(2)]; end
            if p2>ax(2), D.clip = [ax(1) p2]; end
        end
    end
    
    % Events (top-down: listeners)
    methods
        function updateDown(D,~,evnt)
            fn4D_dbstack(['S2D ' evnt.flag])
            switch evnt.flag
                case 'sizes'
                    displaygrid(D)
                case 'slice'
                    displaydata(D,evnt.selflag,evnt.ind,evnt.value)
                case 'grid'
                    displaygrid(D)
                case 'labels'
                    displaylabels(D)
                case 'units'
                    displaylabels(D)
                case 'ij2'
                    displayline(D)
                case 'zoom'
                    displayzoom(D)
                case 'selection'
                    displayselection(D)
                    if strcmp(evnt.selflag,'reset')
                        displaydata(D,'reset',[],[])
                    end
                case 'decoration'
                    displaydecoration(D)
            end
        end
    end
    
    % User
    methods
        function displaysignals(D)
            displaydata(D,'change')
        end
        function access(D) %#ok<MANU>
            keyboard
        end
        function exportsignals(D)
            % get data
            data = {D.SI.slice.data};
            s = fn_map(data(:),@size,'array');
            sslic = size(D.SI.slice); sslic(sslic==1) = [];
            if ~any(row(diff(s,1,1)))
                data = cat(3,data{:});
                if s(1,2)==1
                    data = reshape(data,[s(1,1) sslic 1]);
                else
                    data = reshape(data,[s(1,1:2) sslic]);
                end
            else
                data = reshape(data,[sslic 1 1]);
            end
            % export it
            fn_exportvar(data)
        end
    end
    
end
            

