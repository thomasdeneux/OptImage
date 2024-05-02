classdef activedisplayImage < fn4Dhandle
    % function D = activedisplayImage(SI,options...)
    %---
    % different navigaction actions are:
    % LEFT BUTTON
    % - point                                   -> change cursor
    % - hold                                    -> zoom to region/pan/change cursor depending on 'navigation' property
    % - click outside                           -> zoom reset
    % - double-click                            -> user callback function
    %                          (set as D.usercallback = @(D)dosomething(D))
    % MIDDLE BUTTON
    % - click in region, hold and type a number -> reorder selections
    % - point/hold                              -> add point/area to current selection
    % - click outside                           -> cancel current selection
    % RIGHT BUTTON
    % - point in region                         -> hide/show selection
    % - point/hold                              -> add new selection
    % - click outside                           -> cancel all selections
    % - double-click                            -> make this selection first
    % SCROLL WHEEL
    % - scrolling down/up                       -> zoom in/out
    % 
    % additional mouse actions in 'selection edit' mode are (see 'displayonesel' and 'seleditaction'):
    % LEFT BUTTON
    % - on vertex                               -> move vertex
    % - on edge                                 -> create new vertex and move it
    % MIDDLE BUTTON
    % - on vertex                               -> move region
    % - on edge                                 -> move all regions
    % RIGHT BUTTON
    % - on vertex                               -> remove vertex
 
    properties
        clipmode = 'slice';         % 'slice', 'data', 'link1' or 'link2'
        autoclipmode = 'minmax';    % 'minmax', 'Nstd', 'prcA_B'
        autocliporig = 'curview';   % 'full', 'curview', 'cursel'
        dataflag = 'data';
        
        shapemode = 'ellipse';      % 'poly', 'free', 'rect', 'ellipse', 'ring', 'segment', 'openpoly', 'freeline'
        seldims = 'xy';             % 'x', 'y' or 'xy'
        seleditmode = false;
        selround = 1;               % when making rectangle selections, adjust the edges to borders between pixels
        selmultin = 'none';         % 'none', 'point' or 'grid'
        selshow = 'number+shape';
        selcolor = false;
        
        permutexy = false;
        doratio                     % is set in displayratio
        
        movescale = false;
        
        clip
        cmap = 'gray'; 
        logscale = false;
        channelcolors

        opdef = struct;
        
        menustayvisible = false;
        usercallback
    end
    properties (SetObservable = true)
        dolabels = true;
        doticks = true;
        doxbar = false;
        clipcenter = [];
        navigation = 'zoom';        % 'zoom' or 'pan'
        scrollwheel = true;
    end
    properties (Dependent)
        crossshow
        binning
        xlow
        xhigh
        yfun
    end   
    properties (Dependent, GetAccess='private')
        scaledisplay
    end
    properties (SetAccess='private')
        ha
        hf
        img
        
        % some objects are in public access to allow change of some of
        % their properties (color, linewidth, ...)
        cross
        scalebar % scale bar and text
        hdeco
        
        cmapval = gray(256);
        selectionlabels     
        
        currentselection
        currentdisplay
        autoclipvalue
    end    
    properties (Access='private')
        seldisp = {};
        freeze = true;
        
        txt
        buttons
        slider
        menu
        colorbar
        
        menuitems
        
        curselprev
        
        oldaxis
        oldselectiontype = 'normal';
    end    
    properties (Dependent, SetAccess='private')
        axis
    end   
    properties (SetAccess='private')
        CL
        SI
        C2D
        listenaxpos
    end
    
    % Constructor and Destructor
    methods
        function D = activedisplayImage(varargin)
            % function D = activedisplayImage(SI,options...)
            fn4D_dbstack
            
            % options for initialization
            opt = struct( ...
                'in',                   [], ...
                'clip',                 [], ...
                'seldims',              'xy', ...
                'channelcolors',        []);
            if nargin==0 || ~isobject(varargin{1})
                D.SI = sliceinfo(2);
                D.SI.slice.data = rand(15,10);
                [opt optadd] = fn4D_parseInput(opt,varargin{:});
            else
                D.SI = varargin{1};
                [opt optadd] = fn4D_parseInput(opt,varargin{2:end});
            end
            
            % type check
            if (D.SI.nd+D.SI.ndplus)>3 || ((D.SI.nd+D.SI.ndplus)==3 && D.SI.sizesplus(1)>3)
                error('activedisplayImage class can display only up to 3-dimensions data, and with at most 3 elements in the 3rd dimension')
            end
            
            % dimensions of selections
            if D.SI.nd==1
                if ~isempty(opt.seldims) && ~strcmp(opt.seldims,'x')
                    error 'for uni-dimensional projection, the dimension of selections must be ''x'''
                end
                opt.seldims = 'x';
            end
            if ~isempty(opt.seldims)
                D.seldims = opt.seldims;
            end
            
            % figure and axes
            if ~isempty(opt.in)
                if ~ishandle(opt.in) && mod(opt.in,1)==0 && opt.in>0, figure(opt.in), end
                switch get(opt.in,'type')
                    case 'figure'
                        D.hf = opt.in;
                        figure(opt.in)
                        h = findall(opt.in,'type','axes');
                        if isempty(h)
                            D.ha = axes('parent',opt.in);
                        else
                            D.ha = h(1);
                        end
                    case 'uipanel'
                        delete(get(opt.in,'children'))
                        D.ha = axes('parent',opt.in);
                        D.hf = fn_parentfigure(opt.in);
                    case 'axes'
                        D.ha = opt.in;
                        D.hf = fn_parentfigure(opt.in);
                    otherwise
                        error('bad handle')
                end
            else
                D.ha = gca;
                D.hf = get(D.ha,'parent');
            end
            cla(D.ha,'reset')
            if isempty(get(D.hf,'Tag')), set(D.hf,'Tag','used by fn4D'), end
            
            % image
            D.img = image(0,'Parent',D.ha,'hittest','off','CDataMapping','scaled');
            set(D.ha,'CLimMode','manual')
            
            % set some properties
            D.channelcolors = opt.channelcolors;
            if isempty(opt.clip)
                D.clipmode = 'slice';
            else
                D.clip = opt.clip;
                D.clipmode = 'data';
            end
            D.opdef = fn_imageop('par');
            colormap(D.ha,D.cmapval)
            
            % cross
            D.cross(1) = line('Parent',D.ha,'xdata',[0 0]);
            D.cross(2) = line('Parent',D.ha,'ydata',[0 0]);
            D.cross(3) = line('Parent',D.ha,'xdata',0,'ydata',0,'marker','.','linestyle','none'); % a single point
            set(D.cross,'Color','white')
            for i=1:3, set(D.cross(i),'buttondownfcn',@(u,e)movecross(D,i)), end
            
            % scale bar (and listener for re-positioning upon axes resize)
            D.scalebar(1) = line('Parent',D.ha,'Color','white','visible','off', ...
                'linewidth',3);
            D.scalebar(2) = text('Parent',D.ha,'Color','white','visible','off', ...
                'horizontalalignment','center','verticalalignment','middle');
            D.listenaxpos = fn_pixelsizelistener(D.ha,D,@(h,evnt)displayscalebar(D));
            fn4D_enable('off',D.listenaxpos)
            
            % value and buttons; TODO: don't use fn_coordinates inside fn_controlpositions
            D.txt = uicontrol('Parent',D.hf,'style','text','enable','inactive', ...
                'fontsize',8,'horizontalalignment','left');
            fn_controlpositions(D.txt,D.ha,[0 1], [-1 0 123 10]); % was previously [2 -10 100 10]
            D.buttons(1) = uicontrol('Parent',D.hf);
            %             set(D.buttons(1),'style','frame','enable','off','backgroundcolor',[.5 0 0],'foregroundcolor',[.5 0 0]);
            % display image on button that indicates how image luminance
            % and contrast will vary when using it to adjust clipping
            [ii jj] = meshgrid(-9:0,9:-1:0); x=(0-ii)./(jj-ii)-.5; x(end)=0;
            set(D.buttons(1),'enable','inactive','cdata',fn_clip(sin(pi*x),[-1 1],'gray'))
            fn_controlpositions(D.buttons(1),D.ha,[0 1], [122 0 10 10]); % was previously [102 -10 10 10]
%             % button for projection object (note that D.buttons(2) only
%             % acts as a separator)
%             D.buttons(2) = uicontrol('Parent',D.hf,'hittest','off', ...
%                 'backgroundcolor',[.8 .8 0],'foregroundcolor',[.8 .8 0]);
%             fn_controlpositions(D.buttons(2),D.ha,[0 1], [132 0 10 10]); % was previously [112 -10 10 10]
%             D.buttons(3) = uicontrol('Parent',D.hf, ...
%                 'backgroundcolor',[.8 .8 0],'foregroundcolor',[.8 .8 0]);
%             fn_controlpositions(D.buttons(3),D.ha,[0 1], [138 0 2 10]); % was previously [118 -10 2 10]
%             set(D.buttons(2:3),'style','frame','enable','off')
            
            % callbacks (bottom-up)
            set(D.ha,'ButtonDownFcn',@(ha,evnt)Mouse(D))
            set(D.txt,'buttondownfcn',@(hu,evnt)Mouse(D,'outside'))
            set(D.buttons(1),'buttondownfcn',@(hu,evnt)redbutton(D))
            fn_scrollwheelregister(D.ha,@(n)Scroll(D,n),onoff(D.scrollwheel))
            initlocalmenu(D)
%             try
%                 initlocalmenu(D.SI,D.buttons(3))
%             catch
%                 delete(D.button(2:3))
%                 D.buttons(2:3) = [];
%             end
            
            % sliders (note that scroll wheel registration of slider must
            % occur after scroll wheel registration of axes)
            D.slider = [ ...
                fn_slider('parent',D.hf, ...
                'mode','area','layout','right','visible','off','scrollwheel','on', ... % set scrollwheel after visible so the scroll map does not need to be recalculated
                'callback',@(u,evnt)chgzoom(D,'x',u)) ...
                fn_slider('parent',D.hf, ...
                'mode','area','layout','down','visible','off','scrollwheel','on', ...
                'callback',@(u,evnt)chgzoom(D,'y',u))];
            fn_controlpositions(D.slider(1),D.ha,[0 1 1 0], [142 0 -142 10]); % was previously [122 -10 -122 10]
            fn_controlpositions(D.slider(2),D.ha,[1 0 0 1], [0 -1 10 1]); % was previously [0 -9 10 1]
            
            % communication with parent
            addparent(D,D.SI)
            
            % update display - no need to call displayscalebar(D) since
            % this is done automatically whenever axis zoom is changed
            displaydata(D)
            displayxyperm(D)
            displayratio(D)
            displaygrid(D) % automatic displayzoom(D)
            displaylabels(D)
            displaycross(D)
            displayvalue(D)
            displaydecoration(D)
            updateselection(D,'all')
            
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
            
            % selection action
            % (shape)
            if ~isscalar(D.seldims)
                info.shape.poly = uimenu(m,'label','shape select poly','separator','on', ...
                    'callback',@(u,e)xx(@()set(D,'shapemode','poly')));
                info.shape.free = uimenu(m,'label','shape select free', ...
                    'callback',@(u,e)xx(@()set(D,'shapemode','free')));
                info.shape.rect = uimenu(m,'label','shape select rectangle', ...
                    'callback',@(u,e)xx(@()set(D,'shapemode','rect')));
                info.shape.ellipse = uimenu(m,'label','shape select ellipse', ...
                    'callback',@(u,e)xx(@()set(D,'shapemode','ellipse')));
                m1 = uimenu(m,'label','more selection modes...');
                info.shape.ring = uimenu(m1,'label','ring', ...
                    'callback',@(u,e)xx(@()set(D,'shapemode','ring')));
                info.shape.segment = uimenu(m1,'label','line segment', ...
                    'callback',@(u,e)xx(@()set(D,'shapemode','segment')));
                info.shape.openpoly = uimenu(m1,'label','open poly', ...
                    'callback',@(u,e)xx(@()set(D,'shapemode','openpoly')));
                info.shape.freeline = uimenu(m1,'label','free line', ...
                    'callback',@(u,e)xx(@()set(D,'shapemode','freeline')));
                set(info.shape.(D.shapemode),'checked','on')
            end
            % (rounding)
            m1 = uimenu(m,'label','selection rounding');
            info.selround.round0 = uimenu(m1,'label','no rounding', ...
                'callback',@(u,e)xx(@()set(D,'selround',0)));
            info.selround.round1 = uimenu(m1,'label','round to pixel', ...
                'callback',@(u,e)xx(@()set(D,'selround',1)));
            info.selround.round2 = uimenu(m1,'label','2 pixels', ...
                'callback',@(u,e)xx(@()set(D,'selround',2)));
            info.selround.round3 = uimenu(m1,'label','3 pixels', ...
                'callback',@(u,e)xx(@()set(D,'selround',3)));
            info.selround.round4 = uimenu(m1,'label','4 pixels', ...
                'callback',@(u,e)xx(@()set(D,'selround',4)));
            info.selround.round5 = uimenu(m1,'label','5 pixels', ...
                'callback',@(u,e)xx(@()set(D,'selround',5)));
            info.selround.round6 = uimenu(m1,'label','6 pixels', ...
                'callback',@(u,e)xx(@()set(D,'selround',6)));
            info.selround.round10 = uimenu(m1,'label','10 pixels', ...
                'callback',@(u,e)xx(@()set(D,'selround',10)));
            info.selround.more = uimenu(m1,'label','other...', ...
                'callback',@(u,e)xx(@()set(D,'selround',fn_input('pixels',7,1,100))));
            name = ['round' num2str(D.selround)];
            if ~isfield(info.selround,name), name = 'more'; end
            set(info.selround.(name),'checked','on')
            % (individual pixels or super-pixels)
            m1 = uimenu(m,'label','sub-regions');
            info.selmultin.none = uimenu(m1,'label','none', ...
                'callback',@(u,e)xx(@()set(D,'selmultin','none')));            
            info.selmultin.point = uimenu(m1,'label','use points', ...
                'callback',@(u,e)xx(@()set(D,'selmultin','point')));            
            info.selmultin.grid = uimenu(m1,'label','use grid', ...
                'callback',@(u,e)xx(@()set(D,'selmultin','grid')));            
            set(info.selmultin.(D.selmultin),'checked','on')
            
            % selection mode/display/reset
            info.seledit = uimenu(m,'label','advanced selection','separator','on','checked',onoff(D.seleditmode), ...
                'callback',@(u,e)xx(@()set(D,'seleditmode',~D.seleditmode)));
            m1 = uimenu(m,'label','display selection marks');
            uimenu(m1,'label','none', ...
                'callback',@(u,e)xx(@()set(D,'selshow','')));            
            uimenu(m1,'label','shape+number', ...
                'callback',@(u,e)xx(@()set(D,'selshow','shape+number')));            
            uimenu(m1,'label','cross', ...
                'callback',@(u,e)xx(@()set(D,'selshow','cross')));            
            info.selshow.shape = uimenu(m1,'label','shape','separator','on','checked',onoff(strfind(D.selshow,'shape')), ... 
                'callback',@(u,e)xx(@()set(D,'selshow','toggle shape')));            
            info.selshow.number= uimenu(m1,'label','number','checked',onoff(strfind(D.selshow,'number')), ... 
                'callback',@(u,e)xx(@()set(D,'selshow','toggle number')));            
            info.selshow.cross = uimenu(m1,'label','cross','checked',onoff(strfind(D.selshow,'cross')), ... 
                'callback',@(u,e)xx(@()set(D,'selshow','toggle cross')));                        
            info.selcolor = uimenu(m,'label','color selection marks', ...
                'callback',@(u,e)xx(@()set(D,'selcolor',~D.selcolor)));            
            uimenu(m,'label','reorder selections', ...
                'callback',@(u,e)xx(@()reorderselections(D)))
            info.sel(1) = uimenu(m,'label','reset selection display', ...
                'callback',@(u,e)xx(@()updateselection(D,'all')));
            info.sel(2) = uimenu(m,'label','reset selection', ...
                'callback',@(u,e)xx(@()updateselection(D.SI,'reset')));
            
            % navigation
            fn_propcontrol(D,'scrollwheel', ...
                'menu', ...
                {'parent',m,'label','scrollwheel zooming','separator','on'});
            fn_propcontrol(D,'navigation', ...
                {'menu' {'zoom' 'pan'}}, ...
                {'parent',m});
            
            % display
            % (permute x and y)
            info.permutexy = uimenu(m,'label','permute x and y','separator','on', ...
                'callback',@(u,e)xx(@()set(D,'permutexy',~D.permutexy)));
            % (features)
            m1 = uimenu(m,'label','features');
            info.features.crossshow = uimenu(m1,'label','show cross','checked',get(D.cross(1),'visible'), ...
                'callback',@(u,e)xx(@()set(D,'crossshow',~D.crossshow)));
            fn_propcontrol(D.cross,'Color',{'menu' 'w' 'k' 'r' 'b' 'y'},m1,'label','cross color');
            fn_propcontrol(D,'dolabels','menu',{m1,'label','labels'});
            fn_propcontrol(D,'doticks','menu',{m1,'label','ticks'});
            fn_propcontrol(D,'doxbar','menu',{m1,'label','scale bar'});
            % (color map)
            m1 = uimenu(m,'label','color map');
            info.cmap.gray = uimenu(m1,'label','gray', ...
                'callback',@(u,e)xx(@()set(D,'cmap','gray')));
            info.cmap.jet = uimenu(m1,'label','jet', ...
                'callback',@(u,e)xx(@()set(D,'cmap','jet')));
            info.cmap.mapgeog = uimenu(m1,'label','mapgeog', ...
                'callback',@(u,e)xx(@()set(D,'cmap','mapgeog')));
            info.cmap.mapgeogclip = uimenu(m1,'label','mapgeogclip', ...
                'callback',@(u,e)xx(@()set(D,'cmap','mapgeogclip')));
           info.cmap.mapclip = uimenu(m1,'label','mapclip', ...
                'callback',@(u,e)xx(@()set(D,'cmap','mapclip')));
            info.cmap.mapcliphigh = uimenu(m1,'label','mapcliphigh', ...
                'callback',@(u,e)xx(@()set(D,'cmap','mapcliphigh')));
            info.cmap.mapcliplow = uimenu(m1,'label','mapcliplow', ...
                'callback',@(u,e)xx(@()set(D,'cmap','mapcliplow')));
            info.cmap.vdaq = uimenu(m1,'label','vdaq', ...
                'callback',@(u,e)xx(@()set(D,'cmap','vdaq')));
            info.cmap.green = uimenu(m1,'label','green', ...
                'callback',@(u,e)xx(@()set(D,'cmap','green')));
            info.cmap.red = uimenu(m1,'label','red', ...
                'callback',@(u,e)xx(@()set(D,'cmap','red')));
            info.cmap.signcheck = uimenu(m1,'label','blue-yellow', ...
                'callback',@(u,e)xx(@()set(D,'cmap','signcheck')));            
            info.cmap.bluered = uimenu(m1,'label','blue-red', ...
                'callback',@(u,e)xx(@()set(D,'cmap','bluered')));            
            info.cmap.maporient = uimenu(m1,'label','signcheck', ...
                'callback',@(u,e)xx(@()set(D,'cmap','signcheck')));
            info.cmap.maporient = uimenu(m1,'label','maporient', ...
                'callback',@(u,e)xx(@()set(D,'cmap','maporient')));
            info.cmap.user = uimenu(m1,'label','user...', ...
                'callback',@(u,e)xx(@()set(D,'cmap','user')));
            info.logscale = uimenu(m1,'label','log scale','checked',onoff(D.logscale),'separator','on', ...
                'callback',@(u,e)xx(@()set(D,'logscale',~D.logscale)));
            set(info.cmap.(D.cmap),'checked','on') 
            % (clipping mode)
            m1 = uimenu(m,'label','clipping mode');
            info.clip.slice = uimenu(m1,'label','slice', ...
                'callback',@(u,e)xx(@()set(D,'clipmode','slice')));
            info.clip.data = uimenu(m1,'label','data', ...
                'callback',@(u,e)xx(@()set(D,'clipmode','data')));
            info.clip.link1 = uimenu(m1,'label','link1', ...
                'callback',@(u,e)xx(@()set(D,'clipmode','link1')));
            info.clip.link2 = uimenu(m1,'label','link2', ...
                'callback',@(u,e)xx(@()set(D,'clipmode','link2')));
            set(info.clip.(D.clipmode),'checked','on') 
            % (autoclip)
            uimenu(m,'label','do auto clip', ...
                'callback',@(u,e)xx(@()autoclip(D)));
            m1 = uimenu(m,'label','autoclip mode');
            info.autoclip.minmax = uimenu(m1,'label','minmax', ...
                'callback',@(u,e)xx(@()set(D,'autoclipmode','minmax')));
            info.autoclip.std1 = uimenu(m1,'label','1 STD', ...
                'callback',@(u,e)xx(@()set(D,'autoclipmode','std1')));
            info.autoclip.std2 = uimenu(m1,'label','2 STD', ...
                'callback',@(u,e)xx(@()set(D,'autoclipmode','std2')));
            info.autoclip.std3 = uimenu(m1,'label','3 STD', ...
                'callback',@(u,e)xx(@()set(D,'autoclipmode','std3')));
            info.autoclip.std5 = uimenu(m1,'label','5 STD', ...
                'callback',@(u,e)xx(@()set(D,'autoclipmode','std5')));
            info.autoclip.stdn = uimenu(m1,'label','N STD...', ...
                'callback',@(u,e)xx(@()set(D,'autoclipmode',[num2str(fn_input('N__STD',4,'stepper 1 1 10')) 'std'])));
            info.autoclip.prc01_999 = uimenu(m1,'label','[.1% 99.9%]', ...
                'callback',@(u,e)xx(@()set(D,'autoclipmode','prc.1')));
            info.autoclip.prc1_99 = uimenu(m1,'label','[1% 99%]', ...
                'callback',@(u,e)xx(@()set(D,'autoclipmode','prc1')));
            info.autoclip.prc5_95 = uimenu(m1,'label','[5% 95%]', ...
                'callback',@(u,e)xx(@()set(D,'autoclipmode','prc5')));
            info.autoclip.prct = uimenu(m1,'label','percentiles...', ...
                'callback',@(u,e)xx(@()set(D,'autoclipmode',fn_strcat(fn_input('percentiles',[5 95]),'prc','_'))));
            val = D.autoclipmode; val(val=='.') = [];
            if ~isfield(info.autoclip,val), val = fn_switch(strfind(val,'std'),'stdn','prct'); end
            set(info.autoclip.(val),'checked','on')
            % (autoclip origin)
            info.autoclip.full = uimenu(m1,'label','use full image','separator','on', ...
                'callback',@(u,e)xx(@()set(D,'autocliporig','full')));
            info.autoclip.curview = uimenu(m1,'label','use current view', ...
                'callback',@(u,e)xx(@()set(D,'autocliporig','curview')));
            info.autoclip.cursel = uimenu(m1,'label','use current selection', ...
                'callback',@(u,e)xx(@()set(D,'autocliporig','cursel')));
            info.autoclip.curselmask = uimenu(m1,'label','use current selection and mask', ...
                'callback',@(u,e)xx(@()set(D,'autocliporig','curselmask')));
            set(info.autoclip.(D.autocliporig),'checked','on')
            % (clip centering)
            fn_propcontrol(D,'clipcenter',{'menu' {0 1 []} {'center on 0' 'center on 1'}}, ...
                {'parent',m,'label','clip centering'});
            % (user clip)
            info.usrclip = uimenu(m,'label','user clip...', ...
                'callback',@(u,e)xx(@()set(D,'clip',fn_input('clip',D.clip))));
            % (binning)
            m1 = uimenu(m,'label','binning');
            info.bin.bin1 = uimenu(m1,'label','no binning', ...
                'callback',@(u,e)xx(@()set(D,'binning',0)));
            info.bin.bin2 = uimenu(m1,'label','2x2', ...
                'callback',@(u,e)xx(@()set(D,'binning',2)));
            info.bin.bin3 = uimenu(m1,'label','3x3', ...
                'callback',@(u,e)xx(@()set(D,'binning',3)));
            info.bin.bin4 = uimenu(m1,'label','4x4', ...
                'callback',@(u,e)xx(@()set(D,'binning',4)));
            info.bin.bin8 = uimenu(m1,'label','8x8', ...
                'callback',@(u,e)xx(@()set(D,'binning',8)));
            info.bin.bin16 = uimenu(m1,'label','16x16', ...
                'callback',@(u,e)xx(@()set(D,'binning',16)));
            info.bin.bin32 = uimenu(m1,'label','32x32', ...
                'callback',@(u,e)xx(@()set(D,'binning',32)));
            info.bin.other = uimenu(m1,'label','other...', ...
                'callback',@(u,e)xx(@()set(D,'binning',fn_input('binning',5,1,100))));
            name = ['bin' num2str(max(D.binning,1))];
            if ~isfield(info.bin,name), name = 'other'; end
            set(info.bin.(name),'checked','on')

            % more tools
            info.distline = uimenu(m,'label','distance tool','separator','on', ...
                'callback','fn_imdistline');
            uimenu(m,'label','show color bar', ...
                'callback',@(u,e)xx(@()setcolorbar(D)))

            % activedisplayImage object
            uimenu(m,'label','duplicate in new figure','separator','on', ...
                'callback',@(u,e)xx(@()duplicate(D)));
            uimenu(m,'label','display object ''D'' in base workspace', ...
                'callback',@(u,e)xx(@()assignin('base','D',D)));
            uimenu(m,'label','save image to file', ...
                'callback',@(u,e)xx(@()savepicture(D)));
            uimenu(m,'label','save image to clipboard', ...
                'callback',@(u,e)xx(@()savepicture(D,'clipboard')));
            m1 = uimenu(m,'label','more exports');
            uimenu(m1,'label','duplicate in ...', ...
                'callback',@(u,e)xx(@()duplicate(D,'user')));
            uimenu(m1,'label','save full display to file', ...
                'callback',@(u,e)xx(@()savepicture(D,'full')));
            uimenu(m1,'label','save full display to clipboard', ...
                'callback',@(u,e)xx(@()savepicture(D,'full','clipboard')));
            uimenu(m1,'label','export image to Matlab...', ...
                'callback',@(u,e)xx(@()fn_exportvar(D.currentdisplay)))
            
            
            % repair
            info.fix = uimenu(m,'label','repair scroll','separator','on', ...
                'callback',@(u,e)xx(@()repairscroll(D)));

            D.menuitems = info;
        end       
        function repairscroll(D)
            fn_scrollwheelregister(D.ha,'repair')
        end
        function delete(D)
            fn4D_dbstack
            % invoked when the invisible line is deleted: D.ha necessarily
            % still exists, but D.txt, D.buttons, D.menu might have
            % been deleted already
            if ishandle(D.ha), cla(D.ha,'reset'), end
            if ishandle(D.txt), delete(D.txt), end
            if ishandle(D.buttons), delete(D.buttons), end
            if ishandle(D.menu), delete(D.menu), end
            delete(D.listenaxpos)
            delete(D.slider(isvalid(D.slider)))
        end
    end
    
    % Display - except selection
    methods (Access='private')
        function displayxyperm(D)
            if D.permutexy
                set(D.ha,'ydir','normal')
                view(D.ha,90,90)
            else
                set(D.ha,'ydir','reverse')
                view(D.ha,0,90)
            end
        end
        function displayratio(D)
            % axis image?
            s = D.SI.sizes;
            D.doratio = all(s>1) && D.SI.nd>1 && ~iscell(D.SI.units{1}) && ~iscell(D.SI.units{2}) && strcmp(D.SI.units{1},D.SI.units{2});
            if D.doratio && isempty(D.SI.units{1})
                % do not set the aspect ratio if this would result into a
                % too narrow rectangle
                L = D.SI.sizes .* D.SI.grid(:,1)';
                D.doratio = (max(L)/min(L) < 5);
            end
            if D.doratio
                set(D.ha,'dataAspectRatioMode','manual', ...
                    'dataAspectRatio',[1 1 1])
            else
                set(D.ha,'dataAspectRatioMode','auto')
            end            
        end
        function displaygrid(D)
            fn4D_dbstack
            % changes axis, image xdata and ydata, and cross extremities
            s = D.SI.sizes;
            grid = D.SI.grid;
            if D.SI.nd==1
                s(2) = D.SI.sizesplus(1);
                grid(2,:) = [1 0];
            end
            
            % scaling
            xrange = [1 s(1)]*grid(1,1)+grid(1,2);
            yrange = [1 s(2)]*grid(2,1)+grid(2,2);
            set(D.img,'xdata',xrange,'ydata',yrange)
            
            % axis
            sidefact = grid(:,1)/2; sidefact(s==1)=.5;
            side = sidefact * [-1 1];
            xlim = grid(1,2)+[1 s(1)]*grid(1,1)+side(1,:);
            ylim = grid(2,2)+[1 s(2)]*grid(2,1)+side(2,:);
            D.oldaxis = [xlim; ylim];
            
            %zooming
            if D.SI.nd==1, D.axis(2,:) = ylim; end
            displayzoom(D) % set axis -> automatic display update (+ move buttons, position scale bar)
            
            % change cross!
            set(D.cross(1),'YData',ylim)
            set(D.cross(2),'XData',xlim)
        end       
        function displaylabels(D)
            fn4D_dbstack
            labels = D.SI.labels;
            units = D.SI.units;
            if D.SI.nd==1
                labels{2} = '';
                units{2} = '';
            end
            for i=1:2
                if ~isempty(units{i}) && ~iscell(units{i}), labels{i} = [labels{i} ' (' units{i} ')']; end
            end
            if D.doticks
                set(D.ha,'xtickmode','auto','ytickmode','auto')
            else
                set(D.ha,'xtick',[],'ytick',[])
            end
            if D.dolabels
                xlabel(D.ha,labels{1})
                ylabel(D.ha,labels{2})
            else
                xlabel(D.ha,'')
                ylabel(D.ha,'')
            end
        end
        function displayscalebar(D)
            fn4D_dbstack
            if ~D.doxbar || ~D.doratio
                set(D.scalebar,'visible','off')
                fn4D_enable('off',D.listenaxpos)
                return
            else
                set(D.scalebar,'visible','on')
                fn4D_enable('on',D.listenaxpos)
            end
            % find a nice size for bar: in specialized function
            barsize = BarSize(D.ha);
            % label: use units if any - must be the same for x and
            % y!
            label = num2str(barsize);
            if ~isempty(D.SI.units{1})
                units = D.SI.units;
                %                 if ~strcmp(units{1},units{2})
                %                     error('units are not the same for x and y axis')
                %                 end
                label = [label ' ' units{1}];
            end
            % positions
            barorigin = fn_coordinates(D.ha,'b2a',[20 10]','position');
            barpos = [barorigin barorigin+[barsize 0]'];
            textpos = mean(barpos,2) + ...
                fn_coordinates(D.ha,'b2a',[0 10]','vector');
            % set properties
            set(D.scalebar(1),'xdata',barpos(1,:),'ydata',barpos(2,:))
            set(D.scalebar(2),'position',textpos,'string',label)
            if D.movescale
                set(D.scalebar,'hittest','on','buttondownfcn', ...
                    @(hobj,evnt)fn_moveobject(D.scalebar,'latch'))
            else
                set(D.scalebar,'hittest','off')
            end
        end        
        function displaydata(D)
            fn4D_dbstack
            if isempty(D.img), return, end % this can happen at init
            slice = D.SI.slice;
            
            % no data -> no display
            if isempty(slice) || ~isfield(slice,D.dataflag)
                set(D.img,'cdata',NaN)
                return
            end
            
            % prevent automatic display
            D.freeze = true;
            
            if length(slice)>1
                disp('activedisplayImage can display only one selection at a time')
                slice = slice(1);
            end
            % ignore 'active' flag!
            im = fn_float(slice.(D.dataflag));
            
            % apply image operation
            im = fn_imageop(im,D.opdef);
            if ~isreal(im), im = abs(im); end
            D.currentdisplay = im;
            
            % update the default clipping
            autoclipupdate(D)
            
            % display image
            D.freeze = false;
            displaydata2(D)
        end
        function displaydata2(D)
            % call this function if data has not changed but display needs
            % to be updated (usually because of change in clip)
            im = D.currentdisplay;
            if isempty(im), return, end % this can happen at init
            
            % true color image + mask
            im = ColorImage(im,D.cmapval,D.logscale,D.channelcolors,D.clip,get(D.hf,'color'));
            seldimsnum = D.seldims-'w';
            selectionmarks = D.SI.selection.getselset(seldimsnum).singleset;
            if strcmp(D.autocliporig,'curselmask') && ~isempty(D.currentselection)
                mask = selectionmarks(D.currentselection).mask;
                if D.binning>1, mask = fn_bin(mask,D.binning); end
                im = ApplyMask(im,mask);
            end
            if D.selcolor
                im = HighlightSelection(im,selectionmarks,D.seldims);
            end
            
            % display
            set(D.img,'CData',permute(im,[2 1 3]))
        end       
        function displaycross(D)
            fn4D_dbstack
            ij2 = D.SI.ij2;
            % scaling and translation
            pt = IJ2AX(D.SI,ij2);
            set(D.cross(1),'XData',pt([1 1]))
            if D.SI.nd==2
                set(D.cross(2),'YData',pt([2 2]))
                set(D.cross(3),'XData',pt(1),'YData',pt(2))
            end
        end       
        function displayvalue(D)
            fn4D_dbstack
            im = D.currentdisplay;
            ij = D.SI.ij;
            if D.SI.nd==1
                y = get(D.cross(2),'ydata');
                ij(2) = max(1,min(size(im,2),round(y(1))));
            end
            if isempty(im)
                set(D.txt,'String','')
            else
                if D.binning>1, ij = 1+floor((ij-1)/D.binning); end
                if all(ij<=[size(im,1); size(im,2)])
                    % the opposite can happen if D is displaying another
                    % field than data in D.SI.slice, and its size does not
                    % match that of data
                    set(D.txt,'String', ...
                        ['val(' num2str(ij(1)) ',' num2str(ij(2)) ')=' ...
                        num2str(im(ij(1),ij(2)),'%.3g')])
                end
            end
        end        
        function displayzoom(D)
            fn4D_dbstack
            % must be applied AFTER displaygrid
            oldax = D.oldaxis;
            zoom = IJ2AX(D.SI,D.SI.zoom);
            if D.SI.nd==1, zoom(2,:) = D.axis(2,:); end
            % min and max to stay within range
            zoom(:,1) = max(oldax(:,1),zoom(:,1));
            zoom(:,2) = min(oldax(:,2),zoom(:,2));
            if any(diff(zoom,1,2)<=0)
                disp('new zoom is outside of range - do zoom reset')
                zoom = oldax;
            end
            D.axis = zoom;
        end   
        function setcolorbar(D)
            if D.logscale && D.clip(1)>0, clipshow = log10(D.clip); else clipshow = D.clip; end
            D.colorbar = fn_colorbar(clipshow,D.cmap);
        end
        function displaydecoration(D)
            delete(findall(D.ha,'Tag','fn4D_deco'))
            D.hdeco = [];
            deco = D.SI.decoration;
            if isempty(deco), return, end
            for k=1:length(deco.t)
                if ~isequal(deco.t(k).dims,[1 2])
                    disp 'only 2-dimensional decorations can be displayed'
                    continue
                end
                sel = deco.t(k).set;
                for i=1:length(sel)
                    poly = sel(i).poly;
                    if all(strcmp({poly.type},'point2D'))
                        pp = [poly.points];
                        D.hdeco = line(pp(1,:),pp(2,:),'linestyle','none', ...
                            ... 'color','k','marker','o','markerfacecolor','w','markersize',3, ...
                            'color','y','marker','+','markersize',4, ...
                            'parent',D.ha,'tag','fn4D_deco','hittest','off');
                    else
                        sel = convert(sel,'poly2D');
                        poly = [sel.poly];
                        D.hdeco = gobjects(1,length(poly));
                        for j=1:length(poly)
                            pp = poly(j).points;
                            D.hdeco(j) = line(pp(1,:),pp(2,:),'color','y', ...
                                'parent',D.ha);
                        end
                    end
                end
            end
        end
    end
    
    % Display selection
    methods (Access='private')
        function displayselection(D,flag,ind,value)
            fn4D_dbstack
            if isempty(D.selshow)
                delete(findobj(D.ha,'tag','ActDispIm_Sel'))
                return
            end
            
            % some params
            si = D.SI;
            seldimsnum = D.seldims-'w';
            selectionmarks = si.selection.getselset(seldimsnum).singleset;
            nsel = length(selectionmarks);
            
            % display set...
            if fn_ismemberstr(flag,{'all','reset'})
                % 'findobj' allows a cleanup when some objects were not
                % removed correctly
                delete(findobj(D.ha,'tag','ActDispIm_Sel'))
                D.seldisp = cell(1,nsel);
                isel = 1;
                for k=1:nsel
                    displayonesel(D,k,'new',isel);
                    if selectionmarks(k).active, isel = isel+1; end
                end
                return
            end
            
            % or display update
            if ~isempty(D.curselprev) && ~isempty(strfind(D.selshow,'number'))
                set(D.seldisp{D.curselprev}(1),'color','w')
            end
            switch flag
                case 'new'
                    isel = cumsum([selectionmarks.active]);
                    for idx=ind
                        displayonesel(D,idx,'new',isel(idx)); 
                    end
                case {'add','change','affinity'}
                    % might be several indices
                    for k=ind, displayonesel(D,k,'pos'); end
                case 'remove'
                    delete([D.seldisp{ind}])
                    D.seldisp(ind) = [];
                    nsel = length(D.seldisp);
                    if nsel==0, return, end
                    updateselorderdisplay(D)
                case 'active'
                    % might be several indices
                    for k=ind, displayonesel(D,k,'active'), end
                    updateselorderdisplay(D)
                case 'reorder'
                    perm = value;
                    D.seldisp = D.seldisp(perm);
                    updateselorderdisplay(D)
                case 'indices'
                    % nothing to do
            end
            if ~isempty(D.currentselection) && ~isempty(strfind(D.selshow,'number'))
                set(D.seldisp{D.currentselection}(1),'color','r')
            end
        end        
        function displayonesel(D,k,flag,varargin)
            % function displayonesel(D,k,'new',isel)
            % function displayonesel(D,k,'pos')
            % function displayonesel(D,k,'isel',isel)
            % function displayonesel(D,k,'active')
            % function displayonesel(D,k,'edit')
            
            % flags
            [flagnew flagpos flagisel flagactive flagedit] = ...
                fn_flags('new','pos','isel','active','edit',flag);
                
            % Values
            seldimsnum = D.seldims-'w';
            selectionmarks = D.SI.selection.getselset(seldimsnum).singleset;
            selij = selectionmarks(k);
            if flagnew || flagedit || flagpos
                if isscalar(D.seldims)
                    selij2 = convert(selij,'line1D');
                    if D.SI.nd == 1
                        sel = IJ2AX(D.SI,selij2);
                        poly = [sel.poly];
                        points = {poly.points};
                        orthsiz = D.oldaxis(2,:);
                    else
                        % i can't do better than by hand!!!
                        points = {selij2.poly.points};
                        npart = length(points);
                        for i=1:npart
                            points{i} = points{i}*D.SI.grid(seldimsnum,1) + D.SI.grid(seldimsnum,2);
                        end
                        orthdim = 3-seldimsnum;
                        orthsiz = [.5 D.SI.sizes(orthdim)+.5]*D.SI.grid(orthdim,1) + D.SI.grid(orthdim,2);
                    end
                    npart = length(points);
                    for i=1:npart
                        % line
                        points{i} = [points{i}([1 1 2 2 1]) NaN; orthsiz([1 2 2 1 1]) NaN];
                    end
                    polygon = [points{:}];
                else
                    sel = IJ2AX(D.SI,selij);
                    selij2 = convert(selij,'poly2D');
                    sel2 = IJ2AX(D.SI,selij2);
                    polygon = sel2.poly.points;
                end
                center = [nanmean(polygon(1,:)) nanmean(polygon(2,:))];
            end
            if flagnew || flagedit || flagactive || flagisel
                if selij.active
                    colors = fn_colorset;
                    col = colors(mod(k-1,size(colors,1))+1,:);
                    linestyle = '-';
                    visible = 'on';
                else
                    col = 'k';
                    linestyle = '--';
                    visible = 'off';
                end
            end
            if flagnew || flagisel
                isel = varargin{1};
                str = num2str(isel);
            end
            
            % Create / update objects
            if flagnew
                hl = [];
                if strfind(D.selshow,'number')
                    hl(end+1) = text(center(1),center(2),str, ...
                        'Parent',D.ha,'color','w','visible',visible, ...
                        'horizontalalignment','center','verticalalignment','middle', ...
                        'color',fn_switch(k==D.currentselection,'r','w'));
                end
                if strfind(D.selshow,'shape')
                    hl(end+1) = line(polygon(1,:),polygon(2,:),'Parent',D.ha, ...
                        'Color',col,'LineStyle',linestyle, ...
                        'UserData',k); % set user data because this line will be used when in seledit mode
                end
                if strfind(D.selshow,'cross')
                    hl(end+1) = line(center(1),center(2),'Parent',D.ha, ...
                        'Color',col,'LineStyle','none', ...
                        'Marker','+','MarkerSize',4);
                end
                set(hl,'tag','ActDispIm_Sel','HitTest','off')
                D.seldisp{k} = hl;
            else
                hl = D.seldisp{k};
                i=1; ht=[]; hs=[]; hc=[];
                if strfind(D.selshow,'number'), ht=hl(i); i=i+1; end
                if strfind(D.selshow,'shape'),  hs=hl(i); i=i+1; end
                if strfind(D.selshow,'cross'),  hc=hl(i); i=i+1; end
                he = hl(i:end);
                if flagpos
                    set(ht,'position',center)
                    set(hs,'xdata',polygon(1,:),'ydata',polygon(2,:))
                    set(hc,'xdata',center(1),'ydata',center(2))
                elseif flagisel
                    set(ht,'string',str)
                    set([hs hc he],'color',col)
                elseif flagactive
                    set(hs,'color',col,'linestyle',linestyle)
                    set([ht hc he],'visible',visible)
                end
            end
            
            % Advanced selection mode (in this mode, D.seldisp = [ht hl he]
            % because D.selshow = 'number+shape')
            if ~D.seleditmode || flagisel || flagactive, return, end 
            desc = [];
            switch selectionmarks(k).type
                case {'poly2D','mixed','point2D','line2D'} % TODO: not sure about 'point2D'
                    polymark = polygon;
                case 'rect2D'
                    polymark = polygon(:,1:4); % the 5th point of polygon is a repetition of the 1st one
                    desc = [sel.poly.points' sel.poly.vectors'];
                case {'ellipse2D' 'ring2D'}
                    c = sel.poly.points;
                    u = sel.poly.vectors;
                    e = sel.poly.logic;
                    polymark = [c-u c+u];
                    desc = {c u e};
                otherwise
                    error programming
            end
            if flagnew || flagedit
                % right now, hl has 2 elements: number and shape
                set(hl(2),'hittest','on','buttondownfcn', ...
                    @(h,evnt)seleditaction(D,get(h,'userdata'),'line'))
                hl(3) = line(polymark(1,:),polymark(2,:),'Parent',D.ha, ...
                    'Color',col,'tag','ActDispIm_Sel', ...
                    'LineStyle','none','marker','.', ...
                    'UserData',k,'hittest','on','buttondownfcn',...
                    @(h,evnt)seleditaction(D,get(h,'userdata'),'point'));
                if ~isempty(desc),
                    setappdata(hl(3),'description',desc)
                end
                D.seldisp{k} = hl;
            else
                set(hl(3),'xdata',polymark(1,:),'ydata',polymark(2,:));
                if ~isempty(desc)
                    setappdata(hl(3),'description',desc)
                end
            end
        end  
        function updateselorderdisplay(D)
            isel = 0;
            seldimsnum = D.seldims-'w';
            selectionmarks = D.SI.selection.getselset(seldimsnum).singleset;
            for k=1:length(selectionmarks)
                set(D.seldisp{k}(2:end),'userdata',k)
                if selectionmarks(k).active
                    isel = isel+1;
                    displayonesel(D,k,'isel',isel)
                end
            end
        end     
        function seleditstart(D)
            for k=1:length(D.seldisp), displayonesel(D,k,'edit'), end
        end
        function seleditupdateslice(D,ind,dp)
            seldimsnum = D.seldims-'w';
            if nargin==3
                % specify the shape transformation rather than the new shapes
                mov = affinityND('translate2D',dp);
                ijmov = AX2IJ(D.SI,mov);
                updateselection(D.SI,seldimsnum,'affinity',ind,ijmov)
            else
                seldimsnum = D.seldims-'w';
                selectionmarks = D.SI.selection.getselset(seldimsnum).singleset;
                for k=1:length(ind)
                    i = ind(k);
                    hl = D.seldisp{i}(2:3);
                    seltype = selectionmarks(i).type;
                    switch seltype
                        case {'poly2D' 'point2D' 'line2D'}
                            polygon = [get(hl(2),'xdata'); get(hl(2),'ydata')];
                            selax = selectionND(seltype,polygon);
                        case 'mixed'
                            polygon = [get(hl(2),'xdata'); get(hl(2),'ydata')];
                            selax = selectionND('poly2D',polygon);
                        case 'rect2D'
                            desc = getappdata(hl(2),'description');
                            % if shape was only translated, information on
                            % center is desuete -> update
                            xdata = get(hl(2),'xdata'); ydata = get(hl(2),'ydata');
                            desc(1:2) = [xdata(1) ydata(1)];
                            % twist the rectangle the sign of its
                            % width/height has been changed
                            if desc(3)<0, desc(1)=desc(1)+desc(3); desc(3)=-desc(3); end
                            if desc(4)<0, desc(2)=desc(2)+desc(4); desc(4)=-desc(4); end
                            % update selection
                            selax = selectionND('rect2D',desc);
                        case {'ellipse2D' 'ring2D'}
                            desc = getappdata(hl(2),'description');
                            % if shape was only translated, information on
                            % center is desuete -> update
                            desc{1} = [nanmean(get(hl(2),'xdata')); nanmean(get(hl(2),'ydata'))];
                            selax = selectionND(seltype,desc);
                        otherwise
                            error programming
                    end
                    sel(k) = AX2IJ(D.SI,selax); %#ok<AGROW>
                end
                updateselection(D.SI,seldimsnum,'change',ind,sel)
            end
        end
        function seleditend(D)
            for i=1:length(D.seldisp)
                delete(D.seldisp{i}(3:end))
                D.seldisp{i}(3:end) = [];
                set(D.seldisp{i}(2),'hittest','off','buttondownfcn','')
            end
        end
    end
    
    % Update routines
    methods (Access='private')
        function updateselection(D,flag,ind,value)
            fn4D_dbstack(['updateselection ' flag])
            if nargin<3, ind=[]; end
            if nargin<4, value=[]; end
            % current selection
            seldimsnum = D.seldims-'w';
            selectionmarks = D.SI.selection.getselset(seldimsnum).singleset;
            nsel = length(selectionmarks);
            D.curselprev = D.currentselection;
            switch flag
                case 'reorder'
                    perm = value;
                    if isscalar(ind)
                        D.currentselection = find(perm==ind);
                    else
                        D.currentselection = nsel;
                    end
                case {'remove','reset','all'}
                    D.currentselection = nsel;
                    if D.currentselection==0, D.currentselection=[]; end
                case 'indices'
                    % no change in D.currentselection
                otherwise
                    D.currentselection = ind(end);
            end
            % compute selectionlabels
            si = D.SI;
            if fn_ismemberstr(flag,{'new','add'}) && numel(D.selectionlabels)==prod(si.sizes)
                % just need to add the last selection(s)
                if any(ind>255) && isa(D.selectionlabels,'uint8')
                    D.selectionlabels = uint16(D.selectionlabels);
                end
                for idx=ind
                    D.selectionlabels(selectionmarks(idx).dataind) = idx;
                end
            elseif ~isempty(selectionmarks) ...
                    && ~all(si.sizes(seldimsnum)==selectionmarks(1).datasizes)
                % size mismatch, don't update D.selectionlabels;
                % normally, data should be changed soon too (change
                % is si.sizes) and then sizes will match again
                D.selectionlabels = [];
            else
                % restart from zero
                type = fn_switch(nsel<=255,'uint8','uint16');
                if D.SI.nd==1
                    D.selectionlabels = zeros([si.sizes 1],type);
                else
                    D.selectionlabels = zeros([si.sizes(seldimsnum) 1],type);
                end
                for k=1:nsel
                    D.selectionlabels(selectionmarks(k).dataind) = k;
                end
            end
            % update selection display (sel edit mode)
            if D.selcolor
                displaydata2(D)
            end
            if strcmp(flag,'indices'), return, end
            % update selection display
            displayselection(D,flag,ind,value)
        end
    end
    
    % GET/SET clip
    methods
        function set.clip(D,clip)
            fn4D_dbstack
            clip = full(clip);
            if isequal(clip,D.clip) || isempty(clip), return, end
            if size(clip,2)~=2, error('clip should have 2 columns'), end
            % correct bad values for clip
            d = diff(clip,1,2);
            clip(d==0,:) = fn_add(clip(d==0,1),[-1 1]);
            bad = (d<0 | isnan(d));
            clip(bad,:) = repmat([0 1],sum(bad),1);
            % assign value
            D.clip = clip;
            if fn_ismemberstr(D.clipmode,{'link1','link2'})
                D.CL.clip = clip;
            end
            % also assign it to the axes
            set(D.ha,'CLim',clip(1,:));
            % update color bar if any
            if ~isempty(D.colorbar)
                if ~isvalid(D.colorbar)
                    D.colorbar = [];
                else
                    D.colorbar.clip = clip;
                end
            end
            % update display
            if ~D.freeze, displaydata2(D), end
        end   
        function set.clipmode(D,clipmode)
            fn4D_dbstack
            if strcmp(clipmode,D.clipmode), return, end
            if ~fn_ismemberstr(clipmode,{'slice','data','link1','link2'})
                error('wrong ''%s''',clipmode)
            end
            oldclipmode = D.clipmode;
            D.clipmode = clipmode;
            % check mark in uicontextmenu
            if ~isempty(D.menuitems)
                set(D.menuitems.clip.(oldclipmode),'checked','off')
                set(D.menuitems.clip.(clipmode),'checked','on')
            end
            % specific actions
            if fn_ismemberstr(oldclipmode,{'link1','link2'})
                % cancel previous cliplink and listener
                disconnect(D,D.CL), delete(D.C2D) %#ok<*MCSUP>
            end
            switch clipmode
                case 'slice'
                    autoclip(D)
                case {'link1','link2'}
                    D.CL = cliplink.find(clipmode,D.clip);
                    D.clip = D.CL.clip;
                    D.C2D = connect_listener(D.CL,D,'ChangeClip', ...
                        @(cl,evnt)clipfromlink(D,D.CL));
            end
        end      
        function set.clipcenter(D,newval)
            D.clipcenter = newval;
            % update clip if necessary
            autoclipupdate(D)
        end
        function set.autoclipmode(D,newval)
            fn4D_dbstack
            % update value
            if fn_ismemberstr(newval,{'fit' 'mM' 'minmax'})
                newval = 'minmax';
            else
                token = regexpi(newval,'^([\d.]*)st*d$','tokens');
                if isempty(token), token = regexpi(newval,'^st*d([\d.]*)$','tokens'); end
                if ~isempty(token)
                    % Nstd
                    nstd = token{1}{1};
                    if isempty(nstd), nstd = '1'; end
                    if isnan(str2double(nstd)), error 'autoclipmode flag must be ''minmax'' or ''Nstd''', end
                    newval = ['std' token{1}{1}];
                else
                    % prcA_B
                    tokens = regexpi(newval,'^prc([\d.]*)[_-]{0,1}([\d.]*)$','tokens');
                    if isempty(tokens), error 'autoclipmode flag must be ''minmax'', ''Nstd'' or ''prcA_B''', end
                    low = str2double(tokens{1}{1}); high = str2double(tokens{1}{2});
                    if isnan(high), high = 100-low; end
                    newval = sprintf('prc%g_%g',low,high);
                end
            end
            oldval = D.autoclipmode;
            if strcmp(newval,oldval), return, end
            D.autoclipmode = newval;
            % check mark in uicontextmenu
            items = D.menuitems.autoclip;
            oldval(oldval=='.')=[]; newval(newval=='.')=[];
            if ~isfield(items,oldval), oldval = fn_switch(strfind(oldval,'std'),'stdn','prct'); end
            if ~isfield(items,newval), newval = fn_switch(strfind(newval,'std'),'stdn','prct'); end
            set(D.menuitems.autoclip.(oldval),'checked','off')
            set(D.menuitems.autoclip.(newval),'checked','on')
            % update clip
            D.autoclipvalue = [];
            autoclip(D)
        end
        function set.autocliporig(D,newval)
            oldval = D.autocliporig;
            % update value
            if strcmp(oldval,newval), return, end
            if ~fn_ismemberstr(newval,{'full' 'curview' 'cursel' 'curselmask'})
                error 'wrong value for ''autocliporig'''
            end
            D.autocliporig = newval;
            % update mark in uicontextmenu
            set(D.menuitems.autoclip.(oldval),'checked','off')
            set(D.menuitems.autoclip.(newval),'checked','on')
            % update what is necessary
            autoclipupdate(D)
            if any(strcmp({oldval,newval},'curselmask')), displaydata2(D), end
        end
        function clipfromlink(D,CL)
            D.clip = CL.clip;
        end       
        function clip = get.autoclipvalue(D)
            % auto-compute if necessary
            if isempty(D.currentdisplay)
                clip = [-1 1];
            elseif isempty(D.autoclipvalue)
                
                % prepare data from which clipping range will be computed
                im = D.currentdisplay;
                s = size(im);
                switch D.autocliporig
                    case 'full'
                        % use full image, i.e. nothing to do
                        im = fn_reshapepermute(im,{[1 2] 3});
                    case 'curview'
                        % use current view
                        zoom = D.SI.zoom;
                        if D.SI.nd==1, zoom(2,:) = [-Inf Inf]; end
                        if ~all(isinf(zoom(:)))
                            if D.binning>1, zoom = 1+((zoom-1)/D.binning); end
                            zoom(:,1) = max(round(zoom(:,1)),1);
                            zoom(:,2) = min(round(zoom(:,2)),s(1:2)');
                            if all(diff(zoom,1,2)>0)
                                im = im(zoom(1,1):zoom(1,2),zoom(2,1):zoom(2,2),:);
                            end
                        end
                        im = fn_reshapepermute(im,{[1 2] 3});
                    case {'cursel' 'curselmask'}
                        % use current selection, if any!
                        if ~isempty(D.currentselection) && ~isscalar(D.seldims)
                            seldimsnum = D.seldims-'w';
                            selectionmarks = D.SI.selection.getselset(seldimsnum).singleset;
                            sel = selectionmarks(D.currentselection);
                            im = fn_imvect(im,sel.mask,'vector');
                        end
                end
                bad = isinf(im(:));
                if full(any(bad)), im(bad)=NaN; end
                
                % compute
                clipflag = D.autoclipmode;
                if ~isempty(D.clipcenter), clipflag = [clipflag '[' num2str(D.clipcenter) ']']; end
                clip = fn_clip(im,clipflag,'getrange');
                
                % assign
                D.autoclipvalue = clip;
            else
                clip = D.autoclipvalue;
            end
        end       
        function autoclip(D)
            % uses automatic clipping
            D.clip = D.autoclipvalue;
        end    
        function autoclipupdate(D)
            % reset autoclip value (such that the autoclip function will
            % have to re-compute it)
            % then actually calls autoclip function if 'slice' mode
            D.autoclipvalue = [];
            if strcmp(D.clipmode,'slice'), autoclip(D); end
        end
    end
    
    % GET/SET image operations
    methods
        function yfun = get.yfun(D)
            yfun = D.opdef.user;
        end
        function binning = get.binning(D)
            binning = D.opdef.xbin;
        end
        function xlow = get.xlow(D)
            xlow = D.opdef.xlow;
        end
        function xhigh = get.xhigh(D)
            xhigh = D.opdef.xhigh;
        end
        function set.yfun(D,yfun)
            if ~isempty(yfun) && ~isa(yfun,'function_handle')
                error('''yfun'' should be a function handle with one argument')
            end
            D.opdef.user = yfun;
        end
        function  set.binning(D,binning)
            binning = max(1,binning);
            oldbin = max(1,D.binning);
            if ~isfield(D.opdef,'xbin'), error programming, end
            if oldbin==binning, return, end
            % change value and update display
            D.opdef.xbin = binning;
            % update menu marks
            items = D.menuitems.bin;
            oldname = ['bin' num2str(oldbin)];
            if isfield(items,oldname), set(items.(oldname),'checked','off'), end
            newname = ['bin' num2str(binning)];
            if isfield(items,newname), set(items.(newname),'checked','on'), end
        end
        function  set.xlow(D,xlow)
            D.opdef.xlow = xlow;
        end
        function  set.xhigh(D,xhigh)
            D.opdef.xhigh = xhigh;
        end
        function set.opdef(D,op)
            D.opdef = op;
            if ~D.freeze, displaydata(D), end
        end
        function set.menustayvisible(D,val)
            D.menustayvisible = val;
            set(D.menuitems.pin,'checked',onoff(val))
            if val
                set(D.menu,'pos',get(D.hf,'currentpoint'),'visible','on')
            end
        end
    end
    
    % GET/SET other
    methods
        function set.permutexy(D,b)
            if b == D.permutexy, return, end
            D.permutexy = b;
            displayxyperm(D)
        end
        function set.doratio(D,b)
            D.doratio = b;
        end
        function b = get.crossshow(D)
            b = boolean(get(D.cross(1),'visible'));
        end
        function set.crossshow(D,b)
            if b==D.crossshow, return, end
            on_off = on_off(b);
            set(D.menuitems.features.crossshow,'checked',on_off)
            set(D.cross,'visible',on_off)
        end
        function set.dataflag(D,str)
            D.dataflag = str;
            displaydata(D)
        end
        function set.cmap(D,cm)
            cmold = D.cmap;
            if ischar(cm) && strcmp(cm,'user')
                try
                    answer = inputdlg('define color map','colormap',1,{'jet(256)'});
                    if isempty(answer), return, end
                    cm = evalin('base',answer{1});
                    cmnew = 'user';
                catch %#ok<CTCH>
                    return
                end
            elseif ischar(cm)
                cmnew=cm; 
                if strcmp(cmnew,cmold), return, end
                if strcmp(cmnew,'vdaq')
                    cm = vdaqcolors;
                elseif fn_ismemberstr(cmnew,{'mapgeogclip'})
                    cm = feval(cmnew);
                else
                    try
                        cm = feval(cmnew,256);
                        if ~isequal(size(cm),[256 3]), error('wrong color map ''%s''',cmnew), end
                    catch
                        error('wrong color map ''%s''',cmnew)
                    end                        
                end
            else
                cmnew='user'; 
            end
            D.cmap = cmnew;
            D.cmapval = cm;
            colormap(D.ha,cm) % normally useless, but can be usefull when D.img.cdata is hacked
            
            % update check marks
            set(D.menuitems.cmap.(cmold),'checked','off')
            set(D.menuitems.cmap.(cmnew),'checked','on')            
            
            % update color bar if any
            if ~isempty(D.colorbar)
                if ~isvalid(D.colorbar)
                    D.colorbar = [];
                else
                    D.colorbar.cmap = cm;
                end
            end
            
            % update display
            displaydata2(D)
        end
        function set.logscale(D,b)
            if b==D.logscale, return, end
            D.logscale = b;
            set(D.menuitems.logscale,'checked',onoff(b))
            D.displaydata2
        end
        function set.channelcolors(D,a)
            if ~isempty(a) && (size(a,2)~=3 || size(a,1)>3)
                error('''channelcolors'' is a bad channel to colors conversion matrix')
            end
            if any(sum(a,1)>1)
                error 'bad colors'
            end
            D.channelcolors = a;
            % update display
            if ~D.freeze, displaydata2(D), end
        end
        function set.shapemode(D,shapemodenew)
            if ~fn_ismemberstr(shapemodenew,{'poly','free','rect','ellipse','ring','segment','openpoly','freeline'})
                error('wrong shape selection mode ''%s''',shapemodenew)
            end
            shapemodeold = D.shapemode;
            D.shapemode = shapemodenew;
            % update check marks
            try set(D.menuitems.shape.(shapemodeold),'checked','off'), end %#ok<TRYNC>
            try set(D.menuitems.shape.(shapemodenew),'checked','on'), end %#ok<TRYNC>
        end
        function set.seldims(D,seldimsnew)
            % the 'seldims' property indicates in which dimensions
            % selections are drawn, it can be 'x', 'y' or 'xy'
            if seldimsnew == D.seldims, return, end
            if ~fn_ismemberstr(seldimsnew,{'x' 'y' 'xy'}), error 'wrong selection dimension flag', end
            % reset the selection that became invalid
            oldseldimsnum = D.seldims-'w'; 
            nsel = length(D.SI.selection.getselset(oldseldimsnum).singleset);
            if nsel, updateselection(D.SI,oldseldimsnum,'remove',1:nsel), end
            % update 'seldims' property
            D.seldims = seldimsnew;
            % update local menu (selection shape and selection edit items)
            try initlocalmenu(D), end %#ok<TRYNC>
        end
        function set.seleditmode(D,seleditmodenew)
            if seleditmodenew==D.seleditmode, return, end
            D.selshow = 'number+shape'; % this is mandatory, as in 'edit' mode, the shape is used (and dealing with the possible absence of text is boring)
            D.seleditmode = seleditmodenew;
            % update check marks and display
            if D.seleditmode
                set(D.menuitems.seledit,'checked','on')
                seleditstart(D)
            else
                set(D.menuitems.seledit,'checked','off')
                seleditend(D)
            end
        end
        function set.selround(D,roundnew)
            roundold = D.selround;
            if roundnew==D.selround, return, end
            D.selround = roundnew;
            % update check marks and display
            name = ['round' num2str(roundold)];
            if ~isfield(D.menuitems.selround,name), name = 'more'; end
            set(D.menuitems.selround.(name),'checked','off')
            name = ['round' num2str(roundnew)];
            if ~isfield(D.menuitems.selround,name), name = 'more'; end
            set(D.menuitems.selround.(name),'checked','on')
        end
        function set.selmultin(D,selmultinnew)
            selmultinold = D.selmultin;
            if strcmp(selmultinnew,D.selmultin), return, end
            D.selmultin = selmultinnew;
            % update check marks and display
            set(D.menuitems.selmultin.(selmultinold),'checked','off')
            set(D.menuitems.selmultin.(selmultinnew),'checked','on')
        end
        function set.scaledisplay(D,flag)
            switch flag
                case 'tick'
                    D.doticks = true;
                    D.doxbar = false;
                case 'xbar'
                    D.doticks = false;
                    D.doxbar = true;
                case ''
                    D.doticks = false;
                    D.doxbar = false;
            end
        end
        function set.dolabels(D,b)
            D.dolabels = b;
            displaylabels(D)
        end
        function set.doticks(D,b)
            D.doticks = b;
            displaylabels(D)
        end
        function set.doxbar(D,b)
            D.doxbar = b;
            displayscalebar(D)
        end
        function set.selshow(D,flag)
            % Input -> create accurate 'flag' and 'flagc'
            if islogical(flag) || isnumeric(flag)
                flag = fn_switch(flag,'number+shape','');
            end
            flagc = fn_strcut(flag,'+, ');
            if length(flagc)==2 && strcmp(flagc{1},'toggle')
                flagc0 = fn_strcut(D.selshow,'+, ');
                flagc = setxor(flagc0,flagc{2});
            end
            flag = fn_strcat(flagc,'+'); 
            % Check
            allflags = {'shape' 'number' 'cross'};
            if ~all(ismember(flagc,allflags))
                error 'wrong ''selshow'' flag'
            end
            % Any change?
            if strcmp(flag,D.selshow), return, end
            % Update property, display and menu items
            D.selshow = flag;
            displayselection(D,'all')
            for i=1:length(allflags)
                f = allflags{i};
                set(D.menuitems.selshow.(f),'checked',onoff(ismember(f,flagc)))
            end
        end
        function set.selcolor(D,b)
            if b==D.selcolor, return, end
            D.selcolor = b;
            displaydata2(D)
            if b
                set(D.menuitems.selcolor,'checked','on')
            else
                set(D.menuitems.selcolor,'checked','off')
            end
        end
        function set.movescale(D,b)
            D.movescale = b;
            displayscalebar(D)
        end
        function axis = get.axis(D)
            fn4D_dbstack
            axis = [get(D.ha,'xLim'); get(D.ha,'yLim')];
        end
        function set.axis(D,axis)
            fn4D_dbstack
            oldax=D.axis;
            if all(axis==oldax), return, end
            % set axis
            set(D.ha,'xLim',axis(1,:),'yLim',axis(2,:));
            % re-position scale bar
            displayscalebar(D)
            % change slider parameters
            if any(axis(1,:)~=oldax(1,:)), slideraxis(D,1), end
            if any(axis(2,:)~=oldax(2,:)), slideraxis(D,2), end
            % re-compute clipping fitting the shown area
            autoclipupdate(D)
        end
        function slideraxis(D,k)
            if all(D.axis(k,:)==D.oldaxis(k,:))
                set(D.slider(k),'visible','off')
            else
                set(D.slider(k),'visible','on', ...
                    'min',D.oldaxis(k,1),'max',D.oldaxis(k,2),'value',D.axis(k,:))
            end
        end 
        function set.navigation(D,val)
            fn_ismemberstr(val,{'zoom' 'pan'},'doerror')
            D.navigation = val;
        end
        function set.scrollwheel(D,val)
            val = boolean(val);
            if D.scrollwheel==val, return, end
            D.scrollwheel = val;
            fn_scrollwheelregister(D.ha,onoff(val))
        end
    end
    
    % Events (bottom-up: mouse)
    methods (Access='private')
        function Mouse(D,flag)
            fn4D_dbstack
            % different normal mouse actions are:
            % LEFT BUTTON
            % - point                                   -> change cursor
            % - area                                    -> zoom to region
            % - double-click (or click outside of axis) -> zoom reset
            % MIDDLE BUTTON
            % - click in region, hold and type a number -> reorder selections
            % - point/area                              -> add point/area to current selection
            % - click outside                           -> cancel current selection
            % RIGHT BUTTON
            % - point in region                         -> hide/show selection
            % - point/area                              -> add new selection
            % - click outside                           -> cancel all selections
            %
            % flag is any of 'outside', 'pointonly'
            
            % selection type
            oldselectiontyp = D.oldselectiontype;
            selectiontype = get(D.hf,'selectiontype');
            pointonly = (nargin==2 && strcmp(flag,'pointonly'));
            
            % special case - click outside of axis
            ax = axis(D.ha);
            point =  get(D.ha,'CurrentPoint'); point = point(1,[1 2])';
            if (nargin==2 && strcmp(flag,'outside')) ...
                    || point(1)<ax(1) || point(1)>ax(2) ...
                    || point(2)<ax(3) || point(2)>ax(4)
                oldselectiontyp = selectiontype;
                selectiontype = 'outside';
            end
            
            % store current selection type
            D.oldselectiontype = selectiontype;
            
            % 1D or 2D mode?
            do1d = (D.SI.nd==1);
            
            % open or closed selection / shape selection mode
            mouseselmode = fn_switch(D.shapemode,'openpoly','poly','freeline','free',D.shapemode);
            TYPE = fn_switch(D.shapemode,{'poly','free'},'poly2D','rect','rect2D', ...
                'ellipse','ellipse2D','ring','ring2D', ...
                'segment','line2D',{'openpoly','freeline'},'openpoly2D');
            
            % shortcut
            si = D.SI;
            hb = D.ha;
            seldimsnum = D.seldims-'w';
            selectionmarks = D.SI.selection.getselset(seldimsnum).singleset;
            nsel = length(selectionmarks);
            
            % GO!
            switch selectiontype
                case 'normal'                                       % CHANGE VIEW AND/OR MOVE CURSOR
                    if pointonly
                        rect = [row(point) 0 0];
                    elseif strcmp(D.navigation,'zoom') || all(D.axis(:)==D.oldaxis(:))
                        rect = fn_mouse(hb,'rect-^');
                    else
                        rect = [];
                    end
                    if isempty(rect)
                        moved = pan(D);                             % pan
                        if moved
                            return
                        else
                            rect = [row(point) 0 0];
                        end
                    end
                    if all(abs(rect(3:4))'>diff(D.axis,1,2)/50)     % zoom in
                        if do1d
                            D.axis(2,:) = rect(2)+[0 rect(4)];
                            rect = AX2IJ(si,rect(1)+[0 rect(3)]);
                        else
                            rect = AX2IJ(si,[rect(1)+[0 rect(3)]; rect(2)+[0 rect(4)]]);
                        end
                        si.zoom = rect;
                    elseif abs(rect(3))>diff(D.axis(1,:))/50        % zoom in (x only)
                        if do1d
                            rect = AX2IJ(si,rect(1)+[0 rect(3)]);
                        else
                            rect = AX2IJ(si,[rect(1)+[0 rect(3)]; rect(2)+[0 rect(4)]]);
                        end
                        si.zoom(1,:) = rect(1,:);
                    elseif abs(rect(4))>diff(D.axis(2,:))/50        % zoom in (y only)
                        if do1d
                            D.axis(2,:) = rect(2)+[0 rect(4)];
                        else
                            rect = AX2IJ(si,[rect(1)+[0 rect(3)]; rect(2)+[0 rect(4)]]);
                            si.zoom(2,:) = rect(2,:);
                        end
                    else                                            % change xy
                        if do1d
                            set(D.cross(2),'ydata',rect([2 2]))
                            set(D.cross(3),'ydata',rect(2))
                            point = AX2IJ(si,rect(1));
                        else
                            point = AX2IJ(si,rect([1 2])');
                        end
                        si.ij2 = point;
                        if ~isempty(D.usercallback)                 % user callback
                            feval(D.usercallback,D)
                        end
                    end
                case {'extend' 'alt'}                               % EDIT SELECTION / NEW SELECTION
                    ksel = BelongsToSelection(D,point);

                    if strcmp(selectiontype,'extend') && ksel    	% reorder selections
                        typenumber(D,ksel);
                        return
                    end
                    
                    if strcmp(selectiontype,'extend') && ~strcmp(D.selmultin,'none')
                        % since selections must have specific shapes, we
                        % cannot add pixels to the current selection
                        return
                    end
                    
                    % point selection?
                    if pointonly
                        ispt = true;
                    elseif isscalar(D.seldims)
                        rect = fn_mouse(hb,'rect-');
                        ispt = (rect(2+seldimsnum)<diff(D.axis(seldimsnum,:))/200);
                    else
                        polyax = fn_mouse(hb,[mouseselmode '-']);
                        selax = selectionND(TYPE,polyax);
                        ispt = ispoint(selax,min(diff(D.axis,1,2))/200);
                    end
                    
                    if strcmp(selectiontype,'alt') && ispt && ksel  % show/hide
                        updateselection(si,seldimsnum,'active',ksel, ...
                            ~selectionmarks(ksel).active)
                        return
                    end
                    
                    %                                               % selection add or new 
                    if ispt && ~isempty(D.usercallback)             % user callback
                        feval(D.usercallback,D)
                        return
                    elseif ispt
                        if do1d, point = point(1); end
                        point = AX2IJ(si,point);
                        if isscalar(D.seldims)
                            point = point(seldimsnum);
                            if D.selround
                                % round and transform the point into a line around one or several pixel(s)
                                left = floor((point-.5)/D.selround)*D.selround + .5;
                                line = left + [0 D.selround];
                                sel = selectionND('line1D',line);
                            else
                                sel = selectionND('point1D',point);
                            end
                        else
                            if D.selround
                                % transform point into square around one or several pixel(s)
                                corner = floor((point-.5)/D.selround)*D.selround + .5;
                                poly = fn_add(corner,[0 0 1 1; 0 1 1 0]*D.selround);
                                sel = selectionND('poly2D',poly);
                            else
                                sel = selectionND('point2D',point);
                            end
                        end
                    else
                        roundback = fn_switch(D.selround && strcmp(D.selmultin,'none'),D.selround,1); % in rounding + selmultin mode, we define the selection in the binned image, hence roundback is 1 rather than D.selround
                        if isscalar(D.seldims)
                            if do1d
                                line = AX2IJ(si,rect(1) + [0 rect(3)]);
                            else
                                rect = AX2IJ(si,[rect(1) + [0 rect(3)]; rect(2) + [0 rect(4)]]);
                                line = rect(fn_switch(D.seldims,'x',1,'y',2),:);
                            end
                            if D.selround
                                line = round((line-.5)/D.selround)*roundback + .5;
                                if diff(line)<=0, return, end
                            end
                            sel = selectionND('line1D',line);
                        else
                            if D.selround && strcmp(D.shapemode,'rect')
                                rect = AX2IJ(si,[polyax(1) + [0 polyax(3)]; polyax(2) + [0 polyax(4)]]);
                                rect = round((rect-.5)/D.selround)*roundback + .5;
                                if any(diff(rect,1,2)<=0), return, end
                                sel = selectionND('rect2D',[rect(:,1)' diff(rect,1,2)']);
                            else
                                sel = AX2IJ(si,selax);
                            end
                        end
                    end
                    switch selectiontype
                        case 'extend'
                            updateselection(si,seldimsnum,'add',D.currentselection,sel);
                        case 'alt'
                            if ~strcmp(D.selmultin,'none') && ~ispt
                                % convert selection to individual pixel
                                % selections
                                ComputeInd(sel,floor(si.sizes(seldimsnum)/max(D.selround,1)))
                                mask = sel.mask;
                                if D.selround<=1, mask = mask & ~D.selectionlabels; end % take only pixels that are not selected yet
                                if ~any(mask(:)), return, end
                                if isscalar(D.seldims)
                                    points = find(mask)';
                                else
                                    [ii jj] = find(mask);
                                    points = [ii'; jj'];
                                end
                                npoint = size(points,2);
                                sel(npoint) = sel; % pre-allocation
                                switch D.selmultin
                                    case 'point'
                                        TYPE = fn_switch(isscalar(D.seldims),'point1D','point2D');
                                        if D.selround>1, points = (points-.5)*D.selround+.5; end
                                        for i=1:npoint, sel(i) = selectionND(TYPE,points(:,i)); end
                                    case 'grid'
                                        if isscalar(D.seldims)
                                            left = (points-1)*D.selround + .5;
                                            line = fn_add(permute(left,[1 3 2]),[0 D.selround]);
                                            for i=1:npoint, sel(i) = selectionND('line1D',line(1,:,i)); end
                                        else
                                            corner = (points-1)*D.selround + .5;
                                            rect = fn_add(permute(corner,[1 3 2]),[0 0 1 1; 0 1 1 0]*D.selround);
                                            for i=1:npoint, sel(i) = selectionND('poly2D',rect(:,:,i)); end
                                        end
                                end
                            end
                            updateselection(si,seldimsnum,'new',[],sel);
                    end
                case 'open'                                         % MISC
                    switch oldselectiontyp
                        case 'normal'
                            if ~isempty(D.usercallback)             % user callback
                                feval(D.usercallback,D)
                            else                                    % zoom out
                                if do1d
                                    D.axis(2,:) = D.oldaxis(2,:);
                                    si.zoom = [-Inf Inf];
                                else
                                    si.zoom = [-Inf Inf; -Inf Inf];
                                end
                            end
                        case 'extend'
                            % better not to use it: interferes with poly selection
                        case 'alt'
                            ksel = BelongsToSelection(D,point);
                            if ksel && ~isempty(D.selshow)          % reorder (current selection -> first)
                                perm = [ksel setdiff(1:nsel,ksel)];
                                updateselection(si,seldimsnum,'active',ksel,true);
                                updateselection(si,seldimsnum,'reorder',ksel,perm);
                            end
                    end
                case 'outside'                                      % UNDO ZOOM / SELECTION
                    switch oldselectiontyp
                        case 'normal' 
                            rect = fn_mouse(hb,'rect-');
                            if all(rect(3:4))                       % zoom in
                                if do1d
                                    rect = AX2IJ(si,rect(1)+[0 rect(3)]);
                                else
                                    rect = AX2IJ(si,[rect(1)+[0 rect(3)]; rect(2)+[0 rect(4)]]);
                                end
                                si.zoom = rect;
                            else                                    % zoom out
                                if do1d
                                    si.zoom = [-Inf Inf];
                                else
                                    si.zoom = [-Inf Inf; -Inf Inf];
                                end
                            end
                        case {'extend','alt'}
                            switch oldselectiontyp
                                case 'alt'                          % unselect all regions
                                    if nsel>3
                                        answer = questdlg('Remove all selections?','confirmation','Yes','No','Yes');
                                        if ~strcmp(answer,'Yes'), return, end
                                    end
                                    updateselection(si,seldimsnum,'remove',1:nsel);
                                case 'extend'                       % unselect last region
                                    updateselection(si,seldimsnum,'remove',D.currentselection);
                            end
                    end
                otherwise
                    error programming
            end
        end        
        function seleditaction(D,ind,flag)
            htl = D.seldisp{ind};
            hl = htl(2:end);
            [flagpt flaglin] = fn_flags({'point','line'},flag);
            p = get(D.ha,'currentpoint'); p = p(1,1:2)';
            polymark = [get(hl(2),'xdata'); get(hl(2),'ydata')];
            seldimsnum = D.seldims-'w';
            selectionmarks = D.SI.selection.getselset(seldimsnum).singleset;
            switch get(D.hf,'selectiontype')
                case 'normal'               % MOVE POINT
                    shapetype = selectionmarks(ind).type;
                    switch shapetype
                        case {'poly2D' 'mixed' 'line2D' 'openpoly2D'}
                            % note that first and last point in polygon are
                            % the same!
                            if flagpt
                                % closest point
                                dist = sum(fn_add(polymark,-p).^2);
                                [dum idx] = min(dist); idx = idx(1); %#ok<*ASGLU>
                                if idx==1 && ~fn_ismemberstr(shapetype,{'line2D' 'openpoly2D'})
                                    % need to move both the first and last point (which is a repetition of the first)
                                    idx=[1 size(polymark,2)]; 
                                end 
                            else
                                % closest segment (in fact, closest line)
                                a = polymark(:,1:end-1);
                                b = polymark(:,2:end);
                                ab = b-a;
                                ab2 = sum(ab.^2);
                                ap = fn_add(p,-a);
                                abap = ab(1,:).*ap(2,:)-ab(2,:).*ap(1,:);
                                dist = abs(abap) ./ ab2;
                                [dum idx] = min(dist); idx = idx(1);
                                polymark = [a(:,1:idx) p b(:,idx:end)];
                                set(hl,'xdata',polymark(1,:),'ydata',polymark(2,:))
                                idx = idx+1;
                            end
                            fn_moveobject(hl,'point',idx)
                            seleditupdateslice(D,ind)
                        case 'rect2D'
                            desc = getappdata(hl(2),'description');
                            x = desc(1); y = desc(2);
                            w = desc(3); h = desc(4);
                            x2 = x+w; y2 = y+h;
                            if flagpt
                                % move corner
                                pol = [x x2 x2 x; y y y2 y2]; % anti-clockwise from (x,y)
                                dist = sum(fn_add(pol,-p).^2);
                                [dum idx] = min(dist); idx = idx(1);
                            else
                                % move edge
                                dist = abs([p(2)-y p(1)-x2 p(2)-y2 p(1)-x]);
                                [dum idx] = min(dist);
                            end
                            col = get(hl(1),'color');
                            set(hl,'color',.5*[1 1 1])
                            chgrectangle(D.ha,hl,flagpt,idx,desc)
                            fn_buttonmotion({@chgrectangle,D.ha,hl,flagpt,idx,desc},D.hf);
                            set(hl,'color',col)
                            seleditupdateslice(D,ind)
                        case {'ellipse2D' 'ring2D'}
                            desc = getappdata(hl(2),'description');
                            if flagpt
                                % closest of two anchor points
                                dist = sum(fn_add(polymark,-p).^2);
                                [dum idx] = min(dist); idx = idx(1);
                            elseif strcmp(shapetype,'ellipse2D')
                                % eccentricity
                                idx = 0;
                            else
                                % eccentricity or secondary radius?
                                polygon = [get(hl(1),'xdata'); get(hl(1),'ydata')];
                                dist = sum(fn_add(polygon,-p).^2);
                                [dum idx] = min(dist);
                                if idx<length(polygon)/2
                                    % eccentricity
                                    idx = 0;
                                else
                                    % secondary radius
                                    idx = -1;
                                end
                            end
                            col = get(hl(1),'color');
                            set(hl,'color',.5*[1 1 1])
                            chgellipse(D.ha,hl,idx,desc)
                            fn_buttonmotion({@chgellipse,D.ha,hl,idx,desc},D.hf);
                            set(hl,'color',col)
                            seleditupdateslice(D,ind)
                        otherwise
                            error programming
                    end
                case 'extend'               % MOVE SHAPE
                    if flagpt
                        dp = fn_moveobject(htl);
                        seleditupdateslice(D,ind,dp)
                    elseif flaglin
                        dp = fn_moveobject([D.seldisp{:}]);
                        seleditupdateslice(D,1:length(D.seldisp),dp) % move all shapes
                    end
                case 'alt'                  % REMOVE
                    if fn_ismemberstr(selectionmarks(ind).type, ...
                            {'poly2D','mixed'}) && flagpt
                        % closest point -> remove vertex
                        dist = sum(fn_add(polymark,-p).^2);
                        [dum idx] = min(dist); idx = idx(1);
                        if idx==1, idx=[1 size(polymark,2)]; end
                        polymark(:,idx) = [];
                        selax = selectionND('poly2D',polymark);
                        sel = AX2IJ(D.SI,selax);
                        updateselection(D.SI,'change',ind,sel)
                    else
                        % replace the whole shape
                        set(hl,'visible','off')
                        mouseselmode = fn_switch(D.shapemode,'openpoly','poly','freeline','free',D.shapemode);
                        TYPE = fn_switch(D.shapemode,{'poly','free'},'poly2D','rect','rect2D', ...
                            'ellipse','ellipse2D','ring','ring2D', ...
                            'segment','line2D',{'openpoly','freeline'},'openpoly2D');
                        polyax = fn_mouse(D.ha,mouseselmode,'select new shape');
                        selax = selectionND(TYPE,polyax);
                        sel = AX2IJ(D.SI,selax);
                        updateselection(D.SI,'change',ind,sel)
                        set(hl,'visible','on')
                    end
            end
        end
        function ksel = BelongsToSelection(D,pointax)
            if D.SI.nd==1, pointax = pointax(1); end
            ij = round(AX2IJ(D.SI,pointax));
            if isempty(D.selectionlabels) || any(ij<=0) || any(ij(:)'>D.SI.sizes)
                ksel = 0;
            else
                switch D.seldims
                    case 'x'
                        ksel = D.selectionlabels(ij(1));
                    case 'y'
                        ksel = D.selectionlabels(ij(2));
                    case 'xy'
                        ksel = D.selectionlabels(ij(1),ij(2));
                end
            end
        end        
        function typenumber(D,ksel)
            set(D.hf,'userdata',struct('ksel',ksel,'str',''), ...
                'windowkeypressfcn',@(x,evnt)typenumbertype(D,evnt), ...
                'windowbuttonmotionfcn',@(x,evnt)typenumberexec(D,'move'), ...
                'windowbuttonupfcn',@(x,evnt)typenumberexec(D,'stop'))
        end        
        function typenumbertype(D,evnt)
            disp(evnt)
            info = get(D.hf,'userdata');
            info.str = [info.str evnt.Character];
            set(D.hf,'userdata',info)
        end        
        function typenumberexec(D,flag)
            info = get(D.hf,'userdata');
            seldimsnum = D.seldims-'w';
            selectionmarks = D.SI.selection.getselset(seldimsnum).singleset;
            switch flag
                case 'move'
                    point =  get(D.ha,'CurrentPoint'); point = point(1,[1 2])';
                    ksel = BelongsToSelection(D,point);
                    % if there is no change, return
                    if isequal(ksel,info.ksel), return, end
                    set(D.hf,'userdata',struct('ksel',ksel,'str',''));
                case 'stop'
                    set(D.hf,'userdata',[],'windowkeypressfcn','', ...
                        'windowbuttonmotionfcn','','windowbuttonupfcn','')
                    % if no permutation, at least trick to change active sel
                    if info.ksel && isempty(info.str)
                        updateselection(D.SI,seldimsnum,'active',info.ksel,true)
                        return
                    end
            end
            % no action if no region or no valid string
            ksel = info.ksel;
            k = str2double(info.str);
            if ksel==0 || isnan(k), return, end
            % reorder selections
            % selection 'ksel' will be placed just before the active
            % selection 'k'; however we need first to determine the
            % number of this selection when also considering inactive
            % selections
            nsel = length(selectionmarks);
            if k==0 || k>nsel
                k = nsel;
            else
                activesels = find(cat(1,selectionmarks.active));
                k = activesels(k);
            end
            if k<ksel
                perm = [1:k-1 ksel k:ksel-1 ksel+1:nsel];
            elseif k==ksel
                return
            else
                perm = [1:ksel-1 ksel+1:k ksel k+1:nsel];
            end
            updateselection(D.SI,seldimsnum,'reorder',ksel,perm)
        end        
        function polyedit(D,hl,ksel,kp,flag)
            poly = [get(hl,'xdata'); get(hl,'ydata')];
            np = size(poly,2);
            p = get(D.ha,'currentpoint'); p = p(1,1:2)';
            selectiontype = get(D.hf,'selectiontype');
            seldimsnum = D.seldims-'w';
            switch selectiontype
                case 'normal'               % MOVE POINT
                    % move point
                    if strcmp(flag,'line')
                        % create point on line
                        poly = [poly(:,1:kp) p poly(:,kp+1:end)];
                        set(hl,'xdata',poly(1,:),'ydata',poly(2,:));
                        kp = kp+1;
                    end
                    poly = fn_buttonmotion({@moveedge,D.ha,hl,kp},D.hf);
                    updateselpoly(D,'change',ksel,poly)
                case 'extend'               % MOVE POLY
                    poly = fn_buttonmotion({@movepoly,D.ha,hl,poly,p},D.hf);
                    updateselpoly(D,'change',ksel,poly)
                case 'alt'                  % DELETE
                    f1 = find(isnan(poly(1,1:kp(1)-1)),1,'last');
                    f2 = kp(1)+find(isnan(poly(1,kp(1)+1:end)),1,'first');
                    if isempty(f1), f1=0; end
                    if isempty(f2), f2=np+1; end
                    if f2-f1>5 && strcmp(flag,'point')
                        % remove point
                        if length(kp)==2
                            poly = poly(:,[1:kp(1)-1  kp(1)+1:kp(2)-1 ...
                                kp(1)+1 kp(2)+1:end]);
                        else
                            poly = [poly(:,1:kp-1) poly(:,kp+1:end)];
                        end
                        updateselpoly(D,'change',ksel,poly)
                    elseif f1~=0 || f2~=np+1
                        % remove component
                        poly = [poly(:,1:f1-1) poly(:,f2+1:end)];
                        updateselpoly(D,'change',ksel,poly)
                    else
                        % remove selection
                        updateselection(D.SI,seldimsnum,'remove',ksel)
                    end
            end
        end               
        function updateselpoly(D,flag,ksel,polyax)
            poly = AX2IJ(D.SI,polyax);
            updateselection(D.SI,D.seldims-'w',flag,ksel,selectionND('poly2D',poly));
        end
        function moved = pan(D)
            fn4D_dbstack
            set(D.hf,'pointer','hand')
            p0 = get(D.ha,'currentpoint'); p0 = p0(1,1:2)';
            moved = false;
            fn_buttonmotion(@chgaxis,D.hf)
            function chgaxis
                moved = true;
                p = get(D.ha,'currentpoint'); p = p(1,1:2)';
                movax = p0-p;
                sidedist = D.oldaxis - D.axis;
                movax = max(sidedist(:,1),min(sidedist(:,2),movax));
                D.axis = fn_add(D.axis,movax);
            end
            set(D.hf,'pointer','arrow')
            if moved
                if D.SI.nd==1
                    D.SI.zoom = AX2IJ(D.SI,D.axis(1,:));
                else
                    D.SI.zoom = AX2IJ(D.SI,D.axis);
                end
            end
        end
        function movecross(D,il)
            fn4D_dbstack
            if ~strcmp(get(D.hf,'selectiontype'),'normal')
                % execute callback for axes
                Mouse(D)
                return
            end           
            set(D.hf,'pointer',fn_switch(il,1,'left',2,'top',3,'cross'))
            do1d = (D.SI.nd==1);
            si = D.SI;
            anymove = false;
            fn_buttonmotion(@movecrosssub,D.hf)
            set(D.hf,'pointer','arrow')
            if ~anymove
                % execute callback for axes
                Mouse(D,'pointonly')
                return
            end
            function movecrosssub
                anymove = true;
                p = get(D.ha,'currentpoint'); p = p(1,1:2);
                if do1d
                    if il~=1
                        set(D.cross(2),'ydata',p([2 2]))
                        set(D.cross(3),'ydata',p(2))
                    end
                    if il~=2
                        si.ij2 = AX2IJ(si,p(1));
                    end
                else
                    ij2 = AX2IJ(si,p([1 2])');
                    switch il
                        case 1 % move x only
                            si.ij2(1) = ij2(1);
                        case 2 % move y only
                            si.ij2(2) = ij2(2);
                        case 3 % move x and y
                            si.ij2 = ij2;
                    end
                end
            end
        end
    end

    % Events (bottom-up: scroll wheel, buttons, sliders)
    methods (Access='private')
        function Scroll(D,nscroll)
            si = D.SI;
            p = get(D.ha,'currentpoint');
            origin = AX2IJ(si,p(1,1:2)');
            zoom = si.zoom;
            zoom = [max(.5,zoom(:,1)) min(si.sizes'+.5,zoom(:,2))]; 
            zoomfactor = 1.5^nscroll;
            newzoom = fn_add(origin,zoomfactor*fn_subtract(zoom,origin));
            si.zoom = newzoom;
        end
        function redbutton(D,flag)
            fn4D_dbstack
            if nargin<2
                % pressed red button
                flag = get(D.hf,'selectiontype');
                kchannel = 1;
            else
                if ~strcmp(flag,'yellow'), error programming, end
                flag = 'normal';
                kchannel = fn_switch(get(D.hf,'selectiontype'),'normal',2,'extend',3,0);
                if kchannel==0, return, end
            end
            switch flag
                case 'normal'       % change clip
                    for i=size(D.clip,1)+1:kchannel, D.clip(i,:) = D.clip(1,:); end
                    clip0 = D.clip(kchannel,:);
                    p0 = get(D.hf,'currentpoint');
                    ht = uicontrol('style','text','position',[2 2 200 17],'parent',D.hf);
                    moveclip(D,ht,p0,clip0,kchannel); % this displays the bottom-left numbers
                    % change clip
                    fn_buttonmotion({@moveclip,D,ht,p0,clip0,kchannel},D.hf)
                    delete(ht)
                case 'extend'       % toggle advanced selection
                    D.seleditmode = ~D.seleditmode;
                case 'open'         % use default clipping
                    autoclip(D)
            end
        end        
        function chgzoom(D,flag,U)
            fn4D_dbstack
            ax = D.axis;
            switch flag
                case 'x'
                    ax(1,:) = get(U,'value');
                case 'y'
                    ax(2,:) = get(U,'value');
                otherwise
                    error programming
            end
            if U.sliderscrolling
                set(D.ha,'xLim',ax(1,:),'yLim',ax(2,:));
            else
                % trick to detect a change in axis, and auto-update display
                set(D.ha,'xLim',ax(1,:)+[0 ax(1,2)*100*eps],'yLim',ax(2,:)+[0 ax(2,2)*100*eps]);
                if D.SI.nd==1
                    D.SI.zoom = AX2IJ(D.SI,ax(1,:));
                else
                    D.SI.zoom = AX2IJ(D.SI,ax);
                end
            end
        end       
        function duplicate(D,hobj)
            if nargin<2, hobj=figure; end
            if strcmp(hobj,'user')
                answer = inputdlg({'Figure','Subplot'},'Location of duplicate',1,{'1','111'});
                if isempty(answer), disp interrupted, return, end
                figure(str2double(answer{1}))
                hobj = subplot(str2double(answer{2}));
            end
            activedisplayImage(D.SI,'in',hobj, ...
                'dataflag',D.dataflag, ...
                'clipmode',D.clipmode,'clip',D.clip,'cmap',D.cmap,'logscale',D.logscale, ...
                'selshow',D.selshow, ...
                'seleditmode',D.seleditmode,'shapemode',D.shapemode, ...
                'crossshow',D.crossshow,'doticks',D.doticks,'dolabels',D.dolabels,'doxbar',D.doxbar);
        end
    end
    
    % (public for access by fn_buttonmotion)
    methods
        function moveclip(D,ht,p0,clip0,kchannel)
            % new clip
            p = get(D.hf,'currentpoint');
            dp = p-p0;
            r = clip0(2)-clip0(1);
            if ~isempty(D.clipcenter)
                clip0 = D.clipcenter + [-1 1]*r/2;
                dp = [-1 1]*(dp(2)-dp(1))/2;
            end
            FACT = 1/100;
            clipk = clip0 + dp*(r*FACT);
            c = mean(clipk);
            x = clipk(2)-c;
            thr = r/10;
            if x<=thr % avoid clipmin to become more than clipmax
                clipk = c + [-1 1] * thr;
            end
            
            % display
            set(ht,'string',sprintf('min: %.3f,  max: %.3f',clipk(1),clipk(2)))
            D.clip(kchannel,:) = clipk;
        end
    end
    
    % Events (top-down: listeners)
    methods
        function updateDown(D,~,evnt)
            fn4D_dbstack(['S2D ' evnt.flag])
            switch evnt.flag
                case 'sizes'
                    D.currentdisplay = [];
                    displayratio(D)
                    displaygrid(D) % automatic displayzoom(D)
                    % try to update D.selectionlabels; will work only if
                    % the datasizes in selection are already set correctly
                    updateselection(D,'indices')
                case 'sizesplus'
                    if D.SI.ndplus>1 || (D.SI.ndplus==1 && D.SI.sizesplus(1)>3)
                        error('activedisplayImage class can display only two-dimensional data slices, in at most 3 channels')
                    end
                case 'slice'
                    displaydata(D)
                    displayvalue(D)
                case 'grid'
                    displaygrid(D) % automatic displayzoom(D)
                    % conversion IJ2AX has changed -> update selection
                    % display
                    displayselection(D,'all')
                case 'labels'
                    displaylabels(D)
                case 'units'
                    displayratio(D)
                    displayscalebar(D)
                    displaylabels(D)
                case 'ij2'
                    displaycross(D)
                case 'ij'
                    displayvalue(D)
                case 'zoom'
                    displayzoom(D)
                case 'selection'
                    seldimsnum = D.seldims-'w';
                    if ~strcmp(evnt.selflag,'reset') 
                        if ~any(ismember(seldimsnum,evnt.dims)), return, end % orthogonal selection, for example when image show (space x time) and selections displayed are along space, temporal selection are just ignored
                        if ~all(ismember(seldimsnum,evnt.dims)), disp 'selection change is not compatible with display', return, end
                    end
                    updateselection(D,evnt.selflag,evnt.ind,evnt.value)
                    if D.selcolor, displaydata2(D), end
                    % note: if selflag is 'indices' it might be that the
                    % datasizes has changed in the new selection before
                    % D.SI.sizes has changed (see above)
                case 'decoration'
                    displaydecoration(D)
            end
        end
    end
    
    % Other functions
    methods 
        function reorderselections(D)
            seldimsnum = D.seldims-'w';
            selectionmarks = D.SI.selection.getselset(seldimsnum).singleset;
            if isscalar(seldimsnum)
                selectionmarks = convert(selectionmarks,'line1D');
            else
                selectionmarks = convert(selectionmarks,'poly2D');
            end
            poly = [selectionmarks.poly];
            center = fn_map(@(x)nanmean(x,2),{poly.points},'array')';
            [c ord] = sortrows(fliplr(center)); % reorder first according to y, second to x
            updateselection(D.SI,seldimsnum,'reorder',[],ord);
        end
    end
    
    % Misc
    methods
        function access(D) %#ok<MANU>
            keyboard
        end
        function savepicture(D,varargin)
            [dofull doclipboard] = fn_flags({'full' 'clipboard'},varargin);
            if ~doclipboard
                fname = fn_savefile('*.png','Select image file where to save.');
            end
            if dofull
                % temporarilly hide buttons
                set([D.buttons D.txt],'visible','off')
                c = onCleanup(@()set([D.buttons D.txt],'visible','on'));
                rect = fn_pixelpos(D.ha,'strict');
                if doclipboard
                    x = getframe(D.hf,rect);
                    imclipboard('copy',x.cdata)
                else
                    fn_savefig(D.hf,fname,rect)
                end
            else
                z = D.SI.zoom;
                z = [max([1; 1],floor(z(:,1))) min(D.SI.sizes(:),ceil(z(:,2)))];
                x = D.currentdisplay(z(1,1):z(1,2),z(2,1):z(2,2),:);
                if doclipboard
                    x = fn_clip(x,D.clip,D.cmapval);
                    imclipboard('copy',permute(x,[2 1 3]))
                else
                    fn_saveimg(x,fname,D.clip,1,D.cmapval)
                end
            end
        end
    end
    
end

%----------------
% TOOLS: display
%----------------

function im = ColorImage(im,cmap,logscale,channelcolors,clip,nancolor)

s = size(im);
s3 = size(im,3);
if length(s)>3 || s3>3
    error('activedisplayImage can display only 2D slices, in at most 3 channels')
end
im = reshape(im,s(1)*s(2),s3);

if isempty(channelcolors)
    if s3==1, channelcolors = 1; else channelcolors = eye(s3,3); end
elseif any(size(channelcolors)~=[s3 3])
    channelcolors = channelcolors(1:s3,:);
end
    
% clipping
for i=size(clip,1)+1:s3, clip(i,:) = clip(1,:); end
clip = clip(1:s3,:);
m = clip(:,1)'; 
M = clip(:,2)'; 

% colormap if only one channel
im(isinf(im)) = NaN;
nanpix = any(isnan(im),2);
if s3==1
    ncol = size(cmap,1);
    if logscale && m>0 %#ok<BDSCI>
        im = (log(im)-log(m)) * (ncol/(log(M)-log(m)));
    else
        im = (im-m) * (ncol/(M-m));
    end
    im = 1+uint16(im);
    im = min(im,length(cmap));
    im = cmap(im,:);
else
    % clipping each channel
    im = im - repmat(m,s(1)*s(2),1);
    im = im ./ repmat((M-m),s(1)*s(2),1);
    im = min(1,max(0,im));
    % assign appropriate colors
    if all(ismember(channelcolors(:),[0 1])) && 0
    else
        % column operation
        im = im * channelcolors;
    end
end
im(nanpix,:) = repmat(nancolor,full(sum(nanpix)),1);

% final reshape
im = reshape(im,[s(1) s(2) 3]);
end

%---
function im2 = ApplyMask(im,mask)

s = size(im);
s2 = [s(1)*s(2) s(3:end)];
im = reshape(im,s2);
im2 = nan(s2);
im2(mask,:) = im(mask,:);
im2 = reshape(im2,s);

end

%---
function im = HighlightSelection(im,selectionmarks,seldims) 
% Input im should be a true-color gray image

s = size(im);
if ~isscalar(seldims), im = reshape(im,s(1)*s(2),3); end

% highlight sets with colors
colors = fn_colorset;
ncol = size(colors,1);
for i=1:length(selectionmarks)
    sel = selectionmarks(i);
    if ~sel.active, continue, end
    indsel = sel.dataind;
    % each pixel becomes 67% its original gray color + 33% the
    % highlighting color
    for j=1:3
        switch seldims
            case 'x'
                im(indsel,:,j) = (2*im(indsel,:,j) + colors(mod(i-1,ncol)+1,j))/3;
            case 'y'
                im(:,indsel,j) = (2*im(:,indsel,j) + colors(mod(i-1,ncol)+1,j))/3;
            case 'xy'
                im(indsel,j) = (2*im(indsel,j) + colors(mod(i-1,ncol)+1,j))/3;
        end
    end
end

im = reshape(im,s);
end


%-----------------------
% TOOLS: scale bar size
%-----------------------

function x = BarSize(ha)

% minimal size (25 pix) in axes coordinates
xmin = fn_coordinates(ha,'b2a',[25 0],'vector');
xmin = xmin(1);

% desired size in axes coordinates
xlim = get(ha,'xlim');
x = max(xmin,diff(xlim)/5);

% round to a nice value
x10 = 10^floor(log10(x));
x = x / x10;
vals = [1 2 5];
f = find(x>=vals,1,'last');
x = vals(f) * x10;

end

%---------------------
% TOOLS: change shape
%---------------------

%---
function chgrectangle(ha,hl,flagpt,idx,desc)

% if flagpt
%     % move corner
%     pol = [x x2 x2 x; y y y2 y2]; % anti-clockwise from (x,y)
%     dist = sum(fn_add(pol,-p).^2);
%     [dum idx] = min(dist); idx = idx(1);
% else
%     % move edge
%     dist = abs([p(2)-y p(1)-x2 p(2)-y2 p(1)-x]);
%     [dum idx] = min(dist);
% end

% update coordinates
x = desc(1); y = desc(2);
w = desc(3); h = desc(4);
x2 = x+w;    y2 = y+h;
p = get(ha,'currentpoint'); p = p(1,1:2)';
pol = [x x2 x2 x; y y y2 y2]; % anti-clockwise from (x,y)
% (which coordinates to move?)
if flagpt
    codes = fn_switch(idx,1,[1 1],2,[2 1],3,[2 2],4,[1 2]);
else
    codes = fn_switch(idx,1,[0 1],2,[2 0],3,[0 2],4,[1 0]);
end
% (change x coordinate for relevant corners)
switch codes(1)
    case 1
        pol(1,[1 4]) = p(1);
        desc(1) = p(1);
        desc(3) = x2-p(1);
    case 2
        pol(1,[2 3]) = p(1);
        desc(3) = p(1)-x;
end
% (change y coordinate for relevant corners)
switch codes(2)
    case 1
        pol(2,[1 2]) = p(2);
        desc(2) = p(2);
        desc(4) = y2-p(2);
    case 2
        pol(2,[3 4]) = p(2);
        desc(4) = p(2)-y;
end

% update display
set(hl(1),'xdata',pol(1,[1:4 1]),'ydata',pol(2,[1:4 1]))
set(hl(2),'xdata',pol(1,:),'ydata',pol(2,:))
setappdata(hl(2),'description',desc)
drawnow update

end

%---
function chgellipse(ha,hl,idx,desc)

[x u logic] = deal(desc{:});
p = get(ha,'currentpoint'); p = p(1,1:2)';

switch idx
    case 1
        dp = p-(x-u);
        x = x+dp/2;
        u = u-dp/2;
        v = [u(2); -u(1)];
    case 2
        dp = p-(x+u);
        x = x+dp/2;
        u = u+dp/2;
        v = [u(2); -u(1)];
    case 0
        % eccentricity
        v = [u(2); -u(1)];
        normu2 = sum(u.^2);
        xp = p-x;
        uc = sum(xp.*u)/normu2;
        vc = sum(xp.*v)/normu2;
        logic(1) = abs(vc / (sin(acos(uc))));
    case -1
        % internal radius
        v = [u(2); -u(1)];
        e = logic(1);
        xp = p-x;
        logic(2) = sqrt((xp'*u)^2 + (xp'*v/e)^2) / norm(u)^2;
end
phi = linspace(0,2*pi,20);
udata = cos(phi);
vdata = logic(1)*sin(phi);
if length(logic)==2
    udata = [udata NaN logic(2)*udata];
    vdata = [vdata NaN logic(2)*vdata];
end
xdata = x(1) + u(1)*udata + v(1)*vdata;
ydata = x(2) + u(2)*udata + v(2)*vdata;
polymark = fn_add(x,[u -u]);
set(hl(1),'xdata',xdata,'ydata',ydata)
set(hl(2),'xdata',polymark(1,:),'ydata',polymark(2,:))
if ~isreal(logic), error programming, end
setappdata(hl(2),'description',{x u logic})
drawnow update

end

%---
function y=nanmean(x,varargin)

x(isnan(x))=[];
y = mean(x,varargin{:});

end





    


