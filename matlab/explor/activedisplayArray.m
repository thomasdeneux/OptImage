classdef activedisplayArray < fn4Dhandle
    % function D = activedisplayArray(SI,'prop1',value1,...)
    
    properties
        userclip = [0 1];
        clip
        clipmode = 'data'; % 'data', 'link1' or 'link2'
        dataflag = 'data';
        
        xbin = 1;
        tbin = 1;
        twin;   % time window length
        
        autopos = ''; % '', 'sel', 'cat', 'selcat' or 'catsel'
        linecol = 'sel'; % 'sel', 'cat', '', or custom color (e.g. [1 0 0])
        
        %         scaledisplay = 'tick' % 'tick', 'xbar' or ''
        %         movescale = false;
       
    end
    properties (SetAccess='private')
        ha
        hf
        
        % some objects are in public access to allow change of some of
        % their properties (color, linewidth, ...)
        cross
        %         scalebar % scale bar and text
        
        tokidx      % time indices corresponding to the twin time window
        datasize
    end    
    properties (Access='private')
        databin
        
        buttonpanels
        buttons     % structure
        menu
        menuitems
        hplot   % nxbin * nybin * nsel cell array
        
        oldselectiontype = 'normal';

        ready = false;
    end    
    properties (SetAccess='private')
        CL
        SI
        C2D
        listenaxpos
    end
    
    % Constructor and Destructor
    methods
        function D = activedisplayArray(varargin)
            fn4D_dbstack
            
            % options for initialization
            opt = struct( ...
                'in',   [], ...
                'clip', [], ...
                'xbin', [], ...
                'tbin', [] ...
                );
            if nargin==0 || ~isobject(varargin{1})
                D.SI = sliceinfo(3); 
                D.SI.slice(1) = struct('active',true,'data',rand(3,4,25,3));
                D.SI.slice(2) = struct('active',true,'data',rand(3,4,25,3));
                opt.in = 1;
                [opt optadd] = fn4D_parseInput(opt,varargin{:});
            else
                D.SI = varargin{1};
                [opt optadd] = fn4D_parseInput(opt,varargin{2:end});
            end
            
            % type check
            if D.SI.nd~=3 || D.SI.ndplus>1
                error('activedisplayFrames displays only 3-dimensions data, and maximum 1 additional dimension')
            end
            
            % clipping
            if ~isempty(opt.clip), D.clip = opt.clip; end 
            D.clipmode = 'data';
            
            % figure and axes
            if isempty(opt.in)
                opt.in = gca;
            else
                fn_isfigurehandle(opt.in); % raises new figure if needed
            end
            switch get(opt.in,'type')
                case 'figure'
                    D.hf = opt.in;
                    clf(D.hf)
                    D.ha = axes('parent',D.hf);
                case 'axes'
                    D.ha = opt.in;
                    D.hf = get(opt.in,'parent');
                otherwise
                    error('bad handle')
            end
            cla(D.ha,'reset')
            if isempty(get(D.hf,'Tag')), set(D.hf,'Tag','used by fn4D'), end
            
            % cross
            D.cross = line('Parent',D.ha,'Color','k','hittest','off');
            
            %             % scale bar
            %             D.scalebar(1) = line('Parent',D.ha,'Color','white','visible','off', ...
            %                 'linewidth',3);
            %             D.scalebar(2) = text('Parent',D.ha,'Color','white','visible','off', ...
            %                 'horizontalalignment','center','verticalalignment','middle');
            
            % trick to make reset of axes trigger object deletion
            line('parent',D.ha,'visible','off','deletefcn',@(x,y)delete(D))
            
            % top-left controls: menus and xbin
            % (containing panel)
            hp = uipanel('parent',D.hf,'bordertype','none');
            D.buttonpanels(1) = hp;
            fn_controlpositions(hp,D.ha,[0 1],[1 -11 185 12])
            % (red and yellow menu buttons)
            D.buttons.redmenu = uicontrol('Parent',hp,'pos',[1 1 12 12], ...
                'style','frame','enable','off', ...
                'backgroundcolor',[.5 0 0],'foregroundcolor',[.5 0 0]);
            uicontrol('Parent',hp,'pos',[13 1 12 12], ...
                'hittest','off','style','frame','enable','off', ...
                'backgroundcolor',[.8 .8 0],'foregroundcolor',[.8 .8 0]);
            D.buttons.yellowmenu = uicontrol('Parent',D.hf,'pos',[21 1 2 12], ...
                'style','frame','enable','off', ...
                'backgroundcolor',[.8 .8 0],'foregroundcolor',[.8 .8 0]);
            % (steppers)
            hp1 = uipanel('parent',hp,'units','pixel','pos',[25 1 80 12],'bordertype','none');
            D.buttons.xbin = fn_control(struct('xbin',{D.xbin 'stepper 1 1 100 1'}), ...
                @(s)set(D,'xbin',s.xbin),hp1,'nobutton');
            hp1 = uipanel('parent',hp,'units','pixel','pos',[105 1 80 12],'bordertype','none');
            D.buttons.tbin = fn_control(struct('tbin',{D.xbin 'stepper 1 1 100 1'}), ...
                @(s)set(D,'tbin',s.tbin),hp1,'nobutton');
            
            % top-right controls: time courses position and scale
            % (containing panel)
            hp = uipanel('parent',D.hf,'bordertype','none');
            D.buttonpanels(2) = hp;
            h = 18; H = 7*h; 
            fn_controlpositions(hp,D.ha,[1 1],[-h+1 -H h H])
            uicontrol('parent',hp,'pos',[1 H-1*h+1 h h], ...
                'string','A','callback',@(u,e)chgtime(D,'toggle'));
            uicontrol('parent',hp,'pos',[1 H-2*h+1 h h], ...
                'string','+','callback',@(u,e)chgtime(D,'+'));
            uicontrol('parent',hp,'pos',[1 H-3*h+1 h h], ...
                'string','-','callback',@(u,e)chgtime(D,'-'));
            uicontrol('parent',hp,'pos',[1 H-4*h+1 h h], ...
                'string','e','callback',@(u,e)chgtime(D,'expand'));
            uicontrol('parent',hp,'pos',[1 H-5*h+1 h h], ...
                'string','r','callback',@(u,e)chgtime(D,'reduce'));
            uicontrol('parent',hp,'pos',[1 H-6*h+1 h h], ...
                'string','<','callback',@(u,e)chgtime(D,'<'));
            uicontrol('parent',hp,'pos',[1 H-7*h+1 h h], ...
                'string','>','callback',@(u,e)chgtime(D,'>'));
            
            % xbin, tbin
            if ~isempty(opt.xbin)
                D.xbin = opt.xbin; 
            else
                % start with some binning if there are too many traces
                D.xbin = max(1,round(mean(D.SI.sizes(1:2))/10));
            end
            if ~isempty(opt.tbin), D.tbin = opt.tbin; end
            D.ready = true;
            
            % update display 
            displayzoom(D)
            displaydata(D)
            displaycross(D)
            displaylabels(D)
            %             displayscalebar(D)
            
            % callbacks (bottom-up)
            set(D.ha,'ButtonDownFcn',@(ha,evnt)Mouse(D))
            fn_scrollwheelregister(D.ha,@(n,m)chgtime(D,'scroll',n,m))
            set(D.buttons.redmenu,'buttondownfcn',@(hu,evnt)redbutton(D))
            initlocalmenu(D)
            initlocalmenu(D.SI,D.buttons.yellowmenu)
            
            % communication with parent
            addparent(D,D.SI)
            
            % set more properties
            if ~isempty(opt.clip), D.clip = opt.clip; end
            if ~isempty(optadd)
                set(D,optadd{:})
            end
        end              
        function initlocalmenu(D)
            fn4D_dbstack
            hb = D.buttons.redmenu;
            delete(get(hb,'uicontextmenu'))
            m = uicontextmenu('parent',D.hf);
            D.menu = m;
            set(hb,'uicontextmenu',m)
            
            % clipping
            m1 = uimenu(m,'label','clipping mode','separator','on');
            info.clip.data = uimenu(m1,'label','clip mode data', ...
                'callback',@(hu,evnt)set(D,'clipmode','data'));
            info.clip.link1 = uimenu(m1,'label','clip mode link1', ...
                'callback',@(hu,evnt)set(D,'clipmode','link1'));
            info.clip.link2 = uimenu(m1,'label','clip mode link2', ...
                'callback',@(hu,evnt)set(D,'clipmode','link2'));
            set(info.clip.(D.clipmode),'checked','on') 
            info.usrclip = uimenu(m,'label','user clip', ...
                'callback',@(hu,evnt)set(D,'clip',get(D,'userclip')));
            
            % decoration
            % (coloring scheme)
            m1 = uimenu(m,'label','color scheme','separator','on');
            info.linecol.sel  = uimenu(m1,'label','selection', ...
                'callback',@(hu,evnt)set(D,'linecol','sel'));
            info.linecol.cat  = uimenu(m1,'label','category', ...
                'callback',@(hu,evnt)set(D,'linecol','cat'));
            info.linecol.none = uimenu(m1,'label','none', ...
                'callback',@(hu,evnt)set(D,'linecol','none'));
            info.linecol.custom = uimenu(m1,'label','custom', ...
                'callback',@(hu,evnt)set(D,'linecol',fn_input('line  color',[0 0 1],'color')));
            curflag = fn_switch(isnumeric(D.linecol),'custom',D.linecol);
            set(info.linecol.(curflag),'checked','on')

            % activedisplayArray object
            info.menu(1) = uimenu(m,'label','duplicate in new figure','separator','on', ...
                'callback',@(hu,evnt)duplicate(D));
            info.menu(2) = uimenu(m,'label','display object ''D'' in base workspace', ...
                'callback',@(hu,evnt)assignin('base','D',D));
            
            D.menuitems = info;
        end       
        function delete(D)
            fn4D_dbstack
            % invoked when the invisible line is deleted: D.ha necessarily
            % still exists, but D.txt, D.buttons, D.menu might have
            % been deleted already
            cla(D.ha,'reset')
            delete(D.buttonpanels(ishandle(D.buttonpanels)))
            delete(D.menu(ishandle(D.menu)))
            delete(D.listenaxpos)
        end
    end
    
    % Display
    methods (Access='private')
        function displayratio(D)
            % axis image?
            s = D.SI.sizes;
            if any(s==1) || ~strcmp(D.SI.units{1},D.SI.units{2})
                set(D.ha,'dataAspectRatioMode','auto')
            else
                set(D.ha,'dataAspectRatioMode','manual', ...
                    'dataAspectRatio',[1 1 1])
            end            
        end
        function displayzoom(D)
            fn4D_dbstack
            xytmin = IJ2AX(D.SI,ones(3,1)-.5);
            xytmax = IJ2AX(D.SI,D.SI.sizes'+.5);
            zoomreset = [xytmin(1:2) xytmax(1:2)];
            zoom = IJ2AX(D.SI,D.SI.zoom); zoom = zoom(1:2,:);
            % min and max to stay within range
            zoom(:,1) = max(zoomreset(:,1),zoom(:,1));
            zoom(:,2) = min(zoomreset(:,2),zoom(:,2));
            if any(diff(zoom,1,2)<=0)
                disp('new zoom is outside of range - do zoom reset')
                zoom = zoomreset;
            end
            axis(D.ha,[zoom(1,1) zoom(1,2) zoom(2,1) zoom(2,2)])
            set(D.ha,'ydir','reverse')
        end   
        function displaydata(D)
            if ~D.ready, return, end % happens at init
            
            slice = D.SI.slice;
            slice = slice([slice.active]);
            
            % no data -> no display
            if isempty(slice) || ~isfield(slice,D.dataflag)
                delete(findobj(D.ha,'Tag','fn4D_line'))
                D.hplot = {}; 
                return
            end

            nsel = numel(slice);
            if nsel==0, error programming, end
            
            % bin the data
            D.databin = {slice.data};
            if D.xbin>1
                for i=1:nsel
                    if D.xbin>1
                        D.databin{i} = fn_bin(D.databin{i},[D.xbin D.xbin 1],'smart');
                    end
                    if D.tbin>1
                        D.databin{i} = fn_bin(D.databin{i},[1 1 D.tbin]);
                    end
                end
            end
            
            % display image
            displaydata2(D)
        end
        function displaydata2(D)
            if ~D.ready, return, end % happens at init
            
            % call this function if data and xbin have not changed but
            % display needs to be updated (e.g. because of change in clip
            % or twin)
            
            delete(findobj(D.ha,'Tag','fn4D_line'))
                
            % this happens at init
            if isempty(D.clip)
                chgtime(D,'set') % this will call displaydata2(D) twice
                return
            elseif isempty(D.twin)
                return
            end
            
            % prepare positioning
            % (spatial scale)
            xytmin = IJ2AX(D.SI,ones(3,1)-.5);
            xytmax = IJ2AX(D.SI,D.SI.sizes'+.5);
            nbin = floor(D.SI.sizes(1:2)'/D.xbin);
            xstep = (xytmax(1:2)-xytmin(1:2))./nbin;
            % (time courses stepping and range)
            nsel = length(D.databin);
            ncat = max(fn_map(@(x)size(x,4),D.databin));
            ystep = diff(D.clip)*.5; % default distance between time courses
            seldec = ystep * fn_switch(D.autopos,{'' 'cat'},0,{'sel' 'catsel'},1,'selcat',ncat+.5);
            catdec = ystep * fn_switch(D.autopos,{'' 'sel'},0,{'cat' 'selcat'},1,'catsel',nsel+.5);
            maxdec = (nsel-1)*seldec + (ncat-1)*catdec;
            ylim = [D.clip(1)-maxdec D.clip(2)];
            % (conversion between time courses amplitude and spatial scale)
            tca2space = -xstep(2)/diff(ylim);
            % (conversion between time frames and spatial scale)
            tok = D.tokidx;
            if D.tbin, tok = fn_bin(tok,D.tbin,'or'); end
            nt = sum(tok);
            tt = .05 + .9*(0:nt-1)/(nt-1);
            tt = tt*xstep(1);
            
            % display
            for ksel=1:nsel
                datak = D.databin{ksel}(:,:,tok,:);
                ncatk = size(datak,4);
                for i=1:nbin(1)
                    xoffset(1) = xytmin(1)+(i-1)*xstep(1);
                    for j=1:nbin(2)
                        xoffset(2) = xytmin(2)+(j)*xstep(2);
                        % draw line
                        datai = squeeze(datak(i,j,:,:));
                        datai = datai - ylim(1) - (ksel-1)*seldec;
                        if catdec, datai = fn_subtract(datai,(0:ncatk-1)*catdec); end %#ok<BDLGI>
                        datai = xoffset(2) + datai*tca2space;
                        hl = line(xoffset(1)+tt,datai,'parent',D.ha, ...
                            'hittest','off','tag','fn4D_line');
                        
                        % color
                        if isnumeric(D.linecol)
                            set(hl,'color',D.linecol)
                        elseif strcmp(D.linecol,'sel')
                            set(hl,'color',fn_colorset(ksel))
                        elseif strcmp(D.linecol,'cat')
                            for j=1:length(hl)
                                set(hl(j),'color',fn_colorset(j))
                            end
                        end
                    end
                end
            end
        end     
        function displaylabels(D) 
            fn4D_dbstack
            labels = D.SI.labels;
            units = D.SI.units;
            for i=1:2
                if ~isempty(units{i}), labels{i} = [labels{i} ' (' units{i} ')']; end
            end
            xlabel(D.ha,labels{1});
            ylabel(D.ha,labels{2});
        end        
        function displayscalebar(D) %#ok<MANU>
            fn4D_dbstack
            %             if ~strcmp(D.scaledisplay,'xbar')
            %                 set(D.scalebar,'visible','off')
            %                 set(D.listenaxpos,'Enabled','off')
            %                 return
            %             else
            %                 set(D.scalebar,'visible','on')
            %                 set(D.listenaxpos,'Enabled','on')
            %             end
            %             % find a nice size for bar: in specialized function
            %             barsize = BarSize(D.ha);
            %             % label: use units if any - must be the same for x and
            %             % y!
            %             label = num2str(barsize);
            %             if ~isempty(D.SI.units{1})
            %                 units = D.SI.units;
            %                 %                 if ~strcmp(units{1},units{2})
            %                 %                     error('units are not the same for x and y zoom')
            %                 %                 end
            %                 label = [label ' ' units{1}];
            %             end
            %             % positions
            %             barorigin = fn_coordinates(D.ha,'b2a',[20 10]','position');
            %             barpos = [barorigin barorigin+[barsize 0]'];
            %             textpos = mean(barpos,2) + ...
            %                 fn_coordinates(D.ha,'b2a',[0 10]','vector');
            %             % set properties
            %             set(D.scalebar(1),'xdata',barpos(1,:),'ydata',barpos(2,:))
            %             set(D.scalebar(2),'position',textpos,'string',label)
            %             if D.movescale
            %                 set(D.scalebar,'hittest','on','buttondownfcn', ...
            %                     @(hobj,evnt)fn_moveobject(D.scalebar,'latch'))
            %             else
            %                 set(D.scalebar,'hittest','off')
            %             end
        end        
        function displaycross(D)
            fn4D_dbstack
            nbin = floor(D.SI.sizes(1:2)'/D.xbin);
            actualbinning = D.SI.sizes(1:2)'./nbin;
            square = [0 0 1 1 0; 0 1 1 0 0];
            square = fn_mult(square,D.SI.grid(1:2,1).*actualbinning);
            ij2 = D.SI.ij2(1:2);
            ij2 = floor((ij2-.5)./actualbinning).*actualbinning + .5;
            bottomleft = D.SI.grid(1:2,2) + D.SI.grid(1:2,1).*ij2; 
            square = fn_add(bottomleft,square);
            set(D.cross,'xdata',square(1,:),'ydata',square(2,:))
        end       
        function displayvalue(D)
            fn4D_dbstack
            ij = D.SI.ij;
            x = D.currentdisplay;
            if isempty(x)
                set(D.txt,'String','')
            else
                k = ceil(ij./[D.xbin; D.xbin; D.tbin]);
                set(D.txt,'String', ...
                    ['val(' num2str(k(1)) ',' num2str(k(2))  ',' num2str(k(3)) ')=' ...
                    num2str(x(k(1),k(2),k(3)))])
            end
        end        
    end
        
    % GET/SET
    methods
        function set.clip(D,clip)
            fn4D_dbstack
            D.clip = clip;
            % propagate change in the case of clip link
            if fn_ismemberstr(D.clipmode,{'link1','link2'})
                D.CL.clip = clip;
            end
           % update display
            displaydata2(D)
        end   
        function set.userclip(D,clip)
            fn4D_dbstack
            D.userclip = clip;
        end      
        function set.clipmode(D,clipmode)
            fn4D_dbstack
            if strcmp(clipmode,D.clipmode), return, end
            if ~fn_ismemberstr(clipmode,{'data','link1','link2'})
                error('wrong clip mode ''%s''',clipmode)
            end
            oldclipmode = D.clipmode;
            D.clipmode = clipmode;
            % check mark in uicontextmenu
            if isempty(D.menuitems), return, end % happens at init
            set(D.menuitems.clip.(oldclipmode),'checked','off')
            set(D.menuitems.clip.(clipmode),'checked','on')
            % specific actions
            if fn_ismemberstr(oldclipmode,{'link1','link2'})
                % cancel previous cliplink and listener
                disconnect(D,D.CL), delete(D.C2D)
            end
            switch clipmode
                case {'link1','link2'}
                    D.CL = cliplink.find(clipmode,D.clip); %#ok<*MCSUP>
                    D.clip = D.CL.clip;
                    D.C2D = connect_listener(D.CL,D,'ChangeClip', ...
                        @(cl,evnt)clipfromlink(D,D.CL));
            end
        end       
        function clipfromlink(D,CL)
            fn4D_dbstack
            D.clip = CL.clip;
        end       
        function set.twin(D,twin)
            fn4D_dbstack
            nt = D.SI.sizes(3);
            D.twin = min(twin,nt*D.SI.grid(3,1));
            imin = max(1,D.SI.zoom(3,1));
            imax = min(nt,D.SI.zoom(3,2));
            iwin = min(D.twin / D.SI.grid(3,1),imax-imin);
            imin = fn_coerce(floor(D.SI.ij2(3)-iwin/2),imin,imax-ceil(iwin));
            imax = imin + ceil(iwin);
            D.tokidx = false(1,nt);
            D.tokidx(imin:imax) = true;
            % display update
            displaydata2(D)
        end
        function set.dataflag(D,str)
            D.dataflag = str;
            displaydata(D)
        end
        %         function set.scaledisplay(D,flag)
        %             if ~fn_ismemberstr(flag,{'tick','xbar',''})
        %                 error('wrong value for ''scaledisplay'' property')
        %             end
        %             D.scaledisplay = flag;
        %             displaylabels(D)
        %             displayscalebar(D)
        %         end
        %         function set.movescale(D,b)
        %             D.movescale = b;
        %             displayscalebar(D)
        %         end
        function set.xbin(D,n)
            n = max(1,n);
            if D.xbin==n, return, end
            D.xbin = n;
            % update control
            D.buttons.xbin.xbin = n;
            % update display
            displaydata(D)
            displaycross(D)
        end
        function set.tbin(D,n)
            n = max(1,n);
            if D.tbin==n, return, end
            D.tbin = n;
            % update control
            D.buttons.tbin.tbin = n;
            % update display
            displaydata(D)
        end
        function set.linecol(D,val)
            oldval = D.linecol;
            if isequal(oldval,val), return, end
            if ~(isnumeric(val) && isequal(size(val),[1 3])) && ~fn_ismemberstr(val,{'sel','cat','none'})
                error('wrong property value')
            end
            % update check marks
            oldflag = fn_switch(isnumeric(oldval),'custom',oldval);
            curflag = fn_switch(isnumeric(val),'custom',val);
            set(D.menuitems.linecol.(oldflag),'checked','off')
            set(D.menuitems.linecol.(curflag),'checked','on')
            % set property
            D.linecol = val;
            % update display
            displaydata2(D)
        end
    end
    
    % Events (bottom-up: time buttons, scroll wheel, mouse)
    methods (Access='private')
        function chgtime(D,flag,varargin)
            % change time position or scaling, or data scaling
            % note that display will be updated automatically by calls to
            % set.clip or set.twin
            if strcmp(flag,'scroll')
                [scrollcount modifiers] = deal(varargin{:});
                if fn_ismemberstr('control',modifiers)
                    flag = 'scrollzoom';
                else
                    flag = 'scrollshift';
                end
            else
                scrollcount = 0;
            end
            switch flag
                case 'set'
                    fulldata = [D.databin{:}];
                    D.clip = [min(fulldata(:)) max(fulldata(:))];
                    D.twin = Inf; % automatic clipping and display update
                case 'toggle'
                    nsel = length(D.SI.slice);
                    ncat = max(fn_map(@(x)size(x,4),{D.SI.slice.data}));
                    D.autopos = fn_switch(D.autopos, ...
                        '',         fn_switch(nsel>1,'sel',ncat>1,'cat',''), ...
                        'sel',      fn_switch(ncat==1,'',nsel==1,'cat','selcat'), ...
                        'selcat',   fn_switch(ncat==1,'','cat'), ...
                        'cat',      fn_switch(nsel==1 || ncat==1,'','catsel'), ...
                        'catsel',   '');
                    displaydata2(D) % display update
                case '+'
                    D.clip = mean(D.clip) + [-.5 .5]*diff(D.clip)/sqrt(2);
                case '-'
                    D.clip = mean(D.clip) + [-.5 .5]*diff(D.clip)*sqrt(2);
                case {'expand' 'reduce' 'scrollzoom'}
                    step = 1.1;
                    nstep = fn_switch(flag,'expand',5,'reduce',-5,'scrollzoom',scrollcount);
                    D.twin = D.twin*(step^nstep); % automatic display update
                case {'<' '>' 'scrollshift'}
                    step = D.twin/D.SI.grid(3,1)/20;
                    nstep = fn_switch(flag,'<',-5,'>',5,'scrollshift',-scrollcount);
                    D.SI.ij2(3) = fn_coerce(D.SI.ij2(3) + nstep*step, 1, D.SI.sizes(3));
                otherwise
                    error 'unknown flag'
            end
        end
        function Mouse(D,outsideflag)
            fn4D_dbstack
            % different mouse actions are:
            % - point with left button            -> change spatial point
            % - area with left button             -> spatial zoom
            % - double-click with left button     -> zoom reset
            %   (or click with left button outside of axis)
            
            % selection type
            oldseltype = D.oldselectiontype;
            selectiontype = get(D.hf,'selectiontype');
            
            % special case - click outside of axis
            ax = axis(D.ha); 
            point =  get(D.ha,'CurrentPoint'); point = point(1,[1 2])';
            if (nargin==2 && outsideflag) ...
                    || point(1)<ax(1) || point(1)>ax(2) ...
                    || point(2)<ax(3) || point(2)>ax(4)
                oldseltype = selectiontype;
                selectiontype = 'outside';
            end
            
            % store current selection type
            D.oldselectiontype = selectiontype;
                        
            % GO!
            switch selectiontype
                case 'normal'                           % CHANGE VIEW AND/OR MOVE CURSOR
                    rect = fn_mouse(D.ha,'rect-');
                    rect = rect(:);
                    if all(rect(3:4))                   % zoom in
                        rect = [rect(1:2) rect(1:2)+rect(3:4); 0 0];
                        rect = AX2IJ(D.SI,rect);
                        D.SI.zoom(1:2,:) = rect(1:2,:);
                    else                                % change xy
                        point = AX2IJ(D.SI,[rect(1:2,1); 0]);
                        D.SI.ij2(1:2,:) = point(1:2);
                    end
                case 'extend'                         	% 
                case 'alt'                              % 
                case {'open' 'outside'}                 % SPATIAL ZOOM RESET
                    D.SI.zoom(1:2,:) = [-Inf Inf; -Inf Inf];
                otherwise
                    error programming
            end
        end        
    end

    % Events (bottom-up: buttons, sliders)
    methods (Access='private')
        function redbutton(D)
            fn4D_dbstack
            disp 'no action set yet for red button'
        end        
        function chgscroll(D)
            fn4D_dbstack
            scroll = get(D.slider,'value');
            if D.slider.sliderscrolling
                ylim = get(D.ha,'yLim');
                ylim = ylim + (scroll-ylim(1));
                set(D.ha,'yLim',ylim);
            else
                D.scroll = scroll;
            end
        end       
        function duplicate(D,hobj)
            if nargin<2, hobj=figure; end
            D2 = activedisplayFrames(D.SI,'in',hobj);
            colormap(D2.ha,colormap(D.ha))
            D2.clipmode = D.clipmode;
            D2.clip = D.clip;
        end
    end
    
    % (public for access by fn_buttonmotion)
    methods
        function moveclip(D,ht,p0,clip0)
            % new clip
            p = get(D.hf,'currentpoint');
            r = clip0(2)-clip0(1);
            FACT = 1/100;
            clip = clip0 + (p-p0)*(r*FACT);
            if diff(clip)<=0, clip = mean(clip)+[-1 1]; end
            
            % display
            set(ht,'string',sprintf('min: %.3f,  max: %.3f',clip(1),clip(2))) 
            D.clip = clip;
        end
    end
    
    % Events (top-down: listeners)
    methods
        function updateDown(D,~,evnt)
            fn4D_dbstack(['S2D ' evnt.flag])
            switch evnt.flag
                case 'sizes'
                    displayzoom(D)
                    displaydata(D)
                case 'sizesplus'
                case 'slice'
                    displaydata(D)
                case 'grid'
                    displayzoom(D)
                case 'labels'
                    displaylabels(D)
                case 'units'
                    displaylabels(D)
                case 'ij2'
                    displaycross(D)
                    D.twin = D.twin; % will update tokidx and call displaydata2(D)
                case 'ij'
                case 'zoom'
                    displayzoom(D)
                    D.twin = D.twin; % will update tokidx and call displaydata2(D)
                case 'selection'
            end
        end
    end
    
    % Misc
    methods
        function access(D) %#ok<MANU>
            keyboard
        end
    end
end

%-----------------------
% TOOLS: scale bar size
%-----------------------

function x = BarSize(ha) %#ok<DEFNU>

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


    


