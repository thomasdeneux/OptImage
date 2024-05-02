classdef activedisplayFrames < fn4Dhandle
    % function D = activedisplayFrames(SI,'prop1',value1,...)
    
    properties
        clipmode = 'slice'; % 'slice', 'data', 'link1' or 'link2'
        dataflag = 'data';
        
        scaledisplay = 'tick' % 'tick', 'xbar' or ''

        cmap = 'gray'; 
        
        ncol = 3;
        xbin = 1;
        tbin = 1;
        
        movescale = false;        
        
        yfun
    end
    properties (SetAccess='private')
        ha
        hf
        
        % some objects are in public access to allow change of some of
        % their properties (color, linewidth, ...)
        cross
        scalebar % scale bar and text
        lines
        
        cmapval = gray(256);
        
        currentdisplay
        autoclipvalue
        
        zoom
        datasize
        nrow
        nrowvis
    end    
    properties (Access='private')
        img
        
        posmode     % 'axes' or 'figure'
        txt
        stepper
        buttons
        slider
        menu
        
        menuitems
        
        oldselectiontype = 'normal';
    end    
    properties (Dependent)
        clip
        marker
        markersize
    end   
    properties (Dependent, SetAccess='private')
        scroll
    end   
    properties (SetAccess='private')
        CL
        SI
        C2D
        listenaxpos
    end
    
    % Constructor and Destructor
    methods
        function D = activedisplayFrames(varargin)
            fn4D_dbstack
            
            % options for initialization
            opt = struct( ...
                'in',   [], ...
                'clip', [], ...
                'ncol', 3 ...
                );
            if nargin==0 || ~isobject(varargin{1})
                D.SI = sliceinfo(3); 
                D.SI.slice.data = rand(3,4,25,5);
                opt.in = 1;
                [opt optadd] = fn4D_parseInput(opt,varargin{:});
            else
                D.SI = varargin{1};
                [opt optadd] = fn4D_parseInput(opt,varargin{2:end});
            end
            
            % type check
            if D.SI.nd~=3 || D.SI.ndplus>1
                error('activedisplayFrames displays only 3-dimensions data')
            end
            if D.SI.ndplus>1
                warning 'activedisplayFrames function is known to have bugs when displaying 4-dimensional data'
            end
            
            % optional values
            if isempty(opt.clip)
                D.clipmode = 'slice'; 
            else
                D.clipmode = 'data';
            end
            D.ncol = opt.ncol;
            
            % figure and axes
            if isempty(opt.in)
                opt.in = gcf;
            else
                fn_isfigurehandle(opt.in); % raises new figure if needed
            end
            D.posmode = get(opt.in,'type');
            switch D.posmode
                case 'figure'
                    D.hf = opt.in;
                    clf(D.hf)
                    D.ha = axes('parent',D.hf,'units','pixel');
                case 'axes'
                    D.ha = opt.in;
                    D.hf = get(opt.in,'parent');
                otherwise
                    error('bad handle')
            end
            cla(D.ha,'reset')
            if isempty(get(D.hf,'Tag')), set(D.hf,'Tag','used by fn4D'), end
            
            % image
            D.img = image(0,'Parent',D.ha,'hittest','off','CDataMapping','scaled');
            set(D.ha,'CLimMode','manual','xtick',[],'ytick',[])
            
            % cross
            D.cross(1) = line('Parent',D.ha,'Color','w','hittest','off', ...
                'linestyle','none','marker','+','markersize',7);
            D.cross(2) = line('color','k','linewidth',3,'parent',D.ha);
            
            % scale bar 
            D.scalebar(1) = line('Parent',D.ha,'Color','white','visible','off', ...
                'linewidth',3);
            D.scalebar(2) = text('Parent',D.ha,'Color','white','visible','off', ...
                'horizontalalignment','center','verticalalignment','middle');
            
            % top-left controls all in a single container
            % (containing panel)
            hp = uipanel('parent',D.hf,'bordertype','none');
            fn_controlpositions(hp,D.ha,[0 1],[1 1 365 12])
            % (value)
            D.txt = uicontrol('Parent',hp,'pos',[1 1 100 12], ...
                'style','text','enable','inactive', ...
                'fontsize',8,'horizontalalignment','left');
            % (red and yellow menu buttons)
            D.buttons(1) = uicontrol('Parent',hp,'pos',[101 1 12 12], ...
                'backgroundcolor',[.5 0 0],'foregroundcolor',[.5 0 0]);
            D.buttons(2) = uicontrol('Parent',hp,'pos',[113 1 12 12], ...
                'hittest','off', ...
                'backgroundcolor',[.8 .8 0],'foregroundcolor',[.8 .8 0]);
            D.buttons(3) = uicontrol('Parent',D.hf,'pos',[121 1 2 12], ...
                'backgroundcolor',[.8 .8 0],'foregroundcolor',[.8 .8 0]);
            set(D.buttons,'style','frame','enable','off')
            % (steppers)
            hp1 = uipanel('parent',hp,'units','pixel','pos',[125 1 80 12],'bordertype','none');
            D.stepper = fn_control(struct('xbin',{D.xbin 'stepper 1 1 100 1'}), ...
                @(s)set(D,'xbin',s.xbin),hp1,'nobutton');
            hp1 = uipanel('parent',hp,'units','pixel','pos',[205 1 80 12],'bordertype','none');
            D.stepper(2) = fn_control(struct('tbin',{D.tbin 'stepper 1 1 100 1'}), ...
                @(s)set(D,'tbin',s.tbin),hp1,'nobutton');
            hp1 = uipanel('parent',hp,'units','pixel','pos',[285 1 80 12],'bordertype','none');
            D.stepper(3) = fn_control(struct('ncol',{D.ncol 'stepper 1 1 100 1'}), ...
                @(s)setncol(D,s.ncol),hp1,'nobutton');
            
            % value and buttons; TODO: don't use fn_coordinates inside fn_controlpositions
            
            % sliders
            D.slider = fn_slider('parent',D.hf,'mode','point','layout','down', ...
                'callback',@(u,evnt)chgscroll(D),'visible','off');
            fn_controlpositions(D.slider,D.ha,[1 0 0 1], [-1 -1 12 1]);
            
            % update display 
            displaysize(D)
            displaylabels(D)
            displayvalue(D)
            displayscalebar(D)
            
            % callbacks (bottom-up)
            set(D.ha,'ButtonDownFcn',@(ha,evnt)Mouse(D))
            set(D.txt,'buttondownfcn',@(hu,evnt)Mouse(D,true))
            set(D.buttons(1),'buttondownfcn',@(hu,evnt)redbutton(D))
            initlocalmenu(D)
            initlocalmenu(D.SI,D.buttons(3))
            
            % trick to make reset of axes trigger object deletion
            line('parent',D.ha,'visible','off','deletefcn',@(x,y)delete(D))
            
            % communication with parent
            addparent(D,D.SI)
            
            % listeners (axes position)
            if strcmp(D.posmode,'figure')
                % axes position will best fit the changes of the figure
                % size
                D.listenaxpos = fn_pixelsizelistener(D.hf,D,@(h,evnt)displaysize(D));
            else
                % figure size is ignored
                D.listenaxpos = fn_pixelsizelistener(D.ha,D,@(h,evnt)displaysize(D));
            end
            
            % set more properties
            if ~isempty(opt.clip), D.clip = opt.clip; end
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
            
            m1 = uimenu(m,'label','color map');
            info.cmap.gray = uimenu(m1,'label','gray', ...
                'callback',@(m,evnt)set(D,'cmap','gray'));
            info.cmap.jet = uimenu(m1,'label','jet', ...
                'callback',@(m,evnt)set(D,'cmap','jet'));
            info.cmap.mapgeog = uimenu(m1,'label','mapgeog', ...
                'callback',@(m,evnt)set(D,'cmap','mapgeog'));
            info.cmap.mapclip = uimenu(m1,'label','mapclip', ...
                'callback',@(m,evnt)set(D,'cmap','mapclip'));
            info.cmap.mapcliphigh = uimenu(m1,'label','mapcliphigh', ...
                'callback',@(m,evnt)set(D,'cmap','mapcliphigh'));
            info.cmap.mapcliplow = uimenu(m1,'label','mapcliplow', ...
                'callback',@(m,evnt)set(D,'cmap','mapcliplow'));
            info.cmap.vdaq = uimenu(m1,'label','vdaq', ...
                'callback',@(m,evnt)set(D,'cmap','vdaq'));
            info.cmap.signcheck = uimenu(m1,'label','signcheck', ...
                'callback',@(m,evnt)set(D,'cmap','signcheck'));
            info.cmap.green = uimenu(m1,'label','green', ...
                'callback',@(m,evnt)set(D,'cmap','green'));
            info.cmap.user = uimenu(m1,'label','user', ...
                'enable','off');
            set(info.cmap.(D.cmap),'checked','on') 
            
            m1 = uimenu(m,'label','clipping mode','separator','on');
            info.clip.slice = uimenu(m1,'label','clip mode slice', ...
                'callback',@(hu,evnt)set(D,'clipmode','slice'));
            info.clip.data = uimenu(m1,'label','clip mode data', ...
                'callback',@(hu,evnt)set(D,'clipmode','data'));
            info.clip.link1 = uimenu(m1,'label','clip mode link1', ...
                'callback',@(hu,evnt)set(D,'clipmode','link1'));
            info.clip.link2 = uimenu(m1,'label','clip mode link2', ...
                'callback',@(hu,evnt)set(D,'clipmode','link2'));
            set(info.clip.(D.clipmode),'checked','on') 
            info.usrclip = uimenu(m,'label','user clip...', ...
                'callback',@(hu,evnt)set(D,'clip',fn_input('clip',D.clip)));
            
            info.distline = uimenu(m,'label','distance tool','separator','on', ...
                'callback','fn_imdistline');

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
            delete(D.txt(ishandle(D.txt)))
            delete(D.buttons(ishandle(D.buttons)))
            delete(D.stepper(ishandle(D.stepper)))
            delete(D.menu(ishandle(D.menu)))
            delete(D.slider(ishandle(D.slider)))
            delete(D.listenaxpos)
        end
    end
    
    % Display
    methods (Access='private')
        function displaysize(D)            
            fn4D_dbstack
            
            % set axes position if position mode is 'figure'
            if strcmp(D.posmode,'figure')
                [w h] = fn_pixelsize(D.hf);
                set(D.ha,'pos',[12 12 w-24 h-24]);
            end
            
            % size of an individual frame
            D.datasize(1:2) = floor(D.SI.sizes(1:2)/D.xbin);
            z = .5 + (D.SI.zoom(1:2,:)-.5)/D.xbin;
            z(:,1) = max(1,round(z(:,1)));
            z(:,2) = min(D.datasize(1:2)',round(z(:,2))); %#ok<UDIM>
            D.zoom = z;
            
            % number of rows -> visibility, slider
            D.datasize(3) = floor(D.SI.sizes(3)/D.tbin);
            if D.SI.ndplus>=1
                D.datasize(4) = D.SI.sizesplus(1); 
                D.nrow = D.datasize(3);
                D.ncol = D.datasize(4);
            else
                D.datasize(4) = 1; 
                D.nrow = ceil(D.datasize(3)/D.ncol);
            end
            zratio = (D.zoom(2,2)-D.zoom(2,1)+1)/(D.zoom(1,2)-D.zoom(1,1)+1);
            axessize = fn_pixelsize(D.ha);
            axesratio = axessize(2)/axessize(1);
            nrowvi = D.ncol * axesratio / zratio;
            if ceil(nrowvi)>=D.nrow
                D.nrowvis = D.nrow;
                set(D.slider,'value',1,'visible','off','scrollwheel','off')
            else
                D.nrowvis = max(1,round(nrowvi));
                set(D.slider,'visible','on','scrollwheel','default', ...
                    'min',1,'max',D.nrow+1-D.nrowvis,'sliderstep',[0 D.nrowvis/(D.nrow-D.nrowvis)])
            end
            
            % image display properties
            frsiz = D.zoom(:,2)-D.zoom(:,1)+1;
            set(D.img,'xdata',[1+1/(2*frsiz(1)) D.ncol+1-1/(2*frsiz(1))],'ydata',[1+1/(2*frsiz(2)) D.nrow+1-1/(2*frsiz(2))])
            if ~strcmp(D.SI.units{1},D.SI.units{2})
                set(D.ha,'dataAspectRatioMode','auto')
            else
                set(D.ha,'dataAspectRatioMode','manual', ...
                    'dataAspectRatio',[frsiz(2)*D.SI.grid(2,1) frsiz(1)*D.SI.grid(1,1) 1])
            end
            set(D.ha,'xlim',[1 D.ncol+1],'ylim',[1 D.nrowvis+1])
            
            % lines to separate frames
            delete(D.lines(ishandle(D.lines)))
            D.lines = zeros(1,D.ncol-1+D.nrow-1);
            for i=1:D.ncol-1, D.lines(i) = line([i+1 i+1],[1 D.nrow+1],'color','k'); end
            for i=1:D.nrow-1, D.lines(D.ncol-1+i) = line([1 D.ncol+1],[i+1 i+1],'color','k'); end
            
            % update data
            displaydata(D)
            
            % update cross - automatic scrolling if falls outside of view
            displaycross(D)
       end
        function displaydata(D)
            fn4D_dbstack
            if isempty(D.img), return, end % this can happen at init
            slice = D.SI.slice;
            
            % no data -> no display
            if isempty(slice), return, end
            
            if length(slice)>1
                disp('activedisplayFrames can display only one selection at a time')
                slice = slice(1);
            end
            % ignore 'active' flag!
            x = slice.(D.dataflag);
            
            % apply user function is any
            if ~isempty(D.yfun)
                x = D.yfun(x);
            end
            
            % binning
            x = fn_bin(x,[D.xbin D.xbin D.tbin]);
            
            % store current display
            D.currentdisplay = x;
            
            % update the default clipping (and change display if 'slice' mode)
            currentclip = D.clip;
            autoclipupdate(D)
            
            % display frames (not necessary if clip has changed and
            % automatic updates have taken place)
            if all(D.clip==currentclip), displaydata2(D), end
        end
        function displaydata2(D)
            % call this function if data has not changed but display needs
            % to be updated (usually because of change in clip)
            if isempty(D.img), return, end % this can happen at init
            x = fn_float(D.currentdisplay);
            
            % cut the image according to zoom
            ax = [floor(D.zoom(:,1)) ceil(D.zoom(:,2))];
            x = x(ax(1,1):ax(1,2),ax(2,1):ax(2,2),:,:);
            
            % normalize according to clip
            x = (x-D.clip(1))/diff(D.clip);
            x = max(0,min(1-1e-6,x));
            
            % transform into color indices
            ncolors = size(D.cmapval,1);
            x = 1+uint16(floor(x*ncolors));
            
            % special color index (1) for empty space
            s = size(x);
            if prod(s(3:end))<D.ncol*D.nrow, x(1,1,D.ncol*D.nrow) = 0; end
            x = x+1;
            
            % special reshaping and permute
            if D.SI.ndplus>=1
                x = reshape(x,s(1),s(2),D.nrow,D.ncol);
                x = permute(x,[2 3 1 4]);
            else
                x = reshape(x,s(1),s(2),D.ncol,D.nrow);
                x = permute(x,[2 4 1 3]);
            end
            
            % make color image (already in Matlab convention for images)
            bordercol = [0 0 0];
            cm = single([bordercol; D.cmapval]);
            im = cm(x(:),:);
            im = reshape(im,s(2)*D.nrow,s(1)*D.ncol,3);
            
            % display
            set(D.img,'CData',im)
        end       
        function displaylabels(D) %#ok<MANU>
            fn4D_dbstack
%             if strcmp(D.scaledisplay,'tick')
%                 set(D.ha,'xtickmode','auto','ytickmode','auto')
%                 labels = D.SI.labels;
%                 units = D.SI.units;
%                 for i=1:2
%                     if ~isempty(units{i}), labels{i} = [labels{i} ' (' units{i} ')']; end
%                 end
%                 xlabel(D.ha,labels{1});
%                 ylabel(D.ha,labels{2});
%             else
%                 set(D.ha,'xtick',[],'ytick',[])
%                 xlabel(D.ha,'');
%                 ylabel(D.ha,'');
%             end
        end        
        function displayscalebar(D) %#ok<MANU>
            fn4D_dbstack
%             if ~strcmp(D.scaledisplay,'xbar')
%                 set(D.scalebar,'visible','off')
%                 set(D.listenaxpos,'Enable','off')
%                 return
%             else
%                 set(D.scalebar,'visible','on')
%                 set(D.listenaxpos,'Enable','on')
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
            ij2 = D.SI.ij2;
            % point location in a frame
            [pt pixpos] = ij2pt(D,ij2); %#ok<ASGLU>
            % set cross, but be aware of point outside of display!
            if any(pixpos<0 | pixpos>1)
                set(D.cross(1),'Visible','off')
            else
                frames = 0:prod(D.datasize(3:4))-1;
                set(D.cross(1),'Visible','on', ...
                    'XData',1+mod(frames,D.ncol)+pixpos(1), ...
                    'YData',1+floor(frames/D.ncol)+pixpos(2))
            end
            % automatic scroll if current frame is not seen
            k = ceil(ij2(3)/D.tbin);
            krow = ceil(k/D.ncol);
            krow = max(1,min(D.nrow,krow));
            nshifts = floor((krow-D.scroll)/D.nrowvis);
            D.scroll = D.scroll + nshifts*D.nrowvis;
            % indicate the current frame
            kcol = fn_mod(k,D.ncol);
            set(D.cross(2),'xdata',[kcol kcol kcol+1 kcol+1 kcol],'ydata',[krow krow+1 krow+1 krow krow])
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
        
    % GET/SET clip
    methods
        function clip = get.clip(D)
            fn4D_dbstack
            clip = get(D.ha,'CLim');
        end       
        function set.clip(D,clip)
            fn4D_dbstack
            if all(clip==get(D.ha,'CLim')), return, end
            if ~isequal(size(clip),[1 2]), error('clip should be a 2-element row vector'), end
            if ~(diff(clip)>0) % true if NaN!!
                if diff(clip)==0
                    clip = clip + [-1 1];
                else
                    clip = [0 1];
                end
            end
            set(D.ha,'CLim',clip);
            if fn_ismemberstr(D.clipmode,{'link1','link2'})
                D.CL.clip = clip;
            end
           % update display
            displaydata2(D)
        end   
        function set.clipmode(D,clipmode)
            fn4D_dbstack
            if strcmp(clipmode,D.clipmode), return, end
            if ~fn_ismemberstr(clipmode,{'slice','data','link1','link2'})
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
                case 'slice'
                    autoclip(D)
                case {'link1','link2'}
                    D.CL = cliplink.find(clipmode,D.clip);
                    D.clip = D.CL.clip;
                    D.C2D = connect_listener(D.CL,D,'ChangeClip', ...
                        @(cl,evnt)clipfromlink(D,D.CL));
            end
        end       
        function clipfromlink(D,CL)
            D.clip = CL.clip;
        end       
        function clip = get.autoclipvalue(D)
            % auto-compute if necessary
            if isempty(D.currentdisplay)
                clip = [-1 1];
            elseif isempty(D.autoclipvalue)
                x = D.currentdisplay;
                % select only visible frames
                rowmin = floor(D.scroll);
                rowmax = ceil(D.scroll+D.nrowvis-1);
                frames = 1 + (rowmin-1)*D.ncol : ...
                    min(D.datasize(3),rowmax*D.ncol);
                x = x(:,:,frames);
                % cut the frames according to zoom
                ax = [floor(D.zoom(:,1)) ceil(D.zoom(:,2))];
                x = x(ax(1,1):ax(1,2),ax(2,1):ax(2,2),:);
                % compute and assign clip
                clip = [min(x(:)) max(x(:))];
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
    
    % GET/SET other
    methods
        function set.dataflag(D,str)
            D.dataflag = str;
            displaydata(D)
        end
        function set.cmap(D,cm)
            cmold = D.cmap;
            if ischar(cm)
                cmnew=cm; 
                if strcmp(cmnew,cmold), return, end
                if ~fn_ismemberstr(cmnew,{'gray','jet','mapgeog','mapclip','mapcliphigh','mapcliplow','vdaq','signcheck','green'})
                    error('wrong color map ''%s''',cmnew)
                end
                if strcmp(cmnew,'vdaq')
                    cm = vdaqcolors;
                else
                    cm = feval(cmnew,256);
                end
            else
                cmnew='user'; 
            end
            D.cmap = cmnew;
            D.cmapval = cm;
            
            % update check marks
            set(D.menuitems.cmap.(cmold),'checked','off')
            set(D.menuitems.cmap.(cmnew),'checked','on')            
            
            % update display
            displaydata2(D)
        end
        function m=get.marker(D)
            m = get(D.cross(1),'marker');
        end
        function set.marker(D,m)
            set(D.cross(1),'marker',m)
        end
        function x=get.markersize(D)
            x = get(D.cross(1),'markersize');
        end
        function set.markersize(D,x)
            set(D.cross(1),'markersize',x)
        end
        function set.scaledisplay(D,flag)
            if ~fn_ismemberstr(flag,{'tick','xbar',''})
                error('wrong value for ''scaledisplay'' property')
            end
            D.scaledisplay = flag;
            displaylabels(D)
            displayscalebar(D)
        end
        function set.movescale(D,b)
            D.movescale = b;
            displayscalebar(D)
        end
        function set.yfun(D,yfun)
            if ~isempty(yfun) && ~isa(yfun,'function_handle')
                error('''yfun'' should be a function handle with one argument')
            end
            D.yfun = yfun;
            displaydata(D)
        end
        function scroll = get.scroll(D)
            fn4D_dbstack
            % scroll value indicates which row is displayed on the top (it
            % does not have to be an integer value, e.g. if half a row is
            % displayed)
            ylim = get(D.ha,'yLim');
            scroll = ylim(1);
        end
        function set.scroll(D,scroll)
            fn4D_dbstack
            if scroll==D.scroll, return, end
            % set scroll
            ylim = get(D.ha,'yLim');
            ylim = ylim + (scroll-ylim(1));
            set(D.ha,'yLim',ylim);
            % re-position scale bar
            displayscalebar(D)
            % change slider parameters
            set(D.slider,'value',scroll)
            % re-compute clipping fitting the shown area
            autoclipupdate(D)
        end
        function setncol(D,n)
            n = round(n);
            if D.ncol==n, return, end
            D.ncol = n;
            if isempty(D.stepper), error programming, end % should not happen any more
            D.stepper(3).ncol = n;
            displaysize(D)
        end
        function set.xbin(D,n)
            n = max(1,n);
            if D.xbin==n, return, end
            D.xbin = n;
            D.stepper(1).xbin = n; %#ok<*MCSUP>
            displaysize(D)
        end
        function set.tbin(D,n)
            if D.tbin==n, return, end
            D.tbin = n;
            D.stepper(2).tbin = n;
            displaysize(D)
        end
    end
    
    % Events (bottom-up: mouse)
    methods (Access='private')
        function Mouse(D,outsideflag)
            fn4D_dbstack
            % different mouse actions are:
            % - point with left button            -> change cursor
            % - area with left button             -> zoom to region
            % - double-click with left button     -> zoom reset
            %   (or click with left button outside of axis)
            % - click in region with middle       -> reorder selections
            %   button, hold and type a number
            % - point/area with middle button     -> add point/area to current selection
            % - click with middle button outside  -> cancel current selection
            % - point with right button in region -> hide/show selection
            % - point/area with right button      -> add new selection
            % - click with right button outside   -> cancel all selections
            
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
                        rect = [rect(1:2) rect(1:2)+rect(3:4)];
                        rect = pt2ij(D,rect);
                        if rect(3,1)==rect(3,2);
                            % valid zooming only if within a single frame
                            D.SI.zoom(1:2,:) = rect(1:2,:);
                        end
                    else                                % change xy
                        point = pt2ij(D,rect(1:2));
                        D.SI.ij2 = point;
                    end
                case 'extend'                         	% EDIT SELECTION
                case 'alt'                              % NEW SELECTION
                case 'open'                             % MISC
                    switch oldseltype
                        case 'normal'                   % zoom out
                            D.SI.zoom(1:2,:) = [-Inf Inf; -Inf Inf];
                        case 'extend'
                            % better not to use it: interferes with poly selection
                        case 'alt'
                    end
                case 'outside'                          % REMOVE SELECTION
                    switch oldseltype
                        case 'normal' 
                            D.SI.zoom(1:2,:) = [-Inf Inf; -Inf Inf];
                        case {'extend','alt'}
                    end
                otherwise
                    error programming
            end
        end        
    end

    % Events (bottom-up: buttons, sliders)
    methods (Access='private')
        function redbutton(D)
            fn4D_dbstack
            switch get(D.hf,'selectiontype')
                case 'normal'       % change clip
                    clip0 = D.clip;
                    p0 = get(D.hf,'currentpoint');
                    ht = uicontrol('style','text','position',[2 2 200 17]);
                    moveclip(D,ht,p0,clip0);
                    % change clip
                    fn_buttonmotion({@moveclip,D,ht,p0,clip0})
                    delete(ht)
                case 'extend'       % toggle advanced selection
                case 'open'         % use default clipping
                    autoclip(D)
            end
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
                    displaysize(D)
                case 'sizesplus'
                    if D.SI.ndplus>0
                        error('activedisplayFrames class can display only 3D data')
                    end
                case 'slice'
                    displaydata(D)
                    displayvalue(D)
                case 'grid'
                    displayscalebar(D)
                case 'labels'
                    displaylabels(D)
                case 'units'
                    displayscalebar(D)
                case 'ij2'
                    displaycross(D)
                case 'ij'
                    displayvalue(D)
                case 'zoom'
                    displaysize(D)
                case 'selection'
            end
        end
    end
    
    % Coordinates conversions
    methods
        function ij = pt2ij(D,xy)
            if size(xy,1)==1, xy=xy'; end
            if size(xy,1)~=2, error('set of axes coordinates should have two rows'), end
            % position of frame in image, and of point in frame including
            % border (indices start at 0) 
            frpos = floor(xy)-1;
            xy = mod(xy,1);
            frsiz = D.zoom(:,2)-D.zoom(:,1)+1;
            pixpos = fn_mult(frsiz,xy);
            % position of point in data
            ij = [.5+(fn_add(pixpos,D.zoom(:,1)-1))*D.xbin; ...
                .5+(frpos(1,:)+frpos(2,:)*D.ncol+.5)*D.tbin];
        end
        function [xy xyfr] = ij2pt(D,ij)
            % position of frame in image, and of point in frame including
            % border (indices start at 0) 
            k = ceil(ij(3)/D.tbin);
            frpos = [mod((k-1),D.ncol); floor((k-1)/D.ncol)];
            pixpos = (ij(1:2)-.5)/D.xbin-(D.zoom(:,1)-1);
            % position of point in axes
            frsiz = D.zoom(:,2)-D.zoom(:,1)+1;
            xyfr = pixpos./frsiz;
            xy = (frpos+1) + xyfr;
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


    


