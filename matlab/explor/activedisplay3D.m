classdef activedisplay3D < fn4Dhandle

    properties
        userclip = [0 1];
        clipmode = 'slice'; % 'view', 'slice', 'data', 'link1' or 'link2'
        
        scaledisplay = 'tick' % 'tick', 'xbar' or ''
        dopixelratio = true;
        
        cmap = 'user'; % not used yet % 'gray', 'jet', 'mapclip', 'mapcliphigh', 'signcheck' or 'green' - special case 'user'
        
        yfun
        usercallback
    end
    
    properties (SetAccess='private')
        ha
        hf
        
        % some objects are in public access to allow change of some of
        % their properties (color, linewidth, ...)
        cross
        scalebar % scale bar and text
        
        channelcolors
        %docolorimage
        %selectionlabels     
        
        currentdisplay
        autoclipvalue
    end
    
    properties (SetAccess='private') %(Access='private')
        dims = [1 2; 3 2; 1 3];

        img
        %seldisp = {};
        
        handle
        poslistener
        txt
        buttons
        slider
        menu
        
        menuitems
        
        %curselprev
        %cursel
        
        oldaxis
        oldselectiontype = 'normal';
    end
    
    properties (Dependent)
        clip
    end
    
    properties (Dependent, SetAccess='private')
        axis
        xyz
    end
    
    properties (SetAccess='private')
        CL
        SI
        C2D
        listenaxpos
    end
    
    % Constructor and Destructor
    methods
        function D = activedisplay3D(varargin)
            % function activedisplay3D([P,][options])
            fn4D_dbstack
            
            % options for initialization
            opt = struct( ...
                'in',                   [], ...
                'clip',                 [], ...
                'channelcolors',        [] ...
                );
            if nargin==0 || ~isobject(varargin{1})
                close % XX remove that
                G = geometry;
                D.SI = projection(G,[1 2 3]);
                D.SI.data = rand(4,5,6); % XX remove that
                [opt optadd] = fn4D_parseInput(opt,varargin{:});
            else
                D.SI = varargin{1};
                [opt optadd] = fn4D_parseInput(opt,varargin{2:end});
            end
            
            % type check
            if D.SI.nd~=3 || D.SI.ndplus>1 || (D.SI.ndplus ==1 && D.SI.sizesplus(1)>3)
                error('activedisplay3D class can display only three-dimensional data slices, in at most 3 channels')
            end
            
            % figure and axes
            if ~isempty(opt.in)
                figure(opt.in);
                D.hf = opt.in;
            else
                D.hf = gcf;
            end
            clf(D.hf), delete(get(D.hf,'children'))
            % axes will be re-positionned later when axis will be set
            for i=1:3, D.ha(i) = subplot(2,2,i); end
            set(D.ha,'units','pixel')
            fn_pixelsizelistener(D.hf,D,@(u,e)axespositions(D))
            if isempty(get(D.hf,'Tag')), set(D.hf,'Tag','fn4D'), end
            
            % handle to change the size of the 3 axes
            D.handle = uicontrol('parent',D.hf, ...
                'enable','inactive','buttondownfcn',@(u,e)axespositionsmanual(D), ...
                'visible',onoff(~D.dopixelratio));
            
            % image
            for i=1:3
                D.img(i) = image(0,'Parent',D.ha(i),'hittest','off','CDataMapping','scaled');
            end
            set(D.ha,'CLimMode','manual')
            
            % cross
            for i=1:3
                D.cross(i,1) = line('Parent',D.ha(i),'Color','white', ...
                    'ButtonDownFcn',@(hl,evnt)movebar(D,i,1));
                D.cross(i,2) = line('Parent',D.ha(i),'Color','white', ...
                    'ButtonDownFcn',@(hl,evnt)movebar(D,i,2));
            end
            
            % scale bar (and listener for re-positioning upon axes resize)
            D.scalebar(1) = line('Parent',D.ha(1),'Color','white','visible','off', ...
                'linewidth',3);
            D.scalebar(2) = text('Parent',D.ha(1),'Color','white','visible','off', ...
                'horizontalalignment','center','verticalalignment','middle');
            disp 'warning: position listener ''listenaxpos'' disabled for compatibility with new Matlab version'
            %             D.listenaxpos = connect_listener(D.ha(1),D,'Position','PostSet', ...
            %                 @(h,evnt)displayscalebar(D));
            %             D.listenaxpos.Enable = 'off';
            
            % value and buttons; TODO: don't use fn_controlpositions
            D.txt = uicontrol('Parent',D.hf,'style','text','enable','inactive', ...
                'fontsize',8,'horizontalalignment','left');
            fn_controlpositions(D.txt,       D.hf,[1 0],[-234 92 200 12]);
            D.buttons(1) = uicontrol('Parent',D.hf, ...
                'backgroundcolor',[.5 0 0],'foregroundcolor',[.5 0 0]);
            fn_controlpositions(D.buttons(1),D.hf,[1 0],[ -34 92 12 12]);
            D.buttons(2) = uicontrol('Parent',D.hf,'hittest','off', ...
                'backgroundcolor',[.8 .8 0],'foregroundcolor',[.8 .8 0]);
            fn_controlpositions(D.buttons(2),D.hf,[1 0],[ -22 92 12 12]);
            D.buttons(3) = uicontrol('Parent',D.hf, ...
                'backgroundcolor',[.8 .8 0],'foregroundcolor',[.8 .8 0]);
            fn_controlpositions(D.buttons(3),D.hf,[1 0],[  -14 92 4 12]);
            set(D.buttons,'style','frame','enable','off')
            
            % stepping buttons
            D.buttons(4) = uicontrol('parent',D.hf,'CallBack',@(u,evnt)chgpt(D,'left'),'String','-');
            fn_controlpositions(D.buttons(4),D.hf,[1 0],[-100 50 20 20]);
            D.buttons(5) = uicontrol('parent',D.hf,'CallBack',@(u,evnt)chgpt(D,'right'),'String','-');
            fn_controlpositions(D.buttons(5),D.hf,[1 0],[ -60 50 20 20]);
            D.buttons(6) = uicontrol('parent',D.hf,'CallBack',@(u,evnt)chgpt(D,'up'),'String','|');
            fn_controlpositions(D.buttons(6),D.hf,[1 0],[ -80 70 20 20]);
            D.buttons(7) = uicontrol('parent',D.hf,'CallBack',@(u,evnt)chgpt(D,'down'),'String','|');
            fn_controlpositions(D.buttons(7),D.hf,[1 0],[ -80 30 20 20]);
            D.buttons(8) = uicontrol('parent',D.hf,'CallBack',@(u,evnt)chgpt(D,'forward'),'String','\');
            fn_controlpositions(D.buttons(8),D.hf,[1 0],[ -60 30 20 20]);
            D.buttons(9) = uicontrol('parent',D.hf,'CallBack',@(u,evnt)chgpt(D,'backward'),'String','\');
            fn_controlpositions(D.buttons(9),D.hf,[1 0],[-100 70 20 20]);            

            % sliders
            D.slider = cell(3,2); % trick for initializing the array and not have automatic creation of non-necessary fn_slider objects
            for i=1:3
                D.slider{i,1} = fn_slider('parent',D.hf,'mode','point', ...
                    'callback',@(u,evnt)chgzoom(D,D.dims(i,1),u),'visible','off');
                fn_controlpositions(D.slider{i,1},D.ha(i),[0 1 1 0], [0 0 0 10]);
                D.slider{i,2} = fn_slider('parent',D.hf,'mode','point','layout','down', ...
                    'callback',@(u,evnt)chgzoom(D,D.dims(i,2),u),'visible','off');
                fn_controlpositions(D.slider{i,2},D.ha(i),[1 0 0 1], [0 0 10 0]);
            end
            D.slider = reshape([D.slider{:}],3,2);
            
            % optional values
            D.channelcolors = opt.channelcolors;
            %D.docolorimage = (D.seleditmode || D.SI.ndplus);
            
            % update display - no need to call displayscalebar(D) since
            % this is done automatically whenever axis zoom is changed
            displaydata(D)
            displaygrid(D) % automatic displayzoom(D)
            displaylabels(D)
            displaycross(D)
            displayvalue(D)
            %updateselection(D,'all')
            
            % callbacks (bottom-up)
            for i=1:3, set(D.ha(i),'ButtonDownFcn',@(ha,evnt)Mouse(D,i)), end
            %set(D.txt,'buttondownfcn',@(hu,evnt)Mouse(D,true))
            set(D.buttons(1),'buttondownfcn',@(hu,evnt)redbutton(D))
            initlocalmenu(D)
            initlocalmenu(D.SI,D.buttons(3))
            
            % clipping - automatic updates
            if isempty(opt.clip)
                D.clipmode = 'slice';
            else
                D.clip = opt.clip;
                D.clipmode = 'data';
            end
            
            % communication with parent
            addparent(D,D.SI)
            
            % trick to make reset of axes trigger object deletion
            for i=1:3
                line('parent',D.ha(i),'visible','off','deletefcn',@(x,y)delete(D))
            end
            
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
            
            m1 = uimenu(m,'label','color map');
            info.cmap.gray = uimenu(m1,'label','gray', ...
                'callback',@(m,evnt)set(D,'cmap','gray'));
            info.cmap.jet = uimenu(m1,'label','jet', ...
                'callback',@(m,evnt)set(D,'cmap','jet'));
            info.cmap.mapclip = uimenu(m1,'label','mapclip', ...
                'callback',@(m,evnt)set(D,'cmap','mapclip'));
            info.cmap.mapcliphigh = uimenu(m1,'label','mapcliphigh', ...
                'callback',@(m,evnt)set(D,'cmap','mapcliphigh'));
            info.cmap.signcheck = uimenu(m1,'label','signcheck', ...
                'callback',@(m,evnt)set(D,'cmap','signcheck'));
            info.cmap.green = uimenu(m1,'label','green', ...
                'callback',@(m,evnt)set(D,'cmap','green'));
            info.cmap.user = uimenu(m1,'label','user', ...
                'enable','off');
            set(info.cmap.(D.cmap),'checked','on') 
            
            %m1 = uimenu(m,'label','clipping mode','separator','on');
            m1 = m;
            info.clip.view = uimenu(m1,'label','clip mode view','separator','on', ...
                'callback',@(hu,evnt)set(D,'clipmode','view'));
            info.clip.slice = uimenu(m1,'label','clip mode slice', ...
                'callback',@(hu,evnt)set(D,'clipmode','slice'));
            info.clip.data = uimenu(m1,'label','clip mode data', ...
                'callback',@(hu,evnt)set(D,'clipmode','data'));
            info.clip.link1 = uimenu(m1,'label','clip mode link1', ...
                'callback',@(hu,evnt)set(D,'clipmode','link1'));
            info.clip.link2 = uimenu(m1,'label','clip mode link2', ...
                'callback',@(hu,evnt)set(D,'clipmode','link2'));
            set(info.clip.(D.clipmode),'checked','on') 
            info.usrclip = uimenu(m,'label','user clip', ...
                'callback',@(hu,evnt)set(D,'clip',fn_input('clip',D.clip(1,:))));
            
%             m1 = uimenu(m,'label','shape select mode','separator','on');
%             info.shape.poly = uimenu(m1,'label','poly', ...
%                 'callback',@(hu,evnt)set(D,'shapemode','poly'));
%             info.shape.free = uimenu(m1,'label','free', ...
%                 'callback',@(hu,evnt)set(D,'shapemode','free'));
%             info.shape.ellipse = uimenu(m1,'label','ellipse', ...
%                 'callback',@(hu,evnt)set(D,'shapemode','ellipse'));
%             set(info.shape.(D.shapemode),'checked','on')             
%             info.seledit = uimenu(m,'label','advanced selection','separator','on', ...
%                 'callback',@(hu,evnt)set(D,'seleditmode',~get(D,'seleditmode')));
%             if D.seleditmode, set(info.seledit,'checked','on'), end
%             info.selshow = uimenu(m,'label','display selection marks', ...
%                 'callback',@(hu,evnt)set(D,'selshow',~get(D,'selshow')));            
%             if D.selshow, set(info.seledit,'checked','on'), end
%             info.selopen = uimenu(m,'label','open selections', ...
%                 'callback',@(hu,evnt)set(D,'openselections',~get(D,'openselections')));            
%             if D.selshow, set(info.seledit,'checked','on'), end
%             info.sel(1) = uimenu(m,'label','reset selection display', ...
%                 'callback',@(hu,evnt)updateselection(D,'all'));
%             info.sel(2) = uimenu(m,'label','reset selection', ...
%                 'callback',@(hu,evnt)updateselection(D.SI,'reset'));
            
            m1 = uimenu(m,'label','features','separator','on');
            uimenu(m1,'label','show cross','checked',get(D.cross(1),'visible'), ...
                'callback',@showcross);
            function showcross(m2,evnt) %#ok<INUSD>
                onoff = fn_switch(get(m2,'checked'),'switch');
                set(m2,'checked',onoff)
                set(D.cross,'visible',onoff)
            end
            info.dopixelratio = uimenu(m,'label','constrain pixel ratio','checked',onoff(D.dopixelratio), ...
                'callback',@(hu,evnt)set(D,'dopixelratio',~D.dopixelratio));
            
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
            for i=1:3, cla(D.ha(i),'reset'), end
            delete(D.txt(ishandle(D.txt)))
            delete(D.buttons(ishandle(D.buttons)))
            delete(D.menu(ishandle(D.menu)))
            %if ishandle(D.poslistener), delete(D.poslistener), end % i don't understand how it becomes invalid...
        end
    end
    
    % Display
    methods (Access='private')
        function axespositions(D)
            % distances between axes - side with/without ticks
            D1 = 15;
            if ~strcmp(D.scaledisplay,'tick')
                D2 = D1;
            elseif isempty(get(get(D.ha(1),'ylabel'),'string')) ...
                    && isempty(get(get(D.ha(3),'ylabel'),'string'))
                D2 = 35; 
            else
                D2 = 50;
            end
            siz = diff(D.axis,1,2);
            posf = get(D.hf,'position');
            hspace = posf(3) - 2*(D2+D1);
            vspace = posf(4) - 2*(D2+D1);
            
            % now it's complicate, because a ratio is imposed only if units
            % are the same, and we don't have a singleton dimention 
            % -> basically, 5 different cases, do it by hand
            % it by hand
            units = D.SI.units;
            constraint = zeros(1,3);
            if D.dopixelratio
                for i=1:3
                    dd = [i 1+mod(i,3)];
                    constraint(i) = (~any(D.SI.sizes(dd)==1) ...
                        && strcmp(units{dd(1)},units{dd(2)}));
                end
            end
            switch num2str(constraint,'%i')
                case '000'
                    % no constraint at all
                    xspace  = hspace/2;
                    yspace  = vspace/2;
                    zspace(1) = hspace-xspace;
                    zspace(2) = vspace-yspace;
                case '100'
                    % constraint only between x and y, just leave 50% of
                    % the space in horizontal and vertical for z
                    rxy = siz(1)/siz(2);
                    xspace = (hspace/2+vspace/2*rxy)/2;
                    yspace = (hspace/2/rxy+vspace/2)/2;
                    zspace(1) = hspace-xspace;
                    zspace(2) = vspace-yspace;
                    % this 'guess' might compress too much one of the z
                    if zspace(1)<20
                        zspace(1) = 20;
                        xspace = hspace - zspace(1);
                        yspace = xspace/rxy;
                        zspace(2) = vspace - yspace;
                    elseif zspace(2)<20
                        zspace(2) = 20;
                        yspace = vspace - zspace(2);
                        xspace = yspace*rxy;
                        zspace(1) = hspace - xspace;
                    end
                case '001'
                    % constraint only between x and z
                    xspace = hspace * siz(1)/(siz(1)+siz(3));
                    zspace = hspace-xspace;
                    yspace = vspace-zspace;
                    % not enough space for y?
                    if yspace<vspace/2
                        yspace = vspace/2;
                        zspace = vspace-yspace;
                        xspace = zspace * siz(1)/siz(3);
                    end
                case '010'
                    % the symmetric case (constraint only between y and z)
                    yspace = vspace * siz(2)/(siz(2)+siz(3));
                    zspace = vspace-yspace;
                    xspace = hspace-zspace;
                    % not enough space for y?
                    if xspace<hspace/2
                        xspace = hspace/2;
                        zspace = hspace-xspace;
                        yspace = zspace * siz(2)/siz(3);
                    end
                case '111'
                    % full constraint
                    xspace = hspace * siz(1)/(siz(1)+siz(3));
                    zspace(1) = hspace-xspace;
                    yspace = vspace * siz(2)/(siz(2)+siz(3));
                    zspace(2) = vspace-yspace;
                    % need to reduce one of the two horizontal/vertical
                    if zspace(1)>zspace(2)
                        zspace = zspace(2);
                        xspace = zspace * siz(1)/siz(3);
                    else
                        zspace = zspace(1);
                        yspace = zspace * siz(2)/siz(3);
                    end
            end
            if isscalar(zspace), zspace = [zspace zspace]; end
            % only one of hgap and vgap should be non-zero
            hgap = (hspace - xspace - zspace(1))/2;
            vgap = (vspace - yspace - zspace(2))/2;
            
            % finally, set the axes positions!
            set(D.ha(1),'position',[D2+hgap              D2+vgap+zspace(2)+D1+D2 xspace    yspace])
            set(D.ha(2),'position',[D2+hgap+xspace+D1+D2 D2+vgap+zspace(2)+D1+D2 zspace(1) yspace])
            set(D.ha(3),'position',[D2+hgap              D2+vgap                 xspace    zspace(2)])
            
            % and the handle
            side = 12;
            set(D.handle,'position',[D2+hgap+xspace+(D1+D2)/2-side/2 D2+vgap+zspace(2)+(D1+D2)/2-side/2 side side])
        end
        
        function axespositionsmanual(D)
            set([D.ha],'units','pixel') 
            pos = fn_get(D.ha,'position');
            left = pos{3}(1);
            bottom = pos{3}(2);
            d1 = pos{2}(1) - (pos{3}(1)+pos{3}(3));
            d2 = pos{1}(2) - (pos{3}(2)+pos{3}(4));
            hfirst0 = pos{3}(3); hspace = hfirst0 + pos{2}(3);
            vfirst0 = pos{3}(4); vspace = vfirst0 + pos{2}(4);
            side = 12;

            p0 = get(D.hf,'currentpoint'); 
            fn_buttonmotion(@movehandle)
            function movehandle
                dp = get(D.hf,'currentpoint')-p0;
                hfirst = hfirst0+dp(1);
                vfirst = vfirst0+dp(2);
                set(D.ha(1),'position',[left               bottom+vfirst+d2 hfirst        vspace-vfirst])
                set(D.ha(2),'position',[left+hfirst+d1     bottom+vfirst+d2 hspace-hfirst vspace-vfirst])
                set(D.ha(3),'position',[left               bottom           hfirst        vfirst])
                set(D.handle,'position',[left+hfirst+d1/2-side/2 bottom+vfirst+d2/2-side/2 side side])
            end
        end
        
        function displaygrid(D)
            fn4D_dbstack
            % changes axis, image xdata and ydata, and cross extremities
            s = D.SI.sizes;
            grid = D.SI.grid;
            
            % scaling
            range = zeros(3,2);
            for i=1:3, range(i,:) = [1 s(i)]*grid(i,1)+grid(i,2); end
            for i=1:3, set(D.img(i),'xdata',range(D.dims(i,1),:),'ydata',range(D.dims(i,2),:)), end
            
            % axis
            sidefact = grid(:,1)/2; sidefact(s==1)=.5;
            side = sidefact * [-1 1];
            lims = zeros(3,2);
            for i=1:3
                lims(i,:) = grid(i,2)+[1 s(i)]*grid(i,1)+side(i,:);
            end
            D.oldaxis = lims;
            
            %zooming
            displayzoom(D) % set axis -> automatic display update (+ move buttons, position scale bar)
            
            % change cross!
            for i=1:3
                set(D.cross(i,1),'YData',lims(D.dims(i,2),:))
                set(D.cross(i,2),'XData',lims(D.dims(i,1),:))
            end
        end
        
        function displaylabels(D)
            fn4D_dbstack
            if strcmp(D.scaledisplay,'tick')
                labels = D.SI.labels;
                units = D.SI.units;
                for i=1:3
                    if ~isempty(units{i}), labels{i} = [labels{i} ' (' units{i} ')']; end
                end
                for i=1:3
                    set(D.ha(i),'xtickmode','auto','ytickmode','auto')
                    xlabel(D.ha(i),labels{D.dims(i,1)});
                    ylabel(D.ha(i),labels{D.dims(i,2)});
                end
            else
                set(D.ha,'xtick',[],'ytick',[])
                for i=1:3
                    xlabel(D.ha(i),'');
                    ylabel(D.ha(i),'');
                end
            end
        end
        
        function displayscalebar(D)
            fn4D_dbstack
            if ~strcmp(D.scaledisplay,'xbar')
                set(D.scalebar,'visible','off')
                %                 D.listenaxpos.Enable = 'off';
                return
            else
                set(D.scalebar,'visible','on')
                %                 D.listenaxpos.Enable = 'on';
            end
            % find a nice size for bar: in specialized function
            barsize = BarSize(D.ha(1));
            % label: use units if any - must be the same for x and
            % y!
            label = num2str(barsize);
            if ~isempty(D.SI.units{1})
                units = D.SI.units;
                if ~strcmp(units{1},units{2})
                    error('units are not the same for x and y axis')
                end
                label = [label ' ' units{1}];
            end
            % positions
            barorigin = fn_coordinates(D.ha(1),'b2a',[20 10]','position');
            barpos = [barorigin barorigin+[barsize 0]'];
            textpos = mean(barpos,2) + ...
                fn_coordinates(D.ha(1),'b2a',[0 10]','vector');
            % set properties
            set(D.scalebar(1),'xdata',barpos(1,:),'ydata',barpos(2,:))
            set(D.scalebar(2),'position',textpos,'string',label)
        end
        
        function displaydata(D)
            fn4D_dbstack
            slice = D.SI.slice;
            
            % no data -> no display
            if isempty(slice), return, end
            
            if length(slice)>1
                error('activedisplay3D can display only one selection at a time')
            end
            % ignore 'active' flag!
            im = slice.data;
            
            % apply user function if any
            if ~isempty(D.yfun)
                im = D.yfun(im);
            end
            
            % depending whether the SI property is a projection3D object or
            % not, im can be a cell array of 3 2d-arrays, or a 3d array
            if ~iscell(im)
                ij = D.SI.ij;
                im = {squeeze(im(:,:,ij(3),:)) permute(squeeze(im(ij(1),:,:,:)),[2 1 3]) ...
                    squeeze(im(:,ij(2),:,:))};
            end
            D.currentdisplay = im;
            
            % update the default clipping (and change display if 'slice' mode)
            autoclipupdate(D)
            
            % display
            for i=1:3
                set(D.img(i),'CData',permute(im{i},[2 1 3]))
            end
        end
        
        function displaycross(D)
            fn4D_dbstack
            ij2 = D.SI.ij2;
            % scaling and translation
            pt = IJ2AX(D.SI,ij2);
            for i=1:3
                set(D.cross(i,1),'XData',pt([1 1]*D.dims(i,1)))
                set(D.cross(i,2),'YData',pt([1 1]*D.dims(i,2)))
            end
        end
        
        function displayvalue(D)
            fn4D_dbstack
            ij = D.SI.ij;
            im = D.currentdisplay;
            if isempty(im)
                set(D.txt,'String','')
            else
                set(D.txt,'String', ...
                    ['val(' num2str(ij(1)) ',' num2str(ij(2)) ',' num2str(ij(3)) ')=' ...
                    num2str(im{1}(ij(1),ij(2),:),'%.2g ')])
            end
        end
        
        function displayzoom(D)
            fn4D_dbstack
            % must be applied AFTER displaygrid
            oldax = D.oldaxis;
            zoom = IJ2AX(D.SI,D.SI.zoom);
            % min and max to stay within range
            zoom(:,1) = max(oldax(:,1),zoom(:,1));
            zoom(:,2) = min(oldax(:,2),zoom(:,2));
            if any(diff(zoom,1,2)<=0)
                disp('new zoom is outside of range - do zoom reset')
                zoom = oldax;
            end
            D.axis = zoom;
        end
        
        function displayselection(D,flag,ind,value)
            fn4D_dbstack
            if ~D.selshow
                delete(findobj(D.ha,'tag','ActDispIm_Sel'))
                return
            end
            
            % some params
            nsel = length(D.SI.selectionmarks);
            colors = fn_colorset;
            ncol = length(colors);
            SI = D.SI;
            
            % display set...
            if fn_ismemberstr(flag,{'all','reset'})
                % 'findobj' allows a cleanup when some objects were not
                % removed correctly
                delete(findobj(D.ha,'tag','ActDispIm_Sel'))
                D.seldisp = cell(1,nsel);
                isel = 1;
                selectionmarks = SI.selectionmarks;
                for k=1:nsel
                    displayonesel(D,colors,ncol,k,isel);
                    if selectionmarks(k).active, isel = isel+1; end
                end
                return
            end
            
            % or display update
            if ~isempty(D.curselprev)
                set(D.seldisp{D.curselprev}(2),'fontangle','normal')
            end
            switch flag
                case 'new'
                    isel = sum([SI.selectionmarks.active]);
                    displayonesel(D,colors,ncol,ind,isel);
                case {'add','change'}
                    sel = SI.selectionmarks(ind);
                    polygon = PolyUnion(sel.poly,D.openselections);
                    polygon = IJ2AX(D.SI,polygon);
                    set(D.seldisp{ind}(1),'xdata',polygon(1,:),'ydata',polygon(2,:))
                    set(D.seldisp{ind}(2),'position',[nanmean(polygon(1,:)) nanmean(polygon(2,:))])
                    if D.seleditmode, seleditinit(D,ind), end
                case 'remove'
                    delete([D.seldisp{ind}])
                    D.seldisp(ind) = [];
                    nsel = length(D.seldisp);
                    if nsel==0, return, end
                    updateselorderdisplay(D,colors,ncol)
                case 'active'
                    % might be several indices
                    for k=ind
                        if SI.selectionmarks(k).active
                            col = colors(mod(k-1,ncol)+1,:);
                            linestyle = '-';
                            visible = 'on';
                        else
                            col = 'k';
                            linestyle = '--';
                            visible = 'off';
                        end
                        set(D.seldisp{ind}(1),'color',col,'linestyle',linestyle)
                        set(D.seldisp{ind}(2),'visible',visible)
                    end
                    updateselorderdisplay(D,colors,ncol)
                case 'reorder'
                    perm = value;
                    D.seldisp = D.seldisp(perm);
                    updateselorderdisplay(D,colors,ncol)
                case 'indices'
                    % nothing to do
            end
            if ~isempty(D.cursel)
                set(D.seldisp{D.cursel}(2),'fontangle','oblique')
            end
        end
        
        function displayonesel(D,colors,ncol,k,isel)
            sel = D.SI.selectionmarks(k);
            polygon = PolyUnion(sel.poly,D.openselections);
            polygon = IJ2AX(D.SI,polygon);
            
            if sel.active
                col = colors(mod(k-1,ncol)+1,:);
                linestyle = '-';
                visible = 'on';
            else
                col = 'k';
                linestyle = '--';
                visible = 'off';
            end
            
            str = num2str(isel);
            hl = line(polygon(1,:),polygon(2,:),'Parent',D.ha, ...
                'Color',col,'LineStyle',linestyle);
            ht = text(nanmean(polygon(1,:)),nanmean(polygon(2,:)),str, ...
                'Parent',D.ha,'color','w','visible',visible, ...
                'horizontalalignment','center','verticalalignment','middle');
            if k==D.cursel, set(ht,'fontangle','oblique'), end
            set([hl ht],'tag','ActDispIm_Sel','HitTest','off')
            D.seldisp{k} = [hl ht];
            
            if D.seleditmode, seleditinit(D,k), end           
        end
        
        function updateselorderdisplay(D,colors,ncol)
            isel = 1;
            selectionmarks = D.SI.selectionmarks;
            for k=1:length(selectionmarks)
                set(D.seldisp{k}(2),'string',num2str(isel))
                if selectionmarks(k).active
                    col = colors(mod(k-1,ncol)+1,:);
                    set(D.seldisp{k}(1),'color',col)
                    isel = isel+1;
                end
                if D.seleditmode, seleditinit(D,k), end
            end
        end
        
        function seleditinit(D,k)
            hl = D.seldisp{k}(1);
            set(hl,'Marker','.')
            
            % draw a lot of hit-able things!!
            polygon = [get(hl,'xdata'); get(hl,'ydata')];
            col = get(hl,'color');
            np = size(polygon,2);
            nadd = np-1-2*sum(isnan(polygon(1,:)));
            seldisp = D.seldisp{k};
            delete(seldisp(3:end)); %#ok<PROP>
            seldisp = [seldisp(1:2) zeros(1,2*nadd)]; %#ok<PROP>
            kadd = 3;
            
            % lines
            for i=2:np
                if isnan(polygon(1,i)) || isnan(polygon(1,i-1)), continue, end
                seldisp(kadd) = line('parent',D.ha, ...
                    'xdata',polygon(1,[i-1 i]),'ydata',polygon(2,[i-1 i]), ...
                    'Color',col,'LineStyle','--','Tag','ActDispIm_Sel', ...
                    'ButtonDownFcn',@(l,evnt)polyedit(D,hl,k,i-1,'line'), ...
                    'Interruptible','off');
                kadd = kadd+1;
            end
            
            % points
            for i=1:np
                if isnan(polygon(1,i)), continue, end
                if ~D.openselections && (i==1 || isnan(polygon(1,i-1)))
                    ringstart = i;
                    continue
                end
                if ~D.openselections && (i==np || isnan(polygon(1,i+1)))
                    ip = [ringstart i];
                else
                    ip = i;
                end
                seldisp(kadd) = line(polygon(1,i),polygon(2,i),'parent',D.ha, ...
                    'Color',col,'Marker','.','Tag','ActDispIm_Sel', ...
                    'ButtonDownFcn',@(l,evnt)polyedit(D,hl,k,ip,'point'), ...
                    'Interruptible','off');
                kadd = kadd+1;
            end
            
            % finish
            if kadd-1~=length(seldisp), error programming, end
            D.seldisp{k} = seldisp;
        end
    end
    
    % Update routines
    methods (Access='private')
        function updateselection(D,flag,ind,value)
            fn4D_dbstack(['updateselection ' flag])
            if nargin<3, ind=[]; end
            if nargin<4, value=[]; end
            % current selection
            D.curselprev = D.cursel;
            switch flag
                case 'reorder'
                    perm = value;
                    D.cursel = find(perm==ind);
                case {'remove','reset','all'}
                    D.cursel = length(D.SI.selectionmarks);
                    if D.cursel==0, D.cursel=[]; end
                case 'indices'
                    % no change in D.cursel
                otherwise
                    D.cursel = ind(end);
            end
            % compute selectionlabels
            SI = D.SI;
            switch flag
                case {'new','add'}
                    D.selectionlabels(SI.selectionmarks(ind).dataind) = ind;
                otherwise
                    selectionmarks = SI.selectionmarks;
                    if ~isempty(selectionmarks) ...
                            && ~all(D.SI.sizes==selectionmarks(1).datasizes)
                        % size mismatch, don't update D.selectionlabels;
                        % normally, data should be changed soon too (change
                        % is si.sizes) and then sizes will match again
                        D.selectionlabels = [];
                    else
                        D.selectionlabels = zeros(D.SI.sizes);
                        for k=1:length(selectionmarks)
                            D.selectionlabels(selectionmarks(k).dataind) = k;
                        end
                    end
            end
            % update selection display (sel edit mode)
            if D.seleditmode
                displaydata(D)
            end
            if strcmp(flag,'indices'), return, end
            % update selection display
            displayselection(D,flag,ind,value)
        end
    end
    
    % GET/SET clip
    methods
        function clip = get.clip(D)
            fn4D_dbstack
            clip = zeros(3,2);
            for i=1:3, clip(i,:) = get(D.ha(i),'CLim'); end
        end
        
        function set.clip(D,clip)
            fn4D_dbstack
            if size(clip,1)==1, clip = repmat(clip,3,1); end
            if ~isequal(size(clip),[3 2]), error('clip should be a 3x2 array'), end
            if all(clip==D.clip), return, end
            for i=1:3
                clipi = clip(i,:);
                if diff(clipi)<=0 || any(isnan(clipi)) || any(isinf(clipi))
                    if diff(clipi)==0
                        clipi = clipi + [-1 1];
                    else
                        clipi = [0 1];
                    end
                end
                set(D.ha(i),'CLim',clipi);
            end
            if fn_ismemberstr(D.clipmode,{'link1','link2'})
                D.CL.clip = clip(1,:);
            end
        end
        
        function set.userclip(D,clip)
            fn4D_dbstack
            D.userclip = clip;
        end
        
        function set.clipmode(D,clipmode)
            fn4D_dbstack
            if strcmp(clipmode,D.clipmode), return, end
            if ~fn_ismemberstr(clipmode,{'view','slice','data','link1','link2'})
                error('wrong clip mode ''%s''',clipmode)
            end
            oldclipmode = D.clipmode;
            D.clipmode = clipmode;
            % check mark in uicontextmenu
            set(D.menuitems.clip.(oldclipmode),'checked','off')
            set(D.menuitems.clip.(clipmode),'checked','on')
            % specific actions
            if fn_ismemberstr(oldclipmode,{'link1','link2'})
                % cancel previous cliplink and listener
                disconnect(D,D.CL), delete(D.C2D)
            end
            switch clipmode
                case {'view','slice'}
                    % re-compute autoclip and update display
                    autoclipupdate(D)
                case {'link1','link2'}
                    D.CL = cliplink.find(clipmode,D.clip);
                    D.clip = D.CL.clip;
                    D.C2D = connect_listener(D.CL,D,'ChangeClip', ...
                        @(cl,evnt)clipfromlink(D,D.CL));
            end
        end
        
        function clip = get.autoclipvalue(D)
            % auto-compute if necessary
            if isempty(D.currentdisplay)
                clip = [-1 1];
            elseif isempty(D.autoclipvalue)
                im = D.currentdisplay;
                zoom = D.SI.zoom;
                if ~all(isinf(zoom(:)))
                    zoom(:,1) = max(round(zoom(:,1)),1);
                    zoom(:,2) = min(round(zoom(:,2)),D.SI.sizes');
                    im{1} = im{1}(zoom(1,1):zoom(1,2),zoom(2,1):zoom(2,2));
                    im{2} = im{2}(zoom(3,1):zoom(3,2),zoom(2,1):zoom(2,2));
                    im{3} = im{3}(zoom(1,1):zoom(1,2),zoom(3,1):zoom(3,2));
                end
                clip = zeros(3,2);
                for i=1:3, clip(i,:) = [min(im{i}(:)) max(im{i}(:))]; end
                if ~strcmp(D.clipmode,'view')
                    clip = [min(clip(:,1)) max(clip(:,2))];
                end
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
    end
    methods (Access='private')
        function autoclipupdate(D)
            % reset autoclip value (such that the autoclip function will
            % have to re-compute it)
            % then actually calls autoclip function if 'slice' mode
            D.autoclipvalue = [];
            if fn_ismemberstr(D.clipmode,{'view','slice'}), autoclip(D); end
        end
        
        function clipfromlink(D,CL)
            D.clip = CL.clip;
        end
    end
    
    % GET/SET other
    methods
        function set.cmap(D,cm)
            if ischar(cm)
                cmnew=cm; 
                if ~fn_ismemberstr(cmnew,{'gray','jet','mapclip','mapgeog','mapcliphigh','signcheck','green'})
                    error('wrong color map ''%s''',cmnew)
                end
                cm = feval(cmnew,256);
            else
                cmnew='user'; 
            end
            cmold = D.cmap; %#ok<NASGU>
            D.cmap = cmnew;
            
            % % update check marks
            % set(D.menuitems.cmap.(cmold),'checked','off')
            % set(D.menuitems.cmap.(cmnew),'checked','on')            
            
            % set colormap [TODO: no more figure colormap]
            colormap(cm) 
        end
        function set.scaledisplay(D,flag)
            if ~fn_ismemberstr(flag,{'tick','xbar',''})
                error('wrong value for ''scaledisplay'' property')
            end
            D.scaledisplay = flag;
            displaylabels(D)
            displayscalebar(D)
            axespositions(D)
        end
        function set.dopixelratio(D,val)
            if val==D.dopixelratio, return, end
            D.dopixelratio = val;
            % update marks
            set(D.handle,'visible',onoff(~val))
            set(D.menuitems.dopixelratio,'checked',onoff(val))
            %set(D.poslistener,'Enabled',onoff(val))
            
            % update display
            axespositions(D)
        end
        function set.yfun(D,yfun)
            if ~isempty(yfun) && ~isa(yfun,'function_handle')
                error('''yfun'' should be a function handle with one argument')
            end
            D.yfun = yfun;
            displaydata(D)
        end
        function xyz = get.xyz(D)
            fn4D_dbstack
            xyz = IJ2AX(D.SI,D.SI.ij);
        end
        function axis = get.axis(D)
            fn4D_dbstack
            axis = [get(D.ha(1),'xLim'); get(D.ha(1),'yLim'); get(D.ha(2),'xLim')];
        end
        function set.axis(D,axis)
            fn4D_dbstack
            oldax=D.axis;
            if all(axis==oldax), return, end
            % set axis
            for i=1:3
                set(D.ha(i),'xLim',axis(D.dims(i,1),:),'yLim',axis(D.dims(i,2),:));
            end
            % re-position axes
            axespositions(D)
            % re-position scale bar
            displayscalebar(D)
            % change slider parameters
            for i=1:3
                if any(axis(i,:)~=oldax(i,:)), slideraxis(D,i), end
            end
            % re-compute clipping fitting the shown area
            autoclipupdate(D)
        end
        function slideraxis(D,i)
            idx = find(D.dims==i); % indices of 2 sliders
            if all(D.axis(i,:)==D.oldaxis(i,:))
                set(D.slider(idx),'visible','off')
            else
                step = diff(D.axis(i,:));
                m = D.oldaxis(i,1); M = D.oldaxis(i,2)-step;
                set(D.slider(idx),'visible','on','value',D.axis(i,1), ...
                    'min',m,'max',M,'sliderstep',[0 step/(M-m)])
            end
        end 
    end
    
    % Events (bottom-up: mouse)
    methods (Access='private')
        function Mouse(D,i,outsideflag)
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
            % 
            % i indicates in which axes the event occured
            
            % selection type
            oldselectiontype = D.oldselectiontype;
            selectiontype = get(D.hf,'selectiontype');
            hb = D.ha(i); 
            
            % special case - click outside of axis
            ax = axis(hb); %#ok<CPROP>
            point =  get(hb,'CurrentPoint'); point = point(1,[1 2])';
            if (nargin==3 && outsideflag) ...
                    || point(1)<ax(1) || point(1)>ax(2) ...
                    || point(2)<ax(3) || point(2)>ax(4)
                oldselectiontype = selectiontype;
                selectiontype = 'outside';
            end
            
            % store current selection type
            D.oldselectiontype = selectiontype;
            
%             % open or closed selection
%             if D.openselections
%                 TYPE = 'line2D';
%                 CHG = 'change';
%             else
%                 TYPE = 'poly2D';
%                 CHG = 'add';
%             end
            
            % shortcut
            SI = D.SI;
%             nsel = length(SI.selectionmarks);
            
            % GO!
            dd = D.dims(i,:);
            switch selectiontype
                case 'normal'                           % CHANGE VIEW AND/OR MOVE CURSOR
                    rect = fn_mouse(hb,'rect-');
                    if all(rect(3:4))                   % zoom in
                        ax = D.axis;
                        ax(dd,:) = [rect(1)+[0 rect(3)]; rect(2)+[0 rect(4)]];
                        % change also the zoom in the third dimension,
                        % according to units considerations
                        k = setdiff(1:3,dd);
                        sameunits = strcmp(D.SI.units(dd),D.SI.units{k});
                        if any(sameunits)
                            width = rect(3:4);
                            width = max(width(sameunits));
                            ax(k,:) = D.xyz(k) + [-width/2 width/2];
                        end
                        SI.zoom = AX2IJ(SI,ax);
                    else                                % change xy
                        pt = D.xyz;
                        pt(dd) = rect([1 2]);
                        SI.ij2 = AX2IJ(SI,pt);
                    end
                case 'extend'                         	% EDIT SELECTION
%                     ksel = BelongsToSelection(D,point);
%                     if ksel
%                         typenumber(D,ksel);             % reorder selections
%                     else 
%                         polyax = fn_mouse(hb,[D.shapemode '-'])';
%                         poly = AX2IJ(SI,polyax);
%                         if size(poly,2)==1              % add one point to active region
%                             poly = repmat(round(poly(:,1)),1,4) + ...
%                                 [-.5 -.5 .5 .5; -.5 .5 .5 -.5];
%                             updateselection(SI,CHG,D.cursel, ...
%                                 selectionND(TYPE,poly));
%                         else                            % add poly to active region
%                             updateselection(SI,CHG,D.cursel, ...
%                                 selectionND(TYPE,poly(:,[1:end 1])));
%                         end
%                     end
                case 'alt'                              % NEW SELECTION
%                     polyax = fn_mouse(hb,[D.shapemode '-'])';
%                     ksel = BelongsToSelection(D,polyax(:,1));
%                     poly = AX2IJ(SI,polyax);
%                     if size(poly,2)==1 && ksel          % show/hide
%                         updateselection(SI,'active',ksel, ...
%                             ~SI.selectionmarks(ksel).active)
%                     elseif size(poly,2)==1              % new point
%                         poly = repmat(round(poly(:,1)),1,4) + ...
%                             [-.5 -.5 .5 .5; -.5 .5 .5 -.5];
%                         updateselection(SI,'new',[], ...
%                             selectionND(TYPE,poly));
%                     else                                % new poly
%                         updateselection(SI,'new',[], ...
%                             selectionND(TYPE,poly));
%                     end
                case 'open'                             % MISC
                    switch oldselectiontype
                        case 'normal'                   % zoom out
                            if ~isempty(D.usercallback)             % user callback
                                feval(D.usercallback,D)
                            else                                    % zoom out
                                SI.zoom = [-Inf Inf; -Inf Inf; -Inf Inf];
                            end
                        case 'extend'
                            % better not to use it: interferes with poly selection
                        case 'alt'
%                             ksel = BelongsToSelection(D,point);
%                             if ksel                     % reorder (current selection -> last)
%                                 perm = [setdiff(1:nsel,ksel) ksel];
%                                 updateselection(SI,'active',ksel,true);
%                                 updateselection(SI,'reorder',ksel,perm);
%                             end
                    end
                case 'outside'                          % REMOVE SELECTION
                    switch oldselectiontype
                        case 'normal' 
                            rect = fn_mouse(hb,'rect-');
                            if all(rect(3:4))                   % zoom in
                                ax = D.axis;
                                ax(dd,:) = [rect(1)+[0 rect(3)]; rect(2)+[0 rect(4)]];
                                % change also the zoom in the third dimension,
                                % according to units considerations
                                k = setdiff(1:3,dd);
                                sameunits = strcmp(D.SI.units(dd),D.SI.units{k});
                                if any(sameunits)
                                    width = rect(3:4);
                                    width = max(width(sameunits));
                                    ax(k,:) = D.xyz(k) + [-width/2 width/2];
                                end
                                SI.zoom = AX2IJ(SI,ax);
                            else                                % zoom out
                                SI.zoom = [-Inf Inf; -Inf Inf; -Inf Inf];
                            end
                        case {'extend','alt'}
%                             switch oldselectiontype
%                                 case 'alt'              % unselect all regions
%                                     updateselection(SI,'remove',1:nsel);
%                                 case 'extend'           % unselect last region
%                                     updateselection(SI,'remove',D.cursel);
%                             end
                    end
                otherwise
                    error programming
            end
        end
        
        function movebar(D,i,j)
            fn4D_dbstack
            if ~strcmp(get(D.hf,'selectiontype'),'normal')
                % execute callback for axes
                Mouse(D,i)
                return
            end           
            fn_buttonmotion({@movebarsub,D,i,j});
            % any way, execute also a 'click' in the graph
            p = get(gca,'currentpoint');
            pt = D.xyz;
            pt(D.dims(i,:)) = p(1,1:2);
            D.SI.ij2 = AX2IJ(D.SI,pt);
        end
        
        function ksel = BelongsToSelection(D,pointax)
            
            ij = round(AX2IJ(D.SI,pointax));
            if isempty(D.selectionlabels)
                ksel = 0;
            else
                ksel = D.selectionlabels(ij(1),ij(2));
            end
        end
        
        function typenumber(D,ksel)
            
            set(D.hf,'userdata',struct('ksel',ksel,'str',''), ...
                'keypressfcn',@(x,evnt)typenumbertype(D,evnt), ...
                'windowbuttonmotionfcn',@(x,evnt)typenumberexec(D,'move'), ...
                'windowbuttonupfcn',@(x,evnt)typenumberexec(D,'stop'))
            
        end
        
        function typenumbertype(D,evnt)
            info = get(D.hf,'userdata');
            info.str = [info.str evnt.Key];
            set(D.hf,'userdata',info)
        end
        
        function typenumberexec(D,flag)
            info = get(D.hf,'userdata');
            switch flag
                case 'move'
                    point =  get(D.ha,'CurrentPoint'); point = point(1,[1 2])';
                    ksel = BelongsToSelection(D,point);
                    % if there is no change, return
                    if isequal(ksel,info.ksel), return, end
                    set(D.hf,'userdata',struct('ksel',ksel,'str',''));
                case 'stop'
                    set(D.hf,'userdata',[],'keypressfcn','', ...
                        'windowbuttonmotionfcn','','windowbuttonupfcn','')
                    % if no permutation, at least trick to change active sel
                    if info.ksel && isempty(info.str)
                        updateselection(D.SI,'active',info.ksel,true)
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
            selmarks = D.SI.selectionmarks;
            nsel = length(selmarks);
            if k==0 || k>nsel
                k = nsel+1;
            else
                activesels = find(cat(1,selmarks.active));
                k = activesels(k);
            end
            if k<ksel
                perm = [1:k-1 ksel k:ksel-1 ksel+1:nsel];
            elseif k==ksel
                return
            else
                perm = [1:ksel-1 ksel+1:k-1 ksel k:nsel];
            end
            updateselection(D.SI,'reorder',ksel,perm)
        end
        
        function polyedit(D,hl,ksel,kp,flag)
            poly = [get(hl,'xdata'); get(hl,'ydata')];
            np = size(poly,2);
            p = get(D.ha,'currentpoint'); p = p(1,1:2)';
            selectiontype = get(D.hf,'selectiontype');
            switch selectiontype
                case 'normal'               % MOVE POINT
                    % move point
                    if strcmp(flag,'line')
                        % create point on line
                        poly = [poly(:,1:kp) p poly(:,kp+1:end)];
                        set(hl,'xdata',poly(1,:),'ydata',poly(2,:));
                        kp = kp+1;
                    end
                    poly = fn_buttonmotion(@moveedge,D.ha,hl,kp);
                    updateselpoly(D,'change',ksel,poly)
                case 'extend'               % MOVE POLY
                    poly = fn_buttonmotion(@movepoly,D.ha,hl,poly,p);
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
                        updateselection(D.SI,'remove',ksel)
                    end
                case 'alt'
                    % remove selection
                    polyremove(D,ksel)
            end
        end       
        
        function updateselpoly(D,flag,ksel,polyax)
            poly = AX2IJ(D.SI,polyax);
            if D.openselections, TYPE='line2D'; else TYPE='poly2D'; end
            updateselection(D.SI,flag,ksel,selectionND(TYPE,poly));
        end
    end

    % Events (bottom-up: buttons, sliders)
    methods (Access='private')
        function redbutton(D)
            fn4D_dbstack
            switch get(D.hf,'selectiontype')
                case 'normal'       % change clip
                    clip0 = mean(D.clip,1);
                    p0 = get(D.hf,'currentpoint');
                    ht = uicontrol('style','text','position',[2 2 200 17]);
                    moveclip(D,ht,p0,clip0);
                    % change clip
                    fn_buttonmotion({@moveclip,D,ht,p0,clip0})
                    delete(ht)
                case 'extend'       % toggle advanced selection
                    %D.seleditmode = ~D.seleditmode;
                case 'open'         % use default clipping
                    autoclip(D)
            end
        end        

        function chgpt(D,dirflag)
            SI = D.SI;
            switch dirflag
                case 'left'
                    SI.ij2(1) = SI.ij2(1)-1;
                case 'right'
                    SI.ij2(1) = SI.ij2(1)+1;
                case 'up'
                    SI.ij2(2) = SI.ij2(2)-1;
                case 'down'
                    SI.ij2(2) = SI.ij2(2)+1;
                case 'backward'
                    SI.ij2(3) = SI.ij2(3)-1;
                case 'forward'
                    SI.ij2(3) = SI.ij2(3)+1;
                otherwise
                    error programming
            end
        end
        
        function chgzoom(D,i,U)
            fn4D_dbstack
            % i is the dimension affected
            ax = D.axis;
            ax(i,:) = ax(i,:) + (get(U,'value')-ax(i,1));
            if U.sliderscrolling
                D.axis = ax;
            else
                D.SI.zoom = AX2IJ(D.SI,ax);
            end
        end
        
        function duplicate(D,hf)
            if nargin<2, hf=figure; end
            D1 = activedisplay3D(D.SI,'in',hf);
            colormap(D1.ha(1),colormap(D.ha(1)))
            %D1.seleditmode = D.seleditmode;
            D1.clipmode = D.clipmode;
            D1.clip = D.clip;
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
            set(ht,'string',sprintf('min: %.3f,  max: %.3f',clip(1),clip(2))) %#ok<PROP>
            D.clip = clip;
        end
                
        function movebarsub(D,i,j)
            p = get(gca,'currentpoint');
            pt = D.xyz;
            pt(D.dims(i,j)) = p(1,j);
            D.SI.ij2 = AX2IJ(D.SI,pt);

            %             p = get(D.ha(i),'currentpoint'); p = p(1,j);
            %             k = D.dims(i,j);
            %             [ii jj] = find(D.dims==k);
            %             if jj(1)==1, set(D.cross(ii(1),1),'XData',[p p])
            %             else set(D.cross(ii(1),2),'YData',[p p]), end
            %             if jj(2)==1, set(D.cross(ii(2),1),'XData',[p p])
            %             else set(D.cross(ii(2),2),'YData',[p p]), end
            %             if ~isa(D.SI,'projection3D') %?
            %                 displaydata(D)
            %                 %                 try
            %                 %                     data = D.SI.data;
            %                 %                     xyz = D.xyz;
            %                 %                     xyz(k) = p;
            %                 %                     ij = AX2IJ(D.SI,xyz);
            %                 %                     switch k
            %                 %                         case 1
            %                 %                             set(D.img(2),'cdata',squeeze(data(round(ij(1)),:,:,:)));
            %                 %                         case 2
            %                 %                             set(D.img(3),'cdata',permute(squeeze(data(:,round(ij(2)),:,:)),[2 1 3]));
            %                 %                         case 3
            %                 %                             set(D.img(1),'cdata',squeeze(data(:,:,round(ij(3)),:)));
            %                 %                     end
            %                 %                 end
            %             end
        end
    end
    
    % Events (top-down: listeners)
    methods
        function updateDown(D,~,evnt)
            fn4D_dbstack(['S2D ' evnt.flag])
            switch evnt.flag
                case 'sizes'
                    D.currentdisplay = [];
                    displaygrid(D) % automatic displayzoom(D)
                    % try to update D.selectionlabels; will work only if
                    % the datasizes in selection are already set correctly
                    %updateselection(D,'indices')
                case 'sizesplus'
                    if SI.ndplus>1 || SI.sizesplus(1)>3
                        error('activedisplayImage class can display only two-dimensional data slices, in at most 3 channels')
                    end
                case 'slice' 
                    displaydata(D)
                    displayvalue(D)
                case 'grid'
                    displaygrid(D) % automatic displayzoom(D)
                    % conversion IJ2AX has changed -> update selection
                    % display
                    %displayselection(D,'all')
                case 'labels'
                    displaylabels(D)
                case 'units'
                    displayscalebar(D)
                    displaylabels(D)
                case 'ij2'
                    displaycross(D)
                case 'ij'
                    % if D.SI is a projection3D object, any change in
                    % D.SI.ij should induce a change in D.SI.slice;
                    % therefore we don't update the display now since
                    % anyway it will be updated by the change of the slice
                    if ~isa(D.SI,'projection3D')
                        displaydata(D)
                        displayvalue(D)
                    end
                case 'zoom'
                    displayzoom(D)
                case 'selection'
                    %updateselection(D,evnt.selflag,evnt.ind,evnt.value)
                    % note: if selflag is 'indices' it might be that the
                    % datasizes has changed in the new selection before
                    % D.SI.sizes has changed (see above)
            end
        end
    end
    
end


%-----------------------
% TOOLS: scale bar size
%-----------------------

function x = BarSize(ha)

% minimal size (15 pix) in axes coordinates
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

