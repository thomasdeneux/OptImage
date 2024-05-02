classdef activedisplayPlayer < fn4Dhandle
    % function D = activedisplayPlayer(SI,[options])
    
    properties (SetAccess='private')
        SI
        oktime
        dt
        hp
        playbutton
        label
        slider
        timer
        lastrefresh
    end
    properties 
        speed
    end
    properties (Dependent)
        fps
        playing
    end
    
    % Constructor
    methods
        function D = activedisplayPlayer(varargin)
            fn4D_dbstack
            
            % options for initialization
            opt = struct( ...
                'in',                   [] ...
                );
            if nargin==0 || ~isobject(varargin{1})
                si = sliceinfo(1);
                si.slice = struct('data',(1:10)');
                [opt optadd] = fn4D_parseInput(opt,varargin{:});
            else
                si = varargin{1};
                [opt optadd] = fn4D_parseInput(opt,varargin{2:end});
            end
            if si.nd ~= 1
                error 'Slice info must be of dimension 1'
            end
            D.SI = si;
            
            % graphic objects
            if isempty(opt.in)
                hf = figure(99);
                clf(hf)
                opt.in = hf; %uipanel('parent',hf,'pos',[.1 .1 .4 .2]);
            elseif ~ismember(get(opt.in,'type'),{'figure' 'uipanel'})
                error('''in'' option must be a figure or uipanel handle')
            end
            D.hp = opt.in;
            delete(get(opt.in,'children'))
            D.playbutton = uicontrol('parent',opt.in, ...
                'style','togglebutton', ...
                'callback',@(u,e)playpause(D));
            D.buttonicon('play');
            D.label = uicontrol('parent',opt.in,'style','text');
            D.slider = fn_slider('parent',opt.in,'callback',@(u,e)D.sliderupdate());
            fn_pixelsizelistener(opt.in,D,@(u,e)D.positioncontrols())
            propname = fn_switch(get(opt.in,'type'),'figure','Color','uipanel','BackgroundColor');
            connect_listener(opt.in,D,propname,'PostSet',@(u,e)D.positioncontrols())
            D.positioncontrols()
            
            % delete object upon deletion of any button
            addlistener(D.playbutton,'ObjectBeingDestroyed',@(u,e)delete(D));
            addlistener(D.label,'ObjectBeingDestroyed',@(u,e)delete(D));
            addlistener(D.slider,'ObjectBeingDestroyed',@(u,e)delete(D));
            
            % time unit and speed control (sets dt, controls)
            D.checktimeunit()
            
            % timer
            D.timer = timer('timerfcn',@(u,e)D.changeframe(), ...
                'period',.01,'executionmode','fixedRate'); %#ok<CPROP>
            
            % set more properties
            if ~isempty(optadd)
                set(D,optadd{:})
            end
            
            assignin('base','D',D)
        end
    end
    
    % Speed control
    methods
        function dt = get.dt(D)
            dt = D.SI.grid(1);
        end
        function fps = get.fps(D)
            if D.oktime
                fps = D.speed / D.dt;
            else
                fps = D.speed;
            end
        end
        function set.fps(D,fps)
            if D.oktime
                D.speed = fps * D.dt;
            else
                D.speed = fps;
            end
        end
        function sliderupdate(D)
            % speed expressed as relative to realtime if possible
            x = D.slider.value;
            pow = floor(x);
            u = 10*mod(x,1); % between 0 and 9.9999
            if u<8, u = u/2+1; else u = 2.5*u-15; end % now u is between between 1 and 9.9999
            D.speed = u * 10^pow;
            % update label
            updatelabel(D)
        end
        function set.speed(D,speed)
            D.speed = speed;
            % update label
            updatelabel(D)
            % update slider
            updateslider(D)
        end
        function updateslider(D)
            pow = floor(log10(D.speed));
            u = D.speed/10^pow; % between 1 and 9.9999
            if u<5, u = 2*(u-1); else u = (u+15)/2.5; end % now u is between 0 and 9.9999
            D.slider.value = pow + u/10;
        end
        function updatelabel(D)
            if D.oktime
                set(D.label, 'string', [num2str(D.speed) ' x real time'])
            else
                set(D.label, 'string', [num2str(D.speed) ' fps'])
            end
        end
    end
    
    % Layout
    methods
        function positioncontrols(D)
            % layout
            sz = fn_pixelsize(D.hp);
            horizontallayout = sz(1) > 3*sz(2);
            
            % label font size
            maxfontsize = 10;
            psz = fn_objectsize(D.hp,'points');
            labelheightp = fn_switch(horizontallayout,psz(2),psz(2)/2);
            if labelheightp < maxfontsize
                set(D.label,'fontsize',labelheightp)
            else
                psz(2) = maxfontsize;
                set(D.label,'fontsize',maxfontsize,'units','points','pos',[0 0 psz])
            end
            labelsz = fn_pixelsize(D.label);
            labelheight = labelsz(2);

            % control positions
            sz = fn_pixelsize(D.hp);
            if horizontallayout
                % horizontal layout
                offset = sz(2)/15;
                side = sz(2) - 2*offset;
                buttonwidth = sz(2);
                set(D.playbutton,'units','pixel',...
                    'position',[offset offset buttonwidth-2*offset side])
                labelwidth = min(1.6*sz(2), (sz(1)-buttonwidth)/2);
                set(D.label,'units','pixel', ...
                    'position',[buttonwidth+offset sz(2)/2-labelheight/2 labelwidth-2*offset labelheight])
                set(D.slider,'units','pixel', ...
                    'position',[buttonwidth+labelwidth+offset offset sz(1)-buttonwidth-labelwidth-2*offset side])
            else
                % vertical layout
                offset = sz(2)/25;
                h = sz(2)/2;
                set(D.playbutton,'units','pixel', ...
                    'position',[offset h+offset sz(1)-2*offset h-2*offset])
                labelwidth = min(sz(1)/3, sz(2)*1.5);
                set(D.label,'units','pixel', ...
                    'position',[offset offset labelwidth-2*offset h-2*offset])
                set(D.slider,'units','pixel', ...
                    'position',[labelwidth+offset offset sz(1)-labelwidth-2*offset h-2*offset])
            end
            
            % button font size
            psz = fn_objectsize(D.playbutton, 'points');
            set(D.playbutton,'fontsize',psz(2)*.6)
        end
        function checktimeunit(D)
            % check time unit
            switch D.SI.units{1}
                case {'s' 'second'}
                    D.oktime = true;
                    D.dt = D.SI.grid(1);
                case 'ms'
                    D.oktime = true;
                    D.dt = D.SI.grid(1) * 1e-3;
                case {'' 'frame'}
                    D.oktime = false;
                    D.dt = 1;
                otherwise
                    error('unknown time unit ''%s''', D.SI.units{1})
            end
            
            % update controls
            D.buildcontrols()
        end
        function buildcontrols(D)
            % use a funny log scale
            if D.oktime
                % default 1*realtime
                set(D.slider,'minmax',[-1 2],'step',.02)
                D.speed = 1;
            else
                % default 15fps
                set(D.slider,'minmax',[0 3],'step',.02)
                D.speed = 15;
            end
        end
        function buttonicon(D, mode)
            switch mode
                case 'play'
                    set(D.playbutton,'string',char(9658))
                case 'pause'
                    set(D.playbutton,'string',char([10074 10074]))
            end
        end
    end
    
    % Events
    methods
        function updateDown(D,~,evnt)
            fn4D_dbstack
            switch evnt.flag
                case 'units'
                    checktimeunit(D)
            end
        end
    end
    
    % Play
    methods
        function b = get.playing(D)
            b = logical(get(D.playbutton,'value'));
        end
        function set.playing(D,b)
            if b == D.playing(), return, end
            set(D.playbutton,'value',b)
            D.playpause()
        end
        function playpause(D)
            if get(D.playbutton,'value')
                % change button icon to pause
                D.buttonicon('pause')
                % mark last refresh
                nsecperday = 3600*24;
                D.lastrefresh = now * nsecperday;
                % start the timer
                start(D.timer)
            else
                % stop the timer
                stop(D.timer)
                % change button icon back to play
                D.buttonicon('play')
            end
        end
        function changeframe(D)
            % current frame
            tsel = D.SI.selectionmarks;
            if isempty(tsel)
                kframe = D.SI.ij;
            else
                % kframe is the index in the selection!
                kframe = find(tsel.dataind >= D.SI.ij, 1);
                if isempty(kframe), kframe = 0; end
            end

            % by how many frames to step
            nsecperday = 3600*24;
            t = now * nsecperday;
            if t - D.lastrefresh < 1/D.fps
                % frame has not changed, return
                return
            end
            step = floor((t - D.lastrefresh)*D.fps);
            %disp(['step ' num2str(step)])
            D.lastrefresh = t;
            kframe = kframe + step;

            % which frame to show
            tsel = D.SI.selectionmarks;
            if isempty(tsel)
                D.SI.ij2 = fn_mod(kframe + step, D.SI.sizes);
            else
                % kframe is the index in the selection!
                D.SI.ij2 = tsel.dataind(fn_mod(kframe + step, length(tsel.dataind)));
            end
            
            % if the frame goes out of the zoom, change the zoom!!
            if D.SI.ij2 > D.SI.zoom(2)
                D.SI.zoom = D.SI.zoom + diff(D.SI.zoom);
            elseif D.SI.ij2 < D.SI.zoom(1)
                D.SI.zoom = D.SI.zoom - diff(D.SI.zoom);
            end
            
            % process the event queue
            drawnow
        end
    end
end