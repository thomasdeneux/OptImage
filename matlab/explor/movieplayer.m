classdef movieplayer < fn4Dhandle
    % function D = movieplayer(SI,[options])
    
    properties (SetAccess='private')
        SI
        oktime
        dt
        hp
        playbutton
        label
        slider
        stopcommand
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
        function D = movieplayer(varargin)
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
            D.playbutton = uicontrol('parent',opt.in, ...
                'string',char([9658 10074 10074]),'style','togglebutton', ...
                'callback',@(u,e)playpause(D));
            D.label = uicontrol('parent',opt.in,'style','text');
            D.slider = fn_slider('parent',opt.in,'callback',@(u,e)D.sliderupdate());
            fn_pixelsizelistener(opt.in,D,@(u,e)D.positioncontrols())
            propname = fn_switch(get(opt.in,'type'),'figure','Color','uipanel','BackgroundColor');
            connect_listener(opt.in,D,propname,'PostSet',@(u,e)D.positioncontrols())
            D.positioncontrols()
            
            % delete object upon deletion of any button
            addlistener(D.playbutton,'ObjectBeingDestroyed',@(u,e)delete(D))
            addlistener(D.label,'ObjectBeingDestroyed',@(u,e)delete(D))
            addlistener(D.slider,'ObjectBeingDestroyed',@(u,e)delete(D))
            
            % time unit and speed control (sets dt, controls)
            D.checktimeunit()
            
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
                set(D.label, 'string', ['*realtime: ' num2str(D.speed)])
            else
                set(D.label, 'string', ['fps: ' num2str(D.speed)])
            end
        end
    end
    
    % Layout
    methods
        function positioncontrols(D)
            % control positions
            sz = fn_pixelsize(D.hp);
            if sz(1) > 3*sz(2)
                % horizontal layout
                offset = sz(2)/15;
                side = sz(2) - 2*offset;
                buttonwidth = max(sz(2), min(sz(2)*1.5, sz(1)/3));
                set(D.playbutton,'units','pixel',...
                    'position',[offset offset buttonwidth-2*offset side])
                set(D.label,'units','pixel', ...
                    'position',[buttonwidth+offset offset buttonwidth-2*offset side])
                set(D.slider,'units','pixel', ...
                    'position',[2*buttonwidth+offset offset sz(1)-2*buttonwidth-2*offset side])
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
            
            % icon! format 5/3
            sz = fn_objectsize(D.playbutton, 'points');
            set(D.playbutton,'fontsize',min(sz(1)/2,sz(2)))
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
                % start the movie
                if D.stopcommand
                    disp 'strange, a stop command is active before starting'
                    D.stopcommand = false;
                end
                tsel = D.SI.selectionmarks;
                if isempty(tsel)
                    kframe = D.SI.ij;
                else
                    % kframe is the index in the selection!
                    [~, kframe] = min(abs(D.SI.ij - tsel.dataind));
                end
                nsecperday = 3600*24;
                tprev = now * nsecperday;
                while ~D.stopcommand
                    % by how many frames to step
                    t = now * nsecperday;
                    if t - tprev < 1/D.fps
                        % wait some time before showing the next frame
                        pause(1/D.fps - (t - tprev))
                        step = 1;
                        tprev = tprev + 1/D.fps;
                    else
                        % we missed the time to show the next frame, check
                        % which frame it is most appropriate to show now
                        pause(.001)
                        step = floor((t - tprev)*D.fps);
                        tprev = t;
                    end
                    kframe = kframe + step;
                    
                    % which frame to show
                    tsel = D.SI.selectionmarks;
                    if isempty(tsel)
                        D.SI.ij2 = fn_mod(kframe + step, D.SI.sizes);
                    else
                        % kframe is the index in the selection!
                        D.SI.ij2 = tsel.dataind(fn_mod(kframe + step, length(tsel.dataind)));
                    end
                    disp(D.SI.ij)
                end
                D.stopcommand = false;
            else
                % stop the movie
                if D.stopcommand
                    disp 'strange, a stop command is already activated!'
                else
                    D.stopcommand = true;
                end
            end
        end
    end
end