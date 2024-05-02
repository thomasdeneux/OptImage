classdef activedisplaySlider < fn4Dhandle
% function activedisplaySlider(SI,varargin)
% function activedisplaySlider(G,dim,varargin)

    properties
        hf
        hu
        SI
    end
    
    methods
        function D = activedisplaySlider(varargin)
            % function activedisplaySlider(SI,varargin)
            % function activedisplaySlider(G,dim,varargin)
            %---
            % possible options in varargin: 'in', 'layout'
            fn4D_dbstack
            
            % options for initialization
            opt = struct( ...
                'in',           [], ...
                'layout',       'auto', ...
                'scrollwheel',  [] ...
                );
            if nargin==0 || ~isobject(varargin{1})
                SI = sliceinfo(1);
                SI.slice = struct('data',(1:10)');
                [opt optadd] = fn4D_parseInput(opt,varargin{:});
            elseif isa(varargin{1},'sliceinfo')
                SI = varargin{1};
                [opt optadd] = fn4D_parseInput(opt,varargin{2:end});
            elseif isa(varargin{1},'geometry')
                [G dim] = deal(varargin{1:2});
                SI = projection(G,dim);
                SI.data = shiftdim(zeros(G.sizes(dim),1),-(dim-1));
                [opt optadd] = fn4D_parseInput(opt,varargin{3:end});
            else
                error argument
            end
            D.SI = SI;
            
            % type check
            if SI.nd~=1
                error('activedisplaySlider must rely on a sliceinfo object with nd=1')
            end
            
            % associated graphic object
            if isempty(opt.in)
                opt.in = figure; 
                fn_setfigsize(opt.in,400,30)
            end
            if ~ishandle(opt.in) && mod(opt.in,1)==0 && opt.in>0, figure(opt.in), end
            switch get(opt.in,'type')
                case 'figure'
                    D.hf = opt.in;
                    figure(opt.in)
                    clf(D.hf,'reset')
                    set(D.hf,'menubar','none')
                    hp = uipanel('parent',D.hf);
                case 'uipanel'
                    hp = opt.in;
                    D.hf = get(hp,'parent');
                    while ~strcmp(get(D.hf,'type'),'figure')
                        D.hf = get(D.hf,'parent');
                    end
                otherwise
                    error 'activedisplaySlider object should be defined inside a figure or an uipanel'
            end
            if isempty(get(D.hf,'Tag')), set(D.hf,'Tag','used by fn_4D'), end
            if isempty(opt.scrollwheel), opt.scrollwheel = strcmp(opt.in,'figure'); end
            
            % slider and event (bottom-up)
            D.hu = fn_slider(hp,'mode','point', ...
                'min',1,'max',SI.sizes, ...
                'layout',opt.layout, ...
                'sliderstep',[0 1/(SI.sizes-1)], ...
                'scrollwheel','on', ...
                'callback',@(hu,evnt)event(D));
            if opt.scrollwheel, setscrollwheel(D.hu), end
            set(D.hu,'deletefcn',@(u,e)delete(D))
            
            % communication with parent
            addparent(D,D.SI)
            
            % update display (here, just sets the correct value)
            displayvalue(D);
            
            % set more properties
            if ~isempty(optadd)
                set(D,optadd{:})
            end
        end
        
        function delete(D)
            if ishandle(D.hu), delete(D.hu), end
        end
        
        function event(D)
            fn4D_dbstack
            D.SI.ij2 = get(D.hu,'value');
        end
        
        function updateDown(D,S,evnt)
            fn4D_dbstack
            switch evnt.flag
                case 'ij2'
                    displayvalue(D)
                case 'sizes'
                    set(D.hu,'max',D.SI.sizes,'sliderstep',[0 1/(D.SI.sizes-1)])
            end
        end
        
        function displayvalue(D)
            fn4D_dbstack        
            set(D.hu,'value',D.SI.ij2);
        end
    end
    
end




