classdef activedisplayList < fn4Dhandle
    % function D = activedisplayList(SI,[options])

    properties
        selmultin = false;
        scrollwheel = 'on'; % 'on', 'off' or 'default'
    end
    
    properties (SetAccess='private')
        hu
        hf
        hp
        hlabel
        currentsel
        menu
        menuitems
        tempsel
        SI
    end
    
    % Constructor and Destructor
    methods
        function D = activedisplayList(varargin)
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
            D.SI = si;
            
            % add a fake data so that SI matches the size of its parent
            % geometry object
            if isa(si,'projection') && isempty(si.data) && si.sizes~=si.G.sizes(si.proj)
                s = ones(1,max(2,si.proj));
                s(si.proj) = si.G.sizes(si.proj);
                si.data = zeros(s);
            end
            
            % type check
            if si.nd~=1
                error('activedisplayList must rely on a sliceinfo object with nd=1')
            end
            
            % figure and axes
            if isempty(opt.in), opt.in = gcf; end
            if ~ishandle(opt.in) && mod(opt.in,1)==0 && opt.in>0, figure(opt.in), end
            switch get(opt.in,'type')
                case {'figure' 'uipanel'}
                    D.hp = opt.in;
                    delete(get(opt.in,'children'))
                    if strcmp(get(opt.in,'type'),'figure')
                        figure(opt.in), set(opt.in,'menubar','none')
                    end
                    D.hu = uicontrol('parent',D.hp,'units','normalized','pos',[0 0 1 1]);
                case 'uicontrol'
                    D.hu = opt.in;
                    D.hp = get(opt.in,'parent');
                case 'axes'
                    % replace axes by an uicontrol
                    ha = opt.in;
                    D.hp = get(ha,'parent');
                    D.hu = uicontrol('parent',D.hp,'units',get(ha,'units'),'pos',get(ha,'pos'));
                    delete(ha)
                otherwise
                    error('bad handle')
            end
            D.hf = fn_parentfigure(D.hp);
            if isempty(get(D.hf,'Tag')), set(D.hf,'Tag','used by fn4D'), end
            
            % context menu
            initlocalmenu(D)
            
            % list and event (bottom-up)
            set(D.hu,'style','listbox','min',0,'max',2, ...
                'callback',@(hu,evnt)event(D,'select'))
            if boolean(D.scrollwheel)
                D.scrollwheel = 'on'; % this will automaticall register scroll wheel
            end
            
            % communication with parent
            addparent(D,D.SI)
            
            % auto-delete
            set(D.hu,'deletefcn',@(u,evnt)delete(D))

            % update display (here, just sets the correct value)
            displayselection(D)
            displaycross(D)
            displaylabel(D)
            
            % set more properties
            if ~isempty(optadd)
                set(D,optadd{:})
            end
        end
        function initlocalmenu(D)
            delete(D.menu)
            D.menu = uicontextmenu('parent',D.hf);
            m = D.menu;
            set(D.hu,'UIContextMenu',m)
            
            uimenu(m,'label','new singleton selections','callback',@(u,e)event(D,'newuni'))
            uimenu(m,'label','new group selection','callback',@(u,e)event(D,'newgroup'))
            uimenu(m,'label','add to selection','callback',@(u,e)event(D,'add'))
            
            uimenu(m,'label','remove highlighted group(s)','callback',@(u,e)event(D,'rmgroup'),'separator','on')
            uimenu(m,'label','remove all groups','callback',@(u,e)event(D,'rmgroupall'))
            uimenu(m,'label','remove highlighted individuals','callback',@(u,e)event(D,'rmuni'))
            uimenu(m,'label','remove all individuals','callback',@(u,e)event(D,'rmuniall'))
            uimenu(m,'label','remove highlighted selections','callback',@(u,e)event(D,'rm'))
            uimenu(m,'label','remove all selections','callback',@(u,e)event(D,'rmall'))

            D.menuitems.selmultin = uimenu(m,'separator','on','checked',onoff(D.selmultin), ...
                'label','temporary selection: individuals','callback',@(u,e)set(D,'selmultin',~D.selmultin));
            uimenu(m,'label','select all','accelerator','A','callback',@(u,e)event(D,'selectall'))
            
            m1 = uimenu(m,'label','scroll wheel','separator','on');
            D.menuitems.scrollwheel = uimenu(m1,'label','activated', ...
                'checked',D.scrollwheel,'callback',@(u,e)set(D,'scrollwheel',fn_switch(D.scrollwheel,'toggle')));
            uimenu(m1,'label','make default in figure', ...
                'callback',@(u,e)set(D,'scrollwheel','default'));
        end
        function delete(D)
            if ishandle(D.hu), delete(D.hu), end
            if ishandle(D.hlabel), delete(D.hlabel), end
        end
    end
       
    % Events
    methods
        function event(D,flag,varargin)
            fn4D_dbstack
            
            % modify flag
            switch flag
                case 'select'
                    if strcmp(get(D.hf,'selectiontype'),'open')
                        flag = 'unisel';
                    end
                case 'selectall'
                    set(D.hu,'value',1:D.SI.sizes)
                    flag = 'select';
            end            
            
            % get the current 'temporary' selection
            val = get(D.hu,'value');
            nval = length(val);
            selectionmarks = D.SI.selection.getselset(1).singleset;
            selinds = {selectionmarks.dataind};
            isunisel = fn_map(@isscalar,selinds);
            nsel = length(selectionmarks);
            ntemp = length(D.tempsel);
            if nsel<ntemp || ~isequal(selinds(nsel-ntemp+1:nsel),D.tempsel)
                % selections have been modified outside from the list,
                % loosing the distinction between 'solid' and 'temporary'
                % selections (all become solid)
                D.tempsel = {};
                ntemp = 0;
            end
            nselsolid = nsel-ntemp;
            
            % action
            switch flag
                case 'select'
                    if isscalar(val)
                        D.SI.ij2 = val;
                    end
                    if isempty(val) || (isscalar(val) && nselsolid==0) ...
                            || (isscalar(val) && isequal({val},D.tempsel))
                        % remove temporary selection in the following
                        % cases:
                        % - user unselected all list items
                        % - no solid selection and a single selected item
                        % - repeated selection of temporaray selection item
                        if ntemp
                            D.tempsel = {};
                            updateselection(D.SI,'remove',nselsolid+(1:ntemp))
                        end
                    else
                        % new temporary selection
                        if D.selmultin
                            for i=1:nval, sel(i) = selectionND('point1D',val(i)); end %#ok<AGROW>
                            tempselnew = num2cell(val);
                        else
                            sel = selectionND('point1D',val);
                            tempselnew = {val};
                        end
                        ntempnew = length(tempselnew);
                        % the update happens potentially in several steps;
                        % memory of the temporary selection must be updated
                        % accordingly to preserve the mechanism for marking
                        % differently 'solid' and 'temporary' selections
                        D.downnotify = false;
                        nchg = min(ntemp,ntempnew);
                        if ntemp>nchg
                            D.tempsel = tempselnew(1:ntempnew);
                            updateselection(D.SI,'remove',nselsolid+(nchg+1:ntemp))
                        end
                        if nchg>0
                            D.tempsel = tempselnew(1:nchg);
                            updateselection(D.SI,'change',nselsolid+(1:nchg),sel(1:nchg))
                        end
                        if ntempnew>nchg
                            D.tempsel = tempselnew(1:ntempnew);
                            updateselection(D.SI,'new',[],sel(nchg+1:ntempnew))
                        end
                        % ... and update the display manually
                        D.displayselection()
                        D.downnotify = true;
                    end
                case 'unisel'
                    % double-click -> make new selection with current
                    % index, or remove it
                    if ~isscalar(val), return, end
                    kunisel = find(isunisel(1:nselsolid));
                    f = ([selinds{kunisel}]==val);
                    f = kunisel(f);
                    if ~isempty(f)
                        updateselection(D.SI,'remove',f)
                    elseif isequal(D.tempsel,{val}) && isequal(selinds{nsel},val)
                        % nothing to do: "temporary selection" has been
                        % "solidified"; but de-highlight item in list
                        D.tempsel = {};
                    else
                        updateselection(D.SI,'new',[],selectionND('point1D',val))
                    end
                    displayselection(D)
                case 'newuni'
                    set(D.hu,'value',[]) % let's not get confused!
                    D.tempsel = {};
                    for i=1:nval, sel(i) = selectionND('point1D',val(i)); end %#ok<AGROW>
                    nchg = min(ntemp,nval);
                    if ntemp>nchg
                        updateselection(D.SI,'remove',nchg+1:ntemp)
                    end
                    if nchg>0
                        updateselection(D.SI,'change',nselsolid+(1:nchg),sel(1:nchg))
                    end
                    if nval>nchg
                        updateselection(D.SI,'new',[],sel(nchg+1:nval))
                    end
                case 'newgroup'
                    set(D.hu,'value',[]) % let's not get confused!
                    D.tempsel = {};
                    sel = selectionND('point1D',val);
                    if ntemp
                        if ntemp>1, updateselection(D.SI,'remove',nselsolid+(2:ntemp)), end
                        updateselection(D.SI,'change',nselsolid+1,sel)
                    else
                        updateselection(D.SI,'new',[],sel)
                    end
                case 'add'
                    set(D.hu,'value',[]) % let's not get confused!
                    sel = selectionND('point1D',val);
                    if ntemp, updateselection(D.SI,'remove',nselsolid+(1:ntemp)), end
                    if nselsolid
                        updateselection(D.SI,'add',nselsolid,sel)
                    else
                        updateselection(D.SI,'new',[],sel)
                    end
                case {'rmgroup' 'rmgroupall' 'rmuni' 'rmuniall' 'rm' 'rmall'}
                    set(D.hu,'value',[]) % let's not get confused!
                    if strfind(flag,'all')
                        range = 1:nsel;
                        flag = strrep(flag,'all','');
                    else
                        range = val;
                    end
                    rmmask = false(1,nsel);
                    switch flag
                        case 'rmgroup'
                            rmmask(range) = ~isunisel(range);
                        case 'rmuni'
                            rmmask(range) = isunisel(range);
                        case 'rm'
                            rmmask(range) = true;
                    end
                    rmmask(nsel-ntemp+1:nsel) = true;
                    updateselection(D.SI,'remove',find(rmmask))
                case 'scroll'
                    n = varargin{1};
                    D.SI.ij2 = fn_coerce(D.SI.ij2+n,1,D.SI.sizes);
            end
            if isempty(D.SI.selectionmarks)
                set(D.hu,'value',D.SI.ij)
            end
        end
        function updateDown(D,~,evnt)
            fn4D_dbstack
            switch evnt.flag
                case 'ij2'
                    displaycross(D)
                case 'selection'
                    displayselection(D)
                case 'sizes'
                    displayselection(D)
                    displaycross(D)
                case 'labels'
                    displaylabel(D)
                case 'units'
                    displayselection(D)
            end
        end
    end
        
    % Get/Set - scroll wheel
    methods
        function set.scrollwheel(D,flag)
            switch flag
                case 'on'
                    fn_scrollwheelregister(D.hu,@(n)event(D,'scroll',n)) %#ok<MCSUP>
                    D.scrollwheel = 'on';
                case 'default'
                    fn_scrollwheelregister(D.hu,@(n)event(D,'scroll',n),'default') %#ok<MCSUP>
                    D.scrollwheel = 'on';
                case 'off'
                    fn_scrollwheelregister(D.hu,flag) %#ok<MCSUP>
                    D.scrollwheel = 'off';
                otherwise
                    error 'scrollwheel value must be ''off'', ''on'' or ''default'''
            end
            set(D.menuitems.scrollwheel,'checked',D.scrollwheel) %#ok<MCSUP>
        end
    end
    
    % Get/Set
    methods
        function set.selmultin(D,val)
            if val==D.selmultin, return, end
            D.selmultin = val;
            % update menu item
            set(D.menuitems.selmultin,'checked',onoff(val)) %#ok<MCSUP>
            % update selection
            event(D,'select')
        end
    end
    
    % Display
    methods
        function displayselection(D)
            fn4D_dbstack
            
            % init list with names of items
            n = D.SI.sizes;
            if iscell(D.SI.units{1})
                itemnames = D.SI.units{1};
                % handle mismatch between length of list and number of
                % names
                navail = length(itemnames);
                if n<=navail
                    str = itemnames(1:n);
                else
                    str = [itemnames(1:navail) fn_num2str(navail+1:n,'%i ','cell')];
                end
            else
                grid = D.SI.grid;
                values = (1:n)*grid(1) + grid(2);
                str = fn_num2str(values,'%.4g ','cell');
            end
            
            selectionmarks = D.SI.selection.getselset(1).singleset;
            selinds = {selectionmarks.dataind};
            nsel = length(selinds);
            
            % "temporaray" selections (i.e. currently selected by the user
            % in the list) are marked differently than "solid" selections
            ntemp = length(D.tempsel);
            if nsel<ntemp || ~isequal(selinds(nsel-ntemp+1:nsel),D.tempsel)
                % selections have been modified outside from the list,
                % loosing the distinction between 'solid' and 'temporary'
                % selections (all become solid)
                D.tempsel = {};
                ntemp = 0;
            end
            nselsolid = nsel-ntemp;
            
            % mark solid selections
            isunisel = fn_map(@isscalar,selinds(1:nselsolid));
            for ksel=find(isunisel)
                ind = selinds{ksel};
                str{ind} = [str{ind} '[' num2str(ksel) ']'];
            end
            for ksel=find(~isunisel)
                for ind = selinds{ksel}
                    str{ind} = [str{ind} '[group' num2str(ksel) ']'];
                end
            end
            
            % mark temporary selections
            for i = 1:length(D.tempsel)
                marker = fn_switch(isscalar(D.tempsel{i}),'*','G');
                for ind = D.tempsel{i}
                    str{ind} = [str{ind} marker];
                end
            end
            
            % update display!
            set(D.hu,'string',str)
        end
        function displaycross(D)
            fn4D_dbstack        
            set(D.hu,'value',D.SI.ij);
        end
        function displaylabel(D)
            if isempty(D.SI.labels)
                delete(D.hlabel)
                D.hlabel = [];
            else
                D.hlabel = uicontrol('parent',D.hp,'style','text','string',D.SI.labels{1}, ...
                    'horizontalalignment','center');
                fn_controlpositions(D.hlabel,D.hu,[0 1 1 0],[0 2 0 15])
            end
        end
            
    end
    
end
     
   