function label = tpv_ao_meslist(V) 
% function f(V) 
% function label = f('label') 

if ischar(V) && strcmp(V,'label')
    label = 'List of trials (MES)';
    return
end

% build list of trials
prefix = cell(1,V.ntrial);
list = cell(1,V.ntrial);
for k=1:V.ntrial
    Tk = V.content.trials(k);
    str = num2str(k);
    str = [str repmat(' ',1,2*(3-length(str))) '[' Tk.type];
    if isfield(Tk.fullinfo,'MES') && isfield(Tk.fullinfo.MES(1),'FileSubindex')
        str = [str num2str(Tk.fullinfo.MES(1).FileSubindex)];
    end
    str = [str '] ']; %#ok<*AGROW>
    prefix{k} = str;
    if isfield(Tk.addinfo,'comment')
        str = [str Tk.addinfo.comment]; 
    end
    list{k} = str;
end

% figure
hf = figure(8298);
fn_setfigsize(hf,400,600)
set(hf,'handlevisibility','off','numbertitle','off','name','MES Trial List', ...
    'deletefcn',@(u,e)deletefig())

% figure for image
hfi = [];
imdisplay = [];

% controls
ul = uicontrol('parent',hf,'style','listbox','string',list,'callback',@(u,e)settrial());
um = uicontrol('parent',hf,'style','frame','buttondownfcn',@(u,e)moveframe());
ut = uicontrol('parent',hf,'style','edit','max',2,'horizontalalignment','left','callback',@(u,e)writecomment());
set([ul um ut],'units','normalized')
setframe(.2)
delete(findobj(hf,'type','uimenu'))
m = uimenu(hf,'label','Image');
items.image = uimenu(m,'label','show image','callback',@(u,e)showimage('toggle','image'));
items.points = uimenu(m,'label','show points','callback',@(u,e)showimage('toggle','points'),'checked','on');
% items.numbers = uimenu(m,'show numbers','callback',@(u,e)showimage('toggle','numbers','checked','on'));

% listen to some changes in V
hl_trial = addlistener(V,'ktrial','PostSet',@(u,e)showtrial());
hl_close = addlistener(V,'EventCloseFile',@(u,e)close(hf));
showtrial()


%-- sub-functions --
    function settrial()
        ktrial = get(ul,'value');
        showtrial(ktrial)
        if strcmp(get(hf,'selectiontype'),'open')
            V.ktrial = ktrial;
        end
    end

    function showtrial(ktrial)
        if nargin<1, ktrial=V.ktrial; end
        set(ul,'value',ktrial)
        Tk = V.content.trials(ktrial);
        if isfield(Tk.addinfo,'comment')
            set(ut,'string',Tk.addinfo.comment)
        else
            set(ut,'string','')
        end
        showimage()
    end

    function showimage(flag,item)
        do.image = boolean(get(items.image,'checked'));
        do.points = boolean(get(items.points,'checked'));
        %do.numbers = boolean(get(items.numbers,'checked'));
        if nargin>=1 && strcmp(flag,'toggle')
            do.(item) = ~do.(item);
            set(items.(item),'checked',onoff(do.(item)))
        end
        if ~do.image || ~isfield(Tk.fullinfo,'MES'), return, end

        fn_watch(hf,'startnow')
        if isempty(hfi) || ~ishandle(hfi)
            hfi = figure('integerhandle','off','numbertitle','off','handlevisibility','off', ...
                'name','MES Image');
        end
        if isempty(imdisplay) || ~isvalid(imdisplay.D)
            clf(hfi)
            ha = axes('parent',hfi);
            imdisplay = struct;
            imdisplay.SI = sliceinfo(2,'units',{'um' 'um'});
            imdisplay.D = activedisplayImage(imdisplay.SI,'in',ha, ...
                'channelcolors',[0 1 0; 1 0 0],'scaledisplay','tick','selshow',do.points);
        elseif nargin>=1 && strcmp(item,'points')
            imdisplay.D.selshow = do.points;
            fn_watch(hf,'stop')
            return
        end
        if strcmp(Tk.type,'linescan')
            h1 = Tk.fullinfo.MES(1);
            h3 = Tk.fullinfo.MES(end-1);
            [x y] = fn_loadvar(Tk.file,h3.ImageName,Tk.fullinfo.MES(end).ImageName);
            x = cat(3,x,y);
            imdisplay.SI.slice.data = x;
            imdisplay.SI.grid = [h3.WidthStep h3.WidthOrigin; h3.HeightStep h3.HeightOrigin];
            lines2 = h1.info_Linfo.lines(2);
            anysegment = any(fn_map(@(a)size(a,2),lines2.line1)>1);
            if anysegment
                polyreal = lines2.line2; % actual scanned points
            else
                polyreal = [lines2.line1{:}]; % ROI points (there are no ROI segments)
            end
            polyreal = polyreal(1:2,:);
            polyidx = AX2IJ(imdisplay.SI,polyreal);
            sel = selectionND('point2D',num2cell(polyidx(1:2,:),1),[h3.Width h3.Height]);
            selset = selectionset([h3.Width h3.Height],2,sel);
            imdisplay.SI.setselection(selset)
        else
            imdisplay.SI.slice.data = 0;
        end
        fn_watch(hf,'stop')
    end

    function setframe(hpos)
        d = .01;
        hpos = fn_coerce(hpos,2*d,1-2*d);
        fn_set([ul um ut],'position',{[0 hpos+d 1 1-hpos-d] [0 hpos 1 1e-3] [0 0 1 hpos-d]})
    end

    function moveframe()
        p0 = get(hf,'CurrentPoint'); y0 = p0(2);
        figpos = get(hf,'pos'); H = figpos(4);
        frpos = get(um,'pos'); hpos0 = frpos(2);
        fn_buttonmotion(@moveframesub,hf)
        function moveframesub()
            p = get(hf,'CurrentPoint'); y = p(2);
            setframe(hpos0+(y-y0)/H)
        end
    end
    
    function deletefig()
        delete(hl_trial)
        delete(hl_close)
        if ishandle(hfi), close(hfi), end
    end

    function writecomment()
        ktrial = get(ul,'value');
        if ~isscalar(ktrial), set(ut,'string',''), return, end
        comment = get(ut,'string');
        list{ktrial} = [prefix{ktrial} comment];
        set(ul,'string',list)
        V.content.trials(ktrial).addinfo.comment = comment;
    end

end