function label = script_triallist(V) 
% function f(V) 
% function label = f('label') 

if ischar(V) && strcmp(V,'label')
    label = 'List of trials';
    return
end
 
% build list of trials
prefix = cell(1,V.ntrial);
list = cell(1,V.ntrial);
for k=1:V.ntrial
    Tk = V.content.trials(k);
    num = num2str(k);
    s = Tk.sizes;
    s = s(1:find(s>1,1,'last'));
    str = [num repmat(' ',1,2*(3-length(num))) '[' Tk.type ' ' fn_strcat(s,'x')];
    str = [str '] '];
    if isfield(Tk.addinfo,'description'), str = [str ' "' Tk.addinfo.description '"']; end %#ok<AGROW>
    prefix{k} = str;
    if isfield(Tk.addinfo,'comment')
        str = [str Tk.addinfo.comment]; 
    end
    list{k} = str;
end

% figure
hf = figure(8298);
fn_setfigsize(hf,400,600)
set(hf,'handlevisibility','off','numbertitle','off','name','Trial List', ...
    'color',get(V.hf,'color'),'menubar','none','deletefcn',@(u,e)deletefig())

% controls
ul = uicontrol('parent',hf,'style','listbox','string',list,'callback',@(u,e)settrial());
um = uicontrol('parent',hf,'style','frame','buttondownfcn',@(u,e)moveframe());
ut = uicontrol('parent',hf,'style','edit','max',2,'horizontalalignment','left','callback',@(u,e)writecomment());
set([ul um ut],'units','normalized')
setframe(.2)
delete(findall(hf,'type','uimenu'))

% listen to some changes in V
hl_trial = addlistener(V,'ktrial','PostSet',@(u,e)showtrial());
hl_close = addlistener(V,'EventCloseFile',@(u,e)close(hf));

% display current trial
showtrial()
if isempty(get(ut,'string')), set(ut,'string','Write your comment here'), end

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
    end

    function setframe(hpos)
        dx = 0;
        dy = .01;
        hpos = fn_coerce(hpos,2*dy,1-2*dy);
        fn_set([ul um ut],'position',{[dx hpos+dy 1-2*dx 1-hpos-2*dy] [0 hpos 1 1e-3] [dx dy 1-2*dx hpos-2*dy]})
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