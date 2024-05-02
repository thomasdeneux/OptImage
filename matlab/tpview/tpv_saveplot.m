function tpv_saveplot(varargin)
% function tpv_saveplot([V][,fname|'auto'][,'close'])

% input
V = [];
fname = '';
doclose = false;
for k=1:length(varargin)
    a = varargin{k};
    if isa(a,'tpview')
        V = a;
    elseif ischar(a)
        switch a
            case 'close'
                doclose = true;
            otherwise
                fname = a;
        end
    else
        error argument
    end
end
if isempty(V), V = evalin('base','V'); end

ha = V.grob.time;
pixpos = fn_pixelpos(ha);

% spaces between axes sides and figure sides
left = 25;
right = 10;
bottom = 40;
top = 25;

% new axes
hf = figure('color','w');
fn_setfigsize(hf,left+pixpos(3)+right,bottom+pixpos(4)+top)

hb = copyobj(ha,hf);
set(hb,'units','pixel','pos',[left bottom pixpos(3:4)],'units','normalized')
hc = findobj(hb);
set(hc,'buttondownfcn','','deletefcn','');

% some more properties
set(hb,'box','on')

% % children
% lines = copyobj(findobj(ha,'tag','fn4D_line'),hb);
% set(lines,'buttondownfcn','') % remove callbacks
% copyobj(get(ha,'xlabel'),hb)
% copyobj(get(ha,'title'),hb)

% save 
if isempty(fname)
    fname = fn_savefile('Select file where to save figure');
    if ~fname, return, end
elseif strcmp(fname,'auto')
    fname = fn_autofigname(V.hf);
end
fn_savefig(hf,fname)

% close
if doclose
    close(hf)
end