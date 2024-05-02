function h=tps_displaystim(stim,varargin)
% function h=tps_displaystim(stim[,ha][,'realwidth'][,patch options])
%---
% Display stimulation as a set of vertical lines or rectangles spanning
% vertically the graph.
%
% Input:
% - stim    a vector or a 2xn array - stimulation times, and the second row
%           is the duration of each stim
% - ha      axes handle [default: current axes]
% - 'realwidth'     flag to use to force the display of stimulation as
%           rectangles when stimulation duration are specified (otherwise,
%           is duration are too short, i.e. less than 2 pixels on screen,
%           stimulations are marked with vertical lines rather than thin
%           rectangles) 
% - patch options   pairs of property names / values for the displayed
%           lines or rectangle; note that 'color' name will automatically
%           be replaced by 'facecolor' when rectangles are displayed

if nargin==0, help tps_displaystim, return, end

% input
if isempty(stim), return, end
if size(stim,1)==1, stim(2,:) = 0; end
% (axes handle)
if nargin>1 && isscalar(varargin{1}) && ishandle(varargin{1})
    ha = varargin{1};
    varargin(1) = [];
else
    ha = gca; 
end
% ('realwidth' flag)
if length(varargin)>=1 && ischar(varargin{1}) && strcmp(varargin{1},'realwidth')
    doadjustwidth = false;
    varargin(1) = [];
else
    doadjustwidth = true;
end
ax = axis(ha);

% increase stimulation widths to make them visible
if doadjustwidth
    px = fn_coordinates(ha,'b2a',[1 1],'vector'); pxw = px(1);
    okpatch = (stim(2,:)>=2*pxw);
else
    okpatch = (stim(2,:)>0);
end

% show stim
if any(okpatch)
    xdata = [stim(1,okpatch); stim(1,okpatch); stim(1,okpatch)+stim(2,okpatch); stim(1,okpatch)+stim(2,okpatch)];
    ydata = repmat(ax([3 4 4 3])',1,sum(okpatch));
    icol = find(strcmp(varargin,'color'));
    if ~isempty(icol), varargin{icol} = 'facecolor'; end
    h = patch('xdata',xdata,'ydata',ydata,'facecolor',[1 1 1]*.6,'edgecolor','none','parent',ha,varargin{:});
else
    h = [];
end
if any(~okpatch)
    xdata = [stim(1,~okpatch); stim(1,~okpatch)];
    ydata = repmat(ax([3 4])',1,sum(~okpatch));
    h = [h line(xdata,ydata,'color',[1 1 1]*.6,'parent',ha,varargin{:})];
end
for i=1:length(h), uistack(h(i),'bottom'), end % make the patch appear below any other existing object
if nargout==0, clear h, end
axis(ha,ax) % don't know why, sometimes the axis is changed when drawing the patch

