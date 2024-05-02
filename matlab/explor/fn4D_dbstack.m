function fn4D_dbstack(str,flag)
% function fn4D_dbstack                     -> displays function (mode 'all')
% function fn4D_dbstack(str)                -> displays string ('on' and 'all')
% function fn4D_dbstack('set','off|on|all') -> defines mode (default='off')
%---
% displays current function name, with indent according to stack length

%#ok<*NASGU,*UNRCH>

persistent mode 

if isempty(mode), mode = 0; end
if nargin==2 
    mode = fn_switch(flag,'off',0,'on',1,'all',2);
    return
end

if ~mode || (mode==1 && nargin==0), return, end

ST = dbstack;

n = 0;
for k=2:length(ST)
    if ~any(strfind(ST(k).name,'@'))
        n = n+1;
    end
end
blanks = repmat(' ',1,n);

if nargin<1, str = ST(2).name; end

disp([blanks str])

