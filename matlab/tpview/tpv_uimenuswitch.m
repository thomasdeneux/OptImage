function tpv_uimenuswitch(V,keepvisible)
% function tpv_uimenuswitch(V,keepvisible)
%---
% optional actions that will be performed after pressing any uimenu in the
% tpview window: keep the menu visible, show the command

persistent dovis 
dovis = keepvisible;

hf = V.hf;
menus = findall(hf,'type','uimenu');

for i=1:length(menus)
   
    % get callback
    m = menus(i);
    fun = get(m,'callback');
    if isempty(fun)
        continue
    end
    if isa(fun,'function_handle')
        tokens = regexp(char(fun),'^@\([^\)]*\)([^\(]*)','tokens');
        if ~isempty(tokens) && strcmp(tokens{1}{1},'uimenuaction')
            % callback is already calling uimenuaction
            continue
        end
    end

    % modify it to calling uimenuaction
    set(m,'callback',@(u,e)uimenuaction(u,e,fun))
    
end


    function uimenuaction(m,e,callback)
        if dovis
            mi = m;
            while strcmp(get(mi,'type'),'uimenu')
                set(mi,'visible','on')
                mi = get(mi,'parent');
            end
        end
        fn_evalcallback(callback,m,e)
    end

end


