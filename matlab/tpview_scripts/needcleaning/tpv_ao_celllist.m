function label = tpv_ao_celllist(V) 
% function f(V) 
% function label = f('label') 

if ischar(V) && strcmp(V,'label')
    label = 'List of cells';
    return
end


P = projection(V.a4d.G,1);
hf = fn_figure('tpv_ao_celllist');
L = activedisplayList(P,'in',hf);

hl = addlistener(V,'ktrial','PostSet',@(u,e)closelistifneeded());

    function closelistifneeded()
        if ~isvalid(L)
            delete(hl)
        elseif V.nx~=L.SI.sizes
            close(L.hf)
            delete(hl)
        end
    end

end