function [condarray condfieldname] = tps_readconditionname(condname,nc)
    % ex: if nc=8 and condname = '(C0+C1)/(C3-(C4+C5))'
    % -> condarray = {[1 2],[],[4],[5 6]}
    % if condname = 'C7,C8'
    % -> condarray = {[7],[],[],[]; [8],[],[],[]}
    % confieldname is obtained as a hash number
    
    % TODO: check that condition is valid with nc?
    
    % special case 
    if strcmp(condname,'all conditions')
        condarray = cell(nc,4);
        for i=1:nc, condarray{i,1} = i; end
        return
    end
    
    % separate multiple calculations
    tokens = regexp(condname,'([^,]*)','tokens');
    condnames = [tokens{:}];
    
    % loop on calculations
    condarray = cell(0,4);
    for kcond=1:length(condnames)
        tokens = regexp(condnames{kcond},'([^-/]*)(-[^-/]*)?(/[^-/]*)?(-[^-/]*)?','tokens');
        if ~isscalar(tokens), error('bad condition name ''%s''',condname), end
        tokens = tokens{1};
        for i=1:4
            conds = regexp(tokens{i},'([0-9]*)','match');
            for k=1:length(conds), conds{k}=str2num(conds{k}); end %#ok<ST2NM>
            conds = cell2mat(conds);
            condarray{kcond,i} = conds+1;
        end
    end
    
    % create a compact name (but do not use fn_hash, too slow!)
    if nargout>=2
        tmp = false(max([condarray{:}]),4);
        for i=1:4, tmp(condarray{i},i)=true; end
        namebin = char('0'+tmp);
        condfieldname = char('a' + bin2dec(namebin)');
    end
end
