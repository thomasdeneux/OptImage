classdef activedisplayExt < activedisplay
    % opt.command should be of one of the two forms:
    % {@updatefcn,par1,par2,...}    where updatefcn = function(focus,par1,par2,...)
    % {'command','focusassignname'} where 'focusassignname' is the variable
    %                               name to which to assign the focus object in
    %                               the base workspace, and 'command' should
    %                               use this name
    
    properties
        functionflag = false;
        focusassignname = '';
        command
    end
    
    methods
        function D = activedisplayExt(opt,dttr,dttropt,actdispopt)
            D = D@activedisplay(opt,dttr,dttropt,actdispopt);
        end
        
        function ha = init(D,SI,opt)
            command = opt.command;
            if ~iscell(command), error('command must be a cell array'), end
            if isa(command{1},'function_handle')
                D.functionflag = true;
                D.command = command;
            elseif ischar(command{1})
                D.functionflag = false;
                D.command = command{1};
                D.focusassignname = command{2};
            else
                error('command must be a cell array starting with function handle or string')
            end
            
            if ~isempty(opt.in)
                ha = opt.in;
            else
                ha = figure('integerHandle','off','visible','off');
            end
            setappdata(ha,'fn4Dext',true)
            
            update(D,F)
        end
        
        function update(D,F)
            if D.functionflag
                feval(D.command{1},F,D.command{2:end})
            else
                assignin('base',D.focusassignname,F)
                evalin('base',D.command)
            end
        end
    end
    
    methods (Static)
        function [defaultopt requiredopt] = OptionsList(obj)
            defaultopt = struct( ...
                'command',      '', ...
                'in',           [] ...
                );
            requiredopt = {'command'};
        end
    end
end