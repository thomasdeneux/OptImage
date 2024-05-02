classdef tps_stimtable < hgsetget
    % see also: tps_trial, tps_conditionslist
    
    properties
        table = struct('id',cell(1,0),'name',cell(1,0),'stim',cell(1,0),'type',cell(1,0));
        user
    end
    
    % Creation: no need for a constructor
    
    % Get stim / Add new stim
    methods
        function kstim = getentry(T,id)
            n = length(id);
            if n==1
                kstim = find(id==[T.table.id]);
            else
                % this will generate an error if some entry does not appear
                % in the table
                kstim = zeros(1,n);
                for k=1:n, kstim(k) = find(id(k)==[T.table.id]); end
            end
        end
        function stim = getstim(T,id)
            kstim = getentry(T,id);
            if isempty(kstim)
                disp 'could not find stim entry for this identifier'
                stim = [];
                return
            end
            stim = T.table(kstim).stim;
        end
        function stim = getglobalstim(T,ids)
            kstim = getentry(T,ids);
            stims = {T.table(kstim).stim};
           
            %             % valid output only if the non-empty ones are all equal
            %             stim = [];
            %             for k=1:length(stims)
            %                 stimk = stims{k};
            %                 if isempty(stimk)
            %                     % nothing to do
            %                 elseif isempty(stim)
            %                     stim = stimk;
            %                 elseif ~isequal(stimk,stim)
            %                     stim = [];
            %                     return
            %                 end
            %             end

            % union of all stims
            stims = [stims{:}];
            stim = unique(stims','rows')';
        end
        function s = getdetails(T,id)
            if ~isvector(id), error 'id must be a scalar or vector', end
            kstim = getentry(T,id);
            if isempty(kstim)
                disp 'could not find stim entry for this identifier'
                s = fn_structinit(T.table,length(id));
                [s.id] = dealc(id);
                return
            end
            s = T.table(kstim);
        end
        function id = addstim(T,varargin)
            % function id = addstim(T,id)
            % function id = addstim(T,stim)
            % function id = addstim(T,'name1',value1,...)
            % function id = addstim(T,struct)
            %---
            % note that id can be specified by user; if not, or if it is
            % empty, a hash number is used
            
            % Input
            a = varargin{1};
            if isstruct(a)
                s = a;
            elseif ischar(a)
                s = struct(varargin{:});
            elseif isscalar(a)
                s = struct('id',a);
            else
                s = struct('stim',a);
            end
            if isfield(s,'id') && ~isempty(s.id)
                id = s.id;
            else
                id = fn_hash(s,'num',4);
                s.id = id;
            end
            
            % Get entry
            kstim = getentry(T,id);
            if isempty(kstim)
                % create new entry
                kstim = length(T.table)+1;
                T.table(kstim).id = id;
                % sort table by id
                [ids ord] = sort([T.table.id]); %#ok<ASGLU>
                T.table = T.table(ord);
                kstim = find(ord==kstim);
            end
            
            % Fill-in details
            F = fieldnames(s);
            for k=1:length(F)
                f = F{k};
                if ~isfield(T.table(kstim),f) || isempty(T.table(kstim).(f))
                    T.table(kstim).(f) = s.(f);
                elseif ~isequal(T.table(kstim).(f),s.(f)) && ~isempty(s.(f))
                    prompt = 'Incompatible definitions! Which to use?';
                    a = fn_num2str(T.table(kstim).(f));
                    b = fn_num2str(s.(f));
                    answer = questdlg(prompt,'Stim Def',a,b,a);
                    T.table(kstim).(f) = answer;
                end
            end
            
            % Output?
            if nargout==0, clear id, end
        end
        function updatetable(T,s)
            % s is a structure with updated fields compared to T.table; we
            % will apply the updates, while making sure that important
            % information is not removed
            
            % remove fields
            Fold = fieldnames(T.table);
            base = {'id' 'name' 'stim' 'type'};
            Fs = fieldnames(s);
            Fs = Fs(:)';
            Fnew = [base setdiff(Fs,base)];
            frm = setdiff(Fold,Fnew);
            if ~isempty(frm), T.table = rmfield(T.table,frm); end
            
            % add information
            for k=1:length(s)
                kstim = find([T.table.id]==s(k).id);
                if isempty(kstim), continue, else kstim = kstim(1); end
                for i=1:length(Fs)
                    f = Fs{i};
                    if strcmp(f,'id'), continue, end
                    T.table(kstim).(f) = s(k).(f);
                end
            end
            
            % reorder fields
            T.table = orderfields(T.table,Fnew);
        end
        function T = mergetables(varargin)
            % function T = mergetables(T1,T2,....)
            % function T = mergetables(TT) 
            T = tps_stimtable;
            TT = unique([varargin{:}]);
            for k=1:length(TT)
                tablek = TT(k).table;
                for i=1:length(tablek)
                    addstim(T,tablek(i))
                end
            end
        end
    end
    
    
end