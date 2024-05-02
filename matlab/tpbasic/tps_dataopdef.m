classdef tps_dataopdef
    % function op = tps_dataopdef(name,active,value)
    % function op = tps_dataopdef(oldversion)
    %---
    % Input:
    % - name        string
    % - active      logical
    % - value       cell array
    % - oldversion  structure or cell array of cell arrays
    %
    % one can also use cell arrays for name, active, value: this define
    % a vector object
    
    properties 
        name
        active
        value
    end
    
    % Construction, conversion and correction
    methods
        function op = tps_dataopdef(varargin)
            % convert to tps_dataopdef
            switch nargin
                case 0
                    return
                case 1
                    x = varargin{1};
                    if isa(x,'tps_dataopdef')
                        op = x;
                    elseif isempty(x) || (iscell(x) && (isempty(x{1}) || (iscell(x{1}) && isempty(x{1}{1})))) % some strange old versions have strange dataopdef...
                        op = tps_dataopdef.empty(1,0); % empty opdef
                    elseif isstruct(x)
                        for i=1:numel(x)
                            op(i).name = x(i).name; %#ok<AGROW>
                            op(i).active = x(i).active; %#ok<AGROW>
                            op(i).value = x(i).value; %#ok<AGROW>
                        end
                    elseif iscell(x)
                        if ~iscell(x{1})
                            x = {x};
                        end
                        for i=1:numel(x)
                            vali = x{i};
                            op(i).name = vali{1}; %#ok<AGROW>
                            op(i).active = true; %#ok<AGROW>
                            op(i).value = vali(2:end); %#ok<AGROW>
                        end
                    else
                        error argument
                    end
                case 3
                    [op.name op.active op.value] = deal(varargin{:});
                otherwise
                    error argument
            end
            % correct values if necessary
            op = convert(op);
        end
        function op = convert(op)
            % multiple object
            if ~isscalar(op)
                for i=1:length(op), op(i) = convert(op(i)); end
                return
            end
            
            % attempt to correct
            switch lower(op.name)
                case 'heart'
                    switch length(op.value)
                        case 6
                            % up-to-date version, still some names may have
                            % changed!
                        case 4
                            % old version: missing indices
                            op.value = ['indices' ':' op.value(1:4)];
                        case 3
                            % old version: missing indices and precision
                            op.value = ['indices' ':' op.value(1:2) .6 op.value(3)];
                        otherwise
                            error 'wrong number of values for ''HEART'' operation'
                    end
                    switch op.value{6}
                        case 'remove heart'
                            op.value{6} = 'rm';
                        case 'keep heart'
                            op.value{6} = 'keep';
                    end
                    if ~checkup(op), error 'cannot figure out what is wrong with ''HEART'' operation', end
                case 'threshold'
                    % this correction has been replaced by 'mask'
                    if length(op.value)~=1 || isnan(str2double(op.value{1})), error 'invalid operation specification', end
                    op.name = 'mask';
                    op.value = {'threshold on avg. value' op.value{1} 'min'};
                case 'mask'
                    if strcmp(op.value{1},'threshold')
                        op.value{1} = 'threshold on avg. value';
                    end
                otherwise
                    if ~checkup(op)
                        % just try to put default values into missing
                        % entries
                        specs = tps_dataopdef.operation_specifications;
                        kspec = strcmpi({specs.name},op.name);
                        speck = specs(kspec); 
                        nspeck = length(speck.controls);
                        op.value(nspeck+1:end) = [];
                        for i=length(op.value)+1:nspeck
                            op.value{i} = speck.controls(i).default;
                        end
                        if ~checkup(op), error 'operation specification seems not to be valid', end
                    end
            end
        end
        function ok = checkup(op)
            % automatic checks based on specifications
            ok = fn_supercontrol.checkup(op,tps_dataopdef.operation_specifications);
        end
    end
    
    % Specifications
    methods (Static)
        function specs = operation_specifications(V)
            % base structure for control definition
            s0 = struct('style','','string','','length',1,'default',[],'label','','labellength',0,'callback',[],'more',[]);
            specs = struct('name',{},'controls',{});
            % Space normalization
            s1 = s0; s2 = s0; s3 = s0;
            s1.style = 'popupmenu';
            s1.string = {'divide by','subtract','subtract without mean','DF/F'};
            s1.length = 3;
            s1.default = 'divide by';
            s2.style = 'popupmenu';
            s2.string = {'frames','time (s)','current selection'};
            s2.length = 3;
            s2.default = 'frames';
            s2.callback = @chgdefspace;
            function c = chgdefspace(c)
                switch c{2}
                    case 'frames'
                        c{3} = ':';
                    case 'time (s)'
                        c{3} = '0-1';
                    case 'current selection'
                        c{3} = '';
                end
            end
            s3.style = 'edit';
            s3.length = 2;
            s3.default = ':';
            specs(end+1) = struct('name','SPACE','controls',[s1 s2 s3]);
            % Time normalization
            s1 = s0; s2 = s0; s3 = s0;
            s1.style = 'popupmenu';
            s1.string = {'divide by','subtract','subtract without mean'};
            s1.length = 3;
            s1.default = 'divide by';
            s2.style = 'popupmenu';
            s2.string = {'indices','first selection','last selection','out selections','user'};
            s2.length = 3;
            s2.default = 'indices';
            s3.style = 'edit';
            s3.length = 2;
            s3.default = ':';
            specs(end+1) = struct('name','TIME','controls',[s1 s2 s3]);
            % Bleaching
            % very complex trick: pressing the 'EST' button result in
            % running tps_expfit on the current trial (assumed to be the
            % trial average); the computed parameters are stored through a
            % pointer to avoid duplicating large parameters data over all
            % trials 
            s = struct( ...
                'style',    {'stepper' 'pushbutton' 'popupmenu' 'checkbox' 'checkbox'}, ...
                'string',   {'' 'EST' {'remove bleach' 'keep only bleach'} 'glob.' 'blank'}, ...
                'default',  {2 [] 'remove bleach' false false}, ...
                'callback', {[] @expparams [] [] []}, ... % GOD this is complex!!! 
                'length',   {2 1 2 1 1}, ...
                'label',    {'#exp' '' '' '' ''}, ...
                'labellength',  {1 0 0 0 0} ...
                );
            specs(end+1) = struct('name','BLEACH','controls',s);
            function c = expparams(c)
                set(V.hf,'pointer','watch'), drawnow
                if V.nc>1
                    stimtype = {V.trial.stimdetails.type};
                    blankcond = strcmp(stimtype,'blank');
                    if ~any(blankcond)
                        errordlg 'no condition of type ''blank'', please edit stims'
                        set(V.hf,'pointer','arrow')
                        return
                    end
                    blankavg = mean(V.trial.data(:,:,:,blankcond),4);
                else
                    blankcond = 1;
                    blankavg = V.trial.data;
                end
                s = tps_expfit(blankavg,c{1});
                if ~c{4}, s = rmfield(s,'beta'); end % beta will be estimated independtly for each trial
                s.blankcond = blankcond;
                c{2} = fn_pointer(s); % store the pointer as a value of the 'EST' button
                set(V.hf,'pointer','arrow')
            end
            % Camera noise
            s = struct( ...
                'style',    {'popupmenu' 'edit' 'edit'}, ...
                'string',   {{'indices','first selection','last selection','out selections','user'} '' ''}, ...
                'default',  {'indices' ':' ''}, ...
                'length',   {3 2 3}, ...
                'label',    {'' '' 'high-pass:'}, ...
                'labellength',  {0 0 2} ...
                );
            specs(end+1) = struct('name','CAMERANOISE','controls',s);
            % Bi-directional delay correction
            s1 = s0;
            s1.style = 'edit';
            s1.default = '0';
            s1.label = 'delay';
            s1.labellength = .5;
            specs(end+1) = struct('name','BIDIR SCAN','controls',s1);
            % Heart beat
            s = struct( ...
                'style',    {'popupmenu' 'edit' 'stepper' 'stepper' 'stepper' 'popupmenu'}, ...
                'string',   {{'indices' 'first selection' 'last selection' 'out selection' 'user' 'recording'} '' '' '' '' {'rm' 'keep'}}, ...
                'default',  {'indices' ':' 7 1 .6 'rm'}, ...
                'more',     {[] [] [] [] {'min',0,'max',2,'step',.1} []}, ...
                'length',   {3 2 5 5 5 3}, ...
                'label',    {'' '' '~f' '#har' 'prec.' ''}, ...
                'labellength',  {0 0 2 2 2 0});
            specs(end+1) = struct('name','HEART','controls',s);
            % Motion correction
            s = struct( ...
                'style',    {'popupmenu' 'edit' 'stepper' 'stepper'}, ...
                'string',   {{'indices' 'first selection' 'last selection' 'out selection' 'user'} '' '' ''}, ...
                'default',  {'indices' ':' .5 10}, ...
                'more',     {[] [] {'min',.05,'max',.5,'step',.05} {'min',1,'max',Inf,'step',1}}, ...
                'length',   {5 2 4 3}, ...
                'label',    {'mask' '' 'maxshift' 'ref'}, ...
                'labellength',  {2 0 2 1});
            specs(end+1) = struct('name','MOTION','controls',s);            
            % Masking
            s = struct( ...
                'style',    {'popupmenu' 'edit' 'popupmenu'}, ...
                'string',   {{'threshold on avg. value' 'threshold below' 'indices' 'first selection' 'last selection' 'out selection' 'user'} '' {'NaN' '0' '1' 'min' 'max' 'mean'}}, ...
                'default',  {'threshold on avg. value' '1' 'NaN'}, ...
                'length',   {2 1 2}, ...
                'label',    {'' '' 'outside val.'}, ...
                'labellength',  {0 0 1} ...
                );
            specs(end+1) = struct('name','MASK','controls',s);
            % Filtering
            s1 = s0; 
            s1.style = 'edit';
            s1.length = 3;
            s1.default = '';
            s1.labellength = 2;
            s2 = s1; s3 = s1; s4 = s1;
            s1.label = 'x lp:';
            s2.label = 'x hp:';
            s3.label = 't lp:';
            s4.label = 't hp:';
            s5 = s0;
            s5.style = 'checkbox';
            s5.length = 2;
            s6 = s5;
            s5.string = 'd';
            s5.default = false;
            s6.string = 'z';
            s6.default = true;
            s7 = s5;
            s7.string = 'ph';
            s7.default = false;
            s8 = s5;
            s8.string = 's';
            s8.default = false;
            specs(end+1) = struct('name','FILTER','controls',[s1 s2 s3 s4 s5 s6 s7 s8]);
            % Special detrend
            s1 = s0;
            s1.style = 'popupmenu';
            s1.string = {'detrend & normalize' 'detrend only'};
            s1.default = 'detrend & normalize';
            s1.length = 7;
            s2 = s0;
            s2.style = 'edit';
            s2.default = ':';
            s2.length = 9;
            s2.label = 'based on frames';
            s2.labellength = 5;
            specs(end+1) = struct('name','DETREND','controls',[s1 s2]);
            % Binning & Box-car
            s1 = s0; s2 = s0;
            s1.style = 'edit';
            s1.length = 2;
            s1.default = '2';
            s1.label = 'space:';
            s1.labellength = .7;
            s2.style = 'edit';
            s2.length = 1;
            s2.default = '1';
            s2.label = 'time:';
            s2.labellength = .7;
            specs(end+1) = struct('name','BIN','controls',[s1 s2]);
            s1.default = '1';
            specs(end+1) = struct('name','BOX-CAR','controls',[s1 s2]);
            % FFT
            s1 = s0;
            s1.style = 'checkbox';
            s1.string = 'take norm';
            s1.length = 1;
            s1.default = false;
            specs(end+1) = struct('name','FFT','controls',s1);
            % User action
            s1 = s0;
            s1.style = 'edit';
            s1.length = 3;
            s1.default = 'x./(1+x)';
            s2 = s0;
            s2.style = 'pushbutton';
            s2.string = 'edit';
            s2.length = 1;
            s2.callback = @editfun;
            function c = editfun(c)
                str = c{1};
                if exist(str,'file')
                    edit(str)
                else
                    tokens = regexp(c{1},'^([^\(]*)','tokens');
                    if ~isempty(tokens) && exist(tokens{1}{1},'file')
                        edit(tokens{1}{1})
                    else
                        % prompt user for new file name
                        while true
                            ffun = fn_savefile('*.m','Select file name for writing new data operation function');
                            if isequal(ffun,0), return, end
                            [folder fname] = fileparts(ffun);
                            path_list_cell = regexp(path,pathsep,'Split');
                            if ~ismember(folder,path_list_cell)
                                answer = questdlg(['Folder "' folder '" is not on the path. Add it?'],'', ...
                                    'Add to path', 'Select another file name', 'Cancel', 'Add to path');
                                switch answer
                                    case 'Add to path'
                                        addpath(folder)
                                        break
                                    case 'Select another file name'
                                        continue
                                    case 'Cancel'
                                        disp 'interrupted'
                                        return
                                end
                            else
                                break
                            end
                        end
                        % create the file and edit it
                        template = which('dataop_user_template.txt');
                        txt = fn_readtext(template);
                        txt = strrep(txt,'WRITE_FUNCTION_NAME_HERE',fname);
                        fn_savetext(txt,ffun)
                        c{1} = fname; % write the function name in the text control
                        edit(ffun)
                    end
                end                
            end
            s3 = s0;
            s3.style = 'pushbutton';
            s3.string = 'recompute';
            s3.length = 1;
            s3.callback = @recomputefun;
            function c = recomputefun(c)
                % change the value to force re-computation
                c{3} = rand;
            end
            specs(end+1) = struct('name','USER','controls',[s1 s2 s3]);
        end
    end
    
    % Load (handle previous versions)
    methods (Static)
        function op = loadobj(op)
            op = convert(op);
        end
    end
    
end
    
