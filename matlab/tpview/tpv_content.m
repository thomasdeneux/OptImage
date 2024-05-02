classdef tpv_content < hgsetget
    % tpv_content always make sure that signal data is computed for current
    % trial and for all selections
    properties
        version
        trials
        signals
        ktrial = 1;
    end
    properties (SetAccess='private')
        electrophys % preprocessed electrophysiology recordings
    end
    properties (SetAccess='private')
        nsel = 0;       % number of selection
        ij = [1; 1];    % indicates which pixel to get signal from if there is no selection
    end
    properties (Dependent)
        seldotrial      % allow different selections for different trials
        timeline        % always view second space dimension as time (i.e. not only for the trials of type 'linescan')
    end
    properties (Dependent, SetAccess='private')
        ntrial
        nfr
        nx              % number of time courses: nx=nsel, except when nsel=0 (and nx=1)
        datamodes
        datacond
        trial
        signal
    end
    properties (Transient, GetAccess='private')
        internpar = tpv_internpar;  % internal parameters, shared with tpview and with children trials
    end
    properties (Dependent)
        stimtable       % table of different stimulations (a tps_stimtable object)
    end
    properties
        docfile
        user
    end
    
    % Constructor
    methods
        function C = tpv_content
            C.version = 1.41;
            C.trials = tps_trial(0);
            C.electrophys = tps_electrophys;
            C.signals = tps_signal;
            C.signals.datamode = 'data';
            C.signals.datacond = '';
            C.signals.x = reshape(C.signals.x,1,0);
            C.signals.dataopdef = struct('op',[],'link',1);
            C.internpar = tpv_internpar;
        end
    end
    
    % GET dependent
    methods
        function p = get.internpar(C)
            p = C.internpar;
        end
        function modes = get.datamodes(C)
            modes = regexp(C.signals(1).datamode,'[^_]*','match');
            if isscalar(modes), modes{2}=''; end
        end
        function cond = get.datacond(C)
            cond = C.signals(1).datacond;
        end
        function b = get.seldotrial(C)
            b = C.signals(1).seldotrial;
        end
        function b = get.timeline(C)
            b = C.signals(1).timeline;
        end
        function ntrial = get.ntrial(C)
            ntrial = length(C.trials);
        end
        function nx = get.nx(C)
            nx = size(C.signals(1).x,2);
        end
        function nfr = get.nfr(C)
            nfr = C.trials(C.ktrial).nfr;
        end
        function trial = get.trial(C)
            trial = C.trials(C.ktrial);
        end
        function signal = get.signal(C)
            signal = C.signals(1);
        end
        function x = getrecording(C,sub,ktrial)
            % get simultaneous recordings

            if nargin==1, sub=1:length(C.trial.recording); end

            if nargin>=3
                if ~isscalar(ktrial), error 'getting recordings from multiple trials not handled yet', end
                trial = C.trials(ktrial);
            else
                trial = C.trial;
            end
            
            % we display recordings only for non-calculated conditions 
            kcond = condnum(C.trial.nc,C.datacond);
            if isempty(kcond), x = tps_signalx.empty(1,0); return, end
            
            % recordings should be recalled by their name, not by indices
            if isnumeric(sub) || islogical(sub)
                error 'this syntax is not allowed any more'
            end

            if ischar(sub), sub = {sub}; end
            nrec = length(sub);
            if nrec==0
                x = tps_signalx.empty(1,0);
                return
            end
            
            % common stuff
            % (available recording names)
            recnames = {trial.recording.name};
            % (fake data of the appropriate size)
            if C.timeline || strcmp(trial.type,'linescan')
                data = zeros(trial.ny*trial.nfr,1);
            else
                data = zeros(trial.nfr,1);
            end
            % (other)
            donorm = C.signals(1).opdef.normalize__by__mean;
            
            % get recording signals
            for i=1:length(sub)
                recname = sub{i};
                z = tps_signalx;
                z.data = data;
                % check that recording is available
                switch recname
                    case {'data' 'sfr'}
                        % special case: show the image signals as a recording
                        if strcmp(trial.origin,'ScanImage')
                            disp('don''t know how to read the fill fraction for non ScanImage data')
                            recname = '';
                        end
                    otherwise
                        if strcmp(recname,'electrophysiology') && C.electrophys(C.ktrial).n
                            recname = 'E'; % use C.electrophys
                        else
                            idx = find(strcmpi(recname,recnames));
                            if isempty(idx) || length(trial.recording(idx).signal)<=1
                                recname = '';
                            end
                        end
                end
                % get recording
                switch recname
                    case ''
                        % no recording with this name
                        z.tidx = [0 trial.nfr*trial.dt];
                        z.dataop = [0; 0];
                        z.tag = 'empty';
                    case {'data' 'sfr'}
                        acq = trial.fullinfo.acq;
                        ff = acq.fillFraction;
                        lt = acq.msPerLine/1000;
                        pt = lt*ff/trial.nx;
                        del = acq.acqDelay; % TODO : BUG IN SCANIMAGE, 'del' IS NOT UPDATED!!!
                        bidi = acq.bidirectionalScan;
                        % time (add segments for the flybacks)
                        tline = [(0:trial.nx-1)*pt trial.nx*pt lt-pt] + del;
                        tidx = fn_add(tline',(0:trial.ny*trial.nfr-1)*lt);
                        z.tidx = tidx(:);
                        % the data (add segments for the flybacks)
                        dat = trial.(recname);
                        dat = reshape(dat,trial.nx,trial.ny*trial.nfr);
                        if bidi % re-invert lines if bidirectional scan!
                            dat(:,2:2:end) = dat(trial.nx:-1:1,2:2:end);
                        end
                        m = mean(dat(:));
                        dat(end+[1 2],:) = m;
                        dat = single(dat);
                        if donorm, dat = dat/m; end
                        z.dataop = -dat(:) * C.signals(1).opdef.channelscale;
                        % and the tag
                        z.tag = 'image2rec';
                    otherwise
                        if strcmp(recname,'E')
                            % use C.electrophys
                            r = C.electrophys(C.ktrial);
                            z.spikes = r.spikes;
                            z.tag = 'electrophysiology';
                        else
                            % general case: get recording from trial recordings
                            r = trial.recording(idx);
                            z.tag = r.name;
                        end
                        if r.n==0 || r.n==1
                            % there is actually no recording (e.g. file
                            % exists but has no data, etc.)
                            z.tidx = [0 trial.nfr*trial.dt];
                            z.dataop = zeros(2,1);
                        else
                            z.tidx = (0:r.n-1)*r.dt + r.t0;
                            z.dataop = fn_float(r.signal(:,kcond)) * C.signals(1).opdef.channelscale;
                        end
                end
                if isempty(C.datamodes{2})
                    z.data2op = [];
                else
                    z.data2op = z.dataop;
                end
                x(i) = z; 
            end
        end
        function S = get.stimtable(C)
            if C.ntrial
                S = C.trials(1).stimtable;
            else
                S = [];
            end
        end
    end
    
    % Set content - automatic updates
    methods
        % PARAMETERS
        function set.internpar(C,p)
            C.internpar = p;
            [C.trials.internpar] = deal(p); %#ok<MCSUP>
        end
        function set.seldotrial(C,b)
            if b==C.signals(1).seldotrial, return, end
            C.signals(1).seldotrial = b;
            if ~b
                % make all selections equal to that of the current trial
                ktrials = setdiff(1:C.ntrial,C.ktrial);
                signal = C.signals(1);
                sel = [signal.x(C.ktrial,:).sel];
                if isempty(sel), return, end
                for k = ktrials
                    sel = ComputeInd(sel,C.trials(k).sizes(1:2));
                    selcell = num2cell(sel);
                    [signal.x(k,:).sel] = deal(selcell{:});
                end
                signal.shift(:) = 0;
                erasedata(C)
            end
        end
        function set.timeline(C,b)
            if b==C.timeline, return, end
            C.signals(1).timeline = b;
            trialtypes = {C.trials.type};
            ktrials = find(~strcmp(trialtypes,'linescan'));
            erasedata(C,ktrials)
            X = C.signals(1).x;
            ndnew = 2-b;
            if fn_dodebug, disp('tpv_content.set.timeline: improve code to avoid repetition of the same calculation in different trials'), end
            for k=ktrials
                sizk = C.trials(k).sizes(1:ndnew);
                for i=1:C.nsel, X(k,i).sel = convertdim(X(k,i).sel,ndnew,sizk); end
            end
        end
        % DATA
        function settrials(C,T,datamode,datacond,dataopdef,selection)
            C.trials = T;
            ntr = length(T);
            C.ktrial = min(C.ktrial,ntr);
            % initialize the signals (empty selections)
            C.signals = tps_signal;
            C.signals.shift = zeros(ntr,2);
            C.signals.x = tps_signalx.empty(ntr,0);
            % electrophysiology preprocessings
            C.electrophys = tps_electrophys(T);
            % internal parameters
            [T.internpar] = deal(C.internpar);
            % user
            C.user = struct;
            % initialize more?
            if nargin>=3, C.signals.datamode = datamode; end
            if nargin>=4, C.signals.datacond = datacond; end
            if nargin>=5
                if ~iscell(dataopdef), dataopdef = {dataopdef}; end
                [C.trials.opdef] = deal(dataopdef{:}); 
            end
            C.signals(1).dataopdef = struct('op',cell(1,ntr),'link',true); % avoid duplicating the operation definitions
            if nargin>=6, setselection(C,selection), else C.nsel = 0; end
        end
        function addtrials(C,T)
            if isempty(T), return, end
            ntrold = length(C.trials);
            C.trials = [C.trials T];
            C.electrophys = [C.electrophys tps_electrophys(T)];
            ntr = length(C.trials);
            % internal parameters
            [T.internpar] = deal(C.internpar);
            % saved selections are desuete
            C.signals(2:end) = [];
            % merge the stimulation tables
            mergestimtables(C.trials)
            % update the shift
            if ntrold==0, error programming, end
            C.signals.shift(ntrold+1:ntr,:)    = repmat(C.signals.shift(ntrold,:),ntr-ntrold,1);
            % update the signals - don't forget the bug (better to
            % initialize each array element individually)
            C.signals.x(ntr,end) = tps_signalx;
            X = C.signals.x;
            X(ntr,end) = tps_signalx;
            for i=ntrold+1:ntr
                for j=1:C.nx
                    X(i,j)=tps_signalx; 
                end
            end
            C.signals.x = X;
            % copy the selection of the previous last trial (handle change
            % in trial size and/or scale)
            sel = [C.signals.x(ntrold,:).sel];
            ndold = sel(1).nd;
            sizeold = C.trials(ntrold).sizes(1:ndold);
            scalold = C.trials(ntrold).dx; if ndold==2, scalold(2) = C.trials(ntrold).dy; end
            for k=ntrold+1:ntr
                ndnew = 2 - (C.timeline || strcmp(C.trials(k).type,'linescan'));
                sizenew = C.trials(k).sizes(1:ndnew);
                scalnew = C.trials(k).dx; if ndnew==2, scalnew(2) = C.trials(k).dy; end
                if ndnew~=ndold
                    sel = convertdim(sel,ndnew,sizenew);
                    ndold = ndnew;
                    scalold = scalnew;
                    sizeold = sizenew;
                elseif any(scalnew~=scalold)
                    aff = affinityND(fn_switch(ndold,1,'scale1D',2,'scale2D'),scalold./scalnew);
                    sel = selaffinity(sel,aff,sizenew);
                    scalold = scalnew;
                    sizeold = sizenew;
                elseif any(sizenew~=sizeold)
                    sel = ComputeInd(sel,sizenew);
                    sizeold = sizenew;
                end
                for i=1:C.nx, C.signals.x(k,i).sel = sel(i); end 
            end
            % update operation definition
            klinked = find([C.signals(1).dataopdef.link]==1,1,'first');
            if isempty(klinked), klinked = ntrold; end
            [C.signals(1).dataopdef(ntrold+1:ntr)] = deal(C.signals(1).dataopdef(klinked));            
            [C.trials(ntrold+1:ntr).opdef] = deal(C.trials(klinked).opdef); 
        end
        function rmtrial(C,k)
            % function rmtrial(C,ktrials)
            if nargin<2, k = C.ktrial; end
            if length(k)>=C.ntrial, errordlg('cannot remove specified trial(s)'), return, end
            C.trials(k) = [];
            C.electrophys(k) = [];
            for i=1:length(C.signals)
                C.signals(i).x(k,:) = [];
                C.signals(i).shift(k,:) = [];
                C.signals(i).dataopdef(k) = [];
            end
            % decrease ktrial if necessary
            C.ktrial = min(C.ktrial,C.ntrial);
        end
        function permutetrials(C,perm)
            if ~isempty(setxor(perm,1:C.ntrial)), error 'this is not a permutation', end
            C.trials = C.trials(perm);
            C.electrophys = C.electrophys(perm);
            if C.seldotrial
                disp('making all selections equal') % too difficult otherwise
                C.seldotrial = false;
                C.seldotrial = true;
            end
            for k=1:length(C.signals)
                C.signals.x = C.signals.x(perm,:);
                if any(C.signals(k).shift), error strange, end
            end
        end
        % ELECTROPHYSIOLOGY
        function setelectrophys(C,E)
            C.electrophys = E;
        end
        % DATA MODE
        function setdatamode(C,mode)
            if strcmp(mode,C.signals(1).datamode), return, end
            flagsold = C.datamodes;
            C.signals(1).datamode = mode;
            flagsnew = C.datamodes;
            if ~strcmp(flagsnew{1},flagsold{1})
                [C.signals(1).x.data] = deal([]); 
                [C.signals(1).x.dataop] = deal([]); 
                [C.signals(1).x.validspike] = deal(false);
            end
            if ~strcmp(flagsnew{2},flagsold{2})
                [C.signals(1).x.data2] = deal([]); 
                [C.signals(1).x.data2op] = deal([]); 
                [C.signals(1).x.validspike2] = deal(false);
            end
        end
        function setdatacond(C,cond)
            if strcmp(cond,C.signals(1).datacond), return, end
            C.signals(1).datacond = cond;
            erasedata(C)
        end
        function setdataopdef(C,op,link)
            % function setdataopdef(C,op,link)
            % function setdataopdef(C,'link',link)

            % which trials to change: depends on link
            if ischar(op) && strcmp(op,'link')
                % update link only
                C.signals(1).dataopdef(C.ktrial).link = link;
                if ~link, return, end
                % when linking, get operation definition from other linked
                % trial
                linkedtrials = find([C.signals(1).dataopdef.link]);
                linkedtrials = setdiff(linkedtrials,C.ktrial);
                if isempty(linkedtrials), return, end
                ktrials = C.ktrial;
                op = C.trials(linkedtrials(1)).opdef;
            elseif link
                ktrials = find([C.signals(1).dataopdef.link]);
            else
                ktrials = C.ktrial;
            end
            % update opdef in trial(s)
            [C.trials(ktrials).opdef] = deal(op);
            % update the 'link' flag
            C.signals(1).dataopdef(C.ktrial).link = link;
            % erase signals that became obsolete
            flags = C.datamodes;
            if strfind(flags{1},'op')
                [C.signals(1).x(ktrials,:).data] = deal([]);
                [C.signals(1).x(ktrials,:).dataop] = deal([]);
                [C.signals(1).x.validspike] = deal(false);
            end
            if strfind(flags{2},'op')
                [C.signals(1).x(ktrials,:).data2] = deal([]);
                [C.signals(1).x(ktrials,:).data2op] = deal([]);
                [C.signals(1).x.validspike2] = deal(false);
            end
            % erase spike movies that became obsolete
            for k=ktrials
                C.trials(k).usertransient.dataspikebin = [];
                C.trials(k).usertransient.shotnoisespikebin = [];
            end
        end
        function setsignalopdef(C,op)
            if isequal(op,C.signals(1).opdef), return, end
            C.signals(1).opdef = op;
            if C.signals(1).nx==0, return, end
            [C.signals(1).x.dataop]  = deal([]);
            [C.signals(1).x.data2op] = deal([]);
            [C.signals(1).x.validspike] = deal(false);
            [C.signals(1).x.validspike2] = deal(false);
        end
        % SELECTION
        function setselection(C,selection)
            % function setselection(C,selection)
            %---
            % note that the selectionset object selection can be either
            % scalar (same selection is used for all trials) or a vector of
            % length the number of trials
            switch length(selection(1).t)
                case 0
                    C.nsel = 0;
                    C.signals(1).x = tps_signalx.empty(C.ntrial,0);
                    updateselection(C,'ij')
                case 1
                    updateselection(C,'all',[],selection(1).t.set)
                otherwise
                    error 'input selection has multiple sets'
            end
        end
        function updateselection(C,flag,ind,value)
            % this functions updates the selection for all trials, and
            % interpolate the data and apply operation for the current
            % trial
            signal = C.signals(1);
            scanline = C.timeline || strcmp(C.trials(C.ktrial).type,'linescan');
            seldimsnum = fn_switch(scanline,1,[1 2]);
            if C.timeline
                linetrials = true(1,C.ntrial);
            else
                trialtypes = {C.trials.type};
                linetrials = strcmp(trialtypes,'linescan');
            end
            switch flag
                case 'ij'
                    if nargin==4
                        C.ij = value;
                    end
                    if C.nsel==0
                        % get time course at the current pixel
                        % (reinitialize)
                        % [use of a 'tmp' variable saves time in the case
                        % that C.ntrial is large, because tmp = zzz is
                        % faster than signal.x = zzz]
                        tmp = tps_signalx;
                        % [apparently not necessary to do the following,
                        % but for security...] 
                        for i=C.ntrial:-1:2, tmp(i,1) = tps_signalx; end
                        % (set a point selection)
                        ij = C.ij; if isscalar(ij), ij(2)=1; end
                        sel1 = selectionND('point1D',ij(1),C.trial.sizes(1));
                        sel2 = selectionND('point2D',ij(1:2),C.trial.sizes(1:2));
                        for k=1:C.ntrial
                            sizk = C.trials(k).sizes;
                            if C.timeline || strcmp(C.trials(k).type,'linescan')
                                if sizk(1)~=sel1.datasizes, sel1 = ComputeInd(sel1,sizk(1)); end
                                tmp(k).sel = sel1;
                            else
                                if any(sizk(1:2)~=sel2.datasizes), sel2 = ComputeInd(sel2,sizk(1:2)); end
                                tmp(k).sel = sel2;
                            end
                        end
                        signal.x = tmp;
                    end
                case {'all' 'new'}
                    if strcmp(flag,'all')
                        C.nsel=0; 
                        if isa(value,'selectionset')
                            % this can happen at init, when tpview object
                            % inherits its selections from an already
                            % existing focus object
                            value = value.singleset;
                        end
                        if isempty(value)
                            C.signals(1).x = tps_signalx.empty(C.ntrial,0); 
                            updateselection(C,'ij')
                            return
                        end
                    end
                    for ksel=1:length(value)
                        valuek = value(ksel);
                        % MATLAB BUG: if we create new entries in signal.x,
                        % some of them will point to the same object!
                        % We should initialize each new object by hand...
                        % AND MOREOVER, we need to do it smartly otherwise it
                        % is very slow...
                        x = tps_signalx;
                        % apparently not necessary to do all this, but for
                        % security...
                        for i=C.ntrial:-1:2, x(i,1) = tps_signalx; end
                        
                        % set selection
                        % TODO: it is not efficient to compute indices for
                        % all trials at this stage, especially if there are
                        % many trials with different sizes; indices should
                        % be computed on the fly when needed
                        ktr = C.ktrial;
                        switch valuek.nd
                            case 1
                                % original selection is 1D 
                                % direct set for line scan trials
                                valuemem = valuek;
                                for k=find(linetrials)
                                    sizk = C.trials(k).sizes(1);                                    
                                    if sizk~=valuemem.datasizes, valuemem = ComputeInd(valuemem,sizk); end
                                    x(k).sel = valuemem;
                                end
                                % convert to 2D for image trials
                                if ~all(linetrials)
                                    valuemem = convertdim(valuek,2,C.trial.sizes([1 2]));
                                    for k=find(~linetrials)
                                        sizk = C.trials(k).sizes([1 2]);
                                        if sizk~=valuemem.datasizes, valuemem = ComputeInd(valuemem,sizk); end
                                        x(k).sel = valuemem;
                                    end
                                end
                            case 2
                                % original selection is 2D
                                % convert to 1D for line scan trials
                                if any(linetrials)
                                    valuemem = convertdim(valuek,1,C.trial.sizes(1));
                                    for k=find(linetrials)
                                        sizk = C.trials(k).sizes(1);
                                        if sizk~=valuemem.datasizes, valuemem = ComputeInd(valuemem,sizk); end
                                        x(k).sel = valuemem;
                                    end
                                end
                                % take into account the global shift to set
                                % the new selection for image trials
                                selshift = C.signals(1).shift;
                                shiftktr = selshift(ktr,:);
                                shiftmem = shiftktr;
                                valuemem = valuek;
                                for k=find(~linetrials)
                                    % this way avoids redundant computations
                                    shiftk = selshift(k,:);
                                    if any(shiftk~=shiftmem)
                                        mov = affinityND('translate2D',shiftk-shiftktr);
                                        shiftmem = shiftk;
                                        valuemem = selaffinity(valuek,mov);
                                    end
                                    sizk = C.trials(k).sizes(seldimsnum);
                                    if any(sizk~=valuemem.datasizes), valuemem = ComputeInd(valuemem,sizk); end
                                    x(k).sel = valuemem;
                                end
                        end
                    
                        % add selection to tps_signal object
                        if C.nsel==0
                            signal.x = x;
                        else
                            signal.x = [signal.x x];
                        end
                        C.nsel = size(signal.x,2);
                    end
                case 'change'
                    if C.seldotrial
                        % selections are not updated for all trials, but
                        % only for those after C.ktrial that already have
                        % selections identical to C.ktrial selection
                        ntr = C.ntrial;
                        ktr = C.ktrial;
                        selshift = signal.shift;
                        nind = length(ind);
                        sel = [signal.x(ktr:ntr,ind).sel];
                        id = cat(1,sel.id);
                        
                        for j=1:nind
                            i = ind(j);
                            % for which trials the selection is the same
                            % (use id(2))
                            nsame = find(diff(id(ktr:end,2)),1,'first');
                            if isempty(nsame), ksame = ktr:C.ntrial; else ksame = ktr+(0:nsame-1); end
                            % check that shift value are also the same
                            if ~isscalar(ksame) && any(any(diff(selshift(ksame,:),1)))
                                disp('problem with sel change/shift (error programming)')
                            end
                            % set selection
                            for k=ksame
                                selk = value(j);
                                sizk = C.trials(k).sizes(seldimsnum);
                                if any(sizk~=selk.datasizes), selk = ComputeInd(selk,sizk); end
                                [signal.x(k,i).sel] = selk;
                            end
                            % erase data that has changed
                            erasedata(C,ksame,i)
                        end
                    else
                        for j=1:length(ind)
                            i = ind(j);
                            for k=1:C.ntrial
                                selk = value(j);
                                sizk = C.trials(k).sizes(seldimsnum);
                                if any(sizk~=selk.datasizes), selk = ComputeInd(selk,sizk); end
                                [signal.x(k,i).sel] = selk;
                            end
                        end
                        % erase data that has changed
                        erasedata(C,1:C.ntrial,ind)
                    end
                case 'active'
                    [signal.x(:,ind).active] = deal(value);
                    sel = [signal.x(:,ind).sel];
                    [sel.active] = deal(value);
                case 'reorder'
                    perm = value;
                    signal.x = signal.x(:,perm);
                case 'remove'
                    if xor(C.nsel==0,isempty(ind)), error programming, end
                    signal.x(:,ind) = [];
                    %                     if C.nsel==0 % removed all regions
                    %                         signal.x = reshape(signal.x,C.ntrial,0); % ?
                    %                     end
                    C.nsel = size(signal.x,2);
                    if C.nsel==0
                        signal.shift(:) = 0; % reset the global shift to zero
                        updateselection(C,'ij')
                    end
                case 'reset'
                    signal.x(:) = [];
                    signal.x = reshape(signal.x,C.ntrial,0); 
                    C.nsel = 0;
                    updateselection(C,'ij')
                case 'spikes'
                    signal.x(C.ktrial,ind).spikes = value;
                case 'indices'
                    % happens when creating external objects
                otherwise
                    % there should be no add, affinity
                    error programming
            end            
        end
        function updateselectionshift(C,mov,value)
            if ~C.seldotrial, return, end
            % which trials should be updated (current trial and next ones
            % having the same shift value)
            signal = C.signals(1);
            shiftk = signal.shift(C.ktrial,:);
            nsame = find(diff(signal.shift(C.ktrial:end,:),1),1,'first');
            if isempty(nsame)
                ksame = C.ktrial:C.ntrial; 
            else
                ksame = C.ktrial+(0:nsame-1); 
            end
            % update shift values
            shiftk = shiftk + mov.mat(2:3,1)';
            signal.shift(ksame,:) = repmat(shiftk,length(ksame),1);
            % update regions position - note that this cancels local
            % changes in the regions (i.e. for a given i, all
            % x(ksame,i).sel become the same, which might not be the case
            % before)
            for i=1:C.nsel
                [signal.x(ksame,i).sel] = deal(value(i));
            end
            % erase data that has changed
            erasedata(C,ksame)
        end
        function updateselectionhandy(C,ktrials,flag,varargin)
            % check
            if ~C.seldotrial && ~isempty(ktrials) && ~isequal(ktrials,1:C.ntrial)
                error 'selection change in specific trials contradicts identical selections mode'
            end
            % update selections
            switch flag
                case 'xbinchange'
                    scalefactor = 1/varargin{1};
                    A{1} = [1 0; .5*(1-scalefactor) scalefactor];
                    A{2} = [1 0 0; .5*(1-scalefactor) scalefactor 0; .5*(1-scalefactor) 0 scalefactor];
                    for i=1:length(C.signals)
                        for ktr = ktrials
                            X = C.signals(i).x(ktr,:);
                            for ksel = 1:length(X)
                                sel = X(ksel).sel;
                                sel = selaffinity(sel,A{sel.nd});
                                ComputeInd(sel,C.trials(ktr).sizes(1:sel.nd));
                                X(ksel).sel = sel;
                            end
                        end
                    end
                otherwise
                    error('unknown flag ''%s''',flag)
            end
        end
        % SIGNALS STORAGE
        function storesignals(C,idx,selname)
            % check consistency of idx and name
            if idx>length(C.signals)+1 || idx<2 ...
                    || (idx<=length(C.signals) && ~strcmp(C.signals(idx).name,selname))
                error 'idx or name is not valid for saving signals'
            end
            % create new tps_signal object
            x = copy(C.signals(1));
            % name
            x.name = selname;
            % copy the data operation definitions
            [x.dataopdef.op] = deal(C.trials.opdef);
            % save
            C.signals(idx) = x;
        end
        function loadsignals(C,idx)
            if idx==1, errordlg 'this action has no effect!', return, end
            % load the stored selection
            C.signals(1) = copy(C.signals(idx));
            C.nsel = size(C.signals(1).x,2); % TODO: it might happen that a 'pointer' signal becomes a 'selection' signal...
            % update op. def. in trials, then clear the unwanted duplicate
            % information
            [C.trials.opdef] = deal(C.signals(1).dataopdef.op);
            [C.signals(1).dataopdef.op] = deal([]);
        end
    end
    
    % Data computations
    methods
        function clearmemory(C)
            memorypool.clear()
            if isa(C.trials,'tps_trial'), erasedata(C.trials), end
            for k=1:length(C.signals)
                xx = C.signals(k).x;
                if isempty(xx), break, end
                [xx.delay xx.data xx.data2 xx.dataop xx.data2op] = deal([]);
                [xx.validspike xx.validspike2] = deal(false);
            end
        end
        function erasedata(C,ktrials,ind)
            if nargin<2 || isempty(ktrials), ktrials = 1:C.ntrial; end
            if nargin<3, ind=1:C.nx; end
            if isempty(ind), return, end
            signal = C.signals(1);
            xx = signal.x(ktrials,ind);
            [xx.delay xx.data xx.data2 xx.dataop xx.data2op] = deal([]);
            [xx.validspike xx.validspike2] = deal(false);
            [xx.tag] = deal([]);
        end
        function erasespikes(C,ktrials,ind,doall)
            if nargin<2 || isempty(ktrials), ktrials = 1:C.ntrial; end
            if nargin<3, ind=1:C.nx; end
            if nargin<4, doall=false; end
            if isempty(ind), return, end
            % 'signal' spikes
            signal = C.signals(1);
            xx = signal.x(ktrials,ind);
            if ~doall
                % erase only invalid spikes
                idx = ~[xx.validspike];  if any(idx), [xx(idx).spikes xx(idx).spikefit] = deal([]); end
                idx = ~[xx.validspike2]; if any(idx), [xx(idx).spikes2 xx(idx).spikefit2] = deal([]); end %#ok<NASGU>
            else
                [xx.validspike xx.validspike2] = deal(false);
                [xx.spikes xx.spikefit xx.spikes2 xx.spikefit2] = deal([]);
            end
            % 'movie' spikes
            if doall
                for ktr=ktrials
                    Tk.usertransient.dataspikebin = [];
                    Tk.usertransient.shotnoisespikebin = [];
                end
            end
        end
        function computesignals(C)
            fn_progress('compute signals',C.ntrial)
            for k=1:C.ntrial
                fn_progress(k)
                interpdata(C,k)
                operation(C,k)
            end
        end
        function computespikes(C,ktrials,ind,spikepar,docompute)
            % function computespikes(C,ktrials,ind,spikepar,docompute)
            %---
            % Input:
            % - par     structure - parameter set
            %           or flag: can be 'nospike', 'calcium',
            %           'calciumdrift', 'calciumcustom', 'vsd', 'vsdcustom'
            % - force   logical - do recompute spikes that have already
            %           been estimated?
            
            % input
            if nargin<2 || isequal(ktrials,0), ktrials = 1:C.ntrial; end
            if nargin<3, ind = 1:C.nx; end
            if nargin<4 || isempty(spikepar), spikepar = []; end
            if nargin<5, docompute = true; end
            two = any(C.signals(1).datamode=='_');
            
            % check that default parameters are not empty
            if isempty(C.signals(1).spikepar)
                C.signals(1).spikepar = fn_structmerge(struct('fun',@tps_mlspikes),tps_mlspikes('par'));
            end
            
            % set parameters
            defaultpar = C.signals(1).spikepar;
            if ~isempty(spikepar)
                % build parameters structure
                if isstruct(spikepar)
                    % check that the syntax is correct
                    ok = isfield(spikepar,'fun') && fn_ismemberstr(char(spikepar.fun),{'tps_mlspikes' 'tps_vsdspikes'});
                    if ~ok, error 'incorrect spike parameters structure', end
                    spiketype = fn_switch(char(spikepar.fun),'tps_mlspikes','calcium','tps_vsdspikes','vsd');
                else
                    flag = spikepar;
                    token = regexp(flag,'(calcium|vsd|nospike)','tokens');
                    if isempty(token), error 'invalid spike parameter', end
                    spiketype = token{1}{1};
                    switch spiketype
                        case 'nospike'
                            spikepar = struct('fun',[]);
                        case 'calcium'
                            spikepar = fn_structmerge(struct('fun',@tps_mlspikes),tps_mlspikes('par'));
                            switch strrep(flag,'calcium','')
                                case 'default'
                                    % nothing to do
                                case 'drift'
                                    spikepar.dodetrend = true;
                                case 'custom'
                                    if ~isempty(defaultpar) && isequal(defaultpar.fun,@tps_mlspikes)
                                        spikepar = fn_structmerge(spikepar,defaultpar,'skip');
                                    end
                                    spikepar = fn_structedit(spikepar);
                                    if isempty(spikepar), return, end
                                otherwise
                                    error 'invalid spike parameter'
                            end
                        case 'vsd'
                            spikepar = fn_structmerge(struct('fun',@tps_vsdspikes),tps_vsdspikes('par'));
                            switch strrep(flag,'vsd','')
                                case 'default'
                                    % nothing to do
                                case 'custom'
                                    if ~isempty(defaultpar) && isequal(defaultpar.fun,@tps_vsdspikes)
                                        spikepar = fn_structmerge(spikepar,defaultpar,'skip');
                                    end
                                    spikepar = fn_structedit(spikepar);
                                    if isempty(spikepar), return, end
                                otherwise
                                    error 'invalid spike parameter'
                            end
                    end
                end
                
                % set parameters for each signal
                X = C.signals(1).x;
                for k = ktrials
                    for i = ind
                        if ~isequal(X(k,i).spikepar,spikepar)
                            X(k,i).spikepar = spikepar;
                            X(k,i).validspike = false;
                            X(k,i).validspike2 = false;
                        end
                    end
                end
                
                % set the parameters default value
                if ~strcmp(spiketype,'nospike')
                    C.signals(1).spikepar = spikepar;
                    if ~isequal(defaultpar,spikepar)
                        defaultpar = spikepar;
                        for k = 1:C.ntrial
                            C.trials(k).usertransient.dataspikebin = [];
                            C.trials(k).usertransient.shotnoisespikebin = [];
                        end
                    end
                end
            end
            
            % compute spikes
            if docompute
                if isempty(ind), return, end
                interpdata(C,ktrials,ind)
                operation(C,ktrials,ind)
                xx = C.signals(1).x;
                if all([xx(ktrials,ind).validspike]), return, end
                % loop on trials and cells
                fn_progress('compute spikes',length(ktrials)*length(ind))
                kk = 0;
                for ktr = ktrials
                    for i = ind
                        kk = kk+1;
                        x = xx(ktr,i);
                        if x.validspike, continue, end
                        if isempty(x.spikepar), x.spikepar = defaultpar; end
                        par = x.spikepar;
                        if isempty(par.fun)
                            % no spikes
                            x.spikes = []; x.spikefit = []; x.validspike = true;
                            if two, x.spikes2 = []; x.spikefit2 = []; x.validspike2 = true; end
                            continue
                        end
                        fn_progress(kk)
                        par.dt = Tk.dt_sec;
                        if isempty(par.dt), errordlg 'please define frame rate to estimate spikes', end
                        % mark the times of spiking with the correct
                        % multiplicity
                        [nn fit par dum dum drift] = par.fun(x.dataop,rmfield(par,'fun')); %#ok<ASGLU>
                        if par.drift.parameter, fit = [fit drift]; end %#ok<AGROW>
                        x.spikes = fn_timevector(nn,par.dt,'times') + x.tidx(1);
                        x.spikefit = fit;
                        x.validspike = true;
                        if two
                            [nn2 fit2 par2 dum dum drift2] = par.fun(x.data2op,par); %#ok<ASGLU>
                            if par.drift.parameter, fit2 = [fit2 par2.F0+drift2]; end %#ok<AGROW>
                            x.spikes2 = fn_timevector(nn2,par.dt,'times') + x.tidx(1);
                            x.spikefit2 = fit2;
                            x.validspike2 = true;
                        end
                    end
                end
            end
           
            
        end
        function x = spikemovie(C,ktrial,kcond,spikedataflag,combdata)
            % input
            if nargin<5, combdata = []; end
            % check spike parameters
            par = C.signals(1).spikepar;
            T = C.trials(ktrial);
            deltat = T.dt_sec;
            if isempty(par) || isempty(deltat)
                msg = ['cannot estimate spikes: ' fn_switch(isempty(par), ...
                    'set spike parameters first','frame rate is not defined')]; 
                errordlg(msg)
                if isempty(combdata), x = zeros(T.sizes(1:3)); else x = combdata; end
                return
            end
            if isempty(par.fun), return, end
            par.dt = deltat;
            % compute
            spikeflag = [spikedataflag 'spikebin'];
            if isfield(T.usertransient,spikeflag) && ~isempty(T.usertransient.(spikeflag))
                xbin = T.usertransient.(spikeflag);
            else
                xbin = operation(T,spikedataflag,[],false); % get dataop binned (i.e. before it is re-enlarged to the size of data)
                [nxb nyb nfr nc] = size(xbin);
                if T.spacetimebin(2)~=1, error 'cannot handle temporal binning', end
                [xbin] = par.fun(fn_reshapepermute(xbin,{3 [1 2 4]}),par);
                xbin = fn_reshapepermute(xbin,[nfr nxb nyb nc],[2 3 1 4]);
                T.usertransient.(spikeflag) = xbin;
            end
            % combine with data
            bin = T.spacetimebin(1);
            if isempty(combdata)
                % (show the mask)
                xbin = single(xbin);
                xbin(isnan(T.dataopbin)) = NaN;
                x = fn_enlarge(xbin(:,:,:,kcond),T.sizes(1:3));
            else
                % (combine with data)
                x = fn_enlarge(xbin(:,:,:,kcond),T.sizes(1:3));
                nerode = fn_switch(bin<2,0,bin<6,1,2);
                if nerode
                    x = x(:,:);
                    x = x & ~bwmorph(x,'erode',nerode);
                    x = reshape(x,T.sizes(1:3));
                end
                m = min(combdata(:)); M = max(combdata(:));
                combdata(x) = m + (M-m)*1.1;
                x = combdata;
            end
        end
    end
    methods (Access='private')
        function interpdata(C,ktrials,ind)
            if nargin<2 || isequal(ktrials,0), ktrials = 1:C.ntrial; end
            if nargin<3, ind=1:C.nx; end
            flags = C.datamodes;
            two = ~isempty(flags{2});
            stimids = unique([C.trials.stimid]);
            % loop on trials
            xx = C.signals(1).x;
            for ktr=ktrials
                % check which data does not exist yet
                hasnodata = false(1,length(ind));
                for i=1:length(ind)
                    xi = xx(ktr,ind(i));
                    hasnodata(i) = strcmp(xi.tag,'missing data') ...
                        || isempty(xi.data) || (two && isempty(xi.data2));
                end
                if ~any(hasnodata), continue, end
                
                Tk = C.trials(ktr);
                % data flag
                scanline = C.timeline || strcmp(Tk.type,'linescan');
                % data
                dat = getcondition(Tk,flags{1},C.signals(1).datacond,false);
                ok1 = ~isempty(dat);
                % secondary data
                if two
                    dat2 = getcondition(Tk,flags{2},C.signals(1).datacond,false);
                end
                ok2 = two && ~isempty(dat2);
                % special case: only one time instant -> create a second one to
                % get a correct display (horizontal lines)
                if ~scanline && size(dat,3)==1
                    if ok1, dat  = repmat(dat,[1 1 2]); end
                    if ok2, dat2 = repmat(dat2,[1 1 2]); end
                end
                % is the data a binned version?
                [nx ny nfr ncond] = size(dat);
                if ~ok1 && ok2, [nx ny nfr ncond] = size(dat2); end
                sizeorig = Tk.sizes;
                xbin = floor(sizeorig(1)/nx);
                ybin = floor(sizeorig(2)/ny);
                tbin = floor(sizeorig(3)/nfr);
                % reshape, delay
                if scanline
                    np = nx; nporig = sizeorig(1);
                    nfr = ny*nfr; nfrorig = prod(sizeorig(2:3));
                    delay = zeros(nx,1);
                else                    
                    np = nx*ny; nporig = prod(sizeorig(1:2));
                    nfrorig = max(2,sizeorig(3)); % at least 2 time points to have a correct display of single-frame data
                    linedur_tunit = Tk.linedur;
                    linedur_s = linedur_tunit * ...
                        fn_switch(Tk.tunit,'s',1,'ms',1e-3,{'min' 'minute'},60, ...
                        {'h' 'hour'},3600,'day',3600*24,C.internpar.defaultdt);
                    delay = repmat((0:ny-1)*linedur_s,nx,1);
                end
                if ok1, dat = reshape(dat,np,nfr,ncond); end
                if ok2, dat2 = reshape(dat2,np,nfr,ncond); end
                % condition number
                if Tk.nc==1, kcond = find(Tk.stimid==stimids); else kcond=0; end
                % set signals delay / data / data2
                for i=ind(hasnodata)
                    x = C.signals(1).x(ktr,i); % handle object
                    x.kcond = kcond;
                    if x.sel.nd~=2-scanline, error programming, end
                    dataind = x.sel.dataind;
                    dataind(dataind<1 | dataind>nporig)=[];
                    if (xbin>1 || ybin>1) && (ok1 || ok2)
                        if scanline
                            % the easy case
                            dataind = 1 + floor((dataind-1)/xbin);
                        else
                            % the difficult case
                            [ii jj] = ind2sub(sizeorig(1:2),dataind);
                            if ~all(ii), error programming, end % indices out of bound have already been removed
                            ii = 1 + floor((ii-1)/xbin);
                            jj = 1 + floor((jj-1)/ybin);
                            dataind = sub2ind([nx ny],ii,jj);
                        end
                        dataind(~dataind) = []; % can happen when original sizes are not multiple of xbin/ybin
                    end
                    if ok1
                        subdat = dat(dataind,:,:);
                        bad = any(any(isnan(subdat),2),3);
                        if full(any(bad)), subdat(bad,:,:) = []; end
                        x.data = fn_enlarge(shiftdim(mean(double(subdat),1),1),nfrorig);
                        x.dataop = [];
                    else
                        x.data = zeros(nfrorig,1);
                        x.tag = 'missing data';
                    end
                    if ok2
                        subdat = dat2(dataind,:,:);
                        bad = all(all(isnan(subdat),2),3);
                        if any(bad), subdat(bad,:,:) = []; end
                        x.data2 = fn_enlarge(shiftdim(mean(double(subdat),1),1),nfrorig);
                        x.data2op = [];
                    elseif two
                        x.data2 = zeros(nfrorig,1);
                        x.tag = 'missing data';
                    end
                    if ok1 && (~two || ok2), x.tag = []; end
                    if ok1 || ok2
                        x.delay = mean(delay(dataind));
                    else
                        x.delay = 0;
                    end
                end
            end
        end
        function operation(C,ktrials,ind)
            if nargin<2 || isequal(ktrials,0), ktrials = 1:C.ntrial; end
            if nargin<3, ind = 1:C.nx; end
            if isempty(ind), return, end
            signal = C.signals(1);
            xx = signal.x;
            two = any(signal.datamode=='_');
            op = signal.opdef;
                
            % special - operation on shot noise: for some operations, the
            % parameters should be estimated from the data, not from the
            % fake shot noise data; in such case, a memory of data
            % parameters must be kept; if data is not available, operation
            % cannot be performed
            % code: 0 = no shotnoise; 1 = data - shotnoise; 2 = shotnoise
            %       and data not available
            if ~isempty(strfind(signal.datamode,'shotnoise'))
                shotnoiseflag = fn_switch(signal.datamode, ...
                    {'data_shotnoise','dataop_shotnoiseop'}, 1,2);
            else
                shotnoiseflag = 0;
            end
            
            % loop on trials
            for ktr = ktrials
                    
                % if all data already exist, we assume there is no need to
                % change; however, if one does not exist, we recompute all
                % of them (because of possible neuropil subtraction)
                nodataop = false(1,length(ind));
                for i=1:length(ind)
                    nodataop(i) = (isempty(xx(ktr,ind(i)).dataop) || (two && isempty(xx(ktr,ind(i)).data2op)));
                end
                if ~any(nodataop), continue, end
                
                Tk = C.trials(ktr);
                % time constants
                trial = Tk;           
                deltat = trial.dt_sec;
                if isempty(deltat)
                    deltat = C.internpar.defaultdt;
                    t0 = 0;
                else
                    t0 = trial.t0 * (deltat/trial.dt);
                end
                if C.timeline || strcmp(trial.type,'linescan')
                    deltat = deltat/trial.ny; 
                    nfr = trial.ny*trial.nfr;
                else
                   nfr = trial.nfr;
                end
                
                % tidx, delayshow
                nfr = max(nfr,2); % at least two time points to have a correct display of single-frame data
                tidx = t0+(0:nfr-1)*deltat;
                if ~isempty(op.rmframes)
                    tidx(1:op.rmframes) = [];
                end
                if ~isempty(op.bin) && op.bin~=1
                    tidx = fn_bin(tidx,op.bin);
                end
                nfr2 = length(tidx);
                if op.show__line__scan__delays
                    if isfield(trial.user,'delay')
                        % delay (correction for observed jitter)
                        tidx = tidx+trial.user.delay;
                    end 
                    for i=ind
                        if ~isempty(signal.x(ktr,i).data) % if data is empty, so is delay
                            signal.x(ktr,i).tidx = tidx + signal.x(ktr,i).delay;
                            signal.x(ktr,i).delayshown = signal.x(ktr,i).delay;
                        end
                    end
                else
                    [signal.x(ktr,ind).tidx] = deal(tidx);
                    [signal.x(ktr,ind).delayshown] = deal(0);
                end
                
                % dataop, data2op
                for kdata = 1:1+two
                    
                    % get the data
                    if kdata==1
                        data = [signal.x(ktr,ind).data];
                    else
                        data = [signal.x(ktr,ind).data2];
                    end
                    ncond = size(signal.x(ktr,ind(1)).data,2);
                    data = reshape(data,nfr,ncond*length(ind));
                    
                    % remove first frames
                    if ~isempty(op.rmframes)
                        data(1:op.rmframes,:) = [];
                    end
                    
                    %                     % remove visual stimulation artefact
                    %                     stim = C.trials(C.ktrial).stim;
                    %                     if ~isempty(op.lightartefact) && ~isempty(stim)
                    %                         nstim = size(stim,2);
                    %                         stimstart = stim(1,:);
                    %                         stimstart = repmat(stimstart,nfr2,1);
                    %                         stimlen = stim(2,:);
                    %                         stimend = stimstart + repmat(stimlen,nfr2,1);
                    %                         for i=ind
                    %                             tidxi = repmat(signal.x(ktr,i).tidx(:),1,nstim);
                    %                             instim = (tidxi>=stimstart & tidxi<stimend);
                    %                             idx = any(instim,2);
                    %                             data(idx,i) = data(idx,i)-op.lightartefact;
                    %                         end
                    %                     end
                    
                    % heart pulsation artefact removal
                    if ~isempty(op.heartbeat) && ~isempty(C.trial.heartcycle)
                        % parameters (number of frequency for periodic
                        % signal and modulation)
                        par.K = op.heartbeat(1);
                        par.N = op.heartbeat(2);
                        % heart beat
                        % TODO: add some check that, indeed, the channel
                        % was recording a heart signal!
                        heartcycle = C.trial.heartcycle;
                        if ~isempty(op.rmframes)
                            if C.timeline || strcmp(trial.type,'linescan')
                                heartcycle(1:op.rmframes) = [];
                            else
                                ny = length(heartcycle)/nfr; % ny of the original data, not of the current, possibly cropped, data
                                if mod(ny,1), error programming, end
                                heartcycle(1:op.rmframes*ny) = [];
                            end
                        end
                        % artefact removal
                        data = tps_heartbeat(data,heartcycle,par);
                    end
                    
                    % filter
                    if ~isempty(op.lowpass) || ~isempty(op.highpass)
                        b = op.gauss_d_z_mirr; % 4 booleans
                        fflag = [fn_switch(b(1),'g','s') fn_switch(b(2),'d','') ...
                            fn_switch(b(3),'z','') fn_switch(b(4),'m','')];
                        data = fn_filt(data,{op.lowpass/deltat op.highpass/deltat},fflag);
                    elseif op.gauss_d_z_mirr(2)
                        % detrend
                        data = fn_add(mean(data,1),detrend(data));
                    end
                    
                    % box-car filtering
                    if ~isempty(op.boxcar) && op.boxcar~=1
                        data = imfilter(data,ones(op.boxcar,1));
                    end
                    
                    % bin
                    if ~isempty(op.bin) && op.bin~=1
                        data = fn_bin(data,[op.bin 1]);
                        nfrbin = size(data,1);
                    else
                        nfrbin = nfr2;
                    end
                    
                    % user
                    if ~isempty(op.user)
                        try
                            fun = evalin('base',['@(x)' op.user]);
                            tmp = feval(fun,data(:,1));
                            if size(tmp,2)==size(data,2)
                                % function output is the full data
                                data = tmp;
                            else
                                % function output is one column of the data
                                tmp(1,size(data,2)) = 0; % allocate
                                for i=2:size(data,2)
                                    tmp(:,i) = feval(fun,data(:,i));
                                end
                                data = tmp;
                            end
                        catch ME
                            fprintf 'user-defined signal operation failed: '
                            disp(ME.message)
                        end
                        %data = data ./ exp(-tidx'/20) ./ (1 + 0*exp(-tidx'*1));
                    end
                    
                    %                     % subtract dark noise
                    %                     if ~isempty(op.darknoise)
                    %                         data = data-op.darknoise;
                    %                     end
                    
                    % normalize by average
                    if op.normalize__by__mean
                        avg = mean(data,1);
                        data = fn_mult(data,1./avg);
                        if shotnoiseflag==1 && kdata==1
                            % memorize operation parameters
                            opmem{1} = 1./avg;
                        end
                    end
                    
                    % store the data
                    data = mat2cell(data,nfr2,ncond*ones(1,length(ind)));
                    if kdata==1
                        [signal.x(ktr,ind).dataop]  = deal(data{:});
                    else
                        [signal.x(ktr,ind).data2op] = deal(data{:});
                    end
                end
            end
            
        end
    end   
    methods (Static)
        function s = signalopdef_default
            s = struct( ...
                'show__line__scan__delays',    true, ...
                'rmframes',         [], ...
                ...'lightartefact',    [], ...
                'heartbeat',        [], ...
                'lowpass',          [], ...
                'highpass',         [], ...
                'gauss_d_z_mirr',   [true false true false], ...
                'boxcar',           [], ...
                'bin',              [], ...
                'user',             [], ...
                ...'darknoise',        [], ...
                'normalize__by__mean',  false, ...
                'channelscale',     .05 ...
                );
        end
        function s = signalopdef_spec
            s = struct( ...
                'show__line__scan__delays',    'logical', ...
                'rmframes',         'xslider 0 20 1 %i [1]', ...
                ...'lightartefact',    'xlogslider 0 2 .1 [10]', ...
                'heartbeat',        'xstepper 2 1,0 [2 3]', ...
                'lowpass',          'xlogslider -4 1 .02 %.1e [-2]', ...
                'highpass',         'xlogslider -4 1 .02 %.1e [1]', ...
                'gauss_d_z_mirr',   'multcheck 4', ...
                'boxcar',           'xslider 1 20 1 %i [2]', ...
                'bin',              'xslider 1 20 1 %i [2]', ...
                'user',             'xchar 8 [x./(1+x)]', ...
                ...'darknoise',        'xslider 0 20', ...
                'normalize__by__mean',  'logical', ...
                'channelscale',     'loglogslider -6 2 .1' ...
                );
        end
    end
    
    % Advanced output
    methods
        function slice = getslice(C,ktrials,ind)
            % function slice = getslice(C,ktrials,ind)
            %---
            % ktrial and ksel are vector of indices
            if nargin<2 || isequal(ktrials,0), ktrials = 1:C.ntrial; end
            if nargin<3 || isequal(ktrials,0), ind=1:C.nx; end
            interpdata(C,ktrials,ind)
            operation(C,ktrials,ind)
            if C.internpar.guessspikes, computespikes(C,ktrials,ind), end
            slice = C.signals(1).x(ktrials,ind);
        end
        function map = getselectionmap(C,fillmethod,varargin)
            % function map = getselectionmap(C[,fillmethod,filloptions])
            % returns a map of the active selections
            % additional regions can be created in the remaining empty
            % space according to 'fillmethod' and its arguments:
            % - 'square',k          divides the space more or less in
            %                       square blocks of size k  
            % - 'surround',k,isel   creates concentric regions of width k;
            %                       only from one selection if a second
            %                       argument is specified
            if ~C.ntrial || ~checkconsistency(C.trials) || C.seldotrial
                error('cannot create selection map: trials or selections do not match')
            end
            if nargin<2, fillmethod=''; end
            % selection regions
            nx = C.trial.nx;
            ny = C.trial.ny;
            map = zeros(nx,ny);
            active = false(1,C.nsel);
            for k=C.nsel:-1:1
                sel = C.signals(1).x(1,k).sel;
                active(k) = sel.active;
                if active(k)
                    map(sel.dataind) = k;
                end
            end
            % additional regions
            switch fillmethod
                case ''
                    % nothing to do
                case 'square'
                    % create the background
                    w = varargin{1};
                    nx2 = ceil(nx/w);
                    ny2 = ceil(ny/w);
                    square = reshape(1:nx2*ny2,nx2,ny2);
                    square = kron(square,ones(w));
                    square = square(1:nx,1:ny);
                    % add the map onto this background
                    f = logical(map);
                    square = square + C.nsel;
                    square(f) = map(f);
                    map = square;
                case 'surround'
                    w = varargin{1};
                    if length(varargin)==2 && ~isempty(varargin{2})
                        ksel = varargin{2};
                    else
                        ksel = find(active);
                    end
                    numsel = length(ksel);
                    % let's compute the distances btw. pairs of points
                    [xx yy] = ndgrid(1:nx,1:ny);
                    dx = fn_add(xx(:),-xx(:)');
                    dy = fn_add(yy(:),-yy(:)');
                    dpoints = sqrt(dx.^2+dy.^2);
                    % now the distances to each selection
                    dsel = zeros(nx*ny,numsel);
                    for i=1:numsel
                        ind = C.signals(1).x(1,ksel(i)).sel.dataind;
                        dsel(:,i) = min(dpoints(ind,:),[],1);
                    end
                    % now, which selection is the closest
                    [dclosest iclosest] = min(dsel,[],2); 
                    % now, give a region code according to closest
                    % selection and to the distance
                    dclosest = floor(dclosest/w);
                    code = iclosest + dclosest*numsel;
                    code = reshape(code,nx,ny);
                    % add the map onto this background
                    f = logical(map);
                    code = code + C.nsel;
                    code(f) = map(f);
                    map = code;
            end
        end
        function trig = gettriggers(C)
            % returns a structure with global indices (over frames +
            % trials) of events (trial start, stim start and marked spikes)
            if ~checkconsistency(C.trials) || C.ntrial==0
                warning('trials not consistents')
            end
            % start of trials
            trig.trialstart = 1:C.nfr:C.nfr*C.ntrial; % this one was easy!
            % stimulation
            addinfo = C.trials(1).addinfo;
            stim = C.trials(1).stim;
            if isfield(addinfo,'stim') && ~isempty(stim)
                stimstart = floor(stim(1,1)/C.trials(1).dt);
                trig.stimstart = stimstart:C.nfr:C.nfr*C.ntrial;
            end
            % spikes
            if C.nsel>0
                for k=1:C.nsel
                    spikes = {C.signals(1).x(:,k).spikes};
                    for i=1:C.ntrial
                        if isempty(spikes{i}), continue, end
                        if size(spikes{i},1)==1, spikes{i}(2,:)=0; end % old version
                        spikes{i}(1,:) = round(spikes{i}(1,:)/C.trials(1).dt) + (i-1)*C.nfr;
                    end
                    spikes = [spikes{:}]; % row1=times, row2=flags
                    if isempty(spikes), continue, end
                    allflags = unique(spikes(2,:));
                    nflag = length(allflags);
                    for iflag=1:nflag
                        flag = allflags(iflag);
                        f = ['spike' num2str(k)];
                        if flag>0, f = [f '_' num2str(flag)]; end %#ok<AGROW>
                        idx = (spikes(2,:) == flag);
                        trig.(f) = spikes(1,idx);
                    end
                end
            end
        end
    end
    
    % Load object
    methods (Static)
        function C = loadobj(C)
            % old versions
            if isempty(C.nsel)
                C.nsel = C.signals(1).nx; % TODO: it might happen that a 'pointer' signal becomes a 'selection' signal... 
            end
            if isempty(C.electrophys)
                C.electrophys = tps_electrophys(C.trials);
            end
            % now we have a version field!
            if isempty(C.version), C.version = 1.1; end
            % version 1.2: added the concept of stim table; trials may have
            % a stim table, but do not share the same one
            if C.version < 1.2
                mergestimtables(C.trials)
            end
            % version 1.4: auto-link beween V.a4d.SI and V.a4d.G (note that
            % this has probably no consequence on V.content)
            % version 1.41: masked the options 'set heart' in the trial
            % menu
            C.version = 1.41;
            % check whether comment file has been moved
            if ~isempty(C.docfile)
                try
                    C.docfile = tps_locatefile(C.docfile);
                catch %#ok<CTCH>
                    errordlg 'comments file has been lost!'
                end
            end
        end
    end
    
    % Misc
    methods
        function access(C)
            keyboard
        end
        function autorepair(C)
            erasedata(C)
            if ~isempty(C.electrophys) && any([C.electrophys.n])
                answer = 'No'; %questdlg('Reset also electrophysiology processings','Confirmation','Yes','No','No');
                if strcmp(answer,'Yes'), C.electrophys = tps_electrophys(C.trials); end
            end
        end
    end
end




function cond = condnum(nc,condname)

condarray = tps_readconditionname(condname,nc);
if all(fn_map(@isscalar,condarray(:,1))) && all(all(fn_isemptyc(condarray(:,2:end))))
    % set of single conditions
    cond = [condarray{:,1}];
else
    cond = [];
end

end
    