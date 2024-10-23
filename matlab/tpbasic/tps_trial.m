classdef tps_trial < hgsetget
    % function T = tps_trial
    % function T = tps_trial(f[,flag])
    % function T = tps_trial(data[,sfr][,dx(um)|0,dt(s)][,'linescan'][,'zstack',dz][,'scanning'][,'units',xunit,tunit])
    % function T = tps_trial(data[,sfr],T0)
    % function T = tps_trial(fmat[,dx(um)|0,dt(s)][,'linescan'][,'zstack',dz][,'scanning'][,'units',xunit,tunit])
    % function T = tps_trial(fmat,T0)
    % function T = applyprocessing(T0[,'correctshift',nshift][,'translate',shift][,'crop',crop][,'user',fun])
    %---
    % flag can be 
    % - 'avg'   trial data is obtained by averaging data from several files
    %           (data from each file should all be of the same size)
    % - 'cat'   trial data is obtained by concatenating data from several
    %           files in the time domain (data from each file should all be
    %           of the same size, except for the time dimension)
    % - 'mim'   "movie from images": each file contains a single image, all
    %           of which have the same size and should be assembled to make
    %           a single movie 
    % - 'img'   "single images": even if file might all contain single
    %           images all of the same size, they should be considered as
    %           separate trials (if this is the case but neither 'mim' or
    %           'img' flag is specified, user is prompted to confirm
    %           whether he want a single movie trial or multiple
    %           single-image trials)
    % 
    % fmat can be either a file saved with tps_trial.savedata (i.e. has
    % variables data, sfr and rec), or a matlab file containing a single
    % variable (typically, saved using fn_savevar)
    % 
    % special fields of a tps_trial object:
    % - fullinfo    stores the maximal number of information from the
    %               headers of the data file and/or from other files (for
    %               example, acquisition settings saved by LabView)
    % - addinfo     stores part of this information which might be relevant
    %               for display (for example, the depth of the acquisition)
    % - user        this field can be used by the user to store any
    %               additional information; it must be a structure
    % 
    % special for MES and MESC data: flag can be 'prompt' to let the user
    % select which trials to take from the experiment, or a vector of trial
    % indices
    
    properties (SetAccess='private')
        version
    end
    % Acquisition information
    % (data)
    properties (SetAccess='private')
        file = '';
        preprocess = struct('multifile',[],'op',[]);
        analogfile = {};
        origin = '';
    end
    % (header info)
    properties (Transient, SetAccess='private') % should be Dependent, but we need to keep it non-dependent for backward compatibility reasons
        sizes
    end
    properties (SetAccess='private')
        sizes0  % must appear before xbin and tbin in property declarations; see set.sizes0
        sfrchannel = false;
        fullinfo
    end
    properties
        addinfo
        status = 'n'; % 'n' for normal, 'r' for rejected or 's' for special
        type = ''; % 'movie', 'zstack', 'linescan'
        scanning = [];
        xbin = 1;
        tbin = 1;
        dx = [];
        dy = [];
        dz = [];
        dt = [];
        t0 = [];
        xunit
        tunit
        user = struct;
    end
    properties (Dependent, Transient, SetAccess='private')
        nx
        ny
        nfr
        nc
        
        linedur % tunit
        tidx    % tunit
        dt_sec
    end
    properties (Transient)
        usertransient = struct;
    end
    properties
        stimtable       % table of all possible stimulations
        eventlist = struct('name',cell(1,0),'cond',cell(1,0),'time',cell(1,0));
    end
    properties (SetAccess='private')
        stimid = [];    % entry in table
    end
    properties (Dependent, Transient)
        stim            % time course of stimulation of this trial
    end
    properties (Dependent, Transient, SetAccess='private')
        stimdetails     % structure with additional information on stimulation
    end
    % Operation definition
    properties
        opdef = tps_dataopdef.empty(1,0);
    end
    properties (Transient)
        opmem = tps_dataopdef.empty(1,0);
    end
    % Private access data
    properties (Transient)
        internpar = tpv_internpar; % internal parameter (e.g. whether to save memory), shared with tpview and tpv_content
    end
    properties (Transient, Access='private')
        data_unsaved
        recording_data = {[]}; % this initialization is necessary to avoid error in get.recording when the recording data has not be read
        preventdata0comp = false; % this flag is used to prevent reading data when attempting to save data0 property
    end
    properties (Access='private') % data
        recording_header = struct('name',{},'t0',{},'dt',{},'n',{},'signal',{}); % the 'signal' field is actually stored in recording_data
    end
    properties
        data0 % average data frame
    end
    % Automatic calculations
    properties (Dependent, Transient, SetAccess='private')    
        data
        sfr
        shotnoise
        dataop
        sfrop
        shotnoiseop
        
        dataopmem

        recording
    end
    properties (Transient, SetAccess='private')
        heartcycle
    end
    properties (Transient, Access='private')
        basics = struct('data',[],'sfr',[],'shotnoise',[]);
        conditions = struct('data',struct,'sfr',struct,'shotnoise',struct, ...
            'dataop',struct,'sfrop',struct,'shotnoiseop',struct, ...
            'dataopmem',struct);
        operationsteps = struct('data',struct('op',{},'value',{}), ...
            'sfr',struct('op',{},'value',{}), ...
            'shotnoise',struct('op',{},'value',{}));
    end
    
    % Constructor (calls tps_read)
    methods
        function T = tps_trial(varargin)
            % function T = tps_trial
            % function T = tps_trial(f[,flag])
            % function T = tps_trial(data[,sfr][,dx(um)|0,dt(s)][,'linescan'][,'zstack',dz][,'scanning'][,'units',xunit,tunit])
            % function T = tps_trial(data[,sfr],T0)
            % function T = tps_trial(fmat[,dx(um)|0,dt(s)][,'linescan'][,'zstack',dz][,'scanning'][,'units',xunit,tunit])
            % function T = tps_trial(fmat,T0)
            T.version = 2.0;
            if nargin==0, varargin = {0}; end
            x = varargin{1};
            if ischar(x) || isstring(x) ...
                    || (iscell(x) && (ischar(x{1}) || isstring(x{1})))
                % first argument is a file or list of files
                f = varargin{1};
                if strcmp(fn_fileparts(f, 'ext'), '.mat')
                    T = userdefine(T,varargin{:});
                else
                    T = tps_trial.readheader(varargin{:});
                end
            else
                % make trial from data, use provided header information if
                % any
                T = userdefine(T,varargin{:});
            end
        end
        function T = applyprocessing(T0,varargin)
            % multi-trial
            ntrial = numel(T0);
            if ntrial>1
                T = T0;
                for ktrial=1:ntrial
                    T(ktrial) = applyprocessing(T0(ktrial),varargin{:});
                end
                return
            end
            % checks
            if isempty(T0.file)
                error 'cannot apply preprocessing to tps_trial object whose data is not saved'
            end
            if T0.xbin~=1 || T0.tbin~=1
                error 'cannot apply preprocessing to tps_trial object with binning'
            end
            if ~isempty(T0.preprocess.op)
                error 'cannot apply preprocessing to tps_trial object that already has preprocessing'
            end
            % preprocessing structure
            if nargin==2
                preproc = varargin{1};
            else
                % build structure from arguments, cannot use the 'struct'
                % function directly because of cell arrays that could lead
                % to non-scalar structure
                nop = length(varargin)/2;
                preproc = struct('name',cell(1,nop),'value',cell(1,nop));
                for k=1:nop
                    preproc(k).name = varargin{2*k-1};
                    preproc(k).value = varargin{2*k};
                end
            end
            % copy headers from T0
            T = tps_trial;
            copyheader(T,T0)
            % data file
            erasedata(T)
            T.file = T0.file;
            T.preprocess.multifile = T0.preprocess.multifile;
            % set the preprocessing
            T.preprocess.op = preproc;
            % set the data size
            siz = T0.sizes;
            for kop=1:length(preproc)
                opk = preproc(kop);
                switch opk.name
                    case 'correctshift'
                        nshift = opk.value; if iscell(nshift), nshift = nshift{1}; end
                        siz(1) = siz(1)-abs(nshift);
                    case 'crop'
                        crop = opk.value;
                        siz(1) = length(crop{1});
                        siz(2) = length(crop{2});
                    case 'user'
                        % we need to run the function in order to check the
                        % size of its output
                        fun = opk.value;
                        xtest = fun(zeros(siz));
                        siz = size(xtest); siz(end+1:4) = 1;
                        if length(siz)>4, error 'custom user-defined preprocessing is not valid', end
                end
            end
            T.sizes0 = siz;
            T.sizes = siz;
        end
    end
    methods (Access='private')
        function T = userdefine(T,varargin)
            % data (automatic setting of sizes)
            dat = varargin{1};
            % multiple trials?
            if (isnumeric(dat) || islogical(dat)) && ndims(dat)>4
                % multiple trials -> transform to cell array
                dat = num2cell(dat,1:4);
                dosfr = nargin>=3 && (isnumeric(varargin{2}) || islogical(varargin{2}));
                if dosfr, varargin{2} = num2cell(varargin{2},1:4); end
            end
            if iscell(dat)
                dosfr = nargin>=3 && iscell(varargin{2});
                if dosfr, sf = varargin{2}; end
                % multiple trials
                for k=1:numel(dat)
                    if k>1, T(k) = tps_trial; end
                    if dosfr
                        T(k) = userdefine(T(k),dat{k},sf{k},varargin{3:end});
                    else
                        T(k) = userdefine(T(k),dat{k},varargin{2:end});
                    end
                end
                return
            end
            % set data
            if ischar(dat) || isstring(dat)
                % mat file -> postpone reading
                f = dat;
                assert(strcmp(fn_fileparts(f,'ext'),'.mat'), 'file argument to tps_trial.userdefine must be a mat file')
                T.file = f;
                % format as saved by tps_trial.savedata?
                w = whos('-file', f);
                wnames = {w.name};
                if ismember('data', wnames) && all(ismember(wnames, {'data', 'rec', 'sfr'}))
                    assert(strcmp(wnames{1}, 'data'))
                    % presence of sfr channel?
                    isfr = find(strcmp(wnames, 'sfr'));
                    T.sfrchannel = ~isempty(isfr) && (prod(w(isfr).size) > 0);
                    % presence of recording?
                    if ismember('rec', wnames)
                        rec = load(f,'rec');
                        rec = rec.rec;  % cell array
                        if ~isscalar(rec) || ~isempty(rec{1})
                            warning 'non-empty recordings will not be loaded'
                        end
                    end
                else
                    assert(isscalar(w), 'if mat file was not saved using tps_trial.savedata, it must have only one variable')                    
                    T.sfrchannel = false;
                end                    
                % no need to read data, it is enough to set sizes
                T.sizes0 = w(1).size;
            elseif isnumeric(dat) || islogical(dat)
                dosfr = nargin>=3 && (isnumeric(varargin{2}) || islogical(varargin{2})) ...
                    && ~isvector(varargin{2});
                if dosfr
                    sf = varargin{2};
                    varargin(2) = []; % header info in varargin(2:end) rather than varargin(3:end)
                    T.data_unsaved = {dat sf};
                    setdata(T,dat,sf); % sets sfrchannel property
                else
                    T.data_unsaved = {dat []};
                    setdata(T,dat); % sets sfrchannel property
                end
            else
                error argument
            end
            % header information
            if nargin>2 && isa(varargin{2},'tps_trial')
                T0 = varargin{2};
                copyheader(T,T0(1));
                if isempty(strfind(T0.origin,'header]'))
                    T.origin = ['user [' T0.origin ' header]'];
                end
            else
                T.origin = 'user';
                defaultheadervalues(T)
                % information provided by user
                k = 2;
                while k<=nargin-1
                    a = varargin{k};
                    k = k+1;
                    if ischar(a) || isstring(a)
                        switch a
                            case 'zstack'
                                T.type = 'zstack';
                                T.dz = varargin{k};
                                k = k+1;
                            case 'linescan'
                                T.type = 'linescan';
                                T.scanning = true;
                            case 'scanning'
                                T.scanning = true;
                            case 'units'
                                T.xunit = varargin{k};
                                T.tunit = varargin{k+1};
                                k = k+2;
                            otherwise
                                error('unknown flag ''%s''', a)
                        end
                    else
                        if a
                            [T.dx, T.dy] = deal(a);
                            T.xunit = 'um';
                        end
                        T.dt = varargin{k};
                        k = k+1;
                        T.tunit = 's';
                    end
                end
            end
        end
        function copyheader(T,T0)
            % copy all header information, but not the data, size, and
            % preprocessing
            if ~isscalar(T0) || ~isscalar(T), error('T0 and T must be scalar'), end
            T.origin = T0.origin;
            T.fullinfo = T0.fullinfo;
            T.addinfo = T0.addinfo;
            T.status = T0.status;
            T.user = T0.user;
            T.type = T0.type;
            T.scanning = T0.scanning;
            T.xbin = 1; % he he!!
            T.tbin = 1;
            T.dx = T0.dx;
            T.dy = T0.dy;
            T.dz = T0.dz;
            T.dt = T0.dt;
            T.t0 = T0.t0;
            T.xunit = T0.xunit;
            T.tunit = T0.tunit;
            T.stimtable = T0.stimtable;
            T.stimid = T0.stimid;
            T.opdef = T0.opdef;
            T.opmem = T0.opmem;
        end
    end

    % GET/SET normal
    methods
        function set.status(T,flag)
            if ~ischar(flag) || ~isscalar(flag) || ~ismember(flag,'nrs')
                error '''status'' must be one of ''n'' (normal), ''r'' (rejected) or ''s'' (special)'
            end
            T.status = flag;
        end
    end
    % GET/SET dependent
    methods
        % Sizes
        function nx = get.nx(T)
            nx = T.sizes(1);
        end
        function ny = get.ny(T)
            ny = T.sizes(2);
        end
        function nfr = get.nfr(T)
            nfr = T.sizes(3);
        end
        function nc = get.nc(T)
            nc = T.sizes(4);
        end
        function set.dt(T,dt)
            if dt==T.dt, return, end
            T.dt = dt;
            erasecomp(T) % some computations have a time resolution - dependent behavior
        end
        function set.tunit(T,u)
            if strcmp(u,T.tunit), return, end
            T.tunit = u;
            erasecomp(T) % some computations have a time resolution - dependent behavior
        end
        function settime(T,u,dt)
            T.tunit = u;
            T.dt = dt;
        end
        function setspace(T,u,dx)
            T.xunit = u;
            T.dx = dx;
        end
        % Other
        function linedur = get.linedur(T)
            if isempty(T.scanning)
                % don't know whether it is scanning or not
                linedur = [];
            elseif T.scanning
                linedur = T.dt / T.ny;
            else
                linedur = 0;
            end
        end
        function tidx = get.tidx(T)
            if isempty(T.dt), tidx=[]; return, end
            switch T.type
                case 'zstack'
                    tidx = (0:T.nfr-1)*T.dz;
                case {'movie' 'linescan'}
                    tidx = T.t0+(0:T.nfr-1)*T.dt;
                otherwise
                    error('cannot get tidx: type is not set or is wrong')
            end
        end
        function deltat = get.dt_sec(T)
            fact = fn_switch(T.tunit,'s',1,'ms',1e-3,{'min' 'minute'},60, ...
                {'h' 'hour'},3600,'day',3600*24,[]);
            deltat = T.dt * fact;
        end
        function str = fileshort(T)
            if size(T.file,1)<=1
                str = T.file;
            else
                idxdiff = find(any(diff(T.file)),1,'first');
                if isempty(idxdiff), error 'all files are identical!!??', end
                str = T.file(1,1:idxdiff-1);
            end
        end
    end
    % Stimulation
    methods
        function stim = get.stim(T)
            if isempty(T.stimtable) || isempty(T.stimid), error programming, end
            if T.nc>1
                stim = getstimbycond(T,'all conditions');
            else
                stim = getstim(T.stimtable,T.stimid);
            end
        end
        function s = get.stimdetails(T)
            if isempty(T.stimtable) || isempty(T.stimid), error programming, end
            for k=1:T.nc, s(k) = getdetails(T.stimtable,T.stimid(k)); end %#ok<AGROW>
        end
        function [stim conds] = getstimbycond(T,condname)
            if T.nc==1
                stim = getstim(T.stimtable,T.stimid);
                return
            end
            % get all necessary conditions
            if strcmp(condname,'all conditions')
                conds = 1:T.nc;
            else
                condarray = tps_readconditionname(condname,T.nc);
                conds = [condarray{:}];
            end
            % global stim is valid only if all non-empty stims are the same
            stim = getglobalstim(T.stimtable,T.stimid(conds));
        end
        function set.stim(T,stim)
            if isempty(T.version) || T.version<1.2
                % this can happen when loading an object from the time
                % where tps_trial had a 'stim' field
                T.stimtable = tps_stimtable;
                T.stimid = addstim(T.stimtable,stim);
                return
            end
            if isempty(T.stimtable), error programming, end
            if T.nc>1, error 'cannot set stim directly for multi-condition trial, use setstim method instead', end
            T.stimid = addstim(T.stimtable,stim);
        end
        function setstim(T,x)
            % function setstim(T,id)
            % function setstim(T,struct)
            if ~isscalar(T), mergestimtables(T), end
            if any(fn_isemptyc({T.stimtable})), error programming, end 
            if any(diff([T.nc])) || length(x)~=T(1).nc, error 'number of stims of a trial must fit its number of conditions!', end 
            stimids = zeros(1,T(1).nc);
            for k=1:T(1).nc
                stimids(k) = addstim(T(1).stimtable,x(k));
            end
            [T.stimid] = deal(stimids);
        end
        function mergestimtables(T)
            uniquetable = mergetables([T.stimtable]);
            ok = ismember([uniquetable.table.id],[T.stimid]);
            uniquetable.table(~ok) = [];
            [T.stimtable] = deal(uniquetable);
        end
        function compactstimtable(T,table)
            if nargin<2
                mergestimtables(T)
                table = T(1).stimtable;
            end
            table = table(ismember([table.id],[T.stimid]));
            [T.stimtable] = deal(table);
        end
    end
    % GET/SET dependences between (non dependent) fields
    methods
        function set.sizes0(T,siz)
            siz(end+1:4) = 1;
            T.sizes0 = siz;
            % The line below is extremely dangerous; loading works BECAUSE
            % sizes0 later than xbin and tbin in property declarations.
            if ~isempty(siz), T.sizes = max(1,floor([siz(1:2)/T.xbin siz(3)/T.tbin siz(4)])); end
        end
        function setdata(T,data,sfr)
            % set size and sfrchannel properties
            if ~isempty(data)
                T.sizes0 = size(data); % (normally, setting sizes0 and then data, or data and then sizes0 ends with the same result)
                T.sfrchannel = nargin>=3 && ~isempty(sfr);
            elseif fn_dodebug
                disp 'how can i get here?'
                keyboard
            end
            % initial binning
            data = fn_bin(data,[T.xbin T.xbin T.tbin]);
            if nargin>=3, sfr = fn_bin(sfr,[T.xbin T.xbin T.tbin]); end
            % set data and sfr
            if nargin>=3
                % first sfr, because the last one that we want to set, and
                % keep for sure inside the memory pool, is data
                if isempty(T.basics.sfr), T.basics.sfr = memorypool.item(); end
                recover(T.basics.sfr,sfr,[fn_fileparts(fileshort(T),'base') ' - sfr']);
            end
            if isempty(T.basics.data), T.basics.data = memorypool.item(); end
            recover(T.basics.data,data,[fn_fileparts(fileshort(T),'base') ' - data'],10); % weight of 10 for the memory pool item
            % set frame average
            avg = mean(data,3);
            if ~issparse(avg)
                avg = cast(avg,T.internpar.floattype);
            end
            T.data0 = avg;
        end
        function r = get.recording(T)
            r = T.recording_header;
            if ~isempty(r)
                if r(1).n~=0 && isempty(T.recording_data{1})
                    readrecording(T)
                end
                [r.signal] = deal(T.recording_data{:});
            end
        end
        function set.xbin(T,xbin)
            if xbin==T.xbin, return, end
            % data became invalid
            erasedata(T)
            % modify dx and dy accordingly ATTENTION! THIS IS A VERY
            % DANGEROUS WAY OF PROGRAMMING, AND WORKS ONLY BECAUSE AT
            % LOADING, dx AND dy ARE SET LATER THAN xbin
            T.dx = T.dx * (xbin/T.xbin);
            T.dy = T.dy * (xbin/T.xbin); 
            % set xbin
            T.xbin = xbin;
            % update sizes
            T.sizes(1:2) = max(1,floor(T.sizes0(1:2)/xbin));
        end
        function set.tbin(T,tbin)
            if tbin==T.tbin, return, end
            % data became invalid
            erasedata(T)
            % modify temporal resolution accordingly ATTENTION! THIS IS A VERY
            % DANGEROUS WAY OF PROGRAMMING, AND WORKS ONLY BECAUSE AT
            % LOADING, dt IS SET LATER THAN tbin
            T.dt = T.dt * (tbin/T.tbin);
            % set temporal binnig
            T.tbin = tbin;
            % update sizes
            T.sizes(3) = max(1,floor(T.sizes0(3)/tbin));
        end
    end
    % GET auto-compute at first access attempt
    methods
        function data = get.data(T)
            if isempty(T.basics.data), T.basics.data = memorypool.item(); end
            data = T.basics.data.data;
            if isempty(data)
                try
                    readdata(T) % note that this will call setdata, which itself will also set data0
                catch ME
                    disp(['data reading failed: ' ME.message])
                end
                data = T.basics.data.data;
            end
        end
        function data0 = get.data0(T)
            if T.preventdata0comp
                % the flag was set by function saveobj to prevent data
                % reading and rather save an empty data0; after saving
                % however it must be reset to false
                T.preventdata0comp = false;
            elseif isempty(T.data0)
                T.data0 = mean(T.data(:,:,:),3);
            end
            data0 = T.data0;
        end
        function compute_data0_per_blocks(T,nblock)
            % this methods allows computing data0 even when data does not
            % fit into memory
            % note that we can't check easily that data0 is empty, so
            % calling this method several times will result in repeating
            % the computation
            if nblock==1
                % no need to read data per block
                T.data0;
                return
            end
            dat0 = 0; % class: double
            for kblock = 1:nblock
                idx = floor((kblock-1)*T.nfr/nblock)+1:floor(kblock*T.nfr/nblock);
                dat0 = dat0 + sum(T.getdatablock(idx),3);
            end
            T.data0 = dat0 / T.nfr;
        end
        function sfr = get.sfr(T)
            if isempty(T.basics.sfr), T.basics.sfr = memorypool.item(); end
            sfr = T.basics.sfr.data;
            if isempty(sfr) && ~isempty(T.file) && T.sfrchannel
                try
                    readdata(T)
                catch ME
                    disp(['data reading failed: ' ME.message])
                end
                sfr = T.basics.sfr.data;
            end
        end
        function shotnoise = get.shotnoise(T)
            if isempty(T.basics.shotnoise), T.basics.shotnoise = memorypool.item(); end
            shotnoise = T.basics.shotnoise.data;
            if isempty(shotnoise) && ~isempty(T.data)
                % we will create artificial shot noise based
                % - if PMT gain is known, on average value in each pixel +
                %   calibrations done at the same PMT gain value
                % - otherwise, on the average value and standard deviation
                dormfirstsec = (T.nfr*T.dt > 2);
                % data: remove first second (shutter open)
                if strcmp(T.type,'linescan');
                    nt = T.ny*T.nfr;
                    dat = reshape(T.data,T.nx,1,nt,T.nc);
                    if dormfirstsec
                        dat = dat(:,1,ceil(1/T.linedur):end,:);
                    end
                else
                    dat = T.data;
                    if dormfirstsec
                        dat = T.data(:,:,ceil(1/T.dt):end,:);
                    end
                    nt = T.nfr;
                end
                dat = fn_float(dat);
                moyenne = reshape(repmat(mean(dat,3),[1 1 nt 1]),T.sizes);
                
                disp('compute fake shot noise based on signal fast variations in the present data')
                % compute variance from the differences of two
                % successive frames
                variance = var(diff(dat,1,3),0,3)/2;
                
                stddev = reshape(repmat(sqrt(variance),[1 1 nt 1]),T.sizes);
                shotnoise = single(moyenne + stddev.*randn(T.sizes));
                recover(T.basics.shotnoise,shotnoise,[fn_fileparts(fileshort(T),'base') ' - shotnoise'])
            end
        end
        function x = get.dataop(T)
            x = operation(T,'data');
        end
        function x = get.sfrop(T)
            x = operation(T,'sfr');
        end
        function x = get.shotnoiseop(T)
            x = operation(T,'shotnoise');
        end
        function x = get.dataopmem(T)
            x = operation(T,'data',T.opmem);
        end
        function [dat, binning] = getcondition(T,dataname,condname,doenlarge)
            % function dat = getcondition(T,dataname,condname[,doenlarge])
            %---
            % getcondition does not return an empty array if the data does
            % not exist, but a zero array of the trial size
            
            % input
            if nargin<4, doenlarge = true; end
            
            % get 'all conditions' data
            switch dataname
                case 'data0'
                    dat = mean(T.data0,4);
                    return
                case {'data' 'sfr' 'shotnoise'}
                    dat0 = T.(dataname);
                case {'dataop' 'sfrop' 'shotnoiseop'}
                    dat0 = operation(T,dataname(1:end-2),[],false);
                case 'dataopmem'
                    dat0 = operation(T,'data',T.opmem,false);
                otherwise
                    error programming
            end
            
            % compute the requested condition
            if isempty(dat0)
                dat = [];
            elseif T.nc==1 || isempty(condname) || strcmp(condname,'all conditions')
                dat = dat0;
                if doenlarge, dat = fn_enlarge(dat,T.sizes(1:3)); end
            else
                % get condition(s)
                [condarray condfieldname] = tps_readconditionname(condname,T.nc);
                nconds = size(condarray,1);
                
                % memory pool entry
                if isfield(T.conditions.(dataname),condfieldname)
                    item = T.conditions.(dataname).(condfieldname);
                else
                    item = memorypool.item();
                    T.conditions.(dataname).(condfieldname) = item;
                end
                
                % data (compute only if necessary)
                dat = item.data;
                if isempty(dat)
                    if nconds>1, s = size(dat0); dat = zeros([s(1:3) nconds],class(dat0)); end
                    for k=1:nconds
                        condk = condarray(k,:);
                        if isempty(condk{1}) || (isempty(condk{3}) && ~isempty(condk{4}))
                            error 'condition not valid (no numerator, or negative term without a positive term'
                        end
                        c = cell(1,4);
                        for i=1:4
                            if ~isempty(condk{i})
                                c{i} = dat0(:,:,:,condk{i});
                                if ~isscalar(condk{i}), c{i} = mean(c{i},4); end
                            end
                        end
                        a = c{1};
                        if ~isempty(c{2}), a = a-c{2}; end
                        if ~isempty(c{3})
                            b = c{3};
                            if ~isempty(c{4}), b = b-c{4}; end
                            a = fn_float(a)./b;
                        end
                        if nconds==1, dat = a; else dat(:,:,:,k) = a; end
                    end
                    recover(item,dat,[fn_fileparts(fileshort(T),'base') ' - ' dataname ' - ' condname])
                end
                if doenlarge, dat = fn_enlarge(dat,T.sizes(1:3)); end
            end
        end
        function heartcycl = get.heartcycle(T)
            if isempty(T.heartcycle)
                setheart(T)
            end
            heartcycl = T.heartcycle;
        end            
    end
    
    % Recognized extensions
    methods (Static)
        function [ext, prompt] = knownExtensions()
            ext = {'tptrial','cfd','mpd','xml','mes','mesc','blk','vdq','da', ...
                'pcl', ...
                'avi','bmp','png','jpg','tif','tiff','gif', ...
                'mat'};
            prompt = {'*.tptrial','Trials';'*.tif','ScanImage';'*.mes;*.mesc','Femtonics'; ...
                '*.blk;*.vdq','Optical Imaging';'*.da','NeuroPlex';'*.xml','Prairie';'cfd','CFD';'mpd','MPD'; ...
                '*.pcl','ScienceWares Photon Imaging';...
                '*.avi','Movie';'*.bmp;*.png;*.jpg;*.tif;*.tiff;*.gif','Image';...
                '*.mat','Matlab'};
        end
    end
    % Read header
    methods (Static)
        function T = readfile(f,flag)
            if nargin<2, flag=''; end
            if iscell(f), f=char(f); end
            % Read headers
            T = tps_trial.readheader(f,flag);
            if isempty(T), return, end
            % Read data
            readdata(T);
            % Read recordings
            readrecording(T);
        end
        function T = readheader(f,flag)
            if nargin<2, flag=''; end
            % read format-specific headers for each trial
            if iscell(f), f = char(f); end % convert cell array of file names to a 2D char array
            nf = size(f,1);
            if strcmp(flag,'mim')
                % files contain single images, which are assembled to make 
                % a single movie 
                ntrial = 1;
            else
                ntrial = nf;
            end
            T = cell(1,ntrial);
            if ntrial>1, fn_progress('reading header',ntrial), end
            for k=1:ntrial
                if ntrial>1, fn_progress(k), end
                % new trial
                T{k} = tps_trial; erasedata(T{k}) 
                % read header
                fk = deblank(f(k,:));
                [dum1 dum2 ext] = fileparts(fk); %#ok<ASGLU>
                ext = lower(ext(2:end));
                if length(ext)==1, ext = 'oldvdaq'; elseif strfind(ext,'blk'), ext='blk'; end
                switch ext
                    case 'tptrial'
                        s = load('-MAT',fk);
                        F = fieldnames(s);
                        if length(F)~=1, error('mat file should contain exactly one tps_trial object'), end
                        T{k} = s.(F{1});
                        % check if data file exists
                        locatefile(T{k},f)
                    case 'avi'
                        T{k}.file = fk;
                        VR = VideoReader(fk); %#ok<TNMLP>
                        T{k}.sizes = [VR.Width VR.Height VR.NumberOfFrames 1];
                        T{k}.dt = 1/VR.FrameRate;
                        T{k}.tunit = 's';
                    case 'cfd'
                        T{k}.file = fk;
                        header = cfd_headers(fk);
                        readcfdheader(T{k},header)
                    case 'mpd'
                        T{k}.file = fk;
                        header = mpd_headers(fk);
                        readmpdheader(T{k},header)
                    case {'tif' 'tiff'}
                        T{k}.file = fk;
                        try
                            info=imfinfo(fk); % see ScanImage_genericOpenTif.m
                            header=info(1).ImageDescription;
                            header=ScanImage_parseHeader(header);
                            readscanimageheader(T{k},header)
                        catch %#ok<CTCH>
                            % no header information -> some default values
                            defaultheadervalues(T{k})
                            % we need to read the data to set the proper size
                            readdata(T{k})
                        end
                    case 'xml'
                        T{k}.file = fk;
                        header = fn_readxml(fk);
                        readxmlheader(T{k},header)
                    case '2plsm'
                        T{k}.file = fk;
                        read2plsmheaderanddata(T{k},fk)
                    case {'blk' 'vdq'}
                        T{k}.file = fk;
                        header = oi_headBLK(fk);
                        readblkheader(T{k},header)
                    case 'oldvdaq'
                        T{k}.file = fk;
                        % no header information implemented yet
                        defaultheadervalues(T{k})
                    case 'da'
                        % Neuroplex
                        T{k}.file = fk;
                        header = neuroplex_read(fk);
                        readneuroplexheader(T{k},header);
                    case 'pcl'
                        % ScienceWares Photon Imaging System
                        % read both header and data
                        T{k}.file = fk;
                        readplcheaderanddata(T{k},fk)
                    case {'bmp' 'png' 'jpg' 'tif'} %#ok<MDUPC>
                        T{k}.file = fk;
                        % no header information -> some default values
                        defaultheadervalues(T{k})
                        % we need to read the data to set the proper size
                        readdata(T{k})
                    case 'mes'
                        T{k} = tps_trial.openmesfile(fk,flag);
                    case 'mesc'
                        T{k} = tps_trial.openmescfile(fk,flag);
                    case 'mat'
                        T{k}.file = fk;
                        defaultheadervalues(T{k})
                        F = whos('-file',fk);
                        if isscalar(F)
                            T{k}.sizes0 = F.size;
                        else
                            kdata = find(strcmp({F.name},'data'),1);
                            if isempty(kdata)
                                error 'cannot determine which variable inside the .mat file contains the data'
                            end
                            T{k}.sizes0 = F(kdata).size;
                            ksfr = find(strcmp({F.name},'data'),1);
                            T{k}.sfrchannel = logical(prod(F(ksfr).size));
                        end
                    otherwise
                        error('cannot read data with extension ''%s''',ext)
                end
            end
            T = [T{:}];
            if isempty(T), return, end % happens for example when opening MES/MESc data with empty sub-trials selection
            T.mergestimtables
            % all trials are single images of the same size? 
            if ntrial>1 && all([T.nfr]==1) && ~any(row(diff(cat(1,T.sizes)))) && isempty(flag)
                answer = questdlg('All files appear to contain a single image, all of the same size. Do you want to concatenate them to form a single movie?', ...
                    'tps_trial','Yes','No','Yes');
                if answer == 0
                    % interrupted
                    T = tps_trial.empty(1,0);
                    return
                elseif strcmp(answer,'Yes')
                    flag = 'mim';
                end
            elseif strcmp(flag,'img')
                % flag 'img' only serves to avoid user prompting above, and
                % can be forgotten now
                flag = '';
            end
            % flag: average or concatenate
            if ~fn_ismemberstr(flag,{'avg','cat','mim'}), return, end
            TT = T;
            T = TT(1);
            T.erasedata()  % data from TT(1) might have been read, but is not the full data from all TT
            same = struct('required',true,'nfr',true,'addinfo',true,'fullinfo',true);
            for k=2:ntrial
                same.required = same.required && ...
                    isequal({T.sizes([1 2 4]),T.origin,T.type,T.scanning,T.dx,T.dy,T.dz,T.dt}, ...
                    {TT(k).sizes([1 2 4]),TT(k).origin,TT(k).type,TT(k).scanning,TT(k).dx,TT(k).dy,TT(k).dz,TT(k).dt});
                same.nfr = same.nfr && T.nfr==TT(k).nfr;
                same.addinfo = same.addinfo && isequal(T.addinfo,TT(k).addinfo);
                same.fullinfo = same.fullinfo && isequal(T.fullinfo,TT(k).fullinfo);
            end
            if ~same.required || (strcmp(flag,'avg') && ~same.nfr)
                error('trials not compatible')
            end
            if ~same.addinfo, T.addinfo = struct; end
            if ~same.fullinfo, T.fullinfo = struct; end
            if strcmp(flag,'mim')
                T.file = f;
                T.preprocess.multifile = 'cat';
                T.sizes0(3) = nf;
            else
                T.file = char(TT.file);
                T.preprocess.multifile = flag;
                if strcmp(flag,'cat'), T.sizes0(3) = sum([TT.nfr]); end
            end
        end
    end
    methods (Static, Access='private')
        function T = openmesfile(fmes,flag)
            % Get the headers
            [header sessions] = mes_header(fmes);
            nmes = length(header);
            
            % Load the header info
            if ischar(flag) || isstring(flag)
                % Prompt user for a sub-selection?
                if strcmp(flag,'prompt') && nmes>1
                    subchoice = listdlg('PromptString','Select experiments to open', ...
                        'ListString',sessions,'ListSize',[350 400]);
                    if isempty(subchoice), T = tps_trial.empty(1,0); return, end
                    header = header(subchoice);
                    nmes = length(header);
                end
            else
                % Load only trials specified by user
                subchoice = flag;
                nmes = length(subchoice);
                header = header(subchoice);
            end
            
            % Create the tps_trial object, with a few basic features
            T = tps_trial; erasedata(T)
            uniquestimtable = T.stimtable;
            for k=1:nmes
                T(k).file = fmes;
                erasedata(T(k))
                T(k).stimtable = uniquestimtable;
            end
            
            % Read the headers thoroughly
            for k=1:nmes
                readmesheader(T(k),header{k})
            end
        end
        function T = openmescfile(fmesc,flag)
            % Load the header info
            [h desc] = mesc_header(fmesc);
            nmes = size(desc,1);
            if ischar(flag) || isstring(flag)
                % Prompt user for a sub-selection?
                if strcmp(flag,'prompt')
                    subchoice = listdlg('PromptString','Select experiments to open', ...
                        'ListString',desc,'ListSize',[350 400]);
                    if isempty(subchoice), T = tps_trial.empty(1,0); return, end
                    nmes = length(subchoice);
                else
                    subchoice = 1:nmes;
                end
            else
                subchoice = flag;
                nmes = length(subchoice);
            end
            
            % Create the tps_trial object, with a few basic features
            T = tps_trial; erasedata(T)
            uniquestimtable = T.stimtable;
            for k=1:nmes
                T(k).file = fmesc;
                erasedata(T(k))
                T(k).stimtable = uniquestimtable;
            end
            
            % Read the headers thoroughly
            kall = 0;
            k = 0;
            header = struct('MESC',rmfield(h,'sessions'),'session',[],'unit',[]);
            for i=1:length(h.sessions)
                hi = h.sessions(i);
                header.session = rmfield(hi,'Units');
                for j=1:length(hi.Units)
                    kall = kall+1;
                    if ~any(subchoice==kall), continue, end
                    k = k+1;
                    hij = hi.Units(j);
                    header.unit = hij;
                    try
                        readmescheader(T(k),header)
                    catch
                        fprintf('failed reading header of %s, session %i, unit %i\n',fmesc,i,j)
                    end
                    T(k).fullinfo.ksession = i;
                    T(k).fullinfo.kunit = j;
                end
            end
        end
    end
    methods (Access='private')
        function defaultheadervalues(T)
            T.type = 'movie';
            T.scanning = false;
            T.dx = 1;
            T.dy = 1;
            T.dz = 0;
            T.dt = 1;
            T.t0 = 0;
            T.xunit = 'px';
            T.tunit = 'frame';
            T.stimtable = tps_stimtable;
            for k=1:T.sizes0(4)
                T.stimid(k) = addstim(T.stimtable,k-1);
            end
        end
        function readcfdheader(T,header)
            % origin, fullinfo, type
            T.origin = 'cfd';
            T.fullinfo = header;
            if header.ACV.ZStepSize
                T.type = 'zstack';
            elseif header.ACV.flags.lineScan
                T.type = 'linescan';
            else
                T.type = 'movie';
            end
            T.addinfo.time = header.CFH.time;
            
            %             % Acquisition settings
            %             xmlfile = strrep(T.file,'.CFD','.xml');
            %             if exist(xmlfile,'file')
            %                 AcqSet = tps_readacqsettingsxml(xmlfile);
            %             end
            %             if ~isempty(AcqSet)
            %                 T.fullinfo.AcqSet = AcqSet;
            %                 [T.addinfo objective T.recording_header] = readacqsetheader(T,AcqSet);
            %             else
            %                 disp('attention, objective supposed to be 40x')
            %                 objective = 40;
            %             end
            
            % zoom
            T.addinfo.zoom = header.ASV.zoom_fac;
            tmp = T.fullinfo.CFH.zpos/1000;
            if tmp>0 && tmp<500 % looks like z=0 was on the surface
                T.addinfo.depth = sprintf('%.0fum',tmp);
            end
            nx = header.ACV.pixels.x - 3; % 3 first columns are for the channel data
            ny = header.ACV.pixels.y;
            nfr = header.ACV.num.images;
            nc = 1;
            T.sizes0 =  [nx ny nfr nc];
            % (VERY STRANGE: 1 means 1 channel, and 3 means 2 channels!? -
            % see cfd_read)
            T.sfrchannel = (header.ACV.num.chans==3);
            %             if datenum(header.CFH.date) < datenum('22 July 2008')
            %                 % (OLD SMALL SCANNING MIRRORS)
            %                 % (size: field size with objective 40X, zoom .25, scan
            %                 % amplitude 1V is 460um)
            %                 sizepervolt = 460/(objective/40)/(T.addinfo.zoom/.25);
            %             else
            %                 % (size: field size with objective 40X, zoom .25, scan
            %                 % amplitude 1V is 272um)
            %                 sizepervolt = 272/(objective/40)/(T.addinfo.zoom/.25);
            %             end
            % (scan amplitude in ASV is 2^15/10 times the value in volts)
            T.dx = 1; %sizepervolt/header.ACV.pixels.x * header.ASV.range.x/6553.6;
            T.dy = 1; %sizepervolt/header.ACV.pixels.y * header.ASV.range.y/6553.6;
            T.dz = 1; %header.ACV.ZStepSize/1000;
            T.xunit = 'px';
            T.scanning = true;
            linedur = (header.ACV.scandur + header.ACV.retracedur) / 1000; % sec
            T.dt = ny*linedur;
            %T.addinfo.dwelltime = header.ACV.scandur / nx * 1000; % microsec
            T.t0 = 0;
            T.tunit = 's';
        end
        function readscanimageheader(T,header)
            T.origin = 'ScanImage';
            T.fullinfo = header;
            acq = header.acq;
            if acq.numberOfZSlices>1
                T.type = 'zstack';
                nfr = acq.numberOfZSlices;
            elseif acq.linescan
                T.type = 'linescan';
                nfr = acq.numberOfFrames;
            else
                T.type = 'movie';
                nfr = acq.numberOfFrames;
            end
            T.addinfo.time = header.internal.triggerTimeString;
            
            %             % Acquisition settings
            %             xmlfile = strrep(T.file,'.tif','.xml');
            %             try
            %                 AcqSet = tps_readacqsettingsxml(xmlfile);
            %                 T.fullinfo.AcqSet = AcqSet;
            %                 [T.addinfo objective T.recording_header] = readacqsetheader(T,AcqSet);
            %             catch %#ok<CTCH>
            %                 disp('no xml header found, assume objective 40x')
            %                 objective = 40;
            %             end
            
            disp('update code for new ScanImage')
            %             T.addinfo.xscanamp = acq.scanAmplitudeX / acq.zoomFactor;
            nx = acq.pixelsPerLine;
            ny = acq.linesPerFrame;
            nc = 1;
            T.sizes0 =  [nx ny nfr nc];
            T.sfrchannel = logical(acq.acquiringChannel2);
            %             zoom = 100*acq.zoomhundreds + 10*acq.zoomtens + acq.zoomones;
            %             switch header.software.version
            %                 case 3.6
            %                     % size: field size with objective 40X, zoom 1, scan amplitude
            %                     % 4V is 272um
            %                     sizepervolt = 272/(objective/40)/zoom/4;
            %                     T.dx = sizepervolt/nx * acq.scanAmplitudeX;
            %                     T.dy = sizepervolt/ny * acq.scanAmplitudeY;
            %                 case 3.7
            %                     % this time, measured with higher precision:
            %                     % with objective 40X, 19.5um per
            %                     % X optical degree, 25.9375um per Y optical degree, and
            %                     % in standard.ini, voltsPerOpticalDegree=0.667
            %                     %-
            %                     % note that this value of 0.667 leads in reality to
            %                     % 20 degrees <-> 8V, which is rather 0.4 volts per
            %                     % degree!, strange...
            %                     % note then that 272um for 4V corresponds to
            %                     % 272/4*0.4 = 27.2um per degree, which is coherent
            %                     sizeperdegree = [19.5 25.9375]/(objective/40)/zoom;
            %                     if acq.fastScanningX
            %                         [degx degy] = deal(acq.scanAngularRangeFast,acq.scanAngularRangeSlow);
            %                     else
            %                         [degx degy] = deal(acq.scanAngularRangeSlow,acq.scanAngularRangeFast);
            %                     end
            %                     T.dx = sizeperdegree(1)/nx * degx;
            %                     T.dy = sizeperdegree(1)/ny * degy; disp('measure scale again')
            %                 otherwise
            %                     error('unknown ScanImage version')
            %             end
            T.dx = 1;
            T.dy = 1;
            T.dz = 1; %acq.zStepSize; % unit microns
            T.xunit = 'px'; %'um';
            T.scanning = true;
            linedur = acq.msPerLine/1000; % sec
            T.dt = ny*linedur;
            T.t0 = 0;
            T.tunit = 's';
        end
        function readmpdheader(T,header)
            T.origin = 'mpscope';
            T.fullinfo.MPD = header;
            T.type = 'movie';
            %T.addinfo.time = header.CFH.time;
            
            %             % Acquisition settings
            %             xmlfile = strrep(T.file,'.CFD','.xml');
            %             xmlfile = strrep(xmlfile,'.cfd','.xml');
            %             xmlfile = strrep(xmlfile,'.MPD','.xml');
            %             if exist(xmlfile,'file')
            %                 AcqSet = tps_readacqsettingsxml(xmlfile);
            %                 T.fullinfo.AcqSet = AcqSet;
            %                 [T.addinfo objective T.recording_header] = readacqsetheader(T,AcqSet); %#ok<ASGLU>
            %             else
            %                 %disp('attention, objective supposed to be 40x, setup supposed to be > ')
            %                 %objective = 40;
            %             end
            
            %T.addinfo.zoom = header.ASV.zoom_fac;
            %             tmp = T.fullinfo.CFH.zpos/1000;
            %             if tmp>0 && tmp<500 % looks like z=0 was on the surface
            %                 T.addinfo.depth = sprintf('%.0fum',tmp);
            %             end
            nx = header.nx;
            ny = header.ny;
            nfr = header.nfr;
            nc = 1;
            T.sizes0 =  [nx ny nfr nc];
            T.sfrchannel = true;
            %defaultfieldsize = 2720/objective/T.addinfo.zoom;
            T.dx = 1;
            T.dy = 1;
            T.dz = 0;
            T.scanning = true;
            T.dt = 1;
            T.t0 = 0;
        end
        function readxmlheader(T,header)
            T.fullinfo = header;
            frames = header.Sequence.Frame;
            keys = frames(1).PVStateShard.Key;
            info = [];
            for k=1:length(keys)
                info.(keys(k).key) = keys(k).value;
            end
            if length(frames)>1
                keys2 = frames(2).PVStateShard.Key;
                info2 = [];
                for k=1:length(keys2)
                    info2.(keys2(k).key) = keys2(k).value;
                end
            end
            
            T.origin = 'xml';
            switch header.Sequence.type
                case 'ZSeries'
                    T.type = 'zstack';
                otherwise
                    T.type = 'movie';
            end
            disp(['attention, make sure that objective was really ' info.objectiveLens])
            T.addinfo.objective = info.objectiveLens;
            T.addinfo.zoom = info.opticalZoom;
            nx = info.pixelsPerLine;
            ny = info.linesPerFrame;
            nfr = length(frames);
            nc = 1;
            T.sizes0 =  [nx ny nfr nc];
            T.sfrchannel = (length(frames(1).File)>1);
            T.dx = info.micronsPerPixel_XAxis;
            T.dy = T.dx;
            if strcmp(T.type,'zstack') && T.nfr>1
                T.dz = info2.positionCurrent_ZAxis - info.positionCurrent_ZAxis;
            else
                T.dz = 0;
            end
            T.dt = info.framePeriod;
            T.t0 = 0;
            T.scanning = true;
            if T.linedur~=info.scanlinePeriod, error('problem with scan line period'), end
            T.addinfo.dwelltime = info.dwellTime;
        end
        function read2plsmheaderanddata(T,fname) %#ok<MANU,INUSD>
            error('not handled any more - please revise code')
            %             % read file
            %             [desc dat sf times stims stimunext] = tpm_2plsmread(fname);
            %             header = desc(1);
            %
            %             % data and size - TODO: handle multiple packets
            %             T.addinfo.tmpdata.data = dat;
            %             T.addinfo.tmpdata.sfr  = sf;
            %             siz = size(T.addinfo.tmpdata.data);
            %             siz(end+1:4) = 1;
            %             T.sizes0 =  siz;
            %             if T.nfr~=numel(desc) || T.nc~=1, error('size problem'), end
            %
            %             % header
            %             T.fullinfo = header;
            %             T.origin = '2plsm';
            %             T.type = 'movie'; % TODO: add possibility for z-stack
            %             T.addinfo.objective = header.Objective;
            %             % corresponding number of pixels
            %             pdim1 = round(3*header.Xres/4);
            %             pdim2 = header.Yres;
            %             T.sfrchannel = true; % TODO: add sulforhodamine channel
            %             % real dimensions of image in um without flyback
            %             rdim1 = 100*(header.Right-header.Left);
            %             rdim2 = 100*(header.Top-header.Bottom);
            %             % size of pixel in um
            %             T.dx = rdim1/pdim1;
            %             T.dy = rdim2/pdim2;
            %             if strcmp(T.type,'zstack') && T.nfr>1
            %                 T.dz = '?'; % TODO: add possibility for z-stack
            %             else
            %                 T.dz = 0;
            %             end
            %             T.dt = header.Framerate;
            %             T.t0 = 0;
            %             T.scanning = true;
            %             % additional information to display - TODO: edit this according
            %             % to needs
            %             T.addinfo.PockelV = header.PockelV;
            %             T.addinfo.LaserPower = header.LaserPower;
            %             % stimulation - TODO
            %             T.addinfo.stim1 = struct('times',times,'stims',stims,'stimunex',stimunext);
        end
        function readblkheader(T,header)
            T.origin = 'BLK';
            T.fullinfo = header;
            T.type = 'movie';
            
            % Acquisition settings
            nx = header.framewidth;
            ny = header.frameheight;
            nfr = header.nframesperstim;
            nc = header.nstimuli;
            T.sizes0 =  [nx ny nfr nc];
            T.sfrchannel = false;
            
            % at 4X, image 512x512, 1mm=90pixels
            T.dx = 1; %(1000/90)*(camerabin/2)*header.xbinfactor*(4/objective);
            T.dy = T.dx;
            T.dz = 0;
            T.xunit = 'px';
            %             T.dt = header.nvideoframesperdataframe / camerafs;
            T.dt = 1;
            T.t0 = 0;
            T.scanning = false;
            %             T.tunit = 's';
            %             if T.dt==0, T.dt = 1; T.tunit = 'frame'; end
            T.tunit = 'frame';
            
            % Stim number
            if nc==1
                StimNumFromBLKname(T)
            else
                setstim(T,0:nc-1)
            end
            
            % More info from file name
            MoreInfoFromBLKname(T)
        end
        function readneuroplexheader(T,header)
            T.origin = 'NeuroPlex';
            T.fullinfo = header;
            T.type = 'movie';
            
            % Acquisition settings
            nx = header.NColumns;
            ny = header.NRows;
            nfr = header.NFrames;
            T.sizes0 =  [nx ny nfr 1];
            T.sfrchannel = false;
            
            % No spatial information
            T.dx = 1;
            T.dy = T.dx;
            T.dz = 0;
            T.xunit = 'px';

            % Temporal information
            T.dt = header.FrameIntervalMS;
            T.t0 = 0;
            T.scanning = false;
            T.tunit = 'ms';
            
            % Recordings
            % 8 channels are recorded regardless of whether they are used
            % or not, no header information allows to know which are
            % actually used, so we just read all of them
            rec = T.recording_header;
            for i = 1:8
                rec(i).name = sprintf('channel%i',i);
                rec(i).t0 = 0;
                rec(i).dt = header.FrameIntervalMS/1000/header.AcquisitionRatio; % unit: second
                rec(i).n = nfr*header.AcquisitionRatio;
            end
            T.recording_header = rec;
        end
        function readmesheader(T,header)
            T.origin = 'MES';
            T.fullinfo.MES = header;
            h1 = header(1);
            % Description
            if ~isempty(h1.Comment)
                T.addinfo.comment = h1.Comment;
            end
            % Channels
            for k = 1:length(header)
                if k+1>length(header) || ~isempty(header(k+1).Type), break, end
            end
            nchannel = k;
            channels = {header(1:nchannel).Channel};
            T.fullinfo.channels.nchannel = nchannel;
            T.fullinfo.channels.kdata = find(ismember(channels,{'Camera' 'pmtUG' 'pmtLG' 'pmtLGraw' 'PMTg'}));
            T.fullinfo.channels.ksfr = find(ismember(channels,{'pmtUR' 'pmtLR' 'PMTr'}));
            if isempty(T.fullinfo.channels.ksfr), T.fullinfo.channels.ksfr = find(ismember(channels,{'tIR' 'tIRraw'})); end
            if isempty(T.fullinfo.channels.kdata), error 'no green channel data', end
            usedchannels = [T.fullinfo.channels.kdata T.fullinfo.channels.ksfr];
            if length(usedchannels)<nchannel, disp('some channels were not used:'), disp(channels(setdiff(1:nchannel,usedchannels))), end
            T.sfrchannel = ~isempty(T.fullinfo.channels.ksfr);
            % Type and size
            T.scanning = true;
            switch h1.Context
                case {'Cam' 'Photo'}
                    T.type = 'movie';
                    ny = h1.Height;
                    nfr = 1;
                case 'Zstack';
                    T.type = 'zstack';
                    ny = h1.Height;
                    nfr = length(header)/nchannel;
                case 'Measure'
                    switch h1.Type
                        case 'FF'
                            T.type = 'movie';
                            info = h1.FoldedFrameInfo;
                            ny = info.numFrameLines;
                            %                             if h1.creatingMEScRevision < 1425
                            %                                 nfrskip = floor(h1.FoldedFrameInfo.firstFrameStartTime/Tk.dt_sec);
                            %                             else
                            % nfrskip = floor(info.firstFrameStartTime/info.frameTimeLength*1e3);
                            nfr = info.numFrames; % - nfrskip;
                        otherwise
                            T.type = 'linescan';
                            ny = 1;
                            nfr = h1.Height;
                    end
            end
            % % folded frame: there used to be useful extra data
            % % corresponding to return fly in an older version of MES
            % % but this bug was corrected at least at revision 1425 
            % nx = floor(h1.Width/1.6);
            nx = h1.Width; 
            nc = 1;
            T.sizes0 =  [nx ny nfr nc];
            % Spatial resolution
            T.dx = h1.WidthStep;
            if strcmp(h1.Type,'FF')
                T.dy = h1.FoldedFrameInfo.TransverseStep;
            else
                T.dy = h1.HeightStep;
            end
            T.xunit = 'um';
            if strcmp(T.type,'zstack') 
                if isfield(h1,'info_Posinfo')
                    T.dz = h1.info_Posinfo.zstep;
                else
                    disp 'could not read zstack step, use default of 1 micron'
                    T.dz = 1;
                end
            else
                T.dz = 0; 
            end
            % Temporal resolution: this seems difficult!
            if strcmp(T.type,'linescan')
                protocol = h1.info_Protocol.protocol.d;
                protocol_curves = protocol(strcmp({protocol.name},'Curves'));
                tmin = protocol_curves.data(2,2);
                tmax = protocol_curves.data(2,3);
                T.dt = (tmax-tmin)/nfr;
                % remove first, dark, frames (first 300ms)
                frrem = round(300/T.dt);
                T.sizes0(3) = nfr-frrem;
                T.fullinfo.channels.frrem = frrem;
                T.t0 = tmin + frrem*T.dt;
                T.tunit = 'ms';
            elseif strcmp(h1.Type,'FF')
                T.dt = h1.FoldedFrameInfo.frameTimeLength;
                T.t0 = h1.FoldedFrameInfo.firstFrameStartTime;
                T.tunit = 'ms';
            else
                T.dt = 1;
                T.tunit = 'frame';
            end
            % Stimulation
            if strcmp(T.type,'linescan')
                stimprotocol = h1.info_Protocol.protocol.d(7);
                if stimprotocol.func
                    stimdata = stimprotocol.data;
                    f = find(stimdata(1,:)~=0);
                    if stimdata(1,end)~=0 || ~all(stimdata(1,f+1)==0)
                        disp 'problem with stim definition'
                        f(stimdata(1,f-1)~=0) = [];
                    end
                    d = diff(stimdata(2,:));
                    T.stim = [stimdata(2,f) d(f)]/1e3;
                end
            end
            % Electrophysiology (imported ABF)
            h1names = fieldnames(h1);
            recnames = h1names(~fn_isemptyc(strfind(h1names,'ABF__')));
            nrec = length(recnames);
            if nrec
                rec = T.recording_header;
                for i = 1:nrec
                    ri = h1.(recnames{i});
                    if ~strcmp(ri.xunit,'ms'), error 'xunit must be ms', end
                    rec(i).name = strrep(ri.yname,' ','');
                    rec(i).t0 = ri.x(1)/1e3;
                    rec(i).dt = ri.x(2)/1e3;
                    rec(i).n = length(ri.y);
                end
                T.recording_header = rec;
                if isfield(h1,'FileSubindex')
                    varname = num2str(h1.FileSubindex,'DF%.4i');
                else
                    disp 'FileSubindex not defined in h1, please define it'
                    keyboard
                    varname = num2str(h1.FileSubindex,'Df%.4i');
                end
                T.analogfile = {[T.file varname recnames']}; % information needed to load the recording data
            end
            % Laser amplitude in addinfo
            for k=1:length(header)
                hk = header(k);
                if isfield(hk,'info_A0settings') && ~isempty(hk.info_AOsettings)
                    T.addinfo.laser = hk.info_AOsettings.uAO1y.RAMamps(2);
                    break
                end
            end
        end
        function readmescheader(T,header)
            T.origin = 'MESC';
            T.fullinfo = header;
            h = header.unit;
            % util subfunction
            function str = makestr(x)
                str = char(row(x(1:end-1)));
            end
            % Description
            if ~isempty(h.Comment)
                T.addinfo.comment = makestr(h.Comment);
            end
            % Channels
            if ~strcmp(makestr(h.Channel_0_Name),'UG'), error 'first channel is not UG', end
            if isfield(h,'Channel_1_Name') && ~isempty(h.Channel_1_Name) && ~strcmp(makestr(h.Channel_1_Name),'UR'), error 'second channel is not UR', end
            if isfield(h,'Channel_2_Name') && ~isempty(h.Channel_2_Name), error 'more than 2 channels', end
            T.sfrchannel = isfield(h,'Channel_1_Name') && ~isempty(h.Channel_1_Name);
            % Type and size
            T.scanning = true;
            T.type = 'movie';
            nx = double(h.XDim);
            ny = double(h.YDim);
            nfr = double(h.ZDim);
            nc = 1;
            T.sizes0 =  [nx ny nfr nc];
            % Spatial resolution
            T.dx = h.XAxisConversionConversionLinearScale;
            T.dy = h.YAxisConversionConversionLinearScale;
            T.xunit = 'um';
            % Temporal resolution
            T.dt = h.ZAxisConversionConversionLinearScale;
            T.tunit = 'ms';
            % Stimulation - not defined yet
            % Acquisition parameters
            s = fn_readxml(h.MeasurementParamsXML');
            s = s.Params.param.ATTRIBUTE;
            acq = struct;
            for i=1:length(s)
                acq.(s(i).name) = s(i).value;
            end
            T.fullinfo.acq = acq;
        end
        function readplcheaderanddata(T,fname)
            [data param] = pcl_read(fname);
            T.fullinfo.movie_param = param;
            T.type = 'movie';
            T.scanning = false;
            T.dx = param.xbin;
            T.dy = param.xbin;
            T.dz = 0;
            T.dt = param.dt;
            T.t0 = param.tstart;
            T.xunit = 'px';
            T.tunit = 's';
            T.stimtable = tps_stimtable;
            for k=1:T.sizes0(4)
                T.stimid(k) = addstim(T.stimtable,k-1);
            end
            T.setdata(data)
        end
    end
    methods
        function locatefile(T,fexample)
            % check if data files have been moved
            fil = deblank(T(1).file(1,:));
            if exist(fil,'file') && ~isempty(fileparts(fil)), return, end
            if nargin<2, fexample=''; end
            fileold = {T.file};
            filenew = tps_locatefile(fileold,fexample);
            if isequal(filenew,fileold), return, end % locating file seems to have failed
            [T.file] = deal(filenew{:});
        end
        function StimNumFromBLKname(T)
            ok = true;
            for k=1:length(T)
                fk = fn_fileparts(T(k).file,'base');
                tokens = regexp(fk,'C(\d{1,5})(\d{6})*[^\d]','tokens'); % remove date and time
                if isempty(tokens)
                    ok = false;
                    setstim(T(k),0)
                    continue
                end
                tok = tokens{1}{1};
                setstim(T(k),str2double(tok));
            end
            if ~ok, disp 'could not determine trial condition', end
        end
        function MoreInfoFromBLKname(T)
            for k=1:length(T)
                T(k).fullinfo = fn_structmerge(T(k).fullinfo,oi_infoBLKname(T(k).file));
            end
        end
    end
    
    % Read data
    methods        
        function readdata(T)
            if ~T(1).internpar.readdata
                disp 'data reading is disabled'
                return
            end
            if length(T)>1 || size(T.file,1)<=1
                % trials are 'regular'
                isregular = true;
                nf = length(T);
            else
                % there is a single trial obtained from multiple files
                % through averaging or concatenation (as indicated by flag
                % T.preprocess.multifile)
                isregular = false;
                nf = size(T.file,1);
                switch T.preprocess.multifile
                    case 'avg'
                        datatemp = zeros(T.internpar.floattype);
                        if T.sfrchannel, sfrtemp = zeros(T.internpar.floattype); else sfrtemp = []; end
                    case 'cat'
                        datatemp = [];
                        sfrtemp = [];
                    otherwise
                        error('unexpected multifile flag ''%s''', T.preprocess.multifile)
                end
                if strcmp(T.preprocess.multifile,'cat'), curfr = 0; end 
            end
            % check that files exist
            newpath = [];
            if isregular, files = {T.file}; else files = cellstr(T.file); end
            for k=1:nf
                fk = files{k};
                if ~isempty(fk) && ~exist(fk,'file')
                    error('data file not found')
                end                
            end
            if ~isempty(newpath) % at least some files path changed
                if isregular, [T.file] = deal(files{:}); else T.file = char(files); end
            end
            % read files
            if nf>1, fn_progress('reading file',nf), end
            for k=1:nf
                if nf>1, fn_progress(k), end
                
                % Read data
                kT = fn_switch(isscalar(T),1,k);
                Tk = T(kT);
                fk = files{k};
                if isempty(fk)
                    ext = 'unsaved_data';
                else
                    [dum1, dum2, ext] = fileparts(fk); %#ok<ASGLU>
                    ext = lower(ext(2:end));
                    if length(ext)==1, ext='oldvdaq'; elseif strfind(ext,'blk'), ext='blk'; end
                    if strcmp(ext,'tif') && strcmp(Tk.origin,'ScanImage'), ext='scanimage'; end
                end
                sf = [];
                switch ext
                    case 'unsaved_data'
                        dat = Tk.data_unsaved{1};
                        sf = Tk.data_unsaved{2};
                    case 'mat'
                        v = load(fk);
                        F = fieldnames(v);
                        if isscalar(F)
                            dat = v.(F{1});
                        else
                            dat = v.data;
                            if isfield(v,'sfr'), sf = v.sfr; end
                        end
                    case 'avi'
                        VR = VideoReader(fk); %#ok<TNMLP>
                        dat = VR.read();
                        if ndims(dat)==4, dat = squeeze(dat(:,:,1,:)); end  
                        dat = permute(dat,[2 1 3]);
                    case 'cfd'
                        dat = cfd_read(fk);
                        if Tk.nx<size(dat,1)
                            % remove channel data only if headers from the
                            % new version
                            dat = dat(4:end,:,:,:);
                        end
                        if Tk.sfrchannel, sf = dat(:,:,:,2); end
                        dat = dat(:,:,:,1);
                    case 'mpd'
                        [dat sf] = mpd_read(fk,Tk.fullinfo.MPD);                        
                    case 'scanimage'
                        switch Tk.fullinfo.software.version
                            case 3.6
                                % Thomas: i had modified the saving code to
                                % improve speed, so it is also a different
                                % code to read the data!
                                s = Tk.sizes([2 1 3]);
                                switch Tk.type
                                    case {'movie' 'linescan'}
                                        dat = imread(fk);
                                        dat = reshape(dat,[s numel(dat)/prod(s)]);
                                        dat = permute(dat,[2 1 3 4]);
                                        if Tk.sfrchannel
                                            sf  = dat(:,:,:,2);
                                        end
                                        dat = dat(:,:,:,1);
                                    case 'zstack'
                                        s1 = [s(1:2) (1+Tk.sfrchannel)];
                                        nz = s(3);
                                        dat = zeros(prod(s1),nz,'uint16');
                                        for i=1:nz, dat(:,i) = imread(fk,i); end
                                        dat = reshape(dat,[s1 nz]);
                                        if Tk.sfrchannel
                                            sf  = permute(dat(:,:,2,:),[2 1 4 3]);
                                        end
                                        dat = permute(dat(:,:,1,:),[2 1 4 3]);
                                    otherwise
                                        error('unknown type ''%s''',Tk.type)
                                end
                            case 3.7
                                % In release 3.7, data writing has been
                                % improved, so no need to do the same
                                % modification
                                %                                 switch Tk.type
                                %                                     case 'zstack'
                                %                                         nslices = Tk.nfr;
                                %                                         nfr = 1;
                                %                                     otherwise
                                %                                         nslices = 1;
                                %                                         nfr = Tk.nfr;
                                %                                 end
                                curfr = Tk.nfr;
                                nc = 1+Tk.sfrchannel;
                                s = Tk.sizes([2 1]);
                                dat = zeros([s nc curfr],'uint16');
                                try
                                    for i=1:curfr
                                        for j=1:nc
                                            dat(:,:,j,i) = imread(fk,(i-1)*nc+j);
                                        end
                                    end
                                catch %#ok<CTCH>
                                    if i==1 || j~=1, error('check this'), end
                                    curfr = i-1;
                                    Tk.nfr = curfr;
                                    dat = dat(:,:,:,1:curfr);
                                end
                                dat = permute(dat,[2 1 4 3]);
                                if Tk.sfrchannel
                                    sf = dat(:,:,:,2);
                                end
                                dat = dat(:,:,:,1);
                            otherwise
                                error('unknown ScanImage version')
                        end
                    case '2plsm'
                        % data has already been loaded while reading the
                        % headers and was temporarly stored in T.addinfo
                        dat = Tk.addinfo.tmpdata.data;
                        sf  = Tk.addinfo.tmpdata.sfr;
                        Tk.addinfo = rmfield(Tk.addinfo,'tmpdata');
                    case 'mes'
                        h1 = Tk.fullinfo.MES(1);
                        nmes = fn_switch(Tk.type,'zstack',Tk.nfr,1); % this is the not of MES structures
                        nchannel = Tk.fullinfo.channels.nchannel;
                        kdata = Tk.fullinfo.channels.kdata;
                        kim = fn_add(kdata(:),(0:nmes-1)*nchannel); kim = kim(:);
                        datvar = {Tk.fullinfo.MES(kim).IMAGE};
                        dat = struct2cell(load(Tk.file,datvar{:},'-MAT'));
                        dat = cat(3,dat{:});
                        if ~isscalar(kdata), dat = fn_bin(dat,[1 1 length(kdata)]); end
                        if strcmp(Tk.type,'linescan')
                            dat(:,1:Tk.fullinfo.channels.frrem,:) = [];
                        elseif strcmp(h1.Type,'FF')
                            info = h1.FoldedFrameInfo;
                            %                             if h1.creatingMEScRevision >= 1425
                            % nfrskip = floor(h1.FoldedFrameInfo.firstFrameStartTime/Tk.dt_sec);
                            %                             else
                            %                             nfrskip = h1.FoldedFrameInfo.firstFramePos-1;
                            %                             end
                            fprintf('Is it correct that first "real" frame occurs at time %g?\n',info.firstFrameStartTime)
                            if size(dat,1) ~= Tk.nx
                                error 'horizontal size of data is not as expected'
                            end
                            if floor(size(dat,2)/Tk.ny) == Tk.nfr + 1
                                warning 'folded frame has one more image than expected, the last one will be removed'
                            elseif floor(size(dat,2)/Tk.ny) ~= Tk.nfr
                                error 'vertical size of data is not as expected'
                            end
                            dat = dat(1:Tk.nx,1:Tk.ny*Tk.nfr); % Tk.ny*nfrskip+(1:Tk.ny*Tk.nfr));
                        end
                        dat = reshape(dat,Tk.sizes0);
                        ksfr = Tk.fullinfo.channels.ksfr;
                        if ~isempty(ksfr)
                            kim = fn_add(ksfr(:),(0:nmes-1)*nchannel); kim = kim(:);
                            sfvar = {Tk.fullinfo.MES(kim).IMAGE};
                            sf = struct2cell(load(Tk.file,sfvar{:},'-MAT'));
                            sf = cat(3,sf{:});
                            % if h1.creatingMEScRevision>=1000, sf = fliplr(sf); end % strange, but this bug was corrected at least at revision 1425 
                            if ~isscalar(ksfr), sf = fn_bin(sf,[1 1 length(ksfr)]); end
                            if strcmp(Tk.type,'linescan')
                                sf(:,1:Tk.fullinfo.channels.frrem,:) = [];
                            elseif strcmp(h1.Type,'FF')
                                sf = sf(1:Tk.nx,1:Tk.ny*Tk.nfr); % Tk.ny*nfrskip+(1:Tk.ny*Tk.nfr));
                            end
                            sf = reshape(sf,Tk.sizes0);
                        end
                    case 'mesc'
                        ksession = Tk.fullinfo.ksession;
                        kunit = Tk.fullinfo.kunit;
                        try 
                            dat = mesc_read(Tk.file,ksession-1,kunit-1,0);
                        catch
                            fprintf('failed reading data (session %i,unit %i)\n',ksession,kunit)
                            dat = [];
                        end
                        if Tk.sfrchannel
                            sf = mesc_read(Tk.file,ksession-1,kunit-1,1); 
                        end
                    case 'xml'
                        % get only the green channel
                        dat = zeros(T.nx,T.ny,T.nfr,1);
                        frame = T.fullinfo.Sequence.Frame(1);
                        if length(frame.File)>1
                            disp('assumed green channel was the second one')
                            channel = 2;
                        else
                            channel = 1;
                        end
                        % read
                        for i=1:T.nfr
                            frame = T.fullinfo.Sequence.Frame(i);
                            for j=1
                                fname = frame.File(channel).filename;
                                dat(:,:,i,j) = double(imread(fname));
                            end
                        end
                    case {'blk','vdq','oldvdaq'}
                        dat = oi_loadBLK(fk,'array');
                    case 'da'
                        [~, dat] = neuroplex_read(fk);
                    case 'pcl'
                        dat = pcl_read(fk,Tk.fullinfo.movie_param);
                    case {'bmp' 'png' 'jpg' 'tif' 'tiff'}
                        dat = fn_readimg(fk);
                    otherwise
                        error programming
                end
                                
                % Store data
                if isregular
                    if isempty(dat)
                        % failed reading... we don't need however to
                        % generate an error
                        continue
                    else
                        [dat, sf] = preprocessdata(dat, sf, Tk.preprocess.op);
                        setdata(Tk,dat,sf);
                    end
                else
                    if isempty(dat), error 'failed reading data', end
                    switch T.preprocess.multifile
                        case 'avg'
                            datatemp = datatemp + cast(dat,class(datatemp));
                            sfrtemp = sfrtemp + cast(sf,class(datatemp));
                        case 'cat'
                            nfr1 = size(dat,3);
                            if isempty(datatemp)
                                datatemp = zeros([size(dat,1) size(dat,2) T.nfr size(dat,4)], 'like', dat);
                            end
                            datatemp(:,:,curfr+(1:nfr1),:) = dat;
                            if T.sfrchannel, sfrtemp(:,:,curfr+(1:nfr1),:) = sf; end
                            curfr = curfr+nfr1;
                    end
                end
            end
            if ~isregular && strcmp(T.preprocess.multifile,'avg'), datatemp = datatemp/nf; sfrtemp = sfrtemp/nf; end
            if ~isregular
                [datatemp, sfrtemp] = preprocessdata(datatemp, sfrtemp, T.preprocess.op);
                setdata(T,datatemp,sfrtemp)
            end
        end
    end
    
    % Read data block
    methods
        function [dat, sf] = getdatablock(T,idx)
            % function [dat sf] = getdatablock(T,idx)
            %---
            % Read sub-block of data defined by the vector of frame indices
            % 'idx' (idx can go from 1 to the number of frames).
            
            % some checks
            if ~isscalar(T)
                error 'Method getdatablock applies only to a single tps_trial object'
            end
            if ~T(1).internpar.readdata
                error 'data reading is disabled'
            end
            dosfr = nargout>=2;
            
            % Full data?
            if isequal(idx, 1:T.nfr)
                dat = T.data;
                if dosfr, sf = T.sfr; end
                return
            end
            
            % Cut from full data?
            fk = T.file;
            [~, ~, ext] = fileparts(fk);
            ext = lower(ext(2:end));
            if ~ismember(ext,{'mesc'}) || (~isempty(T.basics.data) && ~isempty(T.basics.data.data))
                dat = T.data(:,:,idx,:);
                if dosfr
                    sf = T.sfr(:,:,idx,:);
                end
                return
            end
            
            % Get the raw data block
            if ~isempty(fk) && ~exist(fk,'file')
                error('data file not found')
            end
            switch ext
                case 'mesc'
                    ksession = T.fullinfo.ksession;
                    kunit = T.fullinfo.kunit;
                    try
                        dat = mesc_read(T.file,ksession-1,kunit-1,0,idx-1);
                    catch
                        fprintf('failed reading data (session %i,unit %i)\n',ksession,kunit)
                        dat = [];
                    end
                    if dosfr && T.sfrchannel
                        sf = mesc_read(T.file,ksession-1,kunit-1,1);
                    else
                        sf = [];
                    end
            end
               
            % Preprocess
            nop = length(T.preprocess.op);
            for kop = 1:nop
                opk = T.preprocess.op(kop);
                f = opk.name;
                switch f
                    case 'correctshift'
                        arg = opk.value; if ~iscell(arg), arg = {arg}; end
                        dat = tps_correctshift(dat,arg{:});
                        if ~isempty(sf), sf = tps_correctshift(sf,arg{:}); end
                    case 'translate'
                        shift = opk.value(:,idx);
                        dat = fn_translate(dat,shift,'full','linear');
                        if ~isempty(sf), sf = fn_translate(sf,shift,'full','linear'); end
                    case 'crop'
                        crop = opk.value;
                        dat = dat(crop{1},crop{2},:);
                        if ~isempty(sf), sf = sf(crop{1},crop{2},:); end
                    case 'user'
                        disp 'user-defined preprocessing might be hazardous when performed on sub-block of data'
                        fun = opk.value;
                        dat = fun(dat);
                        if ~isempty(sf), sf = fun(dat); end
                    otherwise
                        error('unknown preprocessing ''%s''',f)
                end
            end
            
            % Binning (TODO: make it part of the other preprocessings!)
            if any([T.xbin T.tbin]>1)
                dat = fn_bin(dat,[T.xbin T.xbin T.tbin]);
                if ~isempty(sf), sf = fn_bin(sf,[T.xbin T.xbin T.tbin]); end
            end
        end
    end
    
    % Analog recordings
    methods
        function readrecording(T)
            % T.recording_header had been initialized when reading headers, it
            % just remains to read the data
            ntr = length(T);
            if ntr==1 && ~isempty(T.preprocess.multifile)
                switch T.origin
                    case 'BLK'
                        return
                    otherwise                        
                        error('averaging/cat not implemented yet for recording')
                end
            end
            fn_progress('read recording',ntr)
            for k=1:ntr
                fn_progress(k)
                if isempty(T(k).recording_header), continue, end
                % Automatic locating of analog recording file? - TODO: do
                % this step already when reading the header!!!
                if isempty(T(k).analogfile)
                    [d, base, ext] = fileparts(T(k).file);
                    switch lower(ext(2:end))
                        case {'mat' 'da'}
                            T(k).analogfile = {T(k).file};
                        case {'cfd' 'mpd' 'tif'}
                            fdat = [d '/' base '.dat'];
                            if exist(fdat,'file'), T(k).analogfile = {fdat}; end
                    end
                end
                % Read recording(s) if any
                if isempty(T(k).analogfile), continue, end
                [T(k).recording_data T(k).analogfile] = tps_readrecording(T(k).analogfile);
            end
        end
        function attachanalogfile(T,files,varargin)
            % function attachanalogfile(T,files,channels[,delay])
            %---
            % this function tries reading analog data files and guessing to
            % which trials they correspond
            if nargin<2, files = fn_getfile; if isequal(files,0), return, end, end
            TIMETOL = .2; %fn_input('time__tolerance',.2);
            if isscalar(TIMETOL), TIMETOL = [-TIMETOL TIMETOL]; end
            
            if ~iscell(files), files = cellstr(files); end
            if strfind(files{1},'_adaq.mat')
                % apparently data was acquired with VDAQ and analog files
                % have adaq format
                n = length(files);
                % get experiment and block number
                exp = zeros(1,n);
                block = zeros(1,n); 
                for i=1:n
                    tokens = regexp(fn_fileparts(files{i},'name'),'^E(\d)*B(\d)*_adaq.mat$','tokens');
                    if ~isscalar(tokens), error 'could not read analog data file', end
                    tokens = tokens{1};
                    exp(i) = str2double(tokens{1});
                    block(i) = str2double(tokens{2});
                end
                % match to trials
                MoreInfoFromBLKname(T)
                channels = [];
                for k=1:length(T)
                    info = T(k).fullinfo;
                    if ~isfield(T(k).fullinfo,'exp'), continue, end
                    i = find(exp==info.exp & block==info.block);
                    if isempty(i), continue, elseif ~isscalar(i), error 'could not match analog data to trials', end
                    v = load(files{i},'sampling_rate','nSamples','nChannels','pretrig');
                    if isempty(channels) % happens the first time only
                        switch length(varargin)
                            case 0
                                [channels delay] = tps_channelnames(v.nChannels,-v.pretrig/1000);
                                if isempty(channels), return, end
                            case 1
                                channels = varargin{1};
                                delay = -v.pretrig/1000;
                            case 2
                                [channels delay] = deal(varargin{:});
                            otherwise
                                error argument
                        end
                        channelnames = channels(:,2)';
                        channels = cell2mat(channels(:,1)');
                    end
                    T(k).analogfile = {{files{i} channels}};
                    T(k).recording_header = struct('name',channelnames,'t0',delay,'dt',1/v.sampling_rate,'n',v.nSamples,'signal',[]);
                    T(k).recording_data = {[]}; % reinitialize, in case some previous analog attachment needs to be overwritten
                end
            elseif regexp(files{1},'.abf$')
                % Axon ABF data
                
                % ABF headers
                % (read all headers)
                nf = length(files);
                hrec = struct('type',cell(1,nf),'startfile',cell(1,nf),'startsweep',cell(1,nf));
                for i=1:nf
                    try
                        [d si h] = abfload(files{i},'info'); %#ok<ASGLU>
                    catch ME
                        disp 'error reading .abf file:'
                        disp(ME.message)
                        hrec(i).type = -1;
                        continue
                    end
                    hrec(i).startfile = h.lFileStartTime;
                    if isfield(h,'sweepStartInPts')
                        % data were acquired in waveform fixed-length mode
                        ptspers = 1e6/h.si;
                        hrec(i).type = 1;
                        hrec(i).startsweep = h.lFileStartTime + h.sweepStartInPts'/ptspers; % rec start time in seconds
                        hrec(i).n = h.sweepLengthInPts;
                    else
                        % data were probably acquired in gap-free mode, we
                        % are not interested
                        hrec(i).type = 0;
                        hrec(i).startsweep = h.lFileStartTime;
                        hrec(i).n = h.dataPtsPerChan;
                    end
                    hrec(i).channels = h.recChNames;
                    hrec(i).si = h.si;
                end
                % (remove files that could not be read)
                bad = ([hrec.type]==-1);
                files(bad) = [];
                hrec(bad) = [];
                nf = length(files);
                % (reorder according to time)
                [dum ord] = sort([hrec.startfile]); %#ok<ASGLU>
                files = files(ord);
                hrec = hrec(ord);
                % (get sweep times)
                rectimes = [hrec.startsweep]; 
                nrec = length(rectimes);
                % (for each sweep, memorize corresponding file and sweep
                % number within file)
                rectype = zeros(1,nrec);
                ifile = zeros(1,nrec);
                isweep = zeros(1,nrec);
                icur = 0;
                for i=1:nf
                    ni = length(hrec(i).startsweep);
                    rectype(icur+(1:ni)) = hrec(i).type;
                    ifile(icur+(1:ni)) = i;
                    isweep(icur+(1:ni)) = 1:ni;
                    icur = icur+ni;
                end
                disp 'ABF recording times:'
                sperday = 24*3600;
                for i=1:nrec
                    fprintf('%i|%i: %s - %.0f\n',ifile(i),isweep(i),datestr(rectimes(i)/sperday,'HH:MM:SS'),rectimes(i));
                end
                
                % MES headers
                nT = length(T);
                linetimes = zeros(1,nT);
                disp 'MES trial times:'
                for k=1:nT
                    h = T(k).fullinfo;
                    if ~any(T(k).status=='nr') || ~isfield(h,'MES') || ~strcmp(T(k).type,'linescan'), continue, end
                    h1 = h.MES(1);
                    date = strrep(h1.MeasurementDate(13:end),',','.');
                    sperday = 24*3600;
                    linetimes(k) = mod(datenum(date),1)*sperday; % line start time in seconds
                    fprintf('%i: %s - %.0f\n',k,date,linetimes(k))
                end
                
                % Find the correspondances (rec and line times should
                % differ by less than 200ms)
                timediff = fn_subtract(linetimes,rectimes'); % nrec x nT
                if all(rectype==0)
                    % all gap-free: recording always starts before line
                    % scan, and line scan starts before next recording
                    okmatch = (timediff>0);
                    okmatch(1:end-1,:) = okmatch(1:end-1,:) & (timediff(2:end,:)<0);
                    % if several line scans are candidate to match a given
                    % recording, the first one is the good one
                    okmatch(:,2:end) = okmatch(:,2:end) & ~okmatch(:,1:end-1);
                    [matchsweep matchline] = find(okmatch);
                    if any(diff(matchsweep)<=0) || any(diff(matchline)<=0), error programming, end
                elseif all(rectype==1)
                    % all fixed-length
                    [matchsweep matchline] = find(timediff>TIMETOL(1) & timediff<TIMETOL(2));
                    if any(diff(matchsweep)<=0) || any(diff(matchline)<=0)
                        disp 'non-unique correspondance(s) between line scans and recordings!'
                        disp 'please edit manually ''matchsweep'' and ''matchline'''
                        keyboard
                    end
                else
                    % mixture - don't know how to deal with it
                    disp 'some recordings are fixed-length, others are gap-free'
                    disp 'please define manually ''matchsweep'' and ''matchline'''
                    keyboard
                end
                nmatch = length(matchsweep);
                msgbox(sprintf('found %i correspondance(s)',nmatch))
                if nmatch==0, return, end
                
                % Correct for slowing down of the AxoScope?
                answer = questdlg('Correct for slowing down of the AxoScope?','question');
                switch answer
                    case 'Yes'
                        DTCORR = 1.05;
                    case 'No'
                        DTCORR = 1;
                    case 'Cancel'
                        return
                end
                
                % Register the recordings to the appropriate trials
                for k=1:nmatch
                    irec = matchsweep(k);
                    ifil = ifile(irec);
                    ktrial = matchline(k);
                    nchannel = length(hrec(ifil).channels);
                    okchannel = true(1,nchannel);
                    channelnames = cell(1,nchannel);
                    for i=1:nchannel
                        switch hrec(ifil).channels{i}
                            case {'Vm' 'Ipatch' '_Ipatch' 'AI #0' 'IN 0'}
                                channelnames{i} = 'electrophysiology';
                            case 'ECG'
                                channelnames{i} = 'heart';
                            case {'Vhold' 'I_MTest 1'}
                                channelnames{i} = 'Vhold';
                            otherwise
                                fprintf('unknown channel name ''%s'', please edit code\n',hrec(ifil).channels{i})
                                channelnames{i} = strrep(hrec(ifil).channels{i},' ','');
                        end
                    end
                    T(ktrial).analogfile = {{files{ifil} isweep(irec) find(okchannel)}};
                    T(ktrial).recording_header = struct('name',channelnames(okchannel),'t0',0,'dt',hrec(ifil).si/1e6/DTCORR,'n',hrec(ifil).n,'signal',[]);
                    T(ktrial).recording_data = {[]}; % reinitialize, in case some previous analog attachment needs to be overwritten
                end
            elseif regexp(files{1},'.DAT$') % Elphy format
                nf = length(files);
                [rec ifile kfile] = deal(cell(1,nf));
                for i=1:nf 
                    rec{i} = elphy_read(files{i});
                    nreci = size(rec{i},2);
                    ifile{i} = i*ones(1,nreci);
                    kfile{i} = 1:nreci;
                end
                rec = [rec{:}]; ifile = [ifile{:}]; kfile = [kfile{:}];
                trials = find([T.status]~='s');
                ntrial = length(trials);
                ntrec = size(rec,1);
                if size(rec,2)~=ntrial
                    trial2rec = 1:min(size(rec,2),ntrial); trial2rec(end+1:ntrial)=0;
                    prompt = sprintf('Number of recordings (%i) do not match number of trials (%i), please set recording indices for eachtrial',size(rec,2),ntrial);
                    answer = inputdlg(prompt,'input',1,{fn_chardisplay(trial2rec)});
                    trial2rec = evalin('base',['[' answer{1} ']']);
                else
                    trial2rec = 1:ntrial;
                end
                for ktrial=find(trial2rec)
                    krec = trial2rec(ktrial);
                    T(ktrial).analogfile = {files{ifile(krec)} kfile(krec)};
                    T(ktrial).recording_header = struct('name','Elphy','t0',0,'dt',1e-3,'n',ntrec,'signal',[]);
                    T(ktrial).recording_data = {rec(:,krec)};
                end
            else
                error 'could not read analog data file'
            end
        end
        function addrecording(T,name,dt,data,fdat)
            % function addrecording(T,name,dt,data[,fdat|dosave])
            %---
            % this function add some analog signal 
            % they won't be saved correctly if last argument is missing!
            kc = find(strcmp({T.recording_header.name},name));
            if isempty(kc), kc = length(T.recording_header)+1; end
            T.recording_header(kc) = struct('name',name,'t0',0,'dt',dt,'n',length(data),'signal',[]);
            T.recording_data{kc} = data;
            if nargin>=5
                if islogical(fdat)
                    if ~fdat, return, end
                    fdat = fn_savefile;
                end
                fn_savevar(fdat,data);
                T.analogfile = fdat;
            end
        end
        function setheart(T)
            % this is only for old recording - this option has been mask in
            % tpview version 
            % TODO: check that channel data is indeed heart beat
            for k=1:numel(T)
                switch T(k).origin
                    case 'cfd'
                        error 'channel data not implemented any more'
                        T(k).heartcycle = tps_heartcycle(channeldata,T(k).linedur);
                    case 'ScanImage'
                        f = find(ismember({T(k).recording.name},'heart'));
                        if isempty(f), continue, end
                        if T.nc>1, error 'one heart recording for several conditions!?', end
                        r = T(k).recording(f);
                        if r.n<=1, continue, end
                        % TODO: no interpolation?
                        %                         heart = interp1((0:r.n-1)*r.dt,r.signal,(0:T(k).ny*T(k).nfr-1)*T(k).linedur);
                        %                         T(k).heartcycle = tps_heartcycle(heart,T(k).linedur);
                        heartcycl = tps_heartcycle(r.signal,r.dt);
                        T(k).heartcycle = interp1((0:r.n-1)*r.dt,heartcycl,T(k).tidx);
                end
            end
        end
        function setrecordingname(T,i,name)
            if isempty(i), i = 1:length(T.recording_header); end
            if ~iscell(name), name = {name}; end
            [T.recording_header(i).name] = deal(name{:});
        end
    end
    
    % Data operation
    methods
        function set.opdef(T,op)
            if ~isa(op,'tps_dataopdef'), op = tps_dataopdef(op); end
            if isequal(op,T.opdef), return, end
            % set
            activechg = ~isequal(op([op.active]),T.opdef([T.opdef.active]));
            T.opdef = op;       
            if ~activechg, return, end % only inactive part has changed
            T.conditions.dataop = struct; %#ok<MCSUP>
            T.conditions.sfrop = struct; %#ok<MCSUP>
            T.conditions.shotnoiseop = struct; %#ok<MCSUP>
        end
        function set.opmem(T,op)
            if ~isa(op,'tps_dataopdef'), op = tps_dataopdef(op); end
            activechg = ~isequal(op([op.active]),T.opmem([T.opmem.active]));
            T.opmem = op;
            if ~activechg, return, end % only inactive part has changed
            T.conditions.dataopmem = struct; 
        end
        function memorizeop(T)
            [T.opmem] = deal(T.opdef);
        end
        function x = operation(T,dataflag,op,doenlarge)
            % operation definition
            if nargin<3 || isempty(op)
                op = T.opdef;
                isdefaultop = true;
            else
                isdefaultop = false;
            end
            if nargin<4
                doenlarge = true;
            end
            
            % saving politics
            dosavedef = true;
            
            % some parameters
            scanline = strcmp(T.type,'linescan');
            
            % check whether part or all of the computation is already
            % stored! 
            activeop = op([op.active]);
            if isempty(activeop), x = T.(dataflag); return, end
            opsteps = T.operationsteps.(dataflag);
            kopavailable = 0;
            for kop=1:length(activeop)+1
                % is operation at step kop the same in the new definition?
                if length(opsteps)<kop || length(activeop)<kop ...
                        || ~isequal(opsteps(kop).op,activeop(kop))
                    break
                end
                % if so, is previous result available?
                if isempty(opsteps(kop).value) 
                    error 'no memorypool item!'
                elseif ~isempty(opsteps(kop).value.data)
                    kopavailable = kop;
                end
            end
            kopsame = kop-1;
            % is at least a part of the steps already available?
            if kopavailable>0
                % yes
                [x mask deltat] = deal(opsteps(kopavailable).value.data{:});
                [nx ny nfr nc] = size(x);
                % do we need additional steps?
                if length(activeop)>kopavailable
                    % yes
                    if isdefaultop
                        % throw away the divergent steps
                        opsteps(kopsame+1:end) = [];
                    else
                        % we are performing a custom operation, so do not
                        % throw away the steps of the default operation and
                        % do not overwrite them!
                        dosavedef = false;
                    end
                end
            else
                % no, we need to start from zero
                if isdefaultop
                    % throw away divergent operations
                    opsteps(kopsame+1:end) = [];
                elseif ~isempty(opsteps)
                    % we are performing a custom operation, so do not
                    % throw away the steps of the default operation and
                    % do not overwrite them!
                    dosavedef = false;
                end
                x = T.(dataflag);
                if isempty(x), return, end % happens for example with sfr
                % conversion to floating number; note that double-precision
                % is needed when looking at signals of scale 1e-4 DF/F
                if ~issparse(x), x = cast(x,T.internpar.floattype); end
                [nx ny nfr nc] = size(x);
                deltat = T.dt_sec;
                if scanline
                    x = reshape(x,nx,1,T.ny*T.nfr,nc);
                    ny = 1;
                    nfr = T.ny*T.nfr;
                    deltat = deltat/T.ny;
                end
                mask = true(nx,ny);
            end
            
            % continue with the additional operations
            try
                for kop=kopavailable+1:length(activeop)
                    opk = activeop(kop);
                    if ~opk.active, continue, end
                    valk = opk.value;
                    % common to some operations: use a sub-mask
                    ival = fn_switch(opk.name,'time',2,{'mask' 'heart' 'cameranoise' 'motion'},1,[]);
                    if ~isempty(ival) && strcmp(valk{ival},'indices')
                        ind = valk{ival+1};
                        if ischar(ind), ind = fn_strcut(ind,','); else ind = {ind}; end
                        if islogical(ind) && isequal(size(ind),[T.nx T.ny])
                            mask1 = ind;
                        else
                            if ~ismember(length(ind),[1 2]), error 'wrong indices specification', end
                            mask1 = false(T.nx,T.ny);
                            mask1(fn_subsref([T.nx T.ny],ind{:})) = true;
                            if scanline && T.ny~=1
                                error 'not implemented yet'
                            elseif any([nx ny nfr]~=[T.nx T.ny T.nfr])
                                xbin = floor(T.nx/nx);
                                ybin = floor(T.ny/ny);
                                tbin = floor(T.nfr/nfr);
                                mask1 = fn_bin(mask1,[xbin ybin tbin],'or');
                            end
                        end
                        mask1 = mask & mask1;
                    end
                    % different operations
                    tstart = tic;
                    switch opk.name
                        case 'space'
                            act = valk{1};
                            switch valk{2}
                                case 'frames'
                                    ind = valk{3};
                                    if ischar(ind) && ~strcmp(ind,':')
                                        ind = eval(['[' ind ']']);
                                    end
                                case 'time (s)'
                                    if isempty(deltat), error 'frame rate not defined', end
                                    % convert time
                                    % syntax 't1-t2'
                                    tokens = regexp(valk{3},'(\d*(.\d+){0,1}) *- *(\d*(.\d+){0,1})','tokens');
                                    if ~isempty(tokens)
                                        tbeg = str2double(tokens{1}{1});
                                        tend = str2double(tokens{1}{2});
                                    else
                                        tbeg = 0;
                                        tend = str2double(valk{3});
                                        if isnan(tend), error('cannot interpret string'), end
                                    end
                                    ind = max(1,1+round(tbeg/deltat)):min(nfr,1+round(tend/deltat));
                                otherwise
                                    error('cannot average using ''%s''',valk{2})
                            end
                            s = substruct('()',{':' ':' ind ':'});
                            m = mean(subsref(x,s),3);
                            if strcmp(act,'subtract without mean'), m = m-mean(m(~isnan(m))); end
                            switch act
                                case 'divide by'
                                    x=fn_div(x,m); 
                                case {'subtract' 'subtract witout mean'}
                                    x=fn_subtract(x,m);
                                case 'DF/F'
                                    x=fn_div(fn_subtract(x,m),m); % =fn_div(x,m)-1
                            end
                        case 'time'
                            act = valk{1};
                            % if linescan mode, we are already in 'line per
                            % line' mode, in fact
                            dolineperline = ~scanline && ~isempty(strfind(act,'line per line'));
                            if dolineperline
                                x1 = fn_mult(x,mask1);
                                m = fn_div(sum(x1,1),sum(mask1,1));
                                if strfind(act,'div')
                                    x=fn_div(x,m);
                                else
                                    m = fn_normalize(m,[3 4],'-');
                                    x=fn_subtract(x,m);
                                end
                            elseif strcmp(act,'special neuropil')
                                %                             if scanline
                                %                                 disp('''special neuropil'' not implemented yet for line scans')
                                %                             else
                                %                                 x = tps_subtractneuropil(T,valk{3});
                                %                             end
                                error('not implemented any more')
                            else
                                m = mean(fn_imvect(x,mask1),1);
                                m = shiftdim(m,-1);
                                if strcmp(act,'subtract without mean'), m = m-mean(m(~isnan(m))); end
                                if strfind(act,'div'), x=fn_div(x,m); else x=fn_subtract(x,m); end
                            end
                        case 'mask'
                            method = valk{1};
                            switch method
                                case {'threshold' 'threshold on avg. value' 'threshold below'}
                                    if strcmp(method,'threshold') && fn_dodebug, error 'replacement should have occured in tps_dataopdef/convert', end
                                    thr = str2double(valk{2});
                                    if isnan(thr), error 'set a number for the threshold', end
                                    m = mean(mean(x,3),4);
                                    if strcmp(method,'threshold below')
                                        mask = mask & (m<thr);
                                    else
                                        mask = mask & (m>thr);
                                    end
                                case 'indices'
                                    mask = mask1;
                            end
                            outsideval = valk{3};
                            if ismember(outsideval,{'NaN' '0' '1'})
                                outsideval = str2double(outsideval);
                            else
                                x1 = fn_imvect(x,mask);
                                outsideval = feval(outsideval,x1(:));
                            end
                            x = reshape(x,[nx*ny nfr nc]);
                            x(~mask,:,:) = outsideval;
                            x = reshape(x,[nx ny nfr nc]);
                        case 'bleach'
                            [n_exp, par, output_mode, do_global, do_blank] = dealc(valk);
                            % checks
                            errorstr = '';
                            if isempty(par)
                                errorstr = 'bleaching parameters are not defined';
                            elseif n_exp ~= length(par.tau)
                                errorstr = 'number of exponentials cannot be changed';
                            elseif do_global ~= isfield(par,'beta')
                                errorstr = '''global'' parameter cannot be changed';
                            end
                            if ~isempty(errorstr)
                                error([errorstr ', please run ''EST'' on average trial'])
                            end
                            % fit
                            if do_blank
                                % use estimation from the blank
                                blank = mean(x(:,:,:,par.blankcond),4);
                                blankcorr = tps_expfit(blank,par);
                                blankbleach = blank-blankcorr;
                                switch output_mode
                                    case 'remove bleach'
                                        x = fn_subtract(x,blankbleach);
                                    case 'keep only bleach'
                                        x = fn_add(mean(x,3),blankbleach);
                                end
                            else
                                xcorr = tps_expfit(x,par);
                                switch valk{3}
                                    case 'remove bleach'
                                        x = xcorr;
                                    case 'keep only bleach'
                                        x = fn_add(mean(x,3),x - xcorr);
                                end
                            end
                        case 'cameranoise'
                            noisesignal = mean(fn_imvect(x,mask1),1);
                            noisesignal = shiftdim(noisesignal,-1);
                            noisesignal = noisesignal-mean(noisesignal(:));
                            if ~isempty(valk{3})
                                if isempty(deltat), error 'frame rate not defined', end
                                timethr = str2double(valk{3});
                                noisesignal = fn_filt(noisesignal,timethr/deltat,'hm',3);
                            end
                            x = fn_subtract(x,noisesignal);
                        case 'bidir scan'
                            delay = str2double(valk{1});
                            x = tpv_bidircorr(x,delay,'same');
                        case 'heart'
                            if isempty(deltat), error 'frame rate not defined', end
                            freqrange = valk{3} * [.65 1.5];
                            nbaseh = valk{4};
                            prec = valk{5};
                            flag = valk{6};
                            x = fn_imvect(x,mask);
                            % Prepare heart cycle estimation
                            switch valk{1}
                                case 'indices'
                                    % use optics to estimate heart
                                    iheartcycl = mask1(mask); % only a subpart is used to estimate heart cycle
                                    xheartcycl = shiftdim(mean(x(iheartcycl,:,:),1),1);
                                    dtheart = deltat;
                                    ttheart = [];
                                case 'recording'
                                    % use heart recording
                                    r = T.recording;
                                    kheart = find(strcmpi({r.name},'heart'));
                                    if ~isscalar(kheart) || isempty(r(kheart).signal)
                                        error 'no heart recording'
                                    end
                                    dtheart = r(kheart).dt;
                                    ttheart = r(kheart).t0 + (0:r(kheart).n-1)*dtheart;
                                    xheartcycl = r(kheart).signal;
                                    xheartcycl = abs(fn_filt(xheartcycl,.05/dtheart,'h'));
                                otherwise
                                    error('unknown flag ''%s'' for the origin of the heart signal to use',valk{1})
                            end
                            % regressors for all trials 
                            nbase = 2*nbaseh;
                            H = zeros(nfr,nc,nbase);
                            for kcond = 1:nc
                                % (heart cycle)
                                heartcycl = tps_heartcycle(xheartcycl(:,kcond),dtheart,freqrange,prec);
                                if ~isempty(ttheart), heartcycl = interp1(ttheart,heartcycl,(0:nfr-1)*deltat); end
                                % (regressors)
                                for kbase=1:nbaseh
                                    H(:,kcond,2*kbase-1) = cos(heartcycl*(2*kbase*pi));
                                    H(:,kcond,2*kbase)   = sin(heartcycl*(2*kbase*pi));
                                end
                            end
                            H = fn_normalize(H,1,'-'); % remove constants from base functions
                            % regression: on blanks if possible
                            stimtype = {T.stimdetails.type};
                            kblank = strcmp(stimtype,'blank');
                            if ~any(kblank)
                                A = H; 
                                y = x;
                            else
                                A = H(:,kblank,:);
                                y = x(:,:,kblank);
                            end
                            H = fn_reshapepermute(H,{[1 2] 3});
                            A = fn_reshapepermute(A,{[1 2] 3});
                            y = fn_reshapepermute(y,{1 [2 3]});
                            beta = y * pinv(A)'; % operation on columns
                            % apply correction
                            heart = beta * H';
                            heart = reshape(heart,[sum(mask(:)) nfr nc]);
                            dokeep = strcmp(flag,'keep');
                            if dokeep
                                x = fn_add(mean(x,2),heart); % keep the basline
                            else
                                x = x-heart; 
                            end
                            x = fn_imvect(x,mask,NaN);
                        case 'motion'
                            par = fn_register('par');
                            if ~all(mask(:)), par.mask = mask; end
                            par.maxshift = valk{3};
                            par.ref = valk{4};
                            for i=1:nc
                                [shift e x(:,:,:,i)] = fn_register(x(:,:,:,i),par); %#ok<ASGLU>
                            end
                        case 'filter'
                            xl = str2num(valk{1}); %#ok<ST2NM>
                            xh = str2num(valk{2}); %#ok<ST2NM>
                            tl = str2num(valk{3}); %#ok<ST2NM>
                            th = str2num(valk{4}); %#ok<ST2NM>
                            if any([tl th]) && isempty(deltat), error 'frame rate not defined', end
                            tl = tl/deltat; th = th/deltat;
                            dflag = fn_switch(valk{5},'d','');
                            zflag = fn_switch(valk{6},'z','');
                            phaseflags = fn_switch(valk{7},{'phase'},{});
                            sflag = fn_switch(valk{8},'s','');
                            if ~isempty([xl xh])
                                if isempty(xl), xl = 0; end
                                if isempty(xh), xh = 0; end
                                if all(mask(:))
                                    if ~isempty(zflag), m = mean(x(~isnan(x))); x = x-m; end
                                    x = fn_filt(x,[xl xh],[1 2]);
                                    if ~isempty(zflag), x = x+m; end
                                else
                                    if xl && xh, error 'band-pass with a mask not implemented yet', end
                                    xs = fn_switch(xl,xl,xh);
                                    x1 = reshape(x,[nx*ny nfr nc]);
                                    x1(~mask,:,:) = 0;
                                    if xh && ~isempty(zflag)
                                        m = x1(mask,:,:);
                                        m = mean(m(~isnan(m)));
                                    end
                                    x1 = reshape(x1,[nx ny nfr nc]);
                                    x2 = fn_filt(x1,xs,'l',[1 2]);
                                    maskf = fn_filt(mask,xs,'l',[1 2]);
                                    x2 = fn_div(x2,maskf);
                                    if xh, x2 = x1-x2; end
                                    x2 = reshape(x2,[nx*ny nfr nc]);
                                    if xh && ~isempty(zflag), x2 = fn_add(x2,m); end
                                    x2(~mask,:) = NaN;
                                    x = reshape(x2,[nx ny nfr nc]);
                                end
                            end
                            if ~isempty([tl th])
                                x = fn_imvect(x,mask);
                                if isempty(tl), tl = 0; end
                                if isempty(th), th = 0; end
                                x = fn_filt(x,[tl th],[dflag zflag sflag],2,phaseflags{:});
                                x = fn_imvect(x,mask,NaN);
                            end
                            % go back to single precision if needed
                            x = cast(x,T.internpar.floattype);
                        case 'detrend'
                            mode = valk{1};
                            idx = fn_subsref(nfr,valk{2}); % indices on which to base the detrending regression
                            nidx = length(idx);
                            constant = ones(1,nfr) / sqrt(nidx);          % constant(idx) is a normal vector
                            linear = 1:nfr;
                            linear = linear - (linear(idx)*constant(idx)')*constant;  % now constant(idx) and linear(idx) are orthogonal
                            linear = linear / norm(linear(idx));                      % now constant(idx) and linear(idx) are orthonormal
                            x = fn_imvect(x,mask);
                            for kcond=1:nc
                                xk = x(:,:,kcond);
                                % detrend
                                xk = xk - (xk(:,idx)*linear(idx)')*linear; % operation on columns
                                % normalize?
                                if strcmp(mode,'detrend & normalize')
                                    xk = xk ./ ((xk(:,idx)*constant(idx)')*constant);
                                end
                                x(:,:,kcond) = xk;
                            end
                            x = fn_imvect(x,mask,NaN);
                        case 'bin'
                            binx = str2num(valk{1}); %#ok<ST2NM>
                            if scanline
                                if ~isscalar(binx), error 'in scanline mode, spatial binning must be a scalar value', end
                                binx = [binx 1]; %#ok<AGROW>
                            elseif isscalar(binx)
                                binx = [binx binx]; %#ok<AGROW>
                            end
                            if length(binx)~=2, error('wrong binning value'), end
                            bint = str2double(valk{2});
                            x = fn_bin(x,[binx bint]);
                            mask = fn_bin(mask,binx);
                            [nx ny nfr nc] = size(x);
                            deltat = deltat*bint;
                        case 'box-car'
                            binx = str2num(valk{1}); %#ok<ST2NM>
                            if scanline
                                if ~isscalar(binx), error 'in scanline mode, spatial binning must be a scalar value', end
                                binx = [binx 1]; %#ok<AGROW>
                            elseif isscalar(binx)
                                binx = [binx binx]; %#ok<AGROW>
                            end
                            if length(binx)~=2, error('wrong binning value'), end
                            bint = str2double(valk{2});
                            bin = [binx bint];
                            s = size(x);
                            for dim=1:3
                                if bin(dim)==1, continue, end
                                kernel = shiftdim(ones(bin(dim),1),1-dim) / bin(dim);
                                x = convn(x,kernel,'same');
                                ref = convn(shiftdim(ones(s(dim),1),1-dim),kernel,'same');
                                x = fn_div(x,ref);
                            end
                        case 'fft'
                            x = fft(x,[],3);
                            x(:,:,1) = 0;
                            if valk{1}, x = abs(x); end
                        case 'user'
                            if scanline, x = reshape(x,[nx nfr nc]); end
                            fun = valk{1};
                            if exist(fun,'file')
                                x = feval(fun,x);
                            else
                                fun = evalin('base',['@(x)' valk{1}]);
                                x = fun(x);
                            end
                            if scanline, x = reshape(x,[nx 1 nfr nc]); end
                        otherwise
                            error('unknown operation flag ''%s''',opk.name)
                    end
                    
                    % save the step
                    if dosavedef
                        weight = kop + 5*toc(tstart) + 10*(kop==length(activeop)); % weight of the memory pool item
                        if kop<=kopsame
                            if ~isequal(opsteps(kop).op,opk), error programming, end
                            recover(opsteps(kop).value,{x mask deltat}, ...
                                [fn_fileparts(fileshort(T),'base') ' - ' dataflag ' - ' opk.name],weight);
                        else
                            opsteps(kop).op = opk;
                            opsteps(kop).value = memorypool.item({x mask deltat}, ...
                                [fn_fileparts(fileshort(T),'base') ' - ' dataflag ' - ' opk.name],weight); 
                        end
                    end
                end
            catch ME
                disp(['data operation failed: ' ME.message])
                if fn_dodebug
                    disp 'checking empty output upon operation error'
                    x = [];
                else
                    x = zeros(T.sizes);
                end
            end
            
            % final reshape/enlarge
            if scanline && ~isempty(x), x = reshape(fn_enlarge(x,[T.nx 1 T.ny*T.nfr T.nc]),T.nx,T.ny,T.nfr,T.nc); end
            if doenlarge, x = fn_enlarge(x,T.sizes); end
            
            % save the steps
            T.operationsteps.(dataflag) = opsteps;
        end
        function coregister(T,tbin,xindices,yindices,sigma,display) %#ok<MANU,INUSD>
            error('this function is obsolete, use tps_register instead')
        end
    end
    
    % Load/Save
    methods
        function savedata(T,fname)
            % function savedata(T[,fname])
            if nargin<2
                fname = fn_savefile('*.mat','Enter base name for saving trial data'); 
            end
            fname = cellstr(fname);
            if fname{1}(1)~='/' && fname{1}(2)~=':' % absolute path in Linux/Mac and Windows
                error('file name should be provided with the full path')
            end
            if ~ismember(length(fname),[1 numel(T)]), error('number of names does not match number of trials'), end
            n = numel(T);
            if isscalar(fname) && ~isscalar(T) 
                fname = repmat(fname,n,1);
                m = floor(log10(n))+1;
                num = num2str((1:n)',['_%.' num2str(m) 'i']);
                num = cellstr(num);
            else
                num = cell(1,n);
            end
            if n>1, fn_progress('saving trial',n), end
            for k=1:n
                if n>1, fn_progress(k), end
                if ~isempty(T(k).file), error 'it is not permitted to save data that is already attached to a data file', end
                [p b ext] = fileparts(fname{k});
                if isempty(ext), ext = '.mat'; end
                fname{k} = fullfile(p,[b num{k} ext]);
                data = T(k).data_unsaved{1}; sfr = T(k).data_unsaved{2}; %#ok<NASGU>
                rec = T(k).recording_data; %#ok<NASGU>
                save(fname{k},'data','sfr','rec','-v7.3')
                T(k).file = fname{k};
                T(k).data_unsaved = [];
            end
        end
        function T = saveobj(T)
            % save data if necessary
            if isempty(T.file) || isequal(T.file,0)
                fname = fn_savefile('*.mat','Enter base name for saving trial data');
                savedata(T,fname)
            end
            % do not attempt to read data to save data0
            T.preventdata0comp = true;
        end
    end
    methods (Static)
        function T = loadobj(T)
            % problematic versions...
            if isstruct(T)
                s = T;
                if ~isscalar(s), error('cannot handle non-scalar object'), end
                T = tps_trial; erasedata(T)
                T.basics.data = [];
                T.data_unsaved = [];
                F  = fieldnames(s);
                M = metaclass(T);
                F1 = setdiff({M.PropertyList.Name},{'recording','sizes'});
                F2 = {'nx','ny','nfr','nc','linedur','tidx'};
                for i=1:length(F)
                    f = F{i};
                    if ismember(f,F1)
                        if ismember(f,F2), continue, end
                        T.(f) = s.(f);
                    else
                        % try to do something...
                        switch f
                            case 'recording'
                                T.recording_header = s.recording; 
                            case 'sizes'
                                T.sizes0 = s.sizes;
                            case 'selection'
                                T.user.selection = s.selection;
                            case 'stim'
                                T.stim = s.stim;
                            case 'fileflag'
                                T.preprocess.multifile = s.fileflag;
                            otherwise
                                fprintf('tps_trial: field ''%s'' does not exist any more\n',f);
                        end
                    end
                end
            end
            % old version where information on recordings was stored in
            % AcqSet
            if strcmp(T.origin,'cfd') && isempty(T.recording) && isfield(T.fullinfo,'AcqSet') ...
                    && isfield(T.fullinfo.AcqSet,'recordingchannels')
                AcqSet = T.fullinfo.AcqSet;
                r = AcqSet.recordingchannels;
                for i=1:length(r)
                    T.recording_header(i).name = r{i};
                end
                if ~isempty(r), [T.recording_header.dt] = deal(1/AcqSet.samplingrate); end
            end
            % depth might not have been read yet in old versions
            if strcmp(T.origin,'cfd') && ~isfield(T.addinfo,'depth')
                tmp = T.fullinfo.CFH.zpos/1000;
                if tmp>0 && tmp<500 % looks like z=0 was on the surface
                    T.addinfo.depth = sprintf('%.0fum',tmp);
                end
            end
            % more misdefined values from old versions
            doassumptions = false;
            if strcmp(T.type,'Time'), T.type = 'movie'; end
            if isempty(T.scanning), T.scanning = true; doassumptions = true; end
            if doassumptions, disp('Some assumptions were used for misdefined values, please double-check trial header info'), end
            if ~isa(T.opdef,'tps_dataopdef'), T.opdef = tps_dataopdef(T.opdef); end
            if isfield(T.addinfo,'stim'), T.stim = T.addinfo.stim; T.addinfo = rmfield(T.addinfo,'stim'); end
            % check whether data files have been moved
            locatefile(T);
            % now, we have a 'version' property!
            if isempty(T.version), T.version = 1.0; end
            if T.version < 1.1
                T.xunit = 'um';
                T.tunit = 's';
            end
            if T.version < 1.2
                % version 1.2 added the fields stimtable and stimid, while
                % it made field stim dependent
                if isempty(T.stimtable)
                    % note that T.stimtable might already have been defined
                    % when setting the dependent property 'stim'
                    T.stimtable = tps_stimtable;
                    setstim(T,0:T.nc-1); % this will create a 'blank' stim entry
                elseif T.nc>1 && isscalar(T.stimid)
                    T.stimid = repmat(T.stimid,1,T.nc);
                end
            end
            if T.version < 1.3 && isfield(T.addinfo,'createdintpview')
                % the way to mark special trials has changed
                T.addinfo = rmfield(T.addinfo,'createdintpview');
                T.status = 's';
            end
            if T.version < 1.4 
                % new field 't0' to the recording headers
                if ~isempty(T.recording_header), [T.recording_header.t0] = deal(0); end
            end
            if T.version < 1.5 || ~fn_dodebug
                % version 1.5 made some cleanup in the properties and added
                % fields 'sizes0', 'xbin' and 'tbin'; the change is handles
                % above (see line 'T.sizes0 = s.sizes;')
                if isempty(T.sizes0) && isempty(T.sizes)
                    if fn_dodebug, disp 'need to load data to get size!', end
                    T.sizes0 = size(T.data);
                elseif isempty(T.sizes0)
                    if T.xbin~=1 || T.tbin~=1
                        error 'cannot load data saved in previous version of tpview'
                    end
                    T.sizes0 = T.sizes;
                elseif isempty(T.sizes)
                    T.sizes0 = T.sizes0;
                end
            elseif isempty(T.sizes0) || isempty(T.sizes)
                disp 'problem: T.sizes0 or T.sizes is empty!'
                %keyboard
            end
            if T.version < 1.6
                % one single recording structure contains the data for all
                % the conditions
                if ~isempty(T.recording_header), T.recording_header = T.recording_header(1,:); end
            end
            if T.version<1.7 && strcmp(T.origin,'MES')                
                % add AO amplitude to the addinfo
                header = T.fullinfo.MES;
                for k=1:length(header)
                    hk = header(k);
                    if isfield(hk,'info_AOsettings') && ~isempty(hk.info_AOsettings)
                        T.addinfo.laser = hk.info_AOsettings.uAO1y.RAMamps(2);
                        break
                    end
                end
            end
            T.version = 1.7;
        end
    end
    
    % Display
    methods
        function disp(T)
            %str = '<a href = "matlab:help tps_trial">tps_trial</a> object with following header info:';
            str = 'tps_trial object with following header info:';
            if ~isscalar(T)
                siz = size(T);
                str = [fn_strcat(siz,'x') ' ' str];
            end
            fprintf(['  ' str '\n\n'])
            fields = {'preprocess' 'xbin' 'tbin' 'nx' 'ny' 'nfr' 'nc' 'dx' 'dt' 't0' 'opdef' 'recording' 'stim'};
            nf = length(fields);
            ntrial = length(T);
            values = cell(1,nf);
            for i=1:nf
                f = fields{i};
                switch f
                    case {'preprocess' 'opdef' 'opmem' 'recording'}
                        str = [];
                        for k=1:ntrial
                            switch f
                                case 'preprocess'
                                    s = T(k).preprocess.op;
                                otherwise
                                    s = T(k).(f);
                            end
                            if isempty(s), strk={}; else strk = {s.name}; end
                            if ~iscell(str)
                                str = strk;
                            elseif ~isequal(strk,str)
                                str = '(differ between trials)';
                                break
                            end
                        end
                    case 'stim'
                        stimtable1 = T(1).stimtable;
                        oktable = true;
                        for k=2:ntrial
                            if T(k).stimtable~=stimtable1
                                oktable = false;
                                break
                            end
                        end
                        if oktable
                            table = stimtable1.table;
                            ok = ~any(diff([T.nc]));
                            if ok
                                id = cat(1,T.stimid);
                                ok = ~any(row(diff(id)));
                            end
                            if ok
                                id = id(1,:);
                            else
                                id = unique([T.stimid]);
                            end
                            stimname = cell(1,length(id));
                            for j=1:length(id)
                                k = find([table.id]==id(j));
                                if ~isempty(k) && ~isempty(table(k).name)
                                    stimname{j} = table(k).name;
                                else
                                    stimname{j} = num2str(id(j));
                                end
                            end
                            str = fn_strcat(stimname,', ');
                            if ~ok
                                str = ['differ between trials (' str ')']; %#ok<AGROW>
                            end
                        else
                            str = '(stim tables differ between trials)';
                        end
                    case {'dx' 'dt'}
                        x = unique([T.(f)]);
                        unit = unique({T.([f(2) 'unit'])});
                        if isscalar(x) && isscalar(unit)
                            str = [num2str(x) unit{1}];
                        else
                            str = '(differ between trials)';
                        end
                    otherwise
                        x = unique([T.(f)]);
                        if isscalar(x)
                            str = num2str(x);
                        else
                            str = '(differ between trials)';
                        end
                end
                if iscell(str), str = fn_strcat(str,', '); end
                if isempty(str), str = '[]'; end
                values{i} = str;
            end
            names = fliplr(char(fn_map(@fliplr,fields)));
            for i=1:nf
                disp(['    ' names(i,:) ': ' values{i}])
            end
            fprintf('\n')
            str = ['(raw data is in ''data'' property, processed data in ''dataop'' property' ...
                fn_switch(any([T.sfrchannel]),', second channel data in ''sfr'' and ''sfrop'' properties','') ...
                ')'];
            fprintf(['  ' str '\n\n'])
        end
    end
    
    % User
    methods
        function b = checkconsistency(T)
            if isempty(T), b = true; return, end
            if isscalar(T), b = true; return, end
            pars = get(T,{'sizes' 'dx' 'dy' 'dz' 'dt'});
            pars = num2cell(pars,2);
            b = isequal(pars{:});
            if ~b, return; end
        end
        function y = gettimecourses(T,selection)
            if nargin<2, dataname='dataop'; end
            if ~fn_ismemberstr(dataname,{'data','dataop','sfr','sfrop','shotnoise','shotnoiseop'})
                error('wrong data field ''%s''',dataname)
            end
            % non-scalar object
            if isempty(T)
                y = [];
                return
            elseif ~isscalar(T)
                n = length(T); 
                nfr = T(1).nfr; nsel = length(selection.singleset);
                if checkconsistency(T)
                    y = zeros(nfr,n,nsel);
                    for k=1:n
                        yk = gettimecourses(T(k)); %,dataname)
                        y(:,k,:) = reshape(yk,[nfr 1 nsel]); 
                    end
                else
                    y = cell(size(T));
                    for k=1:n, y{k} = gettimecourses(T(k)); end %,dataname)
                end
                return
            end
            % empty selection
            if isempty(selection), y = []; return, end
            % compute indices
            if isempty(selection.datasizes)
                selection = setdatasizes(selection,[T.nx T.ny]);
            end
            % get time courses (loop on selections)
            nsel = length(selection.singleset);
            y =zeros(T.nfr,nsel);
            dat = reshape(T.(dataname),[T.nx*T.ny T.nfr]);
            for k=1:nsel
                ind = selection.singleset(k).dataind;
                y(:,k) = mean(dat(ind,:),1);
            end
        end
        function nosfr(T)
            [T.sfrchannel] = deal(false);
            for k=1:length(T), T(k).basics.sfr = []; end       
        end
        function erasedata(T)
            for k=1:length(T)
                T(k).basics = struct('data',[],'sfr',[],'shotnoise',[]);
                T(k).data0 = [];
            end 
            erasecomp(T)
        end
        function erasecomp(T)
            if isempty(T), return, end
            [T.heartcycle] = deal([]);
            [T.conditions] = deal(struct('data',struct,'sfr',struct,'shotnoise',struct, ...
                'dataop',struct,'dataopmem',struct,'sfrop',struct,'shotnoiseop',struct));
            [T.operationsteps] = deal(struct('data',struct('op',{},'value',{}), ...
                'sfr',struct('op',{},'value',{}), ...
                'shotnoise',struct('op',{},'value',{})));
        end
    end
    
    % Debug
    methods
        function access(T) %#ok<MANU>
            keyboard
        end
    end
end


%---
function [dat, sf] = preprocessdata(dat,sf,preop)

% Preprocess
npreop = length(preop);
for kop = 1:npreop
    preopk = preop(kop);
    f = preopk.name;
    switch f
        case 'correctshift'
            arg = preopk.value; if ~iscell(arg), arg = {arg}; end
            dat = tps_correctshift(dat,arg{:});
            if ~isempty(sf), sf = tps_correctshift(sf,arg{:}); end
        case 'translate'
            shift = preopk.value;
            dat = fn_translate(dat,shift,'full','linear');
            if ~isempty(sf), sf = fn_translate(sf,shift,'full','linear'); end
        case 'crop'
            crop = preopk.value;
            dat = dat(crop{1},crop{2},:,:);
            if ~isempty(sf), sf = sf(crop{1},crop{2},:,:); end
        case 'user'
            fun = preopk.value;
            dat = fun(dat);
            if ~isempty(sf), sf = fun(dat); end
        otherwise
            error('unknown preprocessing ''%s''',f)
    end
end

end

%---
function [addinfo objective recording] = readacqsetheader(T,AcqSet)

addinfo.setup = ['version ' AcqSet.Setupversion ', telescope ' fn_num2str(AcqSet.Telescope)];
str = AcqSet.Objective;
addinfo.objective = str;
str(str<'0' | str>'9')=' ';
objective = sscanf(str,'%i',1);
if isfield(AcqSet,'Laser')
    % setup version 1.31
    addinfo.Laser = fn_str2double(AcqSet.Laser);
else
    % setup version < 1.31
    addinfo.PockelCell = fn_str2double(AcqSet.PockelCell);
end
if isfield(AcqSet,'PMT')
    % old version
    addinfo.PMT = fn_str2double(AcqSet.PMT);
else
    addinfo.PMT1 = fn_str2double(AcqSet.PMT1);
    addinfo.PMT2 = fn_str2double(AcqSet.PMT2);
end

% Stimulation
if isfield(AcqSet,'stimulationparadigmxml')
    switch AcqSet.stimulationVI
        case 'directions_1channel_makestim.vi'
            c = AcqSet.stimulationparadigmxml.Conditions;
            if length(c)>1 
                if ~isfield(AcqSet,'stimnum')
                    error('old version, stimulus number is not defined')
                elseif length(c)>1
                    c = c(AcqSet.stimnum+1);
                end
            end
            if c.amplitudeV==0
                stim = [];
            else
                if c.risetimems || c.falltimems, disp('stimulation has rise/fall time'), end
                npulse = max(1,c.numberofpulses);
                duration = c.risetimems + c.plateautimems + c.falltimems;
                stim = [c.delays + (0:npulse-1)*(c.intervalbtwpulsesms*1e-3); ...
                    ones(1,npulse)*(duration*1e-3)];
            end
        case 'stim_pulses_v1.vi'
            c = AcqSet.stimulationparadigmxml.Conditions(AcqSet.stimnum+1).Pulses;
            if any([c.risetimems c.falltimems]), disp('stimulation has rise/fall time'), end
            stim = [];
            for k=1:length(c)
                ck = c(k);
                if ck.amplitudeV==0, continue, end
                npulsek = max(1,ck.numberofpulses);
                durationk = ck.risetimems + ck.plateautimems + ck.falltimems;
                stimk = [ck.delays + (0:npulsek-1)*(ck.intervalbtwpulsesms*1e-3); ...
                    ones(1,npulsek)*(durationk*1e-3)];
                if ~isempty(stim) && stimk(1,1)<=stim(1,end), disp('stimulation is problematic'), end
                stim = [stim stimk]; %#ok<AGROW>
            end
        case ''
            stim = [];
        otherwise
            disp(['don''t know how to read paradigm for VI ''' ...
                AcqSet.stimulationVI ''''])
            stim = [];
    end
    T.stim = stim;
elseif isfield(AcqSet,'mode') && strcmp(AcqSet.mode,'free2P') % old version!
    % stim is [time instants; durations]
    % cell array if several conditions
    [p fil] = fileparts(AcqSet.parameterfile); %#ok<ASGLU>
    switch fil
        case 'direction_1channel_rest11s-1Hz10s'
            stim = [11:1:20; ones(1,10)*50e-3];
        case 'direction_1channel_rest11s-2Hz5s'
            stim = [11:.5:15.5; ones(1,10)*50e-3];
        case 'direction_1channel_rest21s-2Hz5s'
            stim = [21:.5:25.5; ones(1,10)*50e-3];
        case 'direction_1channel_rest3s-5Hz2s'
            stim = [3:.2:4.8; ones(1,10)*50e-3];
        case 'direction_1channel_rest10s-5Hz2s'
            stim = [10:.2:11.8; ones(1,10)*50e-3];
        case 'direction_1channel_rest11s-5Hz2s'
            stim = [11:.2:12.8; ones(1,10)*50e-3];
        case 'direction_1channel_rest20s-5Hz2s'
            stim = [20:.2:21.8; ones(1,10)*50e-3];
        case 'direction_1channel_rest21s-5Hz2s'
            stim = [21:.2:22.8; ones(1,10)*50e-3];
        case 'direction_1channel_electrodeResistance'
            stim = [0; 50e-3];
        case 'direction_1channel_square'
            stim = []; %{[] [2:.125:3.875; ones(1,16)*50e-3]};
        case 'direction_1channel_intensity'
            stim = [3:.125:4.875; ones(1,16)*50e-3];
        case 'direction_1channel_8backandforth'
            stim = [3:4:15; ones(1,4)*2];
        case 'direction_1channel_lightpulse'
            stim = [1:.125:1.875; ones(1,8)*.05];
        case 'direction_1channel_longrun'
            stim = [5:20:785; ones(1,40)*.05];
        case 'multstim_012345'
            stim = [3:3:18; ones(1,6)*.03];
        case 'direction_1channel_rest2s-8Hz2s'
            stim = [2:.125:3.875; ones(1,16)*.05];
        case 'direction_1channel_rest2s-8Hz8s'
            stim = [2:.125:9.875; ones(1,64)*.05];
        case 'direction_1channel_pulse_delay250ms'
            stim = [.25; .05];
        otherwise
            error('unknown parameter file, help me to define ''stim''')
    end
    T.stim = stim;
else
    % old version or no stimulation
    T.stim = [];
end

% Recordings
recording = T.recording_header;
if isfield(AcqSet,'recordingchannels') && ~isempty(AcqSet.recordingchannels)
    r = AcqSet.recordingchannels;
    for i=1:length(r)
        recording(i).name = r{i};
    end
    [recording.t0] = deal(0);
    if ~isempty(r), [recording.dt] = deal(1/AcqSet.samplingrate); end
    [recording.n] = deal(-1); % this indicates that we do not know yet the size of the recordings
end

end



