classdef tps_signal < hgsetget
% tps_signal is just a container; there is absolutely no checking whether
% the size of the data structure fits the number of trials
    properties (Access='private')
        version
    end
    properties
        % when creating new properties DON'T FORGET to update the COPY
        % function
        name = '';
        datamode = '';
        datacond
        dataopdef = struct('op',cell(1,0),'link',cell(1,0));
        opdef   
        seldotrial = true;
        timeline = false;       % view time courses as if from line scan
        x = tps_signalx.empty;  % ntrials x nsel - fields can be added
        shift = zeros(0,2);     % how did the regions move compared to the first trial
        spikepar       % default parameters for spike estimation
    end
    properties (Dependent, SetAccess='private')
        nx
        nfr
        nfrop
        nexp
        data
        data2
        dataop
        data2op
        sel
        spikes
    end
    
    % Constructor
    methods
        function U = tps_signal
            U.version = 1.1;
            U.opdef = tpv_content.signalopdef_default();
        end
    end
    
    % GET
    methods
        function nx = get.nx(U)
            nx = size(U.x,2);
        end
        function nfr = get.nfr(U)
            if isempty(U.x), nfr=[]; return, end
            nfr = length(U.x(1).data);
            if nfr==0, nfr=[]; end
            for k=2:numel(U.x)
                if length(U.x(k).data)~=nfr, nfr=[]; return, end
            end
        end
        function nfr = get.nfrop(U)
            if isempty(U.x), nfr=[]; return, end
            nfr = length(U.x(1).dataop);
            if nfr==0, nfr=[]; end
            for k=2:numel(U.x)
                if length(U.x(k).dataop)~=nfr, nfr=[]; return, end
            end
        end
        function nexp = get.nexp(U)
            nexp = size(U.x,1);
        end
        function data = get.data(U)
            data = reshape({U.x.data},U.nexp,U.nx);
            nfr = U.nfr;
            if ~isempty(nfr)
                data = cat(3,data{:});
                if ~isempty(data)
                    data = reshape(data,[nfr size(data,2) U.nexp U.nx]);
                end
            end
        end
        function data2 = get.data2(U)
            data2 = reshape({U.x.data2},U.nexp,U.nx);
            nfr = U.nfr;
            if ~isempty(nfr)
                data2 = cat(3,data2{:});
                if ~isempty(data2)
                    data2 = reshape(data2,nfr,size(data2,2),U.nexp,U.nx);
                end
            end
        end
        function dataop = get.dataop(U)
            dataop = reshape({U.x.dataop},U.nexp,U.nx);
            nfr = U.nfrop;
            if ~isempty(nfr)
                dataop = cat(3,dataop{:});
                if ~isempty(dataop)
                    dataop = reshape(dataop,nfr,size(dataop,2),U.nexp,U.nx);
                end
            end
        end
        function data2op = get.data2op(U)
            data2op = reshape({U.x.data2op},U.nexp,U.nx);
            nfr = U.nfrop;
            if ~isempty(nfr)
                data2op = cat(3,data2op{:});
                if ~isempty(data2op)
                    data2op = reshape(data2op,nfr,size(data2op,2),U.nexp,U.nx);
                end
            end
        end
        function sel = get.sel(U)
            sel = reshape([U.x.sel],U.nexp,U.nx);
        end
        function spikes = get.spikes(U)
            spikes = reshape({U.x.spikes},U.nexp,U.nx);
        end
        function U2 = copy(U)
            U2 = tps_signal;
            U2.name = U.name;
            U2.datamode = U.datamode;
            U2.datacond = U.datacond;
            U2.dataopdef = U.dataopdef;
            U2.opdef = fn_structmerge(tpv_content.signalopdef_default, ...
                U.opdef,'skip'); % correct for changes in the syntax
            U2.seldotrial = U.seldotrial;
            U2.shift = U.shift;
            U2.timeline = U.timeline;
            U2.x = copy(U.x);
        end
    end
    
    % Load
    methods (Static)
        function U = loadobj(U)
            % handle objects saved in previous versions
            
            % shift
            if isempty(U.shift)
                U.shift    = zeros(U.nexp,2);
            end
            % time courses
            if ~isa(U.x,'tps_signalx')
                disp('convert structure to tps_signalx')
                x = U.x;
                [ntrial nsel] = size(x);
                U.x = tps_signalx; U.x(ntrial,nsel).sel = [];
                F = fields(x);
                for k=1:length(F)
                    f = F{k};
                    [U.x.(f)] = deal(x.(f));
                end
            end
            % signal operation
            op = U.opdef;
            if isfield(op,'normalization')
                normflag = fn_switch(op.normalization, ...
                    '(none)',       [0 0 0], ...
                    'avg. value',   [0 0 1], ...
                    'neuropil',     [0 1 0], ...
                    'neuropil-avg', [0 1 1], ...
                    'avg-neuropil', [1 1 0]);
                U.opdef.norm_npil_norm = logical(normflag);
            end
            if isfield(op,'norm_npil_norm')
                U.opdef.normalize__by__mean = op.norm_npil_norm(1) | op.norm_npil_norm(3);
            end
            U.opdef = fn_structmerge(tpv_content.signalopdef_default,U.opdef,'skip');
            % signals and spikes
            for k=1:numel(U.x)
                xk = U.x(k);
                % spike parameter
                if ~isempty(xk.spikepar) && ~isstruct(xk.spikepar)
                    xk.spikepar = [];
                end
            end
            % now we have a 'version' property!
            if isempty(U.version), U.version = 0; end
            % changed dataopdef
            if U.version<1
                if ~iscell(U.dataopdef) && ~isa(U.dataopdef,'tps_dataopdef')
                    U.dataopdef = tps_dataopdef(U.dataopdef);
                end
                U.dataopdef = repmat(struct('op',U.dataopdef,'link',true),1,U.nexp);
            elseif U.version<1.1
                U.dataopdef = struct('op',U.dataopdef,'link',true);
            end
            U.version = 1.1;
        end
    end
        
end