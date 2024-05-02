function [recdata frecout] = tps_readrecording(frec,recidx)
% function [recdata frecout] = tps_readrecording(frec[,recidx])
%---
% generic function to read analog recording
%
% Input:
% - frec    string - name of recording file, 
%           or cell array - name of recording file followed by additional
%           information, such as sweep index, or subselection of channels.
%           The recognized formats are:
%           'file_adaq.mat' or {'file_adaq.mat',kchannels}
%           'file.mat'
%           {'file.mes','measfieldname','recname1','recname2',...}
%           'file.dat'
%           {'file.abf',ksweep,kchannels}
% - recidx  vector of indices - subselection of channels
%
% Output:
% - recdata cell array of vectors - signals
% - frecout updated information (in case data file has been moved to
%           another location)


% Input
if ~iscell(frec)
    frec = {frec};
elseif iscell(frec{1})
    if ~isscalar(frec), error 'handling multiple analog data file not implemented yet', end
    frec = frec{1};
end
more = frec(2:end);
frec = frec{1};
frec = tps_locatefile(frec);
if isempty(more), frecout = frec; else frecout = {[frec more]}; end
if nargin<2, recidx=[]; end % no additional subselection of channels

% Read recording
ext = fn_fileparts(frec,'ext');
if ~isempty(ext), ext = lower(ext(2:end)); end
switch ext
    case {'' 'mat'}
        if strfind(frec,'_adaq')
            warnstate = warning('query','MATLAB:load:variablePatternNotFound');
            if isempty(more) && isempty(recidx)
                % in old version, channels were not saved together with the
                % file name
                nrec = 100;
                channels = 1:nrec;
                recdata = cell(1,0);
                warning('off','MATLAB:load:variablePatternNotFound')
            else
                if isempty(more)
                    channels = recidx;
                elseif isempty(recidx)
                    channels = more{1};
                else
                    channels = more{1}(recidx);
                end
                nrec = length(channels);
                recdata = cell(1,nrec);
            end
            i = 0;
            while i<nrec
                i = i+1;
                v = load(frec,['Cond*Trial0Chan' num2str(channels(i)-1)]);
                Fv = fieldnames(v);
                nv = length(Fv);
                if nv==0
                    % the number of channel was (i-1)
                    break
                elseif i==1
                    ncond = nv;
                elseif ncond~=nv
                    error 'problem with the number of conditions'
                end
                recdata{i} = cell2mat(row(struct2cell(v)));
            end
            warning(warnstate.state,'MATLAB:load:variablePatternNotFound')
        else
            try
                recdata = load(frec,'rec');
                recdata = recdata.rec;
            catch %#ok<CTCH>
                recdata = fn_loadvar(frec);
            end
            if ~iscell(recdata), recdata = {recdata}; end
            if ~isempty(recidx), recdata = recdata(recidx); end
        end
    case 'mes'
        measname = more{1};
        recnames = more(2:end);
        if ~isempty(recidx), recnames = recnames(recidx); end
        nrec = length(recnames);
        if nrec~=length(recnames), error programming, end
        recdata = cell(1,nrec);
        header = load(frec,measname,'-MAT');
        h1 = header.(measname)(1);
        for i=1:nrec
            recdata{i} = h1.(recnames{i}).y;
        end
    case 'dat'
        x = fn_readdatlabview(frec);
        if isempty(x)
            disp 'empty recording'
        end
        recdata = num2cell(x,1);
        if ~isempty(recidx) 
            if isempty(x)
                recdata = repmat({zeros(1,0)},1,length(recidx));
            else
                recdata = recdata(recidx); 
            end
        end
    case 'abf'
        if isempty(more)
            disp 'cannot read recording, sweep and channel not specified'
            recdata = {[]};
            return
        end
        [ksweep kchannel] = deal(more{:});
        if ~isempty(recidx), kchannel = kchannel(recidx); end
        x = abfload(frec,'sweeps',ksweep);
        recdata = num2cell(x(:,kchannel),1);
    case 'da'
        channels = neuroplex_read(frec,'channels');
        recdata = num2cell(channels,1);
        if ~isempty(recidx), recdata = recdata(recidx); end
    otherwise
        error('cannot read recording for ''%s'' extension',ext(2:end))
end


