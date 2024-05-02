function label = tpv_elphy_getstim(V) 
% function f(V) 
% function label = f('label') 

if nargin==0, V = evalin('base','V'); end

if ischar(V) && strcmp(V,'label')
    label = 'Get stims from Elphy';
    return
end

T = V.content.trials;
T = T([T.status]~='s'); % normal and rejected trials, but not special trials
if ~strfind(T(1).origin,'MESC'), error 'only MESC data is supported', end

f = fn_getfile({'*.DAT'; '*.txt'},'Select Elphy data file(s)');
f = cellstr(f);
nf = length(f);

trecord = cell(1,nf);
for k=1:nf
    switch lower(fn_fileparts(f{k},'ext'))
        case '.dat'
            [rec vect par] = elphy_read(f{k}); %#ok<NASGU>
            trecord{k} = vect.trecord(1:size(rec,2));
        case '.txt'
            s = importdata(f{k});
            t = s.data(:,2);
            nep = find(t>0,1,'last');
            if any(t(1:nep)==0) || any(t(nep+1:end)>0), error 'could not read stim definition', end
            trecord{k} = t(1:nep);
    end            
end

% Prepare the cell array to write into Excell file
fmes = T(1).file;
if ~strcmp(fn_fileparts(fmes,'ext'),'.mesc')
    fmes = fn_getfile('*.mesc','Select original MESC file');
end
ntrial = length(T);
info = [T.fullinfo];
kunit = [info.kunit]';
if ~isequal(fmes,0)
    [h sessions] = mesc_header(fmes);
    if length(h.sessions)>1, error 'cannot handle multiple sessions', end
    units = h.sessions.Units;
    nunit = length(units);
    comments = {units.Comment}';
    comments = fn_map(@(x)char(x(1:end-1)'),comments);
else
    nunit = max(kunit);
    comments = repmat({''},1,nunit);
    for i=1:ntrial
        if isfield(T(i).addinfo,'comment')
            comments{T(i).fullinfo.kunit} = T(i).addinfo.comment;
        end
    end
end
X = cell(1+nunit,1);
icol=0;
icol=icol+1;
X{1,icol} = 'trial';
X(1+kunit,icol) = num2cell((1:ntrial)');
icol=icol+1;
X{1,icol} = 'unit';
X(2:1+nunit,icol) = num2cell((1:nunit)');
icol=icol+1;
X{1,icol} = 'comment';
X(2:1+nunit,icol) = comments;
icol=icol+3;
C = char('A'+(icol-1));
X{1,icol} = 'write trial type here';

icol=icol+2;
for i=1:nf
    icol=icol+1;
    X{1,icol} = num2str(i,'DAT%i');
    nv = length(trecord{i});
    X(2:1+nv,icol) = num2cell(trecord{i});
end

% Write into Excel file
xls = fullfile(tempdir,'tpv_elphy_getstim.xls');
if exist(xls,'file'), delete(xls), end
xlswrite(xls,X)
% system(xls); % open the document in Excel, returns once user has closed
% [num txt raw] = xlsread(xls,1,[C '1:' C num2str(1+nunit)]); %#ok<ASGLU>
[num txt raw] = xlsread(xls,-1); %#ok<ASGLU>
trialtype = cell2mat(raw(2:end));

% Set stim
trialtype(isnan(trialtype)) = 0;
trialtype(end+1:nunit) = 0;
trialtype = trialtype(kunit);
for i=1:ntrial, T(i).stim = trialtype(i); end

