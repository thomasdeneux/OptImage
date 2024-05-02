function fname = tps_locatefile(fname,fexample)
% function fnamenew = tps_locatefile(fnameold,fexample)
%---
% if fname does not exist, try to locate the data by substituting the
% path in various ways
%
% Input:
% - fnameold    char array or cell array of char arrays
% - fexample    a possibly new containing folder, or a file inside a
%               possibly containing folder
%
% Output:
% - fnamenew    same format as fnamenew

persistent tdisable
if isempty(tdisable), tdisable = uint64(0); end

if nargin<2, fexample = ''; end

% make a cell array
isc = iscell(fname);
if ~isc, fname={fname}; end
ncell = length(fname);

% try to locate a first file
fil = deblank(fname{1}(1,:));
fil = strrep(fil,'\','/');
[p1 name ext] = fn_fileparts(fil,'path','name','ext');
if exist(fil,'file') && ~isempty(p1)
    % file exists, nothing to do
    if ~isc, fname=fname{1}; end
    return
end 

% try several different folders
totry = {'example file','lastdir'};
[lastbaddir lastgooddir] = fn_userconfig('tps_locatefile');
okdir = false;
for k=1:length(totry)
    switch totry{k}
        case 'example file'
            if isempty(fexample), continue, end
            % cut file name into path components
            p1parts = regexp(fil,'([^/])*');
            p1end = p1parts(end)-2;
            fexample = strrep(fexample,'\','/');
            if exist(fexample,'dir')
                p2base = fexample;
            else
                p2base = fileparts(fexample);
            end
            for i=length(p1parts):-1:1
                p2 = fullfile(p2base,fil(p1parts(i):p1end));
                %disp(['try ' p2])
                okdir = exist(fullfile(p2,name),'file');
                if okdir, break, end
                %disp 'not good'
            end
        case 'lastdir'
            % try the last substitution
            n = length(lastbaddir);
            if ischar(lastbaddir) && length(p1)>=n && strcmp(p1(1:n),lastbaddir)
                p2 = [lastgooddir p1(n+1:end)];
            else
                continue
            end
    end
    okdir = ischar(p2) && exist(fullfile(p2,name),'file');
    if okdir, break, end
end

% prompt user
if okdir
    okfile = true;
elseif toc(tdisable)>10 && feature('ShowFigureWindows')
    %fprintf('no file with path ''%s'', please select new location\n',fil)
    fil = fn_getfile({name; ['*' ext]},sprintf('Locate file %s (previously %s)',name,fil));
    okfile = ischar(fil) && exist(fil,'file');
    if okfile
        [p2 newname] = fn_fileparts(fil,'path','name');
        okdir = okfile && strcmp(newname,name); % name has not changed
    else
        usercancel = true;
    end
else
    okfile = false;
    usercancel = false;
end

% if problem solved, update all files
if okdir
    if isempty(p1) && p2(end)~='/', p2 = [p2 '/']; end
    for i=1:ncell
        fi = fname{i};
        fi = strrep(fi,'\','/');
        fname{i} = [repmat(p2,size(fi,1),1) fi(:,length(p1)+1:end)];
    end
elseif okfile
    if ncell>1 || ~isvector(fname{1})
        error 'cannot handle multiple files, please update code'
    end
    fname{1} = fil;
    disp('Updated file name')
elseif usercancel
    uiwait(msgbox({'''Cancel'' has been pressed while searching for missing data file:' ...
        'There will be no additional user prompt for missing files during the next 10 seconds.'}, ...
        'MESSAGE'))
    drawnow % necessary to get the message box cleared before upcoming calculation are performed
    tdisable = tic;
end

% remember the substitution
if okdir || (okfile && strcmp([name ext],fn_fileparts(fil,'name')))
    if ~okdir, p2 = fn_fileparts(fil,'path'); end
    p1 = strrep(p1,'\','/'); p2 = strrep(p2,'\','/');
    n1 = length(p1); n2 = length(p2);
    n = min(n1,n2);
    lastdiff = find(p1(n1-n+1:n1)~=p2(n2-n+1:n2),1,'last');
    if isempty(lastdiff), lastdiff=0; end
    p1 = p1(1:n1-n+lastdiff); p2 = p2(1:n2-n+lastdiff);
    if ~strcmp(p1,lastbaddir) || ~strcmp(p2,lastgooddir)
        fprintf('replaced folder ''%s'' by ''%s''\n',p1,p2)
        fn_userconfig('tps_locatefile',p1,p2)
    end
end

% output
if ~isc, fname = fname{1}; end


