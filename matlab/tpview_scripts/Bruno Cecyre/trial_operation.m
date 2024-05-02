function label = trial_operation(V) 
% function f(V) 
% function label = f('label') 

if nargin==0, V=evalin('base','V'); end
if ischar(V) && strcmp(V,'label')
    label = 'Trial Operation';
    return
end

% Write script code below
persistent smem

% Prompt user for operation definition
prompt = ['Please define trial operation below.' ...
    ' Use X(i) to denote either the ith trial or,' ...
    ' if i is a list of indices, the average trial over these indices.'];
s = struct( ...
    'data',     {'data'     {'data' 'dataop'}   'Perform on:'}, ...
    'label',    {[]         'label 2'           prompt}, ...
    'expr',     {'X(11) - X(1:10)'  'char'      'Operation:'});
if ~isempty(smem)
    s(1).data = smem.data;
    s(1).expr = smem.expr;
end
s = fn_structedit(s,'title','Trial Operation');
if isempty(s), return, end

% Perform the operation
% (analyze expression)
[istart iend matches] = regexp(s.expr,'X\([^\(\)]*(\([^\(\)]*\))*[^\(\)]*\)','start','end','match');
nmatch = length(matches);
% (gather requested trials)
T = cell(1,nmatch);
try
    for i=1:nmatch
        T{i} = V.content.trials(eval(matches{i}(3:end-1)));
    end
catch ME
    errordlg({'Expression is not valid, it produced the following error:' ...
        ME.message})
    return
end
% (replace in expression the matched patterns)
istart = [istart length(s.expr)+1];
expr = s.expr(1:istart(1)-1);
opflag = s.data(5:end);
try
    for i=1:nmatch
        expr = [expr 'fn_means(T{' num2str(i) '}.data' opflag ')' s.expr(iend(i)+1:istart(i+1)-1)]; %#ok<AGROW>
    end
    data = eval(expr);
catch ME
    errordlg({'Expression is not valid, it produced the following error:' ...
        ME.message})
    return
end
allT = [T{:}];
if all([allT.sfrchannel])
    expr = strrep(expr,'data','sfr');
    sfr = eval(expr);
    data = {data sfr};
end

% Create new trial
T = tps_trial(data,allT(1));
% T.opdef = allT(1).opdef;
if isfield(T.addinfo,'comment'), T.addinfo = rmfield(T.addinfo,'comment'); end
T.addinfo.description = s.expr;
T.status = 's'; % special trial
stimids = cat(1,allT.stimid);
if any(diff(stimids))
    disp 'stim definition of newly created trial might be erroneous'
end

% Add new trial to tpview
addtrials(V.content,T)
if ~isempty(V.disppar.trialsublist), V.disppar.trialsublist(end+1) = V.ntrial; end
display_slidertrial(V)
V.ktrial = V.ntrial;

% It worked, memorize parameters
smem = s;

