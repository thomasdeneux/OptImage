function label = newdatasetfromsignals(V) 
% function f(V) 
% function label = f('label') 

if ischar(V) && strcmp(V,'label')
    label = 'New dataset from signals';
    return
end


ntrial = V.ntrial;
trials = 1:ntrial;
nsel = V.content.nsel;

% compute signals
fn_progress('compute signals',ntrial)
for i=1:ntrial, fn_progress(i), getslice(V.content,trials(i)); end

% get signals
x = getslice(V.content,trials);
data0 = cell(ntrial,nsel);
for i=1:ntrial*nsel, data0{i} = permute(x(i).data,[3 4 1 2]); end % 1 x 1 x time x cond
data = cell(1,ntrial);
for i=1:ntrial, data{i} = cat(1,data0{i,:}); end % cell x 1 x time x cond

% make tps_trial object
T0 = V.content.trials(trials);
dt = T0(1).dt;
T = tps_trial(data,1,dt,'linescan','units','cell',T0(1).tunit);
% (more header info)
[T.status] = deal(T0.status);
[T.addinfo] = deal(T0.addinfo);
[T.stimtable] = deal(T0(1).stimtable);
for i=1:ntrial
    setstim(T(i),T0(i).stimid)
end

% open in tpview
tpview(V.skin,T)


