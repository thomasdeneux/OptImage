function label = compare_conditions(V) 
% function f(V) 
% function label = f('label') 

if nargin==0, V=evalin('base','V'); end
if ischar(V) && strcmp(V,'label')
    label = 'compare conditions';
    return
end

% Write script code below

%% Get the data

% which trials are non-rejected individual trials (not average)
T = V.content.trials;
status = [T.status];
idxtrial = find(status=='n');

% get the signals for these trials
x = V.getsignals(idxtrial,[]);
[ntrial nroi] = size(x);

% check whether there are 2 conditions selected
ncond = size(x(1).dataop,2);
if ncond==1
    errordlg('You must select several conditions in ''Time Courses'' to perform comparison.')
    return
elseif ncond==2
    idxcond = [1 2];
elseif ncond>2
    idxcond = [];
    while length(idxcond)~=2
        idxcond = fn_input(['Which conditions out of the ' num2str(ncond) ' displayed should be compared:'], ...
            [1 2]);
        if isempty(idxcond), disp 'interrupted', return, end
    end
end

% ask which ROI if needed
if nroi>1
    idxroi = fn_input('Please select which ROI to use:',1,1,nroi);
    if isempty(idxroi), disp 'interrupted', return, end
    x = x(:,idxroi);
end

% get the signals from the tps_signalx object
x = cat(3,x.dataop); % we create a 3D array: nframe*ncond*ntrial
x = x(:,idxcond,:);  % we subselect the 2 conditions to be compared -> nframe*2*ntrial

%% Get the time period over which to perform averaging

a4d = V.a4d;        % structure with information from all the interactive displays
SIt = a4d.SIt;      % information about what is displayed in the "time courses"
timeselection = SIt.selectionmarks;     % information about the temporal selection if there is one
if ~isempty(timeselection)
    answer = questdlg('Use the temporal window selected in OptImage?','','Yes','No','Yes');
    if strcmp(answer,'Yes')
        idxtime = timeselection.dataind;
    else
        timeselection = []; 
    end
end
if isempty(timeselection)
    dt = V.content.trials(idxtrial(1)).dt_sec;
    if isempty(dt)
        errordlg 'Please set frame duration first'
    end
    timesegment = [];
    while length(timesegment)~=2
        timesegment = fn_input('Please select time window for averaging (in seconds):',[0 1]);
        if isempty(timesegment), disp 'interrupted', return, end
    end
    nt = size(x,1);
    tt = (0:nt-1)*dt;
    idxtime = find(tt>=timesegment(1) & tt<=timesegment(2));
end

%% Perform temporal averaging

x = x(idxtime,:,:); % subselect
x = mean(x,1);      % average, x is now 1*2*ntrial
x = shiftdim(x,1);  % x is now 2*ntrial

%% Perform the statistical comparison

% compute the p-value to compare the median of two empirical distributions
x1 = x(1,:);
x2 = x(2,:);
pvalue = ranksum(x1,x2);

% try to display the result nicely
hf = figure('name','Compare 2 conditions');
plot(ones(1,ntrial),x1,'o')     % individual measures for 1st condition
hold on
plot(ones(1,ntrial)*2,x2,'o')   % for 2nd condition
plot([1 2],[mean(x1) mean(x2)],'color','r')
plot([1 2],[median(x1) median(x2)],'color','k','linewidth',2)
title(sprintf('p-value = %.2g',pvalue))

% get the name of the 2 compared conditions
condname = V.content.datacond;
if strcmp(condname,'all conditions')
    condname1 = ['C' num2str(idxcond(1)-1)];
    condname2 = ['C' num2str(idxcond(2)-1)];
else
    condnames = fn_strcut(condname,','); % separate condition names where there are commas
    condname1 = condnames{idxcond(1)};
    condname2 = condnames{idxcond(2)};
end

% still improve the esthetics
set(gca,'xlim',[.3 2.7],'xtick',[1 2],'xticklabel',{condname1 condname2})















