function [calculations calcname] = tps_definecondcalc(T)

% check
if ~isscalar(T)
    error 'argument must be a single trial'
end
if T.nc<=1
    errordlg('trial has only one condition, cannot define an operation on conditions')
    calculations = '';
    calcname = '(none)';
    return
end

% a pointer to all the relevant information
info = fn_pointer;
info.T = T;

% make a nice stim table (only the conditions present in the trial, names,
% etc.)
makestimtable(info);

% init controls
initcontrols(info);

% additional information needed to build up condition name
initmem(info)

% wait!
waitfor(info.ctrl.ok)

% output
if ~ishandle(info.ctrl.hf)
    % figure has been closed
    calculations = '';
    calcname = '(none)';
else
    calculations = get(info.ctrl.calcs,'string');
    calcname = get(info.ctrl.calcname,'string');
    close(info.ctrl.hf)
end

%---
function makestimtable(info)

% trial
T = info.T;

% select only conditions present in the trial
ids = T.stimid;
if length(ids)~=T.nc, error 'T.stimid is not properly defined', end
ncond = T.nc;
stimtable = T.stimtable.table;
tableids = [stimtable.id];
entries = zeros(1,ncond); 
for k=1:ncond, entries(k) = find(tableids==ids(k)); end
stimtable = stimtable(entries);

% handle names
if isfield(stimtable,'name')
    names = {stimtable.name};
else
    names = cell(1,ncond);
end
for k=1:ncond
    if isempty(names{k}), names{k} = num2str(k-1,'C%i'); end
end

% list of relevant F
F = fieldnames(stimtable);
F(ismember(F,{'id' 'stim' 'name'})) = [];

% make everything text
for k=1:ncond
    for i=1:length(F)
        f = F{i};
        if ~ischar(stimtable(k).(f)), stimtable(k).(f) = num2str(stimtable(k).(f)); end
    end
end

% set info
info.ncond = ncond;
info.nfilt = length(F);
info.stimtable = stimtable;
info.names = names;
info.F = F;


%---
function initcontrols(info)

% info on conditions
stimtable = info.stimtable;
names = info.names;
F = info.F;

% sizes
ncond = length(stimtable);
nfilt = length(F);

% Size parameters
% (general sizes)
htext = 18;
hctrl = 22;
dtext = (hctrl-htext)/2;
hlist = 15;
D = 5;          % basic separation
DD = 15;
dfilt = 2;
dinner = 4;
% (specific heights)
hpmain = dinner + hctrl + dinner;
hfiltu = hctrl+dfilt;
hfilt  = nfilt*hfiltu - dfilt;
hcond  = min(500,ncond*hlist);
hpcond = dinner + hfilt + D + hcond + dinner;
hfullcalc = 3*htext;
H = DD + hpmain + D + hpcond + D + hfullcalc + D + hfullcalc + DD;
% (specific widths)
wfiltlabel = 100;
wfilt = 150;
wcond = wfiltlabel + D + wfilt;
wop = 40;
wpanel = D + wcond + D + wop + D + wcond + D;
wcondnum = 150;
wbut = 60;
W = D + wpanel + D;

% Figure
hf = 436;
ctrl.hf = hf;
if ishandle(hf), close(hf), end
figure(hf), clf(hf)
set(hf,'numberTitle','off','name','Define calculated condition(s)', ...
    'menubar','none','Resize','off')

fn_setfigsize(hf,W,H)

% Main menu panel
pmain = uipanel('parent',hf,'units','pixel','pos',[D H-DD-hpmain wpanel hpmain]);
ctrl.condnum = uicontrol('parent',pmain,'style','popupmenu', ...
    'pos', [D dinner wcondnum hctrl], ...
    'string',{'computation 1'},'value',1, ...
    'callback',@(u,e)chgcondnum(info));
ctr.condnumadd = uicontrol('parent',pmain,'string','add', ...
    'pos', [D+wcondnum+D dinner wbut hctrl], ...
    'callback',@(u,e)chgcondnum(info,'add'));
ctr.condnumrm = uicontrol('parent',pmain,'string','delete', ...
    'pos', [D+wcondnum+D+wbut+D dinner wbut hctrl], ...
    'callback',@(u,e)chgcondnum(info,'delete'));

% Condition definition panel
pcond = uipanel('parent',hf,'units','pixel','pos',[D DD+hfullcalc+D+hfullcalc+D wpanel hpcond]);
% (filtering)
ctrl.filt = zeros(2,nfilt);
for k=1:2
    xstart = D + fn_switch(k,1,0,2,wcond+D+wop+D);
    for i=1:nfilt
        f = F{i};
        choices = unique({stimtable.(f)});
        ystart = hpcond-dinner-i*hfiltu+dfilt;
        uicontrol('parent',pcond,'style','text','string',f, ...
            'horizontalalignment','left', ...
            'pos', [xstart ystart+dtext wfiltlabel htext]);
        ctrl.filt(k,i) = uicontrol('parent',pcond,'style','popupmenu', ...
            'pos', [xstart+wfiltlabel+D ystart wfilt hctrl], ...
            'string',['no filtering' choices],'value',1, ...
            'callback',@(u,e)readfilters(info,k));
    end
end
% (conditions)
for k=1:2
    xstart = D + fn_switch(k,1,0,2,wcond+D+wop+D);
    ctrl.condlist(k) = uicontrol('parent',pcond,'style','listbox', ...
        'pos', [xstart dinner wcond hcond], ...
        'string', names, 'userdata', 0:info.ncond-1, ...
        'max', 2, 'value', [], ...
        'callback', @(u,e)readconds(info));
end
% (operation)
ctrl.operation = uicontrol('parent',pcond,'style','popupmenu', ...
    'pos', [D+wcond+D (hpcond-hctrl)/2 wop hctrl], ...
    'string', {'/','-'}, 'value', 1, ...
    'callback', @(u,e)readconds(info));

% Condition calculation and name
ctrl.calcs = uicontrol('parent',hf,'style','text', ...
    'pos', [D DD+hfullcalc+D wpanel hfullcalc], ...
    'max', 2, 'horizontalalignment','left');
ctrl.calcname = uicontrol('parent',hf,'style','edit', ...
    'backgroundcolor','w', ...
    'pos', [D DD wpanel hfullcalc], ...
    'max', 2, 'horizontalalignment','left');

% OK button
ctrl.ok = uicontrol('parent',hf,'string','OK', ...
    'pos', [W-D-wbut DD wbut hctrl], ...
    'callback',@(u,e)terminate(info));

% set info
info.ctrl = ctrl;

%---
function initmem(info)

info.ncalc = 1;
info.kcalc = 1;
% store the initial (empty) calculation
info.calculations = struct('ctrlvalues',struct,'filtering',{cell(1,2)},'calc','','name','');
readfilters(info)           % fills 'filtering'
readconds(info)             % fills 'calc' and 'name'
memorizectrlvalues(info)    % fills 'ctrlvalues'
info.defcalc = info.calculations;

%---
function memorizectrlvalues(info)

ctrl = info.ctrl;
info.calculations(info.kcalc).ctrlvalues = struct( ...
    'filt',     {fn_get(ctrl.filt,'value')}, ...
    'cond',     {fn_get(ctrl.condlist,{'string' 'userdata' 'value'})}, ...
    'op',       get(ctrl.operation,'value') ...
    );

%---
function readfilters(info,kside)

if nargin<2, kside = 1:2; end
ctrl = info.ctrl;
filters = ctrl.filt;
for k=kside
    okcond = true(1,info.ncond);
    filtname = cell(1,info.nfilt);
    for i=1:info.nfilt
        f = info.F{i};
        val = get(filters(k,i),'value');
        if val==1, continue, end
        choices = get(filters(k,i),'string');
        value = choices{val};
        okcond = okcond & strcmp({info.stimtable.(f)},value);
        filtname{i} = [value(1:min(5,end))];
    end
    nofilter = all(okcond);
    set(ctrl.condlist(k), ...
        'string',   info.names(okcond), ...
        'userdata', find(okcond)-1, ...
        'value',    fn_switch(nofilter,[],1:sum(okcond)))
    if nofilter, filtname = 'all'; else filtname = fn_strcat(filtname,'_'); end
    info.calculations(info.kcalc).filtering{k} = filtname;
end
readconds(info)

%---
function readconds(info)

ctrl = info.ctrl;
calcs = cell(1,2); % first row is calculation, second row is smart name
calcnames = cell(1,2);
for k=1:2
    nums = get(ctrl.condlist(k),'userdata');
    sel = get(ctrl.condlist(k),'value');
    val = nums(sel);
    if isempty(val)
        str = '';
    else
        str = sprintf('C%i+',val);
        str(end) = [];
    end
    calcs{k} = str;
    if length(val)==length(nums)
        % all filtered conditions are selected: use filtering name
        calcnames{k} = info.calculations(info.kcalc).filtering{k};
    else
        % user selection: print the full selection
        calcnames{k} = str;
    end   
end
if isempty(calcs{2})
    calcs = calcs{1};
    calcnames = calcnames{1};
else
    ops = get(ctrl.operation,'string');
    op = ops{get(ctrl.operation,'value')};
    calcs = ['(' calcs{1} ')' op '(' calcs{2} ')'];
    calcnames = ['(' calcnames{1} ')' op '(' calcnames{2} ')'];
end
info.calculations(info.kcalc).calc = calcs; 
info.calculations(info.kcalc).name = calcnames; 
set(ctrl.calcs,'string',fn_strcat({info.calculations.calc},','))
set(ctrl.calcname,'string',fn_strcat({info.calculations.name},','))

%---
function chgcondnum(info,flag)

% memorize the current calculation
memorizectrlvalues(info)

% update set of calculations
ctrl = info.ctrl;
if nargin==2
    % update set of calculations
    switch flag
        case 'add'
            info.ncalc = info.ncalc+1;
            info.kcalc = info.ncalc;
            info.calculations(info.ncalc) = info.defcalc;
        case 'delete'
            if info.ncalc==1
                % reset the unique calculation
                info.calculations(1) = info.defcalc;
            else
                % delete the current calculation
                info.calculations(info.kcalc) = [];
                info.ncalc = info.ncalc-1;
                info.kcalc = min(info.kcalc,info.ncalc);
            end
    end
    % update popup menu
    set(ctrl.condnum, ...
        'string', fn_map(@(x)sprintf('computation %i',x),1:info.ncalc,'cell'), ...
        'value', info.kcalc)
else
    info.kcalc = get(ctrl.condnum,'value');
end

% update controls
memk = info.calculations(info.kcalc).ctrlvalues;
%fn_set(ctrl.filt,'value',memk.filt)
fn_set(ctrl.condlist,{'string' 'userdata' 'value'},memk.cond)
fn_set(ctrl.operation,'value',memk.op)
set(ctrl.calcs,'string',fn_strcat({info.calculations.calc},','))
set(ctrl.calcname,'string',fn_strcat({info.calculations.name},','))
readfilters(info)

%---
function terminate(info)

delete(info.ctrl.ok)

