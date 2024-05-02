function Treg = tps_register(T,varargin)
% function Treg = tps_register(T[,fname])
%---
% Input:
% - T       tps_trial object
% - fname   dirname or {dirname suffix} or {dirname prefix suffix} - model
%           for naming the individual files where to save the resampled
%           data [default: if fname is not specified, resampled data is not
%           saved; default prefix is ''; default suffix is .reg.mat]
%
% Output:
% - Treg    tps_trial object with all data now aligned (registering and
%           resampling has been applied)


if nargin==0, help tps_register, return, end

% Input
dosave = false; cursel = [];
for k=1:length(varargin)
    a = varargin{k};
    if ischar(a)
        dosave = true;
        fname = a;
        if ischar(fname), fname = {fname}; end
        fbase = fn_fileparts(fname{1},'noext');
        savedir = [fbase '_data'];
        if ~exist(savedir,'dir'), mkdir(savedir), end
        if length(fname)>=2, suffix = fname{end}; else suffix = '.reg'; end
        if isempty(regexp(suffix,'.mat$', 'once')), suffix = [suffix '.mat']; end
        if length(fname)==3, prefix = fname{2}; else prefix = ''; end
    elseif iscell(a)
        cursel = a;
    end
end

% Default output
Treg = []; % empty output will signal an error
    
% Parameters
useprevious = false;
if dosave
    parfile = fullfile(savedir,'regpar.mat');
    if exist(parfile,'file')
        answer = questdlg('Previous registration results exist, use them?');
        switch answer
            case 'Cancel'
                Treg = [];
                return
            case 'Yes'
                useprevious = true;
                s = fn_loadvar(parfile);
        end
    end
end
if ~useprevious
    if isempty(cursel) || isempty(cursel{2})
        reftrial = 1;
        refrange = 1:10;
    else
        [reftrial refrange] = deal(cursel{:});
    end    
    s = struct( ...
        'a',            {[]     'label'     'Define registration action'}, ...
        'data',         {'data' {'data' 'dataop'}   'Perform on'}, ...
        'trials',       {1:length(T)    'double'    'Trials'}, ...
        'regtype',      {'every frame to a unique reference' ...
        {'trial averages only' 'every frame to a unique reference' 'every frame to a trial-specific reference'  'every frame to a trial-specific reference, then register trials'} 'Register:'}, ...
        'iterativeref', {false  'logical'   'Recursive reference'}, ...
        'reftrial',     {reftrial 'double'  'Reference trial (can be several trials)'}, ...
        'refrange',     {refrange 'double'  'Range of frames for the reference (0 = all)'}, ...
        'doroi',        {false  'logical'   'Use a sub-region of interest'}, ...
        ... 'dorepeat',     {false  'logical'   'repeat with updated ref'}, ...
        'b',            {[]     'label'     'Algorithm parameters'}, ...
        'xlowpass',     {0      'double'    'Spatial low-pass (px)'}, ...
        'xhighpass',    {0      'double'    'Spatial high-pass (px)'}, ...
        'tsmoothing',   {0      'double'    'Temporal smoothing of drifts (frame)'}, ...
        'maxshift',     {20     'double'    'Maximum shift distance (% of full size)'}, ...
        'doxreg',       {false  'logical'   'Start each frame with a global registration'}, ...
        'tolx',         {1e-3   'double'    'Tolx'}, ...
        'tolfun',       {1e-4   'double'    'Tolfun'}, ...
        'c',            {[]     'label'     'Output'},...
        'dosavereg',    {false  'logical'   'Save registered data in .mat files (if not, will be recomputed on the fly when needed)'}, ...
        'docut',        {true   'logical'   'Cut to valid pixels'}, ...
        'suffix',       {suffix 'char'      'File suffix for coregistered trials'}, ...
        'datatype',     {'single'   {'double' 'single' 'uint16' 'uint8'}    'Data type'});
    if any(fn_isemptyc({T.file}))
        % Data of the original trials are not saved -> we need to save the
        % data of the corrected trials
        s(1).dosavereg = true;
    end
    s = fn_structedit(s);
    if isempty(s), return, end
    fn_savevar(parfile,s)
end
% (algo parameter)
regtype = fn_switch(s.regtype, ...
    'trial averages only',                                              'trial', ...
    'every frame to a unique reference',                                'frameglobal', ...
    'every frame to a trial-specific reference',                        'framelocal', ...
    'every frame to a trial-specific reference, then register trials',  'framelocalglobal');
regpar = fn_register('par');
regpar.xlowpass = s.xlowpass;
regpar.xhighpass = s.xhighpass;
regpar.tsmoothing = s.tsmoothing;
regpar.maxshift = s.maxshift/100;
regpar.tolx = s.tolx;
regpar.tolfun = s.tolfun;
regpar.doxreg = s.doxreg;
% (file)
if dosave
    filesorig = {T.file};
    suffix = s.suffix;
    % use an additional number suffix for trials?
    donum = any(fn_isemptyc(filesorig)) ...        % at least one of the trial does not have a file name
        || fn_find(filesorig,@ismatrix,'any') ...  % at leat one of the trial is a "multi-file" trial
        || length(unique(filesorig))<length(T);    % at least two trials have the same file name
end

% Checks
T = T(s.trials);
%if s.dorepeat && ~strcmp(regtype,'framelocal'), errordlg 'ref updating applies only to ''frame local'' registration type', return, end
if ~strcmp(regtype,'framelocal')
    if any(any(diff([T.nx; T.ny],1,2))), errordlg 'trial dimensions must match', return, end
    if any(any(diff([T.dx],1,2))), errordlg 'trial spatial resolutions must match', return, end
end
if any([T.sfrchannel])
    answer = questdlg('Secondary channel will not be coregistered. Proceed?','','Yes','Cancel','Yes');
    if ~strcmp(answer,'Yes'), return, end
end

% Size parameters
ntrial = length(T);
nx1 = T(1).nx;
ny1 = T(1).ny;
ntall = [T.nc].*[T.nfr];
ntmax = max(ntall);

% deprecated 'dormline'
dormline = strcmp(T(1).origin,'ScanImage'); % last line in ScanImage acquisitions is bad 
if dormline && fn_dodebug
    error 'this case should be programmed more nicely'
    ny1 = ny1-1;  %#ok<UNRCH>
end

% Get the reference frame if needed and the mask
reffile = fullfile(savedir,'reference.mat');
if useprevious && exist(reffile,'file')
    [ref doinitialcut oki okj mask] = fn_loadvar(reffile);
else
    if ismember(regtype,{'trial' 'global'}) || s.doroi
        if s.iterativeref
            REF = iterativeReference(T,s,regpar);
            if isempty(REF), return, end
            nreftrial = 1;
        else
            nreftrial = length(s.reftrial);
            REF = zeros(T(1).nx,T(1).ny,nreftrial,fn_switch(s.datatype,'double','double','single'));
            fn_progress('get reference frame',nreftrial)
            for i=1:nreftrial
                fn_progress(i)
                datak = T(s.reftrial(i)).(s.data);
                datak = datak(:,:,:,1); % first condition
                if ~isequal(s.refrange,0) && ~strcmp(regtype,'trial')
                    datak = datak(:,:,s.refrange);
                end
                REF(:,:,i) = mean(datak,3);
            end
        end
        % (need to cut a sub-frame from the images?)
        if any(isnan(REF(:)))
            oki = ~all(any(isnan(REF),3),2);
            okj = ~all(any(isnan(REF),3),1);
            doinitialcut = true;
            REF = REF(oki,okj,:);
        else
            doinitialcut = false;
            [oki okj] = deal([]);
        end
        % (register several reference frames?)
        if nreftrial>1
            regpar.ref = mean(REF,3);
            [dum dum XREF] = fn_register(REF,regpar); %#ok<ASGLU>
            ref = nmean(XREF,3);
            clear REF XREF
        else
            ref = REF;
        end
    else
        % no need for a global reference at this stage for 'local' and
        % 'localglobal' registration types
        ref = [];
        doinitialcut = false;
        [oki okj] = deal([]);
    end
    % (region of interest)
    if s.doroi
        mask = fn_maskselect(ref,'poly',false,'gray');
    else
        mask = [];
    end
    % (save)
    if dosave
        fn_savevar(reffile,ref,doinitialcut,oki,okj,mask)
        if ~isempty(ref), fn_saveimg(ref,fullfile(savedir,'reference.png')), end
    end
end
if fn_ismemberstr(regtype,{'trial' 'frameglobal'}), regpar.ref = ref; end
regpar.mask = mask;


% File for saving drift estimation
driftfile = fullfile(savedir,'estimateddrifts.mat');
if useprevious && exist(driftfile,'file')
    [SHIFT AVREG] = fn_loadvar(driftfile);
    kstart = find(isnan(SHIFT(1,1,:)),1,'first');
    if isempty(kstart) && ~strcmp(regtype,'framelocalglobal'), clear AVREG, end % corrected average frames no more needed
else
    kstart = 1;
    SHIFT = nan(2,ntmax,ntrial); % registration parameters
    AVREG = cell(1,ntrial); % corrected average frames
end

% Prepare display
if ~isempty(kstart) || strcmp(regtype,'framelocalglobal')
    hf = 207;
    figure(hf), colormap gray
    set(hf,'numbertitle','off','name','tps_register','tag','tps_register')
    clf(hf)
    ha = subplot(411);
    hlall = plot(.5+(.5:ntmax*ntrial)/ntmax,SHIFT(:,:)','hittest','off');
    hl = fn_spikedisplay(1.5:ntrial,0,2*(nx1+ny1),'color',[1 1 1]*.5,'hittest','off');
    hl(hl==0)=[];
    uistack(hl,'bottom')
    set(ha,'xlim',[.5 ntrial+.5])
    mM = [min(SHIFT(:)) max(SHIFT(:))];
    if ~any(isnan(mM)) && diff(mM)>0, set(ha,'ylim',mM), end
    xlabel trial
    ylabel 'shift (pixel)'
    title 'estimated shift for all trials'
    hb = subplot(412);
    hlcur = plot((0:ntmax-1)*T(1).dt,SHIFT(:,:,1)');
    xlabel(['time (' T(1).tunit ')'])
    ylabel 'shift (pixel)'
    legend(hlcur,'x','y')
    set(ha,'buttondownfcn',@(u,e)trialdisplay(get(ha,'currentpoint')))
    subplot(234), im0 = imagesc(fn_clip(ref','prc1')); axis image, grid on, set(gca,'xticklabel',[],'yticklabel',[]), title 'ref'
    subplot(235), im1 = imagesc(fn_clip(ones(nx1,ny1)','prc1')); axis image, grid on, set(gca,'xticklabel',[],'yticklabel',[]), title 'aligned trial'
    cm = mapgeog;
    subplot(236), im2 = imagesc(fn_clip(ones(nx1,ny1)','prc1',cm)); axis image, grid on, set(gca,'xticklabel',[],'yticklabel',[]), title 'diff'
    drawnow
end
function trialdisplay(ktrial)
    fn_set(hlall,'ydata',num2cell(SHIFT(:,:)',1))
    set(ha,'ylim',[-1 1]*max(abs(SHIFT(:))))
    ref1 = fn_normalize(ref,[1 2],'std');
    if ktrial==0
        % registration result over all trials
        hlcur = plot(1:ntrial,shifttrials,'parent',hb);
        xlabel(hb,'trial')
        fn_axis(hb,'tight',[1 1.1])
        avreg1 = fn_normalize(ref,[1 2],'std');
    else
        % registration result for a specific trial
        ktrial = fn_coerce(round(ktrial(1)),1,ntrial);
        fn_set(hlcur,'ydata',num2cell(SHIFT(:,:,ktrial),2))
        fn_axis(hb,'tight',[1 1.1])
        avreg1 = fn_normalize(AVREG{ktrial},[1 2],'std');
    end
    set(im0,'cdata',fn_clip(ref1','prc1'));
    set(im1,'cdata',fn_clip(avreg1','prc1'));
    set(im2,'cdata',fn_clip(avreg1'-ref','fit[0]',cm));
    drawnow
end

% Register Loop
curshift = [0 0]; % this value is actually never updated
for k=kstart:ntrial    
    fn_progress(['register trial ' num2str(k) '/' num2str(ntrial)])   
  
    % Get data
    datak = threed(T(k).(s.data));
    
    % Cut sub-frames?
    if doinitialcut
        datak = datak(oki,okj,:);
    end
    
    % Coregister
    switch regtype
        case 'frameglobal'
            regpar.shift0 = curshift;
            shiftk = fn_register(datak,regpar);
        case 'trial'
            regpar.shift0 = curshift;
            shiftk = fn_register(mean(datak,3),regpar);
            shiftk = repmat(shiftk,1,2);
        case {'framelocal' 'framelocalglobal'}
            if s.refrange, ref = mean(datak(:,:,s.refrange),3); else ref = mean(datak,3); end            
            regpar.ref = ref;
            regpar.shift0 = [0 0];
            shiftk = fn_register(datak,regpar);
            %             if s.dorepeat
            %                 % (display)
            %                 SHIFT(:,1:ntall(k),k) = shiftk;
            %                 fn_progress 'resample'
            %                 avreg = nmean(fn_translate(datak,-shiftk),3);
            %                 nans = isnan(avreg);
            %                 avreg(nans) = ref(nans);
            %                 AVREG{k} = avreg;
            %                 trialdisplay(k)
            %                 % (second estimation with updated ref)
            %                 ref = avreg; regpar.ref = avreg;
            %                 shiftk = fn_register(datak,regpar);
            %             end
    end
    SHIFT(:,1:ntall(k),k) = shiftk;
    fn_progress 'resample'
    avreg = nmean(fn_translate(datak,-shiftk),3);
    nans = isnan(avreg);
    avreg(nans) = ref(nans);
    AVREG{k} = avreg;
    
    % (display)
    trialdisplay(k)
    
    % save
    if dosave, fn_savevar(driftfile,SHIFT,AVREG), end
end

% Register trials together for 'framelocalglobal' registration type
if strcmp(regtype,'framelocalglobal')
    driftgfile = fullfile(savedir,'estimateddriftsglobal.mat');
    if useprevious && exist(driftgfile,'file')
        shifttrials = fn_loadvar(driftgfile);
    else
        % Register a first time to the middle trial average
        AVREG = cat(3,AVREG{:});
        ref = AVREG(:,:,ceil(ntrial/2));
        regpar.ref = ref;
        shifttrials = fn_register(AVREG,regpar);
        % Display
        fn_progress resample
        avreg = nmean(fn_translate(AVREG,-shifttrials),3);
        trialdisplay(0)
        % Register a second time to the average over all trials
        nans = isnan(avreg);
        avreg(nans) = ref(nans);
        ref = avreg; regpar.ref = avreg;
        shifttrials = fn_register(AVREG,regpar);
        % Display
        fn_progress resample
        avreg = nmean(fn_translate(AVREG,-shifttrials),3);
        trialdisplay(0)
        % Save
        fn_savevar(driftgfile,shifttrials)
    end
    % Combine drifts within trials and between trials
    shifttrials = permute(shifttrials,[1 3 2]);
    SHIFT = fn_add(SHIFT,shifttrials);
end

% Cut to have all pixel values defined
if s.docut && ~strcmp(regtype,'framelocal')
    negshift = -SHIFT(:,:); % corrections that will be applied
    negshift = negshift(:,~isnan(negshift(1,:)));
    mi = 1+max(0,floor(max(negshift(1,:)))+1); % e.g. if the right-most x correction is +0.5, pixels are defined from 2 only
    mj = 1+max(0,floor(max(negshift(2,:)))+1);
    Mi = nx1+min(0,ceil(min(negshift(1,:)))-1); % e.g. if the left-most x correction is -0.5, pixels are defined only until nx-1
    Mj = ny1+min(0,ceil(min(negshift(2,:)))-1);
    if doinitialcut
        mi = mi+(find(oki,1,'first')-1);
        mj = mj+(find(okj,1,'first')-1);
        Mi = Mi-(nx1-find(oki,1,'last'));
        Mj = Mj-(ny1-find(okj,1,'last'));
    end
end

% Resample loop
Treg = T;
if s.dosavereg
    for k=1:ntrial
        fn_progress(['trial ' num2str(k) '/' num2str(ntrial) ' - resample'])
        
        % Name for saving file
        fk = T(k).file;
        if ischar(fk) && ~ismatrix(fk)
            fk = fn_fileparts(T(k).file,'base');
        else
            fk = '';
        end
        if donum
            num = num2str(k,'trial%.3i');
            if isempty(fk), fk = num; else fk = [fk '_' num]; end %#ok<AGROW>
        end
        fk = [savedir '/' prefix fk suffix]; %#ok<AGROW>
        
        % Shift
        shiftk = reshape(SHIFT(:,1:ntall(k),k),[2 T(k).nfr T(k).nc]);
        
        
        if useprevious && exist(fk,'file')
            % Form new trial
            Treg(k) = tps_trial(fk,T(k));
            Treg(k).user.registration = shiftk;
        else
            % Get data
            datak = T(k).(s.data);
            
            % Resample (drift estimation was performed using bicubic interpolation
            % because of differentiability; the final interpolation however is
            % linear to limit the 'smoothing' effect of resampling)            
            datak = fn_translate(datak,-shiftk,'linear',s.datatype);
            
            % Cut
            if s.docut
                if strcmp(regtype,'framelocal')
                    negshift = -shiftk(:,:); % corrections applied
%                     negshift = negshift(:,~isnan(negshift(1,:)));
                    mi = 1+max(0,floor(max(negshift(1,:)))+1); % e.g. if the right-most x correction is +0.5, pixels are defined from 2 only
                    mj = 1+max(0,floor(max(negshift(2,:)))+1);
                    [nx, ny, ~] = size(datak);
                    Mi = nx+min(0,ceil(min(negshift(1,:)))-1); % e.g. if the left-most x correction is -0.5, pixels are defined only until nx-1
                    Mj = ny+min(0,ceil(min(negshift(2,:)))-1);
                end
                datak = datak(mi:Mi,mj:Mj,:,:);
            end
            
            % Form new trial
            Treg(k) = tps_trial(datak,T(k));
            Treg(k).user.registration = shiftk;
            
            % Save
            if dosave
                fn_progress(['trial ' num2str(k) '/' num2str(ntrial) ' - save'])
                fk = fn_fileparts(T(k).file,'base');
                if donum
                    num = num2str(k,'trial%.3i');
                    if isempty(fk), fk = num; else fk = [fk '_' num]; end %#ok<AGROW>
                end
                fk = [savedir '/' prefix fk suffix]; %#ok<AGROW>
                savedata(Treg(k),fk)
            end
        end
    end
else
    if s.docut
        cropping = {'crop',{mi:Mi,mj:Mj}};
    else
        cropping = {};
    end
    for k=1:ntrial
        shiftk = reshape(SHIFT(:,1:ntall(k),k),[2 T(k).nfr T(k).nc]);
        Treg(k) = applyprocessing(T(k),'translate',-shiftk,cropping{:});
    end
end

end

% % Old version could handle differences in spatial resolution between
% trials!!

% % Some constants
% cx = (1+nx)/2; % center of image (x)
% cy = (1+ny)/2; % center of image (y)

%     data0k = mean(datak(:,:,1:10),3);
%     shiftk = fn_register(datak,struct('I0',data0k));
%
%     % Coregister frame0 to reference frame
%     dx = T(k).dx;
%     dy = T(k).dy;
%     % (adapt the reference to match the zoom of current trial)
%     if dx==dx0 && dy==dy0
%         I0k = I0;
%     else
%         % indices in datak to which correspond points in I0
%         inow = cx + ((1:nx)-cx)*(dx0/dx);
%         jnow = cy + ((1:ny)-cy)*(dy0/dy);
%         inew = 1:nx;
%         jnew = 1:ny;
%         iok = (inew>=inow(1) & inew<=inow(end));
%         jok = (jnew>=jnow(1) & jnew<=jnow(end));
%         [inow jnow] = ndgrid(inow,jnow);
%         [inew jnew] = ndgrid(inew(iok),jnew(jok));
%         I0k = interpn(inow,jnow,I0,inew,jnew,'cubic');
%         if ~all(iok) || ~all(jok), data0k = data0k(iok,jok); end
%     end
%     % (register)
%     if k<kstart && docheck
%         shift0k = fn_alignimage(data0k,I0k,hf);
%         [shift0k e] = fn_register(data0k,struct('ref',I0k,'shift0',shift0k));
%         fprintf('error: %.2f\n',e);
%     else
%         [shift0k e] = fn_register(data0k,struct('ref',I0k));
%         fprintf('error: %.2f\n',e);
%     end

%     % Resample
%     % (first, transformation to bring I0k onto I0 pixel)
%     if dx==dx0 && dy==dy0
%         M1 = eye(3);
%     else
%         a = diag([dx/dx0 dy/dy0]);                % one original pixel has size dx; in the new resampling at dx0, it will cover dx/dx0 pixel(s)
%         b = [cx - dx/dx0*cx, cy - dy/dy0*cy];     % b satisfies a*c + b = c
%         M1 = [[a; b] [0; 0; 1]];
%     end
%     % (second, transformation to bring data0k onto I0k)
%     M2 = [[eye(2); -shift0k'] [0; 0; 1]];
%     % (third, transformation to bring each frame of data on data0k)
%     for i=1:nt
%         M3 = [[eye(2); -shiftk(:,i)'] [0; 0; 1]];
%         M = M1*M2*M3;
%         M = [[M([2 1],[2 1]); M(3,[2 1])] [0; 0; 1]]; % invert x and y to fit Matlab conventions
%         % (do the resampling now)
%         tf = maketform('affine',M);
%         datak(:,:,i) = imtransform(datak(:,:,i),tf,'xdata',[1 ny],'ydata',[1 nx]); % bilinear interpolation: ok, but creates problem on the sides
%         if ~isempty(sfrk)
%             sfrk(:,:,i) = imtransform(sfrk(:,:,i),tf,'xdata',[1 ny],'ydata',[1 nx]);
%         end
%     end

%     % Form the new trial and save
%     disp('saving')
%     datak = reshape(datak,[nx*ny nt]);
%     datak(any(isnan(datak),2),:) = NaN; % this will be a bit cleaner
%     datak = reshape(datak,[nx ny nt]);
%     Treg(k) = tps_trial({datak,sfrk},T(k));
%     set(Treg(k),'dx',dx0,'dy',dy0)
%     [p b] = fileparts(T(k).file);
%     savedata(Treg(k),[savedir b '.mat']);
%     erasedata(Treg(k)) % save memory
%     if ~isempty(fname), save(fname,'Treg','-mat'), end


%---
function m = nmean(x,dim)
    nans = isnan(x);
    x(nans) = 0;
    n = sum(~nans,dim);
    n(n==0) = NaN; % prevent divideByZero warnings
    % Sum up non-NaNs, and divide by the number of non-NaNs.
    m = sum(x,dim) ./ n;
end

%---
function ref = iterativeReference(T,s,regpar)
% Align trial on a first reference, then compute a new, better-quality,
% reference, align trials again on this new reference, repeat until user
% satisfaction.

% Data
data = cat(4,T(s.reftrial).data);
data = fn_float(data(:,:,:)); % convert to float to have NaNs
p = fn_pointer;
p.data = data;
p.xdata = data;
nt = size(data,3);
p.shift = zeros(2,nt);
v = fn_imvect(p.data);
v = fn_normalize(v,1,'zscore');
v0 = fn_normalize(mean(v,2),1,'zscore');
p.rms = 1-column(mean(fn_mult(v,v0),1));

% Init display
hf = fn_figure('Iterative Reference');
figure(hf), clf
a(1) = fourd(p.data,'2d','in',subplot(121),'shapemode','ellipse');
a(2) = fourd(p.rms,'plot','in',subplot(122),'mat',{[] [] 3});
uicontrol('string','REINIT','pos',[5 75 100 30],'callback',@(u,e)reinit_ref(p,s,regpar,a))
uicontrol('string','ALIGN->NEW REF','pos',[5 40 100 30],'callback',@(u,e)iterate_ref(p,s,regpar,a))
% p.txt = uicontrol('style','text','pos',[5 40 100 30]);
ok=uicontrol('string','DONE','pos',[5 5 100 30],'callback',@(u,e)delete(u));

waitfor(ok)
if ishandle(hf)
    ref = a(1).SI.slice.data;
    close(hf)
else
    ref=[]; 
end


end

%---
function reinit_ref(p,s,regpar,a)
    p.xdata = p.data;
    a(1).SI.data = p.xdata;
    nt = size(p.data,3);
    p.shift = zeros(2,nt);
    
    v = fn_imvect(p.data);
    v = fn_normalize(v,1,'zscore');
    v0 = fn_normalize(mean(v,2),1,'zscore');
    p.rms = 1-column(mean(fn_mult(v,v0),1));
    a(2).SI.data = p.rms;
end

%---
function iterate_ref(p,s,regpar,a)

% sel = a(2).SI.selection.getselset(1).singleset;
% if isempty(sel)
%     idx = a(2).SI.ij;
% else
%     idx = sel.indices;
% end
% regpar.ref = mean(p.xdata(:,:,idx),3);
regpar.ref = a(1).SI.slice.data;
[p.shift e p.xdata] = fn_register(p.data,regpar); %#ok<ASGLU>
a(1).SI.data = p.xdata;

v = fn_imvect(p.data);
okpix = ~isnan(regpar.ref);
v = fn_normalize(v(okpix,:),1,'zscore');
v0 = fn_normalize(regpar.ref(okpix),1,'zscore');
p.rms = 1-column(mean(fn_mult(v,v0),1));
a(2).SI.data = p.rms;

end


