function label = retinotopymapping(V) 
% Perform analysis on current trial

persistent smem

% Label
% function label = f('label') 
if nargin==0, V=evalin('base','V'); end
if ischar(V) && strcmp(V,'label')
    label = 'Retinotopy mapping';
    return
end

%%
% Checks
T = V.trial;
nc = T.nc;
if ~ismember(nc,[1 2]), errordlg 'Trial needs to have 2 conditions.', return, end
if isempty(T.dt_sec), errordlg 'Frame rate is not defined.', return, end

% Parameters
s = struct( ...
    'stimdelay',    {5          'double'    'Delay of first stimulation'}, ...
    'T',            {10         'double'    'Stimulation period (s)'}, ...
    'trestrict',    {[15 105]   'double'    'Restrict projection to time segment'}, ...
    'nharm',        {5          'stepper 1 1 Inf'   'Numer of harmonics'}, ...
    'hemodelay',    {'auto'     'char'      'Hemodynamic delay (s)'});
if ~isempty(smem)
    s(1) = smem;
end
s = fn_structedit(s);
smem = s;
if isempty(s), return, end
if strcmp(s.hemodelay,'auto')
    s.hemodelay = [];
else
    s.hemodelay = str2double(s.hemodelay);
end

% Further checks
if mod(diff(s.trestrict),s.T)
    errordlg 'Length of restricting time segment must be a multiple of period T.'
    return
end
ncycle = diff(s.trestrict)/s.T;

% Timer
c = fn_watch(V.hf); %#ok<NASGU>

% Some constants
dt = T.dt_sec;
nt = T.nfr;
tt = (0:nt-1)*dt;
frames = find(tt>=s.trestrict(1) & tt<=s.trestrict(2));
nframe = length(frames);

% Regressors
nbase = 2*s.nharm;
stimcycle = (2*pi/s.T)*(tt(frames)-s.stimdelay); % phase 0 at stimulation start times
if isempty(s.hemodelay)
    % hemodynamic delay estimated in each pixel
    X = zeros(nframe,nbase);
    for kbase=1:s.nharm
        X(:,2*kbase-1) = cos(kbase*stimcycle);
        X(:,2*kbase)   = sin(kbase*stimcycle);
    end
else
    X = zeros(nframe,2,nbase);
    for kbase=1:s.nharm
        X(:,1,2*kbase-1) = cos(kbase*stimcycle);
        X(:,1,2*kbase)   = sin(kbase*stimcycle);
        X(:,2,2*kbase-1) = cos(s.hemodelay-(kbase*stimcycle-s.hemodelay));
        X(:,2,2*kbase)   = sin(s.hemodelay-(kbase*stimcycle-s.hemodelay));
    end
    X = reshape(X,nframe*2,nbase);
end

%%
% Projection

switch nc
    case 1
        ktrial = V.ktrial;
        data1 = fn_reshapepermute(V.content.trials(ktrial).dataop,{3 4 [1 2]});
        data2 = fn_reshapepermute(V.content.trials(ktrial+1).dataop,{3 4 [1 2]});
        data = [data1 data2]; clear data1 data2
    case 2
        data = T.dataop;
        data = fn_reshapepermute(data,{3 4 [1 2]}); % nt * nc * (nx*ny)
end
datac = data(frames,:,:); avgc = mean(datac,1); 
data = fn_subtract(data,avgc); datac = fn_subtract(datac,avgc);
if isempty(s.hemodelay)
    beta = X\reshape(datac,[nframe 2*T.nx*T.ny]);
    beta = reshape(beta,[nbase 2 T.nx T.ny]);
else
    beta = X\reshape(datac,[nframe*2 T.nx*T.ny]);
    beta = reshape(beta,[nbase 1 T.nx T.ny]);
end
clear datac

% Display the result of the projection on selected ROIs signals
rois=V.content.signal.x(V.ktrial,:);
nroi = length(rois);
[x xproj1 xproj] = deal(NaN(nt,2,nroi));
for i=1:nroi
    idx = rois(i).sel.dataind;
    x(:,:,i) = mean(data(:,:,idx),3);
    if isempty(s.hemodelay)
        betai = X\x(frames,:,i);
        xproj(frames,:,i) = X*betai;  
        xproj1(frames,:,i) = X(:,1:2)*(X(:,1:2)\x(frames,:,i));
    else
        betai = X\column(x(frames,:,i));
        xproj(frames,:,i) = reshape(X*betai,[nframe 2]);  
        xproj1(frames,:,i) = reshape(X(:,1:2)*(X(:,1:2)\column(x(frames,:,i))),[nframe 2]);
    end

    x(:,:,i) = detrend(x(:,:,i)); % for better display
    
    phasei = atan2(-betai(2,:),-betai(1,:)); % phase of the negative peak, between -pi and pi, size nx*ny*nc
    phasei = (s.T/(2*pi))*phasei; % between -T/2 and T/2
    if ~isempty(s.hemodelay), phasei = s.hemodelay+[1 -1]*(phasei-s.hemodelay); end
    hemodelayi = mod( (phasei(1)+phasei(2))/2 , s.T/2 );
    % hemodelay = nmean(hemodelay(:));
    tuningi = mod( phasei(1)-hemodelayi , s.T );
    fprintf('region %i: phase = %.1f/%.1f, hemodelay = %.1f, tuning = %.1f\n',i,phasei,hemodelayi,tuningi)
end

hf = fn_figure('retinotopy mapping - regions','color','w');
figure(hf)
tstep = (ceil(nt*dt/s.T)+1)*s.T;
hl = fn_gridplot(tt,cat(4,x,xproj1,xproj),'rowcol-',[tstep nstd(x(:))*5]);
hl = reshape(hl,[2 nroi 3]);
set(hl(:,:,2),'color','k','linewidth',2)
xtick = s.trestrict(1) + row(fn_add((0:ncycle)'*s.T,(0:nroi-1)*tstep));
set(gca,'xtick',xtick,'xgrid','on')

%%
% Estimate hemodynamic delay (d) and tuning (x)
% we have
%   phase_cond1 = x + d [T]
%   phase_cond2 = (T-x) + d [T] = -x +d [T]
% hence
%   d = (phase_cond1+phase_cond2) [T/2]    (-> we take the value btw. 0 and T/2)
%   x = phase_cond1 - d [T]                (-> we take the value btw. 0 and T)

phase = squeeze(atan2(-beta(2,:,:,:),-beta(1,:,:,:))); % phase of the negative peak, between -pi and pi, size nx*ny*nc
phase = (s.T/(2*pi))*phase; % between -T/2 and T/2 % nc * nx * ny
if isempty(s.hemodelay)
    phase = permute(phase,[2 3 1]); % nx * ny * nc
else
    phase = cat(3,phase,s.hemodelay-(phase-s.hemodelay));
end
hemodelay = mod( (phase(:,:,1)+phase(:,:,2))/2 , s.T/2 );
% hemodelay = nmean(hemodelay(:));
tuning = mod( phase(:,:,1)-hemodelay , s.T );

sensitivity = squeeze(mean(sqrt(sum(beta(1:2,:,:,:).^2,1)),2));
sensitivity = min(1,sensitivity/max(sensitivity(:))*2);

hf = fn_figure('retinotopy mapping','color','w');
figure(hf)
colormap(gray(256))

dx = 3.14; % pixel size in microns for configuration low - 1280x960
xx = [0 T.nx-1]*(dx*16); % binning 16x16
yy = [0 T.ny-1]*(dx*16);

% subplot(221)
% imagesc(xx,yy,phase(:,:,1)',[-1 1]*s.T/2), axis image
% subplot(222)
% imagesc(xx,yy,phase(:,:,2)',[-1 1]*s.T/2), axis image
G = V.a4d.G;
subplot(221)
fourd(hemodelay,G,'clip',[0 s.T/2])
% imagesc(xx,yy,hemodelay',[0 s.T/2]), axis image
% set(gca,'xtick',0:1e3:3e3,'ytick',0:1e3:3e3,'xticklabel',[],'yticklabel',[]), grid on
title 'hemodynamic delay'
subplot(222)
tuningcol = fn_clip(tuning,[0 s.T],maporient(256),[1 1 1]);
tune_sens = fn_add( fn_mult(sensitivity,tuningcol), fn_mult(1-sensitivity,third([1 1 1]*.6)) );
fourd(permute(tune_sens,[1 2 4 5 6 7 3]),G,'2d','dimsplus',7)
% imagesc(xx,yy,permute(tune_sens,[2 1 3])), axis image, grid on
% set(gca,'xtick',0:1e3:3e3,'ytick',0:1e3:3e3,'xticklabel',[],'yticklabel',[]), grid on
title 'tuning'


%% green image overlay

if ~isfield(V.content.user,'roiimage')
    froi = fn_getfile('*.png;*.bmp','Select ROI image');
    if isequal(froi,0)
        V.content.user.roiimage = T.data0;
    else
        V.content.user.roiimage = fn_readimg(froi);
    end
    V.file_markchange()
end
roiimage = V.content.user.roiimage;

try %#ok<TRYNC>
    try
        acq = V.content.trials(end-1).user.acquisition;
    catch
        acq = V.content.trials(end-1).user.acquistion;
    end
    ROIpos = acq.ROIPosition;
    roiimage = roiimage(ROIpos(1)+(1:ROIpos(3)),ROIpos(2)+(1:ROIpos(4)));
end
if all(size(roiimage)>=2*[T.nx T.ny]), roiimage = fn_bin(roiimage,2); end
roiimage = fn_clip(roiimage);

dx = 3.14; % pixel size in microns for configuration low - 1280x960
xx = [0 size(roiimage,1)-1]*(dx*2); % binning 2x2
yy = [0 size(roiimage,2)-1]*(dx*2);

% transformation between ROI image and (binned) raw data
nbin = round(size(roiimage,1)/V.nx);
if nbin==1
    G1 = G;
else
    G1 = rotation(V.a4d.F,'mat',[V.dx/nbin V.dy/nbin]);
end

figure(hf)
subplot(223)
fourd(roiimage,G1)
% imagesc(xx,yy,roiimage'), axis image
% set(gca,'xtick',0:1e3:3e3,'ytick',0:1e3:3e3,'xticklabel',[],'yticklabel',[]), grid on
title 'calcium'

% tuninglarge = fn_enlarge(tuningcol,size(green));
% senslarge = fn_enlarge(sensitivity,size(green));
% im = fn_add( fn_mult(senslarge/2,tuninglarge), fn_mult(1-senslarge/2,green) );
tunesenselarge = fn_enlarge(tune_sens,size(roiimage));
im = fn_add(tunesenselarge/2,roiimage/2);

figure(hf)
subplot(224)
fourd(permute(im,[1 2 4 5 6 7 3]),G1,'2d','dimsplus',7)
% imagesc(xx,yy,permute(im,[2 1 3])), axis image
% set(gca,'xtick',0:1e3:3e3,'ytick',0:1e3:3e3,'xticklabel',[],'yticklabel',[]), grid on
title 'overlay'


%% Save results in base workspace

assignin('base','roiimage',roiimage)
assignin('base','tune_sens',tune_sens)
assignin('base','overlay',im)
disp 'results saved in variables roiimage, tune_sens and overlay'

