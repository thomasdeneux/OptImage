function label = av_popcorrelation(V) 
% function f(V) 
% function label = f('label') 

if nargin==0, V = evalin('base','V'); end

if ischar(V) && strcmp(V,'label')
    label = 'Population analysis';
    return
end

trials = find([V.content.trials.status]=='n');
% trials = 1:50;
ntrial = length(trials);
T1 = V.content.trials(trials(1));
T = V.content.trials(trials);

% get signals
if strcmp(T1.type,'linescan')
    nt = T1.ny*T1.nfr;
    sels = 1:T1.nx;
    nsel = length(sels);
    data = squeeze(cat(4,V.content.trials(trials).data));
    data = permute(data,[2 1 3]); % nt x nsel x ntrial
    data = data(:,sels,:);
else
    sels = 1:V.content.nsel;
    % sels = [1 9 17 18 26];
    nsel = length(sels);
    nt = T1.nfr;
    fn_progress('compute signals',ntrial)
    for i=1:ntrial, fn_progress(i), getslice(V.content,trials(i),sels); end
    x = getslice(V.content,trials)';
    data = zeros(nt,nsel,ntrial);
    for i=1:ntrial*nsel
        data(:,i) = x(i).data;
    end
end


dt = T1.dt_sec;
tidx = (0:nt-1)*dt;

% population vectors
stim = T1.stim;
normtime = [0.5 1]; %[0 stim(1)];
avgtime = [1.5 3]; %stim(1)+[1 2];
normtype = 'trial'; % 'trial', 'global' or 'none'
inorm = (tidx>=normtime(1) & tidx<=normtime(2));
iavg = (tidx>=avgtime(1) & tidx<=avgtime(2));
switch normtype
    case 'trial'
        x = squeeze(mean(data(iavg,:,:),1)./mean(data(inorm,:,:),1));
    case 'global'
        x = squeeze(fn_div(mean(data(iavg,:,:),1),mean(mean(data,1),3)));
    case 'none'
        x = squeeze(mean(data(iavg,:,:),1));
end

% reorder according to condition
conds = [T.stimid];
ucond = unique(conds);
ncond = length(ucond);
npercond = zeros(1,ncond);
x1 = zeros(nsel,ntrial);
trialspercond = cell(1,ncond);
for i=1:ncond
    trialsi = find(conds==ucond(i));
    trialspercond{i} = trialsi;
    npercond(i) = length(trialsi);
    x1(:,sum(npercond(1:i-1))+(1:npercond(i))) = x(:,trialsi);
end

% prepare figures
figure(1), clf, set(1,'color','w')
figure(2), clf, set(2,'color','w')

% show population vectors
figure(2), subplot(121)
imagesc(x,fn_clip(x,'prc0-98','getrange'))
set(gca,'fontsize',12)
xlabel('trial','fontsize',13)
ylabel('neuron','fontsize',13)

figure(1)
imagesc(x1,fn_clip(x1,'prc0-98','getrange'))
for i=1:(ncond-1)
    z = sum(npercond(1:i))+.5;
    line([z z],.5+[0 nsel],'color','k','linewidth',1)
end
ticks = cumsum(npercond)-npercond/2;
conddetails = T1.stimtable.getdetails(ucond);
condnames = {conddetails.name};
set(gca,'xtick',ticks,'xticklabel',condnames, ...
    'fontsize',12)
xlabel('condition','fontsize',13)
ylabel('neuron','fontsize',13)



% cross-correlations
method = 'correlation';
% method = 'euclidean';
C = 1-squareform(pdist(x',method));
C1 = 1-squareform(pdist(x1',method));

clip = [0 max(1-pdist(x',method))];

figure(2), subplot(122)
imagesc(C,clip)
axis image
set(gca,'fontsize',12)
xlabel('trial','fontsize',13)
ylabel('trial','fontsize',13)

figure(1), subplot(122)
figure(3), set(3,'color','w')
imagesc(C1,clip)
axis image
for i=1:(ncond-1)
    z = sum(npercond(1:i))+.5;
    line(.5+[0 ntrial],[z z],'color','w')
    line([z z],.5+[0 ntrial],'color','w')
end
ticks = cumsum(npercond)-npercond/2;
conddetails = T1.stimtable.getdetails(ucond);
condnames = {conddetails.name};
set(gca,'xtick',ticks,'ytick',ticks,'xticklabel',condnames,'yticklabel',condnames, ...
    'fontsize',12)
xlabel('condition','fontsize',13)
rotateXLabels(gca,90)
ylabel('condition','fontsize',13)

% average inter-trial correlations for every conditions pair
C2 = zeros(ncond,ncond);
for i=1:ncond
    for j=1:ncond
        C2(i,j) = mean(mean(C(trialspercond{i},trialspercond{j})));
    end
end
figure(4)
imagesc(flipud(C2))
axis image
set(gca,'xtick',1:ncond,'xticklabel',condnames, ...
    'ydir','normal','ytick',1:ncond,'yticklabel',condnames(end:-1:1), ...
    'fontsize',12)
xlabel('condition','fontsize',13)
rotateXLabels(gca,90)
ylabel('condition','fontsize',13)
colorbar



