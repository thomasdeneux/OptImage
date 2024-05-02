function [channels delay] = tps_channelnames(nchannel,delay)
% function [channels delay] = tps_channelnames(nchannel[,delay])


s = struct('channel1',cell(1,2));
for i=1:nchannel
    [s.(['channel' num2str(i)])] = deal(['name of channel' num2str(i)],'xchar');
end
spec = s(2);
s = s(1);
if nargout==2
    if nargin<2, s.delay = 0; else s.delay = delay; end
    spec(1).delay = 'double';
    spec(2).delay = 'delay (in seconds)';
end

ok = false;
while ~ok
    s = fn_structedit(s,spec);
    if isempty(s)
        channels = {};
        delay = [];
        return
    end
    channels = {};
    ok = true;
    for i=1:nchannel
        namei = s.(['channel' num2str(i)]);
        if ~regexp(namei,'[a-zA-Z]\w*')
            errordlg 'not a valid variable name'
            ok = false;
            break
        end
        if ~isempty(namei), channels = [channels; {i namei}]; end %#ok<AGROW>
    end
end

if nargout==2
    delay = s.delay;
end