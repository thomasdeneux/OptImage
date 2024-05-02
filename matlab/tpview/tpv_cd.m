function f = tpv_cd(flag)
%---
% this function can be re-defined in order to get different file names for
% tpv_cd('data') and tpv_cd('calibration')

if fn_dodebug
    disp 'avoid using tpv_cd!!!'
    f = fullfile(fn_cd('2ps'),flag);
else
    f = '';
end

if nargout==0
    cd(f)
    clear f
end