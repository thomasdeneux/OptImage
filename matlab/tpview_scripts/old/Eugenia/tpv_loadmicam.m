function label = tpv_loadmicam(V) 
% function f(V) 
% function label = f('label') 

if ischar(V) && strcmp(V,'label')
    label = 'Manage MiCam data';
    return
end


T = V.content.trials;
T.mergestimtables
for i=1:V.ntrial
    fk = fn_fileparts(T(i).file,'base');
    tokens = regexp(fk,'^C\d*_(\d)*_.*$','tokens');
    cond = tokens{1}{1};
    condnum = str2double(cond);
    setstim(T(i),condnum)
end