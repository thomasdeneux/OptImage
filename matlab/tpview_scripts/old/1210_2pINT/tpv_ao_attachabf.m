function label = tpv_ao_attachabf(V) 
% function f(V) 
% function label = f('label') 

if ischar(V) && strcmp(V,'label')
    label = 'Attach ABF recordings';
    return
end

% look for files
recdir = fn_cd('spikes','data/recordings');

date = V.content.trials(1).fullinfo.MES(1).MeasurementDate(1:10);
date = strrep(date,'.','_');

d = dir([recdir '/' date '*.abf']);
files = fn_map(@(f)fullfile(recdir,f),{d.name});

% attach recordings
data_attachanalog(V,files)

% fix recording permutations
tpv_ao_fixabfmix(V)
