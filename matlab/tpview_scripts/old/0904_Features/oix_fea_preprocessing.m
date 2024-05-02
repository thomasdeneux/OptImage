function label = oix_fea_preprocessing(V) 
% function oix_fea_preprocessing(V) 
% function label = oix_fea_preprocessing('label') 

if ischar(V) && strcmp(V,'label')
    label = 'Correct data';
    return
end

fea_preprocessing(V)
