function optimage(varargin)

% add appropriate folders to the Matlab path
basedir = fileparts(which('optimage'));
folders = {'' 'brick' 'colormaps' 'explor' ...
    'matlab exchange' 'matlab exchange/mask2poly' 'tpbasic' 'tpview'};
for i=1:length(folders), folders{i}=fullfile(basedir,folders{i}); end
addpath(folders{:})

% run OptImage
tpview('OptImage',varargin{:});


