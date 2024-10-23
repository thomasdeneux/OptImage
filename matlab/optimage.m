function optimage(varargin)

% add appropriate folders to the Matlab path
basedir = fileparts(which('optimage'));
folders = {'' 'brick' 'colormaps' 'explor' ...
    'matlab exchange' 'matlab exchange/mask2poly' 'tpbasic' 'tpview'};
for i=1:length(folders), folders{i}=fullfile(basedir,folders{i}); end
addpath(folders{:})

% attempt also to add the 'io' folder
res = which('oi_loadBLK');
if isempty(res)
    basematlab = fileparts(fileparts(basedir)); % 2 levels up
    iodir = fullfile(basematlab, 'io');
    if exist(iodir, 'dir')
        addpath(iodir)
    else
        error 'Please add folder ''io'' to path to run OptImage'
    end
end

% run OptImage
tpview('OptImage',varargin{:});


