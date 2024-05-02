function out = oix_registerimage(V) 
% function oix_registerimage(V|C) 
% function label = oix_registerimage('label') 

if ischar(V) && strcmp(V,'label')
    out = 'Registered Image';
    return
end

% Get the tpv_content object
C = V.content;

% Existing registered images
if isfield(C.user,'registeredimages')
    images = C.user.registeredimages;
else
    images = struct('name',{},'image',{},'mat',{});
end

% Prompt for existing image
if ~isempty(images)
    answer = fn_input('choose image','new...',['new...' {images.name}]);
    if isempty(answer), return, end
    donew = strcmp(answer,'new...');
else
    donew = true;
end

% Controls specification
typeflag = 'stepper 1';
s = struct( ...
    'name',     {''         'char'}, ...
    'scale',    {1          'double'}, ...
    'rotation', {0          [typeflag ' -20 20 .5 %g']}, ...
    'xtrans',   {0          [typeflag ' -50 50 .1 %g']}, ...
    'ytrans',   {0          [typeflag ' -50 50 .1 %g']} ...
    );
spec = s(2); s = s(1);

% Read new image or load from content
if donew
    % New image
    fname = fn_getfile('','Select image');
    if ~fname, return, end
    a = fn_readimg(fname);
    
    % Image structure
    im = struct('name',fn_fileparts(fname,'base'),'image',a,'mat',eye(3));
    kim = length(images)+1;
    images(kim) = im;
    C.user.registeredimages = images;
else
    kim = strcmp(answer,{images.name});
    im = images(kim);
    [nx ny ncol] = size(im.image); %#ok<NASGU>
    s.rotation = atan2(im.mat(2,3),im.mat(2,2))*(180/pi);
    rot =  [cosd(s.rotation) sind(s.rotation); -sind(s.rotation) cosd(s.rotation)];
    s.scale = im.mat(2,2) / rot(2,2);
    trans = im.mat(2:3,1) - (-rot+eye(2))*s.scale*[nx/2; ny/2];
    s.xtrans = trans(1)/s.scale;
    s.ytrans = trans(2)/s.scale;
end
s.name = im.name;


% Display image
R = showimage(V,im);

% Controls
str1 = num2str(floor(size(im.image,1)));
str2 = num2str(floor(size(im.image,2)));
spec.xtrans = [typeflag ' -' str1 ' ' str1 ' .1 %g'];
spec.ytrans = [typeflag ' -' str2 ' ' str2 ' .1 %g'];
fn_control(s,spec,@(x)controlimage(V,R,im,kim,x))
    

%---
function R = showimage(V,im)

R = rotation(V.a4d.F,'mat',im.mat);
SI = projection(R,[1 2],'data',im.image);
if size(im.image)==3, SI.dimsplus = 3; end
activedisplayImage(SI,'in',figure);

if nargout==0, clear R, end


%---
function controlimage(V,R,im,kim,x)

% name
im.name = x.name;

% scaling
scale = eye(2)*x.scale;

% rotation (to make it centered on the image center, add a specific
% translation)
rot = [cosd(x.rotation) sind(x.rotation); -sind(x.rotation) cosd(x.rotation)];
[nx ny ncol] = size(im.image); %#ok<NASGU>
trans = (-rot+eye(2))*scale*[nx/2; ny/2];

% translation
trans = trans + [x.xtrans; x.ytrans]*x.scale;

% final transformation
im.mat = [1 0 0; trans rot*scale];
R.mat = im.mat;

% save
V.content.user.registeredimages(kim) = im;









    
    