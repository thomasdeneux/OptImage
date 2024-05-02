function label = map_orientation_colors(V) 
% function f(V) 
% function label = f('label') 

if nargin==0, V=evalin('base','V'); end
if ischar(V) && strcmp(V,'label')
    label = 'Orientation maps -> colors';
    return
end

%% Make sure that the 'Orientation maps' script was executed

if ~(strcmp(V.disppar.im1cond,'C1/C3') && strcmp(V.disppar.im2cond,'C2/C4'))
    errordlg('Execute script ''Orientation maps'' before executing ''Orientation maps -> color''.')
end

%% Add a small spatial high-pass to get nice maps centered on zero

T.opdef = tps_dataopdef({{'space','divide by','frames','1'} ...
    {'filter' '' '100' '' '' false true false false}});

% Update display
V.display_changeview('ktrial')

%% Get images currently displayed in first and second image displays

% images
map1 = fn_clip(V.a4d.SI1.slice.data, V.a4d.im1.clip, [-1 1]); % values between -1 and 1
map2 = fn_clip(V.a4d.SI2.slice.data, V.a4d.im1.clip, [-1 1]); % values between -1 and 1

%% Display it with a color code for orientation preference

angle = atan2(map2,map1); % value between -pi and pi; maybe the correct formula is atan2(map2,map1) instead
saturation = min(sqrt(map1.^2 + map2.^2) / M, 1); % value between 0 and 1

color = fn_clip(angle, [-pi pi], maporient(256));
color = fn_add(.5 * (1-saturation), fn_mult(color, saturation));

figure(1), subplot(211)
imagesc(permute(color,[2 1 3]))
axis image

% display also with blood vessels superimposed
vessels = fn_clip(V.trial.data0(:,:,1),'prc1-1'); % values between 0 and 1
mix = 0.1;
color_vessels = fn_add(color * mix, vessels * (1-mix));
subplot(212)
imagesc(permute(color_vessels,[2 1 3]))
axis image

%% Show the color map

fn_showcolormap(maporient(256), [-180 180])
ylabel({'orientation preference (°)' 'please check that the values are correct'})
