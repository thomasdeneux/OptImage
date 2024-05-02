function label = map_orientation(V) 
% function f(V) 
% function label = f('label') 

if nargin==0, V=evalin('base','V'); end
if ischar(V) && strcmp(V,'label')
    label = 'Orientation maps';
    return
end

% Write script code below

%% Everything can be performed in the program manually rather than by the code below!

%% Select block files

f = fn_getfile('*.BLK','Please select block files for orientation map.');

%% Load and average trials

% below is a slightly complicated syntax for reading blocks and averaging
% them on the fly
T = tps_trial.readheader(f,'avg'); % a "tps_trial" object

% Set tpview to display this tps_trial object
V.setdata(T)

%% Improve the display

% Divide by first frame to remove noise present before stimulation
T.opdef = tps_dataopdef({'space','divide by','frames','1'});

% Update display
V.display_changeview('ktrial')

% In first image display look at condition operation C1/C3 (I don't know what are
% conditions C1, C2, C3, C4, but if they represent 0°, 45°, 90°, 135°, then
% C1/C3 represents 0°/90°)
V.setpar('im1cond','C1/C3')

% And in second image display look at C2/C4
V.setpar('im2','dataop')
V.setpar('im2cond','C2/C4')



