function label = map_dominance(V) 
% function f(V) 
% function label = f('label') 

if nargin==0, V=evalin('base','V'); end
if ischar(V) && strcmp(V,'label')
    label = 'Ocular dominance';
    return
end

% Write script code below

%% Select block files

fleft = fn_getfile('*.BLK','Please select block files for LEFT eye.');
fright = fn_getfile('*.BLK','Please select block files for RIGHT eye.');

%% Load and average left and right eye data

% below is a slightly complicated syntax for reading blocks and averaging
% them on the fly
Tleft = tps_trial.readheader(fleft,'avg'); % a "tps_trial" object
Tright = tps_trial.readheader(fright,'avg');

dataleft = Tleft.data; % its data
dataright = Tright.data;

%% Compute the division right / left and load in tpview

data = dataright ./ dataleft; % the new data
T = tps_trial(data, Tleft); % build a new tps_trial object with this data by using the header info inside Tleft 
V.setdata(T)

%% All the steps below can be performed in the program rather than by command lines

% Divide by first frame to remove noise present before stimulation
T.opdef = tps_dataopdef({'space','divide by','frames','1'});

% Update display
V.display_changeview('ktrial')

% Look at image display of condition C1+C2+C3+C4 (i.e. average the 4
% directions together)
V.setpar('im1cond','C1+C2+C3+C4')

