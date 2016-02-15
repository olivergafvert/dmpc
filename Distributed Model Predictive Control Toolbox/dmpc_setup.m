mpath = mfilename('fullpath');
disp(['Found path: ' mpath]);
if ispc,
    fs = '\'; 
else
    fs = '/'; 
end
temp = strfind( mpath, fs );
mpath = mpath( 1 : temp(end) - 1 );
disp(['Setting project root: ' mpath]);
addpath(genpath(mpath))
savepath

