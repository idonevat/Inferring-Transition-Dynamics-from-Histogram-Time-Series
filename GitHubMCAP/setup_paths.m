function setup_paths()
%SETUP_PATHS Add the shared simulation functions to the MATLAB path.
root = fileparts(mfilename('fullpath'));
addpath(fullfile(root,'common'));
end
