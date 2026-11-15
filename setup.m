function setup()
%SETUP Add all project folders to the MATLAB path.
%   Run once per MATLAB session from the repository root:
%       >> setup
root = fileparts(mfilename('fullpath'));
addpath(fullfile(root, 'model'), fullfile(root, 'data'), ...
        fullfile(root, 'derivation'), fullfile(root, 'identification'), ...
        fullfile(root, 'validation'), fullfile(root, 'simulation'), ...
        fullfile(root, 'lab_scripts'));
fprintf('Rotary double pendulum: paths added (%s)\n', root);
end
