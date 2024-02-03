%% program to study the protein turnover kinetics; version 2.0.0  
%( last upadted on 02-03-2024 by Dr. Abhijit Dasgupta)


close all
clearvars 
warning('off','all')
%% 
% Determine where your m-file's folder is.
folder = fileparts(which(mfilename)); 
% Add that folder plus all subfolders to the path.
addpath(genpath(folder));

%read parameter file
fid     = fopen('JUMPt.params');
param   = textscan(fid,'%s','delimiter','\n');
cellfun(@eval,param{1});
fclose(fid);

params.setting      = setting;
params.input_file   = input_file;
params.bin_size     = bin_size;
params.opti_algo    = optimization_algorithm;
params.purity       = purity_of_SILAC_food;
params.timepoints   = number_of_timepoints;
params.apparent_T50 = apparent_T50_calculation;

fprintf('\n\n  Welcome to JUMPt version 2.0.0 program \n\n\n ')

%Read the input data from input file and distribute the protein different bins 
data = binning(params);
fprintf('\n\n  Completed reading input file \n\n\n  Now fitting the protein data and calculate half-lives\n')

% Calculate the half-lives with different settings
if params.timepoints == 3
    calc_half_lives_3TP(data,params); %Fitting/optimizing proteins to calculate the half-lives
elseif params.timepoints == 4
    calc_half_lives_4TP(data,params); %Fitting/optimizing proteins to calculate the half-lives
elseif params.timepoints == 5
    calc_half_lives_5TP(data,params); %Fitting/optimizing proteins to calculate the half-lives
end


fprintf('\n\n\n  Completed exporting corrected half-lives to the out_file \n\n\n')

if params.apparent_T50 == 1
    Apt_HL_Calculation(params);
    fprintf('\n\n\n  Completed exporting apparent half-lives to the out_file \n\n\n *******  JUMPt program is complete *******\n\n')
end



