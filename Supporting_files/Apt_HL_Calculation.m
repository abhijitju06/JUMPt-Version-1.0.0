function data = Apt_HL_Calculation(params)
%%%%% Apparent T50 calculation
%%%%% Abhijit Dasgupta
%% 
%read parameter file
fid     = fopen('JUMPt.params');
param   = textscan(fid,'%s','delimiter','\n');
cellfun(@eval,param{1});
fclose(fid);
outfile	= join(['results_Apparent_T50', params.input_file], "_");

params.setting      = setting;
params.input_file   = input_file;
params.bin_size     = bin_size;
params.opti_algo    = optimization_algorithm;
params.purity       = purity_of_SILAC_food;

All_data    = readtable(params.input_file);
time_points_ind 	= startsWith(All_data.Properties.VariableNames, 'time_point');
conc_ind            = startsWith(All_data.Properties.VariableNames, 'concentration');
data.t              = table2array(All_data(1,time_points_ind));  
data.t_long         = linspace(0, 32,321);
    
data.TotalLysConc	= table2array(All_data(2,conc_ind));
data.LysRatio    	= table2array(All_data(3,time_points_ind));% Lys data
data.freeLysConc    = table2array(All_data(3,conc_ind));
data.ProtInfo       = (All_data(4:end,1:find(time_points_ind,1)-1)); % protein IDs
data.ProtConc       = table2array(All_data(4:end,conc_ind))';% Protein concentration
SILAC_data          = table2cell(All_data(4:end,time_points_ind));% protein data
SILAC_data_NaN      = cellfun(@(x) ~isa(x,'double'),SILAC_data); % check for possible 0x0 char reads
SILAC_data(SILAC_data_NaN)   = {NaN};% set those to NaN
data.SILAC_data              = cell2mat(SILAC_data)'; % heavy/light ratios,
data.timepoints = All_data.Properties.VariableNames(time_points_ind);

% calculate apparent half-lives to sort the proteins
size_Prot_data  = size(data.SILAC_data,2);
fun_deg = @(b,tspan) exp(-b*tspan); 

for i =1:size_Prot_data 
    Pdata_temp  = data.SILAC_data(:,i);
	non_NaN_index_of_Pdata_temp = find(~isnan(Pdata_temp));
	Pdata_temp  = Pdata_temp(non_NaN_index_of_Pdata_temp); 
    rng default % For reproducibility
	ExpDegRate_temp         = fminsearch(@(b) norm(Pdata_temp - fun_deg(b,data.t(non_NaN_index_of_Pdata_temp))), rand(1,1)) ; 
	apparent_half_lives(i)  = ExpDegRate_temp;
end

apparent_half_lives = transpose(log(2)./apparent_half_lives);
[~, idx]            = sort(apparent_half_lives);

data.SILAC_data     = data.SILAC_data(:,idx); 
data.ProtConc       = data.ProtConc(:,idx);
data.ProtInfo       = data.ProtInfo(idx, :); 

SORT_SILAC = transpose (data.SILAC_data);
sorted_apparent_half_lives = num2cell(sort(apparent_half_lives));
protein_name	= table2array(data.ProtInfo(:,1));
gene_name = table2array(data.ProtInfo(:,2));

file_content = cat(2,protein_name,gene_name,sorted_apparent_half_lives);
col_header={'Protein','Gene','Apparent T50'}; 

output_matrix=[col_header; file_content]; 

xlswrite(outfile,output_matrix);
