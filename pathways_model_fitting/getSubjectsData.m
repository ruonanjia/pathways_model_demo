function [id_unique, data_raw] = getSubjectsData(file_name)
%Extracts subjects data and ids
% Input:
%      - file_name: name of the data sheet
% Output:
%      - id_unique: ids of all subjects
%      - data_raw: data sheet, matlab table

% read .csv
data_raw = readtable(file_name, 'Delimiter', ',','TreatAsEmpty','N/A');

% get rid of NaN lines
% toDelete = isnan(data_raw.is_med);
% data_raw(toDelete,:) = [];


% take our all subject ID
id = data_raw.id;

% get rid of repetition
id_unique = unique(id);


end