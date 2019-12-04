% A.Crimi 2019
% 
% Code related to Kara et al. 2020
%
% Main script loading all images and calling the functions to compute features and plot results
list = dir('*.tif');

results = zeros(length(list)/2,5);
result_count = 1;
 
for kk = 1 : 2:  length(list)
    
    % [n, mean_int, n_pixel, intensities] = 
    [n, n_norm, mean_int, n_pixel, intensities] =  red_vs_green( list(kk+1).name, list(kk).name);
    
    results(result_count,1)  = n;
    results(result_count,2)  = n_norm;
    results(result_count,3)  = mean_int;
    results(result_count,4)  = n_pixel;
    results(result_count,5)  = n_green;
    results(result_count,6)  = intensities;
    
    result_count = result_count +1
    
end

% Save the original features
csvwrite('results_feautes.csv',results);
  
% Average per well field
results_averaged = zeros(length(list)/4,5);  % Individual intensities are not used here, so there are only 5 files
result_count = 1;
for jj = 1 : 2:  length(list)/2
    results_averaged(result_count,:)  = mean( results(jj:jj+1,:)  );
    result_count = result_count +1;
end

count_letters = 0;
plate_list = strings(24*16,1);
for array = char(double('A'):double('P')) % create character array from 'A' to 'Z'
    for num = 1 : 24
        plate_list(count_letters+num) =  strcat(array,num2str(num)); 
    end
    count_letters = count_letters +24;
end 
results_averaged = num2cell(results_averaged);

% Match sampleID&GeneSymbol 
% Find file control
dinfo = dir('.');
pattern = 'Samples'; 
false_match = cellfun(@isempty, regexp({dinfo.name}, pattern) );
dinfo(false_match) = [];   %delete entries that accidentally matched
filename_gene = dinfo.name;
%filename_gene = 'Samples Dest_P_65-72.csv';
fid = fopen(filename_gene, 'rt');  
genelist = textscan(fid,'%s%s%s%s%s%s','HeaderLines',1,'Delimiter',',','EndOfLine','\r\n','ReturnOnError',false);
fclose(fid);
dest_well = genelist{2};
dest_well = dest_well(1:2:end);
sampleID = genelist{5};
sampleID = sampleID(1:2:end);
geneSymbol  = genelist{6};
geneSymbol  = geneSymbol(1:2:end);
matched_sampleID = strings(size(plate_list));
matched_gene = strings(size(plate_list));
for kk = 1 : length(plate_list)
    
    for jj = 1 : length(sampleID)
        tf = strcmp(plate_list(kk),string(dest_well(jj))); 

        if(tf) 
           matched_sampleID(kk) = sampleID(jj);
           matched_gene(kk) = geneSymbol(jj);
        end
    end
end

% Match sourcewell
% Find file control
dinfo = dir('.');
pattern = 'controls'; 
false_match = cellfun(@isempty, regexp({dinfo.name}, pattern) );
dinfo(false_match) = [];   %delete entries that accidentally matched
filename_control = dinfo.name;
%filename_gene = 'Picking List controls DD 65 - 72.csv';
fid = fopen(filename_control, 'rt');  
sourcelist = textscan(fid,'%s%s%s%s','HeaderLines',1,'Delimiter',',','EndOfLine','\r\n','ReturnOnError',false);
fclose(fid);
dest_well = sourcelist{2};
dest_well = dest_well(1:2:end);
source_well= sourcelist{4};
source_well = source_well(1:2:end); 

matched_source = strings(size(plate_list)); 
for kk = 1 : length(plate_list)
    
    for jj = 1 : length(dest_well)
        tf = strcmp(plate_list(kk),string(dest_well(jj))); 

        if(tf) 
           matched_source(kk) = source_well(jj); 
        end
    end
end

%Save everything
tot_feat = [ results_averaged plate_list matched_sampleID matched_gene matched_source];
%csvwrite('avg_results_features.csv',results_averaged);
fid = fopen('avg_results_features.csv','wt');
if fid>0
     fprintf(fid,'CellCount,NormCellCount,MeanInt,Area,GreenCellCount,Well,Sampleid,GeneSymbol,SourceWell\n');
     for k=1:size(tot_feat,1)
         fprintf(fid,'%s,%s,%s,%s,%s,%s,%s,%s,%s\n',tot_feat{k,:}); %%
     end
     fclose(fid);
end

% Plot and save averaged heatmaps
plot_data('avg_results_features.csv')
