% A.Crimi 2019
% 
% Code related to Kara et al. 2020
%
function plot_data(filename)
 
%data = csvread(filename);
fid = fopen(filename, 'rt');  %the 't' is important!
stream = textscan(fid,'%f%f%f%f%f%s%s%s%s','HeaderLines',1,'Delimiter',',','EndOfLine','\r\n','ReturnOnError',false);
fclose(fid); 
data_averaged = cell2mat(stream(:,1:5));
 
% The convetion is to have 
% in the first column Cell count
% in the second column Norm Cell Count
% in the thid column Mean intensity 
% in the fourth column Area
% in the fifth column Green Cell count

cell_count = data_averaged(:,1);
n_cell_count = data_averaged(:,2);
mean_int = data_averaged(:,3);
area = data_averaged(:,4);
GFP_cell_count = data_averaged(:,5);

plot_heatmap(cell_count, 'CellCount');
plot_heatmap(n_cell_count, 'NCellCount');
plot_heatmap(mean_int, 'MeanInt');
plot_heatmap(area, 'Area');
plot_heatmap(GFP_cell_count, 'GFPCellCount');
