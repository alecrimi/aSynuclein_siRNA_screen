# -*- coding: utf-8 -*-
"""
Created on Fri Dec 22 09:55:35 2017

@author: achin
"""
import numpy as np
import pandas as pd

# Create the plate map according to the indexing for selection by position
def create_plate_map (Data_Frame, index_start,index_end, temp_repeat, previous_val):
    Data_Frame_temp = Data_Frame[index_start:index_end]
    del Data_Frame_temp['Gene ID']
    del Data_Frame_temp['siRNA ID']
    del Data_Frame_temp['Row']
    del Data_Frame_temp['Col']
    gene_id = [] 
    if temp_repeat == 0:
        initial_val = index_start  
    else:
        initial_val = previous_val
    index = 0  
    for x in range (index_start,index_end):
        row_col = Data_Frame_temp['Location (Row-Col)'].iloc[index]
        ref_acc = Data_Frame_temp['RefSeq Accession'].iloc[index]
        index = index + 1
        if (row_col[1:len(row_col)] =='23' or row_col[1:len(row_col)] =='24' 
        or  ref_acc != ref_acc):
            gene_id.append('')
        else:
            initial_val = initial_val + 1
            gene_id.append(str(initial_val))
    return Data_Frame_temp, library, gene_id, initial_val

# Load data according to plate name
def get_library (full_library_name):
    #print (full_library_name)
    separator_pos = full_library_name.find('-')
    library_name_no = full_library_name[separator_pos:len(full_library_name)] 
    library_name = library_name_no[1]
    return library_name

def process_file (ind):    
    file = "plate_P_" + str(ind)+ ".csv"
    df = pd.read_csv (file)
    return df
frames = [ process_file(i) for i in range(1, no_plates+1)]
result = pd.concat(frames)
path = 'plate_all.csv' 
result.to_csv (path, index = False)  
    
""" Read the csv sheet to pandas dataframe """
index_start = 0
index_end = 0
gene_id = 1
count_A = 0
count_B = 0
count_C = 0
temp_repeat = 0
Data_Frame = pd.read_csv('Silencer Select hGenome List.csv') # E:\\files\\Silencer Select hGenome List.csv

indexes_A = []
indexes_B = []
indexes_C = []
count_all_rows_A = 0
count_all_rows_B = 0
count_all_rows_C = 0
# Load data according to plate name
for index, row in Data_Frame.iterrows():
     library = get_library(Data_Frame['Plate Name'].iloc[index])
     if library == 'A':
         indexes_A.append (index)
         count_all_rows_A = count_all_rows_A + 1
     if library == 'B':
         indexes_B.append (index)
         count_all_rows_B = count_all_rows_B + 1
     if library == 'C':
         indexes_C.append (index)
         count_all_rows_C = count_all_rows_C + 1

# Purely integer-location based indexing for selection by position
data_frame_A = Data_Frame.iloc[indexes_A]
data_frame_B = Data_Frame.iloc[indexes_B]
data_frame_C = Data_Frame.iloc[indexes_C]

# Save the indexed data into csv files
data_frame_A.to_csv ('Library_A.csv', index = False)
data_frame_B.to_csv ('Library_B.csv', index = False)
data_frame_C.to_csv ('Library_C.csv', index = False)

count_A = 0
temp_repeat = 0
data_frame_A['Plate Name'] = 'Hm Drgbl Gen siRNA Lib-ABC'
no_plates = 0
 
# Map source and destination according to plate picklists reported    
for ind, row in data_frame_A.iterrows():    
    gene = data_frame_A['RefSeq Accession'].iloc[count_A] 
    well =  data_frame_A['Location (Row-Col)'].iloc[count_A] 
    if well == 'A1':
        index_start = count_A
 
    if well == 'P24':
        index_end = count_A + 1
        print (index_start, index_end)
        if temp_repeat == 0:
            previous_val = 0
             
            Data_Frame_temp, library, gene_id_val, initial_val = create_plate_map (data_frame_A, 
                                                                                   index_start,                                                                       
                                                                                   index_end, 
                                                                                   temp_repeat, 
                                                                                   previous_val)

        else:
            previous_val = initial_val
            Data_Frame_temp, library, gene_id_val, initial_val = create_plate_map (data_frame_A, 
                                                                                   index_start, 
                                                                                   index_end, 
                                                                                   temp_repeat, 
                                                                                   previous_val)
        Data_Frame_temp = Data_Frame_temp.assign(geneid=pd.Series(np.array(gene_id_val)).values)   
        Data_Frame_temp.geneid[Data_Frame_temp.geneid != ''] = 'ID' + Data_Frame_temp.geneid
        source_plate_id = [str(no_plates+1) for i in range(0,len(Data_Frame_temp))]
        Data_Frame_temp = Data_Frame_temp.assign(source_plate_id=pd.Series(np.array(source_plate_id)).values) 
        Data_Frame_temp.source_plate_id[Data_Frame_temp.source_plate_id != ''] = 'Source_P' + Data_Frame_temp.source_plate_id
        #Data_Frame_temp['source_plate_id']= 'Source_P' + Data_Frame_temp['source_plate_id'] 
        no_plates = no_plates + 1
        
        # Save the mapping
        path = 'plates\\plate' + '_' + 'P' + '_' + str(no_plates) + '.csv'
        Data_Frame_temp.to_csv (path, index = False)  
        
        temp_repeat = temp_repeat + 1

    """ Count index data frame """
    count_A = count_A + 1
    
