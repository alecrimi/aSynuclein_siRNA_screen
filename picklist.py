# -*- coding: utf-8 -*-
"""
Created on Thu Jan  4 13:15:39 2018

@author: achin
"""
import pandas as pd

def read_plate_layout (path):
    """
    Read plate layout
    """
    plate_layout = "plate_layout.csv"
    df = pd.read_csv (path + plate_layout)
    return df
    
def read_source_plates (path):
    """
    Read source plates from 1 to 62
    """
    range_files = range(1,63)
    file_names = []
    for i in range_files:
        path_file = path +  "plate_P_" + str(i) +  ".csv"
        file_names.append (path_file)
    return file_names
   
def process_file (path, ind):
    file = path + 'plate_P_' + str(ind)+ '.csv'
    df = pd.read_csv (file)
    return df
    
def create_duplicates_source_plates (path, no_plates):   
    range_plates = []
    for i in range(1, no_plates+1):
        repeat_ind = 0
        while repeat_ind < 2:
            range_plates.append(i)
            repeat_ind = repeat_ind + 1
    frames = [ process_file(path, i) for i in range_plates]
    result = pd.concat(frames)
    path_file_duplicates =  (path[0:len(path)-7] 
                             + 'plate_all_duplicates.csv') 
    result.to_csv (path_file_duplicates, index = False)  
     

def row_mapping (value):
    if value == 0:
        row_value = 'A'
    elif value == 1:
        row_value = 'B'
    elif value == 2:
        row_value = 'C'
    elif value == 3:
        row_value = 'D'
    elif value == 4:
        row_value = 'E'
    elif value == 5:
        row_value = 'F'
    elif value == 6:
        row_value = 'G'
    elif value == 7:
        row_value = 'H'
    elif value == 8:
        row_value = 'I'
    elif value == 9:
        row_value = 'J'
    elif value == 10:
        row_value = 'K'
    elif value == 11:
        row_value = 'L'
    elif value == 12:
        row_value = 'M'
    elif value == 13:
        row_value = 'N'
    elif value == 14:
        row_value = 'O'
    elif value == 15:
        row_value = 'P'
    return row_value
   
def write_dest_plate_list (source_well_col, dest_well_col, dest_plate_col, 
                       source_plate_bar_code_col, sample_id_col, gene_symbol_col,
                       transfer_vol_col, path): 
    d = {'Source well': source_well_col,
     'Destination well' : dest_well_col,
     'Destination plate bar code': dest_plate_col,
     'Source plate bar code': source_plate_bar_code_col,
     'Sample id': sample_id_col,
     'Gene symbol': gene_symbol_col,
     'Transfer volume' : transfer_vol_col
     }
    df_picklist= pd.DataFrame(d)
    path_file = path[0:len(path)-7] + 'picklist_all_test.csv' 
    df_picklist.to_csv (path_file, 
                        columns=['Source well', 'Destination well', 
                                 'Destination plate bar code',
                                 'Source plate bar code', 'Sample id', 
                                 'Gene symbol', 'Transfer volume'], index=False)
    
def destination_plate_well (df, sample_position): 
    """ df is a data frame containing the layout of the destination plate 
    / dest_layout """
    list_cols = list(df.columns.values)
    for index, row in df.iterrows():
        for col in list_cols:
            if df[col].iloc[index].isdigit() == True:
                if df[col].iloc[index] == sample_position:
                    print ('yes', index, col)
                    row_value = row_mapping (index)
                    well = row_value + col[3:len(col)]
                    return well
 
path_layout = "E:\\destination plates\\"
path_source_files = 'W:\\Neurolab\\AndraChincisan\\4.sRNA screen\\new_plates\\plates\\'
    
df = read_plate_layout (path_layout)
file_names = read_source_plates (path_source_files) 
file_1 = 'W:\\Neurolab\\AndraChincisan\\4.sRNA screen\\new_plates\\plates\\plate_P_1.csv' 
create_duplicates_source_plates (path_source_files, len(file_names))   
file_1 = 'W:\\Neurolab\\AndraChincisan\\4.sRNA screen\\new_plates\\plate_all_duplicates.csv' 
file_1 = 'W:\\Neurolab\\AndraChincisan\\4.sRNA screen\\new_plates\\plate_1_duplicates.csv'
df_source = pd.read_csv(file_1)
# Fields
count_df1_ind = 0
list_cols = list(df_source.columns.values)
for index, row in df_source.iterrows():
    #for col in list_cols:
    count_df1_ind = count_df1_ind  + 1 

""" Columns for destination plate data frame """
source_well_col = []
source_plate_bar_code_col = []
sample_id_col = []
gene_symbol_col = []
transfer_vol_col = []
source_samples_col = []
dest_plate_col = []
dest_well_col = []

max_dest_samples = 264 # define maximum no of samples in a destination plate
dest_no = 1 # number of destination plates (used as destination plate id)
source_samples = 0 # counter of samples from source plates / up to 264 an then reset
non_emp = 0
for index, row in df_source.iterrows():
    """ Get source plate values """
    dest_well = ''
    source_well = df_source['Location (Row-Col)'].iloc[index] # source well (row, col)
    source_plate_bar_code = df_source['source_plate_id'].iloc[index] # source plate id
    sample_id = df_source['geneid'].iloc[index] # sample id
    gene_symbol = df_source['Gene Symbol'].iloc[index] # unique gene symbol
    transfer_vol = 30 # transfer volume quantity
    #if gene_symbol == gene_symbol:
    source_samples = source_samples + 1
    dest_well = destination_plate_well (df, str(source_samples))
    #print (source_samples, dest_well)
    if source_samples == max_dest_samples:
        dest_no = dest_no + 1
        source_samples = 1
            
    dest_plate = 'Dest_P_' + str(dest_no)
    if not ((source_well[1:len(source_well)] == '23') or 
             (source_well[1:len(source_well)] == '24')):  
        non_emp =non_emp + 1
        #if gene_symbol == gene_symbol and sample_id == sample_id:
        source_well_col.append (source_well)
        source_plate_bar_code_col.append (source_plate_bar_code)
        sample_id_col.append (sample_id)
        gene_symbol_col.append (gene_symbol)
        transfer_vol_col.append (transfer_vol)
        source_samples_col.append (source_samples)
        dest_plate_col.append (dest_plate)
        dest_well_col.append (dest_well)
#            print (source_well, source_plate_bar_code, sample_id, gene_symbol, 
#                   transfer_vol, source_samples, dest_plate, dest_well)

write_dest_plate_list (source_well_col, dest_well_col, dest_plate_col, 
                       source_plate_bar_code_col, sample_id_col, gene_symbol_col,
                       transfer_vol_col, path_source_files)    

    


#dest_no = 0
#
#dest_plate = 'D' + str(dest_no)
#for index, row in df.iterrows():
#    for col in list_cols:
#        if df[col].iloc[index] == 'blank':
#            print ('exclude: ', df[col].iloc[index])
#        elif df[col].iloc[index].isdigit() == True:
#            print ('Samples: ', df[col].iloc[index])
#        else:
#             if df[col].iloc[index] == 'neg_1':
#                 print ('control neg_1: ', df[col].iloc[index])
#             elif df[col].iloc[index] == 'pos_Prnp':
#                 print ('control pos_Prnp: ', df[col].iloc[index])
#             elif  df[col].iloc[index] == 'pos_Hscd': 
#                 print ('control pos_Hscd: ', df[col].iloc[index])
#        count_df_ind = count_df_ind + 1
#
#
#                
#            print ('exclude: ', df[col].iloc[index])
#        elif df[col].iloc[index].isdigit() == True:
#            print ('Samples: ', df[col].iloc[index])
#        else:
#             if df[col].iloc[index] == 'neg_1':
#                 print ('control neg_1: ', df[col].iloc[index])
#             elif df[col].iloc[index] == 'pos_Prnp':
#                 print ('control pos_Prnp: ', df[col].iloc[index])
#             elif  df[col].iloc[index] == 'pos_Hscd': 
#                 print ('control pos_Hscd: ', df[col].iloc[index])
#        count_df_ind = count_df_ind + 1
            
    
#Source Well = Location (Row-Col)
#Destination Well = ?
#Destination Plate Barcode = ?
#Source Plate Barcode = source_plate_id
#Sample ID = geneid (max ID21584 )
#Gene Symbol = Gene Symbol 
#Transfer Volume = 30
#(Left: Mimic Picklist All Plates.xlsx, (Daniel), Right: plate_all.csv(Andra))
