This is the code repository related to the manuscript 
@article{Kara2020Identification, author = {Eleanna Kara 1, Alessandro Crimi 1, Marc Emmenegger 1, Claudia Manzoni 2, 3, Sara Bandres Ciga4, Juan Botía 3, 5, Annie Wiedmer1, Manfredi Carta 1, Daniel Heinzer 1, Merve Avar 1, Andra Chincisan1, Manuela Tan3, Edwin Jabbari3, Cornelis Blauwendraat 4, Sonia Garcia Ruiz 3, Daniel Pease 1, Jiang-An Yin 1, Lorene Mottier 1, Alessandra Carrella 1, Dezirae Schneider 1, Caroline Aemisegger 6, Zhanyun  Fan 7, Jordan Marks 7, Sarah Hopp 7, 8, Huw Morris3, Patrick Lewis 2, 3, Mike Nalls 4, Mina Ryten 3, John Hardy 3, Bradley Hyman 7, Adriano Aguzzi 1  }, title = {IDENTIFICATION OF GENETIC MODIFIERS OF Α-SYNUCLEIN PROPAGATION THROUGH A HIGH-CONTENT, GENOME-WIDE RNAI SCREEN}, booktitle = {  }, year = {  }, month = {  } }  

In this project, we are quantifying the degree of a-synuclein cell to cell transfer and if/how it is modulated by each gene within the genome. To this end, we are using a construct encoding GFP-2a-aSynuclein-RFP and we are quantifying the number of recipient cells (RFP+GFP-) normalized to the total number of transfected cells (RFP+GFP+).  

This repository comprises the scripts for the image analysis and graph generations.
The used images are available on http://figshare.com/...

the script run_all_RG.m calls iterates through the entire dataset computing features and plotting results,

the script red_vs_green.m is the function actually computing the features,

while plot_data.m is the function plotting and saving the heatmaps.



 

