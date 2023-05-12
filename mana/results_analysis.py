import glob
import matplotlib.pyplot as plt
import scipy.cluster.hierarchy as hc 
import pandas as pd
import numpy as np
import networkx as nx
from pyvis import network as net

from .modelling import fullname_equation

def get_node_list(gml_file):
	"""get_node_list.

	Parameters
	----------
	gml_file : str
		the path to the gml file

	Returns
	-------
		a list containing nodes ids from the gml file

	"""
	with open(gml_file,'r') as gml_handler:
		parsed_gml = nx.parse_gml(gml_handler)
	return(list(parsed_gml.nodes))

def visualize_gml(gml_file,window_size=['1000px','1000px'],notebook=True):
	"""visualize_gml.

	Parameters
	----------
	gml_file : str
		the path to the gml file
	window_size : list
		the list of x and y pixel sizes for the output window
	notebook : boolean
		if true enable the pyvis optimized visualisation for notebooks

	Returns
	-------
		an interactive visualisation window of the graph

	"""
	g = net.Network(window_size[0],window_size[1],notebook=notebook)
	nxg = nx.readwrite.gml.read_gml(gml_file)
	g.from_nx(nxg)
	return g.show(gml_file.replace('.gml','.html'))


def dendro_reactions(matrix,title = "No title set"):
	"""dendro_reactions.

	Parameters
	----------
	matrix : pandas dataframe
		a distance matrix
	title : str
		the title of the dendrogram

	Returns
	-------
		a scipy linkage object and show plot the dendrogram

	"""
	plt.figure(figsize=(20, 7))
	plt.title(title)
	# Create dendrogram
	linkage = hc.linkage(matrix, method='ward')
	hc.dendrogram(linkage,labels=list(matrix.index))

	plt.xlabel('Reaction ID')
	plt.ylabel('Euclidean distance')
	plt.show()
	return linkage

def generate_annotation_table(reaction_file,model,hgnc_data,scores_file):
	"""generate_annotation_table.

	Parameters
	----------
	reaction_file : pandas dataframe
		a dataframe containing the partial enumeration results
	model : cobra model
		A cobra object, loaded with the cobra library
	hgnc_data : pandas dataframe
		a hgnc database annotation dataframe
	scores_file : str
		the path to the computed scores file
	Returns
	-------
		Write an excel file with one reaction per line and corresponding annotations in columns

	"""
	comp = {"c":"Cytoplasm","m":"Mitochondrion","x":"Peroxisome","l":"Lysosome","g":"Golgi appartus",\
	 "e":"Extracellular space","r":"Endoplasic reticulum","n":"Nucleus","i":"Mitochondrial intermembrane space"}
	r_direction_hdler = open(reaction_file.replace('.tab','_direction.tab'), 'w')
	computed_scores = pd.read_csv(scores_file,sep='\t')
	computed_scores.index = computed_scores.data_id
	annot_df = pd.read_csv(reaction_file,sep=',')
	#generate columns
	annot_df["Reaction name"] = np.nan
	annot_df["Direction"] = np.nan
	annot_df["Equation"] = np.nan
	annot_df["Associated genes"] = np.nan
	annot_df["Associated genes"] = annot_df["Associated genes"].astype(object)
	annot_df["Pathway in model"] = np.nan
	annot_df["Localisation"] = np.nan
	annot_df["Associated genes links"] = ""
	annot_df["Associated genes links"] = annot_df["Associated genes links"].astype(str)
	for reaction in annot_df.iloc[:,0]:
		rid = reaction.replace('R_','')
		r = model.reactions.get_by_id(reaction.replace('R_',''))
		annot_df.loc[annot_df.iloc[:,0] == reaction,"Reaction name"] = r.name
		if ((computed_scores.loc[rid,'f_treatment']-computed_scores.loc[rid,'f_ctrl']) > 0):
			annot_df.loc[annot_df.iloc[:,0] == reaction,"Direction"] = "UP"
			r_direction_hdler.write(reaction+'\t'+str(1)+"\n")
		elif(computed_scores.loc[rid,'f_treatment']-computed_scores.loc[rid,'f_ctrl'] < 0):
			annot_df.loc[annot_df.iloc[:,0] == reaction,"Direction"] = "DOWN"
			r_direction_hdler.write(reaction+'\t'+str(-1)+"\n")
		else:
			annot_df.loc[annot_df.iloc[:,0] == reaction,"Direction"] = "UNDETERMINED"
			r_direction_hdler.write(reaction+'\t'+str("")+"\n")
		annot_df.loc[annot_df.iloc[:,0] == reaction,"Equation"] = fullname_equation(r)
		#convert HGNC IDs to Hugo symbols
		list_genes = []
		for elem in r.genes:
			list_genes.append(hgnc_data.loc[hgnc_data["HGNC ID"] == elem.id,'Approved symbol'].to_list()[0])
		if((len(list_genes) == 0) & (len(r.genes) != 0)):
			list_genes = r.genes
		annot_df.at[annot_df.loc[annot_df.iloc[:,0] == reaction,"Associated genes"].index[0],"Associated genes"] = list_genes
		#create genecards url with approved symbols
		list_urls = []
		for gene in list_genes:
			list_urls.append("https://www.genecards.org/cgi-bin/carddisp.pl?gene="+str(gene))
		if(len(list_urls) == 0):
			annot_df.at[annot_df.loc[annot_df.iloc[:,0] == reaction,"Associated genes links"].index[0],"Associated genes links"] = ""
		elif(len(list_urls) == 1):
			annot_df.at[annot_df.loc[annot_df.iloc[:,0] == reaction,"Associated genes links"].index[0],"Associated genes links"] = list_urls[0]
		else:
			for i in range(0,len(list_urls)):
				if(8+i> annot_df.shape[1]): #8 is the default number of columns
					annot_df["Associated genes links"+str(i)] = ""
					annot_df["Associated genes links"+str(i)] = annot_df["Associated genes links"+str(i)].astype(str)
				if(i == 0):
					annot_df.at[annot_df.loc[annot_df.iloc[:,0] == reaction,"Associated genes links"].index[0],"Associated genes links"] \
						= annot_df.at[annot_df.loc[annot_df.iloc[:,0] == reaction,"Associated genes links"].index[0],"Associated genes links"] + \
						list_urls[i]
				else:
					annot_df.at[annot_df.loc[annot_df.iloc[:,0] == reaction,"Associated genes links"+str(i)].index[0],"Associated genes links"+str(i)]\
						  = annot_df.at[annot_df.loc[annot_df.iloc[:,0] == reaction,"Associated genes links"+str(i)].index[0],"Associated genes links"+str(i)] + \
						list_urls[i]
		if "array([]," in r.subsystem: #if no subsystem, leave cell empty
			annot_df.loc[annot_df.iloc[:,0] == reaction,"Pathway in model"] = ""
		else:
			annot_df.loc[annot_df.iloc[:,0] == reaction,"Pathway in model"] = r.subsystem.split("'")[1]
		#get compartment, if more than one compartment annotated as transport
		comps = set()
		compsR = set()
		compsP = set()
		for elem in r.reactants:
			comps.add(elem.compartment)
			compsR.add(elem.compartment)
		for elem in r.products:
			comps.add(elem.compartment)
			compsP.add(elem.compartment)
		
		if(len(comps) == 2):
			for elem in annot_df.loc[annot_df.iloc[:,0] == reaction,"Equation"].to_list()[0].split(' '):
				if elem == '-->' or elem == '<--' or elem == '<=>':
					if(len(compsR)>1 and len(compsP)>1):
						#its a co-transport
						annot_df.loc[annot_df.iloc[:,0] == reaction,"Localisation"] = "CoTransport"+'('+comp[list(compsR)[0]]+\
							'/'+comp[list(compsR)[1]]+elem+comp[list(compsP)[0]]+'/'+comp[list(compsP)[1]]+')'
					elif(len(compsR)>1 and len(compsP)==1):
						annot_df.loc[annot_df.iloc[:,0] == reaction,"Localisation"] = "CoTransport"+'('+comp[list(compsR)[0]]+\
							'/'+comp[list(compsR)[1]]+elem+comp[compsP.pop()]+')'
					elif(len(compsR)==1 and len(compsP)>1):
						annot_df.loc[annot_df.iloc[:,0] == reaction,"Localisation"] = "CoTransport"+'('+comp[compsR.pop()]+\
							elem+comp[list(compsP)[0]]+'/'+comp[list(compsP)[1]]+')'
					else:
						annot_df.loc[annot_df.iloc[:,0] == reaction,"Localisation"] = "Transport"+'('+comp[compsR.pop()]+\
							elem+comp[compsP.pop()]+')'
					break
		else:
			annot_df.loc[annot_df.iloc[:,0] == reaction,"Localisation"] = comp[comps.pop()]
	r_direction_hdler.close()
	#Excel writer (jmcnamara on stackoverflow)
	writer = pd.ExcelWriter('working_files/clusters_tables/'+reaction_file.replace('.tab','')+'_table.xlsx', engine='xlsxwriter')
	# Convert the dataframe to an XlsxWriter Excel object.
	annot_df.to_excel(writer, sheet_name="DAR_enriched_table", index=False)
	# Get the xlsxwriter workbook and worksheet objects.
	workbook  = writer.book
	worksheet = writer.sheets["DAR_enriched_table"]
	# Add a text wrap format.
	text_wrap_format = workbook.add_format({'text_wrap': True})
	# Set columns format.
	worksheet.set_column(0,0,15) #Reaction ID
	worksheet.set_column(1,1,30) #Reaction ID, Reaction name
	worksheet.set_column(2,2,10) #Direction
	worksheet.set_column(3,4,30) #Equation, Associated genes
	worksheet.set_column(5,6,45) #Pathway in model, Localisation
	worksheet.set_column(7, annot_df.shape[1], 60, text_wrap_format) # genes columns
	#add an autofilter
	worksheet.autofilter(0, 0, annot_df.shape[0], annot_df.shape[1] - 1)
	# Close the Pandas Excel writer and output the Excel file.
	writer.save()