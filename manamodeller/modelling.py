
import os
import re
import numpy as np
import pandas as pd
from cobra import core
from ast import And, BoolOp, Expression, Name, Or
from tqdm.auto import tqdm
from .utils import make_pickle, make_csvs

def get_GPR_reactions(model):
	"""get_GPR_reactions.

	Parameters
	----------
	model : cobra model
		A cobra object, loaded with the cobra library

	Returns
	-------
	pandas dataframe
		a dataframe with all the reactions in the model having a GPR

	"""
	gpr = pd.DataFrame([], columns=['Reaction ID','Reaction Name','GPR'])
	for r in range(len(model.reactions)):
		reaction = model.reactions[r]
		if reaction.gene_reaction_rule != '':
			line = [reaction.id,reaction.name,reaction.gene_reaction_rule]
			gpr.loc[len(gpr)] = line
	return gpr

def eval_gpr_activity(expr, gh, gl):
	"""eval_gpr_activity.

	This is an adaptation of the eval_gpr function available in cobrapy.
	Instead of evaluating if a GPR is active according to a list of knockout
	genes, it will evaluate if the GPR is regulated by a list of genes

	Exemple of usage :
	provide the list of Highly expressed genes and Lowly expressed genes
	Return all expressions that are True, thus Highly expressed

	evaluate compiled ast of gene_reaction_rule with list of active genes

	Parameters
	----------
	expr : Expression
		The ast of the gene reaction rule
	gh : list
		list of highly expressed genes
	gl : list
		list of lowly expressed genes
	Returns
	-------
	bool
		True if the reaction is active with the given gh and gl lists
		otherwise false
	"""
	if isinstance(expr, Expression):
		return eval_gpr_activity(expr.body, gh, gl)
	elif isinstance(expr, Name):
		if expr.id in gh:
			return 1
		elif expr.id in gl:
			return -1
		else:
			return 0
	elif isinstance(expr, BoolOp):
		op = expr.op
		if isinstance(op, Or):
			return max(eval_gpr_activity(i, gh, gl) for i in expr.values)
		elif isinstance(op, And):
			return min(eval_gpr_activity(i, gh, gl) for i in expr.values)
		else:
			raise TypeError("unsupported operation " + op.__class__.__name__)
	elif expr is None:
		return False
	else:
		raise TypeError("unsupported operation  " + repr(expr))
		
def find_reactions_expression_levels(gprs,gh,gl):
	"""find_reactions_expression_levels.

	Parameters
	----------
	gprs : cobra gpr
		A cobra gpr object, obtained from a cobra model
	gh : list
		list of highly expressed genes
	gl : list
		list of lowly expressed genes

	Returns
	-------
	tuple
		a tuple of three list.
		rh, is the list of reactions identified as active according to gene expression data
		rl, is the list of reactions identified as inactive according to gene expression data
		rn, is the list of reactions not constrained by gene expression data

	"""
	rh = []
	rl = []
	rn = []
	for i in range(len(gprs)):
		res = eval_gpr_activity(core.gene.GPR().from_string(gprs.iloc[i].GPR).body,gh,gl)
		if res > 0:
			rh.append(gprs.iloc[i][0])
		elif res < 0:
			rl.append(gprs.iloc[i][0])
		else:
			rn.append(gprs.iloc[i][0])
	return rh,rl,rn

def get_reactions_ids(model):
	"""get_reactions_ids.

	Parameters
	----------
	model : cobra model
		A cobra object, loaded with the cobra library

	Returns
	-------
	pandas dataframe
		a dataframe with all the reactions and eventually their GPR

	"""
	all_reactions_ids = pd.DataFrame([],columns=['Reaction ID','Reaction Name','GPR'])
	for r in range(len(model.reactions)):
		reaction = model.reactions[r]
		if reaction.gene_reaction_rule == '':
			line = [reaction.id,reaction.name,'No GPR']
			all_reactions_ids.loc[len(all_reactions_ids)] = line
		else:
			line = [reaction.id,reaction.name,reaction.gene_reaction_rule]
			all_reactions_ids.loc[len(all_reactions_ids)] = line
	return all_reactions_ids

def get_gene_list(model):
	"""get_gene_list.

	Parameters
	----------
	model : cobra model
		A cobra object, loaded with the cobra library

	Returns
	-------
	list
		a list containing all the genes in the model

	"""
	gene_list = []
	for gene in model.genes:
		gene_list.append(gene.id)
	return gene_list

def identify_model_gene_ids(model):
	"""identify_model_gene_ids.

	Parameters
	----------
	model : cobra model
		A cobra object, loaded with the cobra library

	Returns
	-------
	str
		The type of gene identifiers used in the model

	"""
	model_genes = get_gene_list(model)
	if 'HGNC:' in model_genes[0]:
		return 'HGNC ID'
	elif 'ENSG' in model_genes[0]:
		return 'Ensembl gene ID'
	elif '_AT' in model_genes[0]:
		return 'NCBI Gene ID'
	else:
		return 'model not implemented'
	
def fullname_equation(reaction):
	"""fullname_equation.

	Parameters
	----------
	reaction : cobra reaction
		A cobra reaction object, obtained from a cobra model

	Returns
	-------
	str
		the reaction with complete metabolites names

	"""
	reaction_metabolites = list(reaction.metabolites)
	#generate metabolites dict
	metabolites_dict = {}
	for met in reaction_metabolites:
		metabolites_dict[met.id] = met.name 
	splitted_equation = reaction.reaction.split(' ')
	fullname_equation = ""
	for elem in splitted_equation:
		if re.match(".*\[.\].*",elem):
			fullname_equation = fullname_equation+metabolites_dict[elem]+" "
		else:
			fullname_equation = fullname_equation+elem+" "
	return fullname_equation
	
def map_single_column(data,hgnc_data,col_to_add):
	"""map_single_column.

	Parameters
	----------
	data : pandas dataframe
		a gene expression dataset
	hgnc_data : pandas dataframe
		a hgnc database annotation dataframe
	col_to_add : str
		the name of the identifier column to map

	Returns
	-------
	pandas dataframe
		the gene expressed dataframe with the mapped identifier added
	"""
	mapped_ids = []
	genes = data['ENTREZID']
	for i in range(len(genes)):
		id = hgnc_data.loc[hgnc_data['NCBI Gene ID'] == str(genes[i])][col_to_add]
		if len(id) == 0:
			mapped_ids.append('NA')
		else:
			mapped_ids.append(id.iloc[0])
	if len(mapped_ids) == len(data):
		data.insert(2,col_to_add,mapped_ids)
	return data

def find_high_low_exprs(uarray_data,threshold_dw_perc=25,threshold_up_perc=75):
	"""find_high_low_exprs.

	Parameters
	----------
	uarray_data : pandas dataframe
		Description of parameter `uarray_data`.
	threshold_dw_perc : int
		The percentile below which we consider that genes are not expresed.
	threshold_up_perc : int
		The percentile above which we consider that genes are highly expressed.

	Returns
	-------
	list
		Return list of highly/lowly expressed gene .

	"""
	# threshold_dw_perc : Low expression percentile
	# threshold_up_perc : High expression percentile
	if pd.DataFrame(uarray_data).shape[1] > 1:
		raise ValueError("More than one column provided, check duplicates")
	if (0 <= threshold_up_perc <= 100) & (0 <= threshold_dw_perc <= 100) & \
	 (threshold_dw_perc < threshold_up_perc):
		threshold_dw = np.percentile(uarray_data,threshold_dw_perc)
		threshold_up = np.percentile(uarray_data,threshold_up_perc)
		low_exprs = uarray_data.loc[uarray_data.iloc[:,] <= threshold_dw]
		high_exprs = uarray_data.loc[uarray_data.iloc[:,] >= threshold_up]
	return [high_exprs,low_exprs]

def preprocess_data(data,gene_id_col,model,pickle_path="",csvs_path="", threshold_dw_perc = 25, threshold_up_perc = 75):
	"""preprocess_data.

	Parameters
	----------
	data : pandas dataframe
		a gene expression dataset
	gene_id_col : str
		the name of the gene identifier column
	model : cobra model
		A cobra object, loaded with the cobra library
	pickle : str
		if not empty, write pkl files containing categorized reactions activity at given location
	csvs : str
		if not empty, write csv files containing categorized reactions activity at given location
	threshold_dw_perc : int
		The percentile below which we consider that genes are not expresed.
	threshold_up_perc : int
		The percentile above which we consider that genes are highly expressed.

	Returns
	-------
		Write a csv or a pickle with for categorized reactions activity

	"""
	pbar = tqdm(total=len(data.filter(like=".CEL",axis=1).columns))
	suffix = '_rh_rl_zscores_'+str(threshold_dw_perc)+'_'+str(threshold_up_perc)
	gprs = get_GPR_reactions(model)
	for i in range(len(data.columns)):
		uarray = data.columns[i]
		if ('.CEL' in uarray):
			uarray_data = data[uarray]
			uarray_data.index = data[gene_id_col]
			gh,gl = find_high_low_exprs(uarray_data,threshold_dw_perc,threshold_up_perc)
			rh,rl,rn = find_reactions_expression_levels(gprs,gh,gl)
			if pickle_path != "":
				picklef = uarray+suffix
				if os.path.exists(pickle_path):
					make_pickle([rh,rl,rn],pickle_path+'/'+picklef+'.pkl')
				else:
					os.mkdir(pickle_path)
					make_pickle([rh,rl,rn],pickle_path+'/'+picklef+'.pkl')
			if csvs_path != "":
				csvsf = uarray+suffix
				if os.path.exists(csvs_path):
					make_csvs([rh,rl,rn],csvs_path+"/",csvsf)
				else:
					os.mkdir(csvs_path)
					make_csvs([rh,rl,rn],csvs_path+"/",csvsf)
			pbar.update(1)
