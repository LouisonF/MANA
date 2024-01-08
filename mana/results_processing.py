import os
import re
import glob
import pandas as pd

from multiprocessing import JoinableQueue
from .utils import launch_multi_proc

def concatenate_reaction_div_enum(path_concat_rxn_enum,path_concat_div_enum, out_dir,col_index="",single_csv=False,ncpus=1):
	"""concatenate_reaction_div_enum.

	Parameters
	----------
	path_concat_rxn_enum : str
		path to the concatenated reaction enum directory
	path_concat_div_enum : str
		path to the concatenated diversity enum directory
	out_dir : str
		path to the csvs output directory
	col_index : str
		column name of a column to be used as index (optional)
	single_csv : boolean
		option for the concatenate_csv function, if True all solutions will be stored in a single csv file
	ncpus : int
		the number of cpus allocated, will enable parallel processing
	Returns
	-------
		a JoinableQueue Object 

	"""
	#source : https://www.freecodecamp.org/news/how-to-combine-multiple-
	#csv-files-with-8-lines-of-code-265183e0854/
	q = JoinableQueue()
	list_dir = glob.glob(path_concat_rxn_enum+"/*.csv")
	for rxn_enum_file in list_dir:
			files = glob.glob(path_concat_div_enum+'/'+os.path.basename(rxn_enum_file).split('_')[0]+'*_solutions.csv')
			files.append(rxn_enum_file)
			if ncpus == 1:
				concatenate_csv(files,out_dir,col_index,single_csv)
			else:
				q.put((concatenate_csv,(files,out_dir,col_index,single_csv)))
	if ncpus == 1:
		pass
	else:
		launch_multi_proc(ncpus,q)
	return q

def concatenate_solutions(csv_dir,out_dir,col_index="",single_csv=False,ncpus=1,
 restart=False):
	"""concatenate_solutions.

	Parameters
	----------
	csv_dir : str
		path to the csvs to concatenate directory
	out_dir : str
		path to the csvs output directory
	col_index : str
		column name of a column to be used as index (optional)
	single_csv : boolean
		option for the concatenate_csv function, if True all solutions will be stored in a single csv file
	ncpus : int
		the number of cpus allocated, will enable parallel processing
	Returns
	-------
		a JoinableQueue Object 

	"""
	#source : https://www.freecodecamp.org/news/how-to-combine-multiple-
	#csv-files-with-8-lines-of-code-265183e0854/
	if (os.path.basename(out_dir) == "full_rxn_enum_set"):
		index_suffix = '_renum'
	elif (os.path.basename(out_dir) == "full_div_enum_set"):
		index_suffix = '_rdivers'
	else:
		return "Wrong outdir name"
	q = JoinableQueue()
	list_dir = list(set([i.split('_')[0] for i in os.listdir(csv_dir)]))
	for i in range(0,len(list_dir)):
		if restart:
			#search if a concatenated csv exist for this file id in the full rxn set
			rgx = re.compile(csv_dir+list_dir[i].split('_')[0]+'.*')
			match_res = list(filter(rgx.match,os.path.basename(csv_dir)))
			if len(match_res) == 0:
				files = glob.glob(csv_dir+list_dir[i].split('_')[0]+'*_solutions.csv')
				if ncpus == 1:
					concatenate_csv(files,out_dir,col_index,single_csv,index_suffix)
				else:
					q.put((concatenate_csv,(files,out_dir,col_index,single_csv,index_suffix)))
		else:
			files = glob.glob(csv_dir+list_dir[i].split('_')[0]+'*_solutions.csv')
			if ncpus == 1:
				concatenate_csv(files,out_dir,col_index,single_csv,index_suffix)
			else:
				q.put((concatenate_csv,(files,out_dir,col_index,single_csv,index_suffix)))
	if ncpus == 1:
		pass
	else:
		launch_multi_proc(ncpus,q)
	return q


def concatenate_csv(filenames,out_dir,col_index,single_csv,index_suffix=""):
	"""concatenate_csv.

	Parameters
	----------
	filenames : str
		list of csv files to concatenate into one csv file
	out_dir : str
		path to the concatenated csv output directory
	col_index : str
		column name of a column to be used as index (optional)
	single_csv : boolean
		option for the concatenate_csv function, if True all solutions will be stored in a single csv file
	index_suffix : str
		suffix to add to csv's row index
	Returns
	-------
		write the concatenated csv in the 

	"""
	list_csvs = []
	index = []
	nrenum = 0
	if col_index == "":
		col_index = list(pd.read_csv(filenames[0]).columns)
	for i in range(len(filenames)):
		tmp = pd.read_csv(filenames[i])
		#check that the number of columns match and get colnames of the first file
		if len(col_index) != len(tmp.columns):
			print("Error")
		else:
			tmp.columns = col_index
		if len(os.path.basename(filenames[i]).split('_')) == 2:
			file_id = os.path.basename(filenames[i]).split('_')[0]+'_full_rxn_enum'
			nrenum = tmp.shape[0]
		else:
			file_id = os.path.basename(filenames[i]).split('_')[0]+'_'+os.path.basename(filenames[i]).split('_')[3]
		index.append(file_id)
		list_csvs.append(tmp)
	combined_csv = pd.concat(list_csvs,ignore_index=False)
	if nrenum > 0:
		#Modify index after reaction_enum solutions
		line_count = pd.RangeIndex(0,combined_csv.shape[0],1)
		index_list = list(os.path.basename(filenames[0]).split('_')[0]+'_' + line_count.astype(str) + str(index_suffix))
		index_list[0:nrenum] = list(combined_csv[0:nrenum]['Solutions_IDS'])
		combined_csv.index = index_list
	else:
		line_count = pd.RangeIndex(0,combined_csv.shape[0],1)
		combined_csv.index = os.path.basename(filenames[0]).split('_')[0]+'_' + line_count.astype(str) + str(index_suffix)
	combined_csv.drop(combined_csv.columns[0],axis=1,inplace=True)
	combined_csv.drop_duplicates(inplace=True) #remove identical solutions
	if single_csv:
		combined_csv.to_csv(out_dir+'/all_solutions.csv', mode='a', encoding='utf-8-sig')
	else:
		combined_csv.to_csv(out_dir+'/'+os.path.basename(filenames[0]).split('_')[0]+'_solutions.csv', encoding='utf-8-sig')
	return


def remove_done_batchs(batch_dir,result_dir,launch_undone = True,relax_param = False,enum_type="reaction_enum", para_batch=False):
	"""remove_done_batchs.

	Parameters
	----------
	batch_dir : str
		path to the batchs directory
	result_dir : str
		path to the modelling result directory
	launch_undone : boolean
		If True, write the master bash file to launch all failed batchs
	relax_param : boolean
		If True, relax the mipgap tolerance parameter
	enum_type : str
		string indicating which type of enumeration is being processed (optional)
	para_batch : boolean
		if True, launch each batch file independantly (instead of parallel on conditions, parallel on batch)

	Returns
	-------
		a list with failed batch names

	"""
	removed_batchs = []
	results = glob.glob(result_dir+'/*solutions.csv')
	for file in results:
		#reconstruct the batch name, with a regex to be able to use it for reaction enum and diversity enum
		cleanfile = os.path.basename(file)
		item = str(cleanfile.split('_')[0])+'_'+str(cleanfile.split('_')[3])+'_.*_enum.sh'
		#look if item match with a batch file,(meaning that the batch is done)
		batchs = os.listdir(batch_dir)
		for batch in batchs:
			if re.search(item,batch):
				os.remove(batch_dir+'/'+batch)
				removed_batchs.append(batch)
	if relax_param == True:
		#listdir again because we do not want to iterate over removed batchs
		batchs = os.listdir(batch_dir)
		for batch in batchs:
			#read current batch
			with open(batch_dir+batch,'r') as f:
				content = f.read()
			#replace mipgap
			with open(batch_dir+batch,'w') as f:
				f.write(re.sub(r'--mipgap .*','--mipgap 0.01',content))
	if para_batch == True:
		batchs = os.listdir(batch_dir)
		for batch in batchs:
			#read current batch
			with open(batch_dir+batch,'r') as f:
				content = f.read()
			if '#!/bin/bash' in content:
				continue
			with open(batch_dir+batch,'w') as f:
				f.write('#!/bin/bash\n#SBATCH -p workq\n#SBATCH --mem=12G\n#SBATCH --cpus-per-task=12\n#SBATCH -t 48:00:00\n#SBATCH -J '+enum_type+'\n#SBATCH -o log_dir/runout_relaunch.out\n#SBATCH '
				'-e log_dir/runerr_relaunch.out\nsource activate cobrapy\n'+content)
		if launch_undone == True:
			with open(batch_dir.split('/')[0]+"/launch_failed_batch_"+enum_type+".sh", "w+") as f:
				f.write('#!/bin/bash\n#SBATCH -p workq\n#SBATCH --mem=12G\n#SBATCH --cpus-per-task=12\n#SBATCH -t 48:00:00\n#SBATCH -J '+enum_type+'\n#SBATCH -o log_dir/runout_relaunch.out\n#SBATCH '
				'-e log_dir/runerr_relaunch.out\nsource activate cobrapy\n ls '+batch_dir+'*enum.sh|xargs -n 1 -P 1 bash')
	return removed_batchs

def remove_zerobiomass_solutions(enum_dir,reaction_list,separator=','):
	"""remove_zerobiomass_solutions.

	Parameters
	----------
	enum_dir : str
		the path to the enumeration directory
	reaction_list : str
		the path to the reaction_list directory
	separator : str
		the character used to separate columns in the file

	Returns
	-------
		overwrite the csv file without solutions with 0 flux in biomass reaction

	"""
	col_list = ['Solutions_IDS'] + list(pd.read_csv(reaction_list).iloc[:,0])
	for file in os.listdir(enum_dir):
		tmp_file = pd.read_csv(enum_dir+'/'+file,sep=separator)
		tmp_file.columns = col_list
		#drop rows where biomass_reaction equals 0
		tmp_file.drop(tmp_file[tmp_file['biomass_reaction'] == 0].index,axis=0,inplace=True)
		#write modified solution file
		tmp_file.to_csv(enum_dir+'/'+file,sep=separator,index=False)
		