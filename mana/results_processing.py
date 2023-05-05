import os
import re
import glob
import progressbar
import pandas as pd
import numpy as np

from multiprocessing import JoinableQueue
from .utils import launch_multi_proc

def concatenate_solutions(csv_dir,out_dir,col_index="",single_csv=False,ncpus=1,
 restart=False,combine_r_d_enum=False):
    #source : https://www.freecodecamp.org/news/how-to-combine-multiple-
    #csv-files-with-8-lines-of-code-265183e0854/
    if combine_r_d_enum:
        csv_dir = 'working_divers'
    os.chdir(csv_dir)
    if csv_dir == 'working_full' or csv_dir == 'working_full/':
        list_done_set = os.listdir('../full_rxn_enum_set')
        index_suffix = '_renum'
    elif csv_dir == 'working_divers' or csv_dir == 'working_divers/':
        index_suffix = '_rdivers'
        if combine_r_d_enum:
            list_done_set = os.listdir('../full_rxn_enum_set')
        else:
            list_done_set = os.listdir('../full_div_enum_set')
    else:
        return "Wrong csv_dir name"
    q = JoinableQueue()
    pbar = progressbar.ProgressBar()
    list_dir = list(set([i.split('_')[0] for i in os.listdir()]))
    for i in range(0,len(list_dir)):
        if restart:
            #search if a concatenated csv exist for this file id in the full rxn set
            rgx = re.compile('.*'+list_dir[i].split('_')[0]+'.*')
            match_res = list(filter(rgx.match,list_done_set))
            if len(match_res) == 0:
                files = glob.glob('*'+list_dir[i].split('_')[0]+'*_solutions.csv')
                if ncpus == 1:
                    concatenate_csv(files,out_dir,col_index,single_csv,index_suffix)
                else:
                    q.put((concatenate_csv,(files,out_dir,col_index,single_csv,index_suffix)))
        else:
            files = glob.glob('*'+list_dir[i].split('_')[0]+'*_solutions.csv')
            if combine_r_d_enum:
                rgx = re.compile('.*'+list_dir[i].split('_')[0]+'.*')
                match_res = list(filter(rgx.match,list_done_set))
                match_res[0] = '../full_rxn_enum_set/'+match_res[0]
                files = match_res + files
            if ncpus == 1:
                concatenate_csv(files,out_dir,col_index,single_csv,index_suffix)
            else:
                q.put((concatenate_csv,(files,out_dir,col_index,single_csv,index_suffix)))
    if ncpus == 1:
        pass
    else:
        launch_multi_proc(ncpus,q)
    os.chdir('..')
    return q

def concatenate_csv(filenames,out_dir,col_index,single_csv,index_suffix=""):
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
        index_list = list(os.path.basename(filenames[0]).split('_')[0]+'_' + combined_csv.index.astype(str) + str(index_suffix))
        index_list[0:nrenum] = list(combined_csv[0:nrenum]['Solutions_IDS'])
        combined_csv.index = index_list
    else:
        combined_csv.index = os.path.basename(filenames[0]).split('_')[0]+'_' + combined_csv.index.astype(str) + str(index_suffix)
    combined_csv.drop(combined_csv.columns[0],axis=1,inplace=True)
    combined_csv.drop_duplicates(inplace=True) #remove identical solutions
    # if len(col_index) > 0:
    # 	combined_csv.columns = list(col_index)
    # else:
    # 	#If no col_index, use col_index of the last file
    # 	combined_csv.columns = tmp.drop('Unnamed: 0',axis=1).columns
    if single_csv:
        combined_csv.to_csv('../'+out_dir+'/all_solutions.csv', mode='a', encoding='utf-8-sig')
    else:
        combined_csv.to_csv('../'+out_dir+'/'+os.path.basename(filenames[0]).split('_')[0]+'_solutions.csv', encoding='utf-8-sig')
    return

def remove_done_batchs(batch_dir,result_dir,launch_undone = True,relax_param = False,enum_type="reaction_enum", para_batch=False):
    removed_batchs = []
    batchs = os.listdir(batch_dir)
    results = glob.glob(result_dir+'/*solutions.csv')
    for file in results:
        #reconstruct the batch name, with a regex to be able to use it for reaction enum and diversity enum
        cleanfile = os.path.basename(file)
        item = str(cleanfile.split('_')[0])+'_'+str(cleanfile.split('_')[3])+'_.*_enum.sh'
        #look if item match with a batch file,(meaning that the batch is done)
        for batch in batchs:
            if re.search(item,batch):
                os.remove(batch_dir+'/'+batch)
                removed_batchs.append(batch)
#             else:
#                 print("Missing result file:\t"+batch)
#                 print("To relaunch:\t"+batch)
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
	col_list = ['Solutions_IDS'] + list(pd.read_csv(reaction_list).iloc[:,0])
	for file in os.listdir(enum_dir):
		tmp_file = pd.read_csv(enum_dir+'/'+file,sep=separator)
		tmp_file.columns = col_list
		#drop rows where biomass_reaction equals 0
		tmp_file.drop(tmp_file[tmp_file['biomass_reaction'] == 0].index,axis=0,inplace=True)
		#write modified solution file
		tmp_file.to_csv(enum_dir+'/'+file,sep=separator,index=False)
        