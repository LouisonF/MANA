
import random
import os
import pandas as pd

def write_div_enum_script(script_path,batch_directory, rxn_enum_set_dir,output_directory, modelfile, weightfile,\
                           reactionFile, prev_sol_dir ='prev_sol_dir/', log_dir='log_dir',dist_anneal=0.9, obj_tol=0.01,\
                              iters=100,para_batchs=False):
    """write_div_enum_script.

	Parameters
	----------
	script_path : str
		path to the diversity_enum.py dexom python script
	batch_directory : str
		path to the directory were batch files should be written
	rxn_enum_set_dir : str
		path to the directory of processed reaction-enum results
	output_directory : str
		path to the directory were diversity-enum modelling results should be written
    modelfile : str
		path to the model's json file
	weightfile : str
		path to the csvs file that contains binarized reactions activity (according to transcriptomic data)
	reactionFile : str
		path to the file that contains the list of reactions in the model
    prev_sol_dir : str
		path to the directory were reaction-enum solutions used as starting point for the diversity enumeration
        process should be saved
	log_dir : str
		path to the directory were log files should be stored
	dist_anneal : float
		dexom-python parameter, 0<=a<=1 controls the distance between each successive solution
    obj_tol : float
		dexom-python parameter, objective value tolerance, as a fraction of the original value
	iters : int
		dexom-python parameter, maximal number of iterations
	para_batchs : boolean
		if True, launch each batch file independantly (instead of parallel on conditions, parallel on batch)
	Returns
	-------
		write batch files ready to launch on a adequatly prepared slurm computing platform

	"""
    #generate as many batch as ranges step for reaction enum
    barcode = os.path.basename(weightfile).split('_')[2]
    with open(reactionFile, "r") as file:
        rxns = file.read().split("\n")
    split_val = (len(rxns) // iters) + 1
    #load the reaction enum set file for this barcode
    enum_set = pd.read_csv(rxn_enum_set_dir+'/'+barcode+'_solutions.csv')
    enum_set.index = enum_set.iloc[:,0]
    enum_set.index.name = 'ids'
    enum_set.drop(enum_set.columns[0],axis=1,inplace=True)
    #Generate stratified random sampling
    nbatch = 0
    prev_lb = 0 #previous number where i % split_val == 0
    for i in range(enum_set.shape[0]):
        if i % split_val == 0:
            #random pick a solution in the range
            tmp_sol = enum_set.iloc[random.randint(prev_lb,i),:]
            prev_lb = i
            prevsol_file = prev_sol_dir+barcode+'_'+str(i)+'.csv'
            pd.DataFrame(tmp_sol).transpose().to_csv(prevsol_file)
            #Now we generate the batch script:
            if para_batchs:
                with open(batch_directory+'/batch/'+barcode+ '_' + str(i) + "_diversity_enum.sh", "w+") as f:
                    f.write('#!/bin/bash\n#SBATCH -p workq\n#SBATCH --mem=12G\n#SBATCH --cpus-per-task=12\n#SBATCH -t 72:00:00\n#SBATCH -J div_enum\n#SBATCH -o %s/runout%s_div.out\n#SBATCH '
                        '-e %s/runerr%s_div.out\nsource activate cobrapy \n'
                        % (str(log_dir),str(barcode),str(log_dir),str(barcode)))
                with open(batch_directory+'/batch/'+barcode+ '_' + str(i) + "_diversity_enum.sh", "dist_anneal") as f:
                    f.write('python %s -o %s/%s_div_enum_%i -m %s -r %s -p %s -a %.5f -i %i --obj_tol %.4f'
                        % (script_path,output_directory, barcode, i, modelfile, weightfile, prevsol_file, dist_anneal, iters, obj_tol))
            else:
                with open(batch_directory+'/batch/'+barcode+ '_' + str(i) + "_diversity_enum.sh", "w+") as f:
                    f.write('python %s -o %s/%s_div_enum_%i -m %s -r %s -p %s -a %.5f -i %i --obj_tol %.4f'
                        % (script_path,output_directory, barcode, i, modelfile, weightfile, prevsol_file, dist_anneal, iters, obj_tol))
            nbatch=nbatch+1
    if para_batchs == False:
        with open(batch_directory+"/runfiles_"+barcode+"_diversity_enum.sh", "w+") as f:
            f.write('#!/bin/bash\n#SBATCH -p workq\n#SBATCH --mem=12G\n#SBATCH --cpus-per-task=12\n#SBATCH -t 72:00:00\n#SBATCH -J div_enum\n#SBATCH -o %s/runout%s_div.out\n#SBATCH '
                    '-e %s/runerr%s_div.out\nsource activate cobrapy\nls %s/batch/%s_*_diversity_enum.sh|xargs -n 1 -P 1 bash'
                    % (str(log_dir),str(barcode),str(log_dir),str(barcode),str(batch_directory),str(barcode)))

def write_rxn_enum_script(script_path,batch_directory,output_directory, modelfile, weightfile,\
                           reactionFile="", log_dir='log_dir',obj_tol=0.001, iters=100,para_batchs=False):
    """write_rxn_enum_script.

	Parameters
	----------
	script_path : str
		path to the diversity_enum.py dexom python script
	batch_directory : str
		path to the directory were batch files should be written
	output_directory : str
		path to the directory were diversity-enum modelling results should be written
    modelfile : str
		path to the model's json file
	weightfile : str
		path to the csvs file that contains binarized reactions activity (according to transcriptomic data)
	reactionFile : str
		path to the file that contains the list of reactions in the model
	log_dir : str
		path to the directory were log files should be stored
    obj_tol : float
		dexom-python parameter, objective value tolerance, as a fraction of the original value
	iters : int
		dexom-python parameter, maximal number of iterations
	para_batchs : boolean
		if True, launch each batch file independantly (instead of parallel on conditions, parallel on batch)
	Returns
	-------
		write batch files ready to launch on a adequatly prepared slurm computing platform

	"""
    barcode = os.path.basename(weightfile).split('_')[2]
    with open(reactionFile, "r") as file:
        rxns = file.read().split("\n")
    rxn_num = (len(rxns) // iters) + 1
    if para_batchs:
        for i in range(rxn_num):
            with open(batch_directory+'/batch/'+barcode+ '_' + str(i) + "_reaction_enum.sh", "w+") as f:
                f.write('#!/bin/bash\n#SBATCH -p workq\n#SBATCH --mem=12G\n#SBATCH --cpus-per-task=12\n#SBATCH -t 24:00:00\n#SBATCH -J rxn_enum\n#SBATCH -o %s/runout%s_div.out\n#SBATCH '
                        '-e %s/runerr%s_div.out\nsource activate cobrapy \n'
                        % (str(log_dir),str(barcode),str(log_dir),str(barcode)))
            with open(batch_directory+'/batch/'+barcode+ '_' + str(i) + "_reaction_enum.sh", "a") as f:
                f.write('python %s -o %s/%s_rxn_enum_%i --range %i_%i -m %s -r %s -l %s '
                        '-t 600 --mipgap %f \n' % (script_path,output_directory,barcode, i, i*iters, i*iters+iters, modelfile, weightfile, reactionFile, obj_tol))
    else:
        for i in range(rxn_num):
            with open(batch_directory+'/batch/'+barcode+ '_' + str(i) + "_reaction_enum.sh", "w+") as f:
                f.write('python %s -o %s/%s_rxn_enum_%i --range %i_%i -m %s -r %s -l %s '
                        '-t 600 --mipgap %f \n' % (script_path,output_directory,barcode, i, i*iters, i*iters+iters, modelfile, weightfile, reactionFile, obj_tol))
            with open(batch_directory+"/runfiles_"+barcode+"_reaction_enum.sh", "w+") as f:
                f.write('#!/bin/bash\n#SBATCH -p workq\n#SBATCH --mem=12G\n#SBATCH --cpus-per-task=12\n#SBATCH -t 24:00:00\n#SBATCH -J rxn_enum\n#SBATCH -o %s/runout%s_div.out\n#SBATCH '
                        '-e %s/runerr%s_div.out\nsource activate cobrapy\nls %s/batch/%s_{0..%i}_reaction_enum.sh|xargs -n 1 -P 1 bash'
                         % (str(log_dir),str(barcode),str(log_dir),str(barcode),str(batch_directory),str(barcode), int(rxn_num-1)))