
#imports
import time
import pickle
import pandas as pd
from tqdm.auto import tqdm
from multiprocessing import Process

def worker(q,_finish):
	"""worker.

	Parameters
	----------
	q : JoinableQueue
		a JoinableQueue object filled with tasks to perform
	_finish : boolean
		a boolean value indicating if the worker should process elements from the JoinableQueue

	Returns
	-------
		None

	"""
	while _finish == False:
		job,args = q.get()
		job(*args)
		q.task_done()

def launch_multi_proc(num_workers,q):
	"""launch_multi_proc.

	Parameters
	----------
	num_workers : int
		the number of workers which will define the number of allowed parallel threads
	q : JoinableQueue
		a JoinableQueue object filled with tasks to perform

	Returns
	-------
		None

	"""
	procs = []
	q_size_tot = q.qsize()
	pbar = tqdm(total=q_size_tot)
	up_n_prev = q_size_tot-q._unfinished_tasks.get_value()
	up_n = 0
	_finish = False
	for i in range(num_workers):
		p = Process(target=worker, args=(q,_finish))
		p.daemon = True
		p.start()
		procs.append(p)
		time.sleep(0.1)
	_finish = True
	while q_size_tot-up_n != 0:
		up_n = q_size_tot-q._unfinished_tasks.get_value()
		if up_n != up_n_prev:
			pbar.update(up_n-up_n_prev)
		up_n_prev = up_n
		time.sleep(2)
	q.join()       # block until all tasks are done
	print ("End of Queue")
	#When queue is empty, terminate all alive workers
	for p in procs:
		p.terminate()
		p.join()

def make_pickle(object,filename):
	"""make_pickle.

	Parameters
	----------
	object : pkl 
		the pkl object to save 
	filename : str
		the path and filename for the pickle file
	notebook : boolean
		if true enable the pyvis optimized visualisation for notebooks

	Returns
	-------
		write a .pkl file at the designated location

	"""
	with open(filename, 'wb') as f:
		pickle.dump(object, f)

def make_csvs(data,out_folder,celfilename):
	"""make_csvs.

	Parameters
	----------
	data : list
		list of lists with data to be transformed into a pandas dataframe
	out_folder : str
		the path where csvs files will be saved
	celfilename : str
		the initial CEL filename (used as identifier)

	Returns
	-------
		an interactive visualisation window of the graph

	"""
	rh_df = pd.DataFrame(data[0],columns=['reactions'])
	rh_df['weights'] = 1
	rl_df = pd.DataFrame(data[1],columns=['reactions'])
	rl_df['weights'] = -1
	rh_rl_df = pd.concat([rh_df,rl_df])
	rh_rl_df.to_csv(out_folder+'bin_reactions_'+celfilename+'.csv', index=False)     