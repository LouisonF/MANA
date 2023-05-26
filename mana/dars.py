import math
import os
import pandas as pd
import numpy as np
import scipy.cluster.hierarchy as hc 

from scipy.spatial import distance_matrix

def calculate_frequencies(data, name):
	"""calculate_frequencies.

	Parameters
	----------
	data : pandas dataframe
		a dataframe containing the partial enumeration results
	name : str
		the name of the pandas Series

	Returns
	-------
		a pandas series with the sum of partial enumeration results for each reaction.

	"""
	prob_vec = [sum(data.iloc[:, col])/data.shape[0] for col in range(data.shape[1])]
	return pd.Series(prob_vec, name=name)

def calculate_frequencies_for_dir(full_enum_path,rList,output_file= ""):
	"""calculate_frequencies_for_dir.

	Parameters
	----------
	full_enum_path : str
		the path to the full_enum directory
	rList : list
		the list of reactions in the model
	output_file : str
		if not empty, the frequencies table is written in csv to the path provided in this param

	Returns
	-------
		a pandas dataframe containing the frequencies table, of each csv file found in the full enum path.

	"""
	csv_files = os.listdir(full_enum_path)
	freq_table = pd.DataFrame()
	freq_table = pd.concat([calculate_frequencies(pd.read_csv(full_enum_path+"/"+csv_file,index_col=0), csv_file.split('_')[0])
							for csv_file in csv_files], axis=1).transpose()
	freq_table.insert(0,"Barcode",freq_table.index)
	if "Barcode" not in rList:
		rList.insert(0,"Barcode")
	freq_table.columns = rList
	if len(output_file) > 0:
		freq_table.to_csv(output_file)
	return freq_table

def calculate_freq_ctrls(rListFile, all_cpds, time, pheno, working_path):
	"""calculate_freq_ctrls.

	Parameters
	----------
	rListFile : str
		the path to the model's reactions list file
	all_cpds : str
		the list of compounds to process, from the props.properties file
	time : str 
		the exposure time to consider
	pheno : pandas dataframe
		the pandas data frame containing Open TG-Gates metadata.
	working_path : str
		the root path of your working directory.

	Returns
	-------
		return a pandas dataframe with the calculated frequencies 
		for all controls partial enumeration results corresponding to input parameters

	"""
	freq_ctrls = pd.DataFrame()
	rList = list(pd.read_csv(rListFile).iloc[:,0])
	cond_list = []
	for cpd in str(all_cpds).split('/'):
		# pool enumerated solutions for each compound's control at specified time
		tmp_fullenum = pd.DataFrame()
		barcodes = list(pheno[(pheno["compound_name"] == cpd) & (pheno["sacri_period"] == time) & \
			(pheno["dose_level"] == "Control")].index)
		for barcode in barcodes:
			try:
				tmp_fullenum = pd.concat([tmp_fullenum, \
					pd.read_csv(working_path+"/control_"+str(time).replace(" ","_")+"/full_enum/"+barcode+\
						"_solutions.csv",index_col=0)])
			except FileNotFoundError:
				print("File not found for "+barcode)
		freq_ctrls = pd.concat([freq_ctrls, calculate_frequencies(tmp_fullenum, cpd)], axis=1)
		cond_list.append(cpd+"_control_"+str(time).replace(" ","_"))
		freq_ctrls.columns = cond_list
	freq_ctrls.index = rList          
	freq_ctrls = freq_ctrls.transpose()
	return freq_ctrls

def rotate(vector, theta, rotation_around=None):
	"""rotate.

	reference: https://en.wikipedia.org/wiki/Rotation_matrix#In_two_dimensions

	Parameters
	----------
	vector : pandas dataframe
		the activation frequencies dataframe to rotate
	theta : float
		rotation angle in radians
	rotation_around : np.array
		A point around which vector will be rotated around. Can be None

	Returns
	-------
		The rotated dataframe

	"""
	vector = np.array(vector)

	if vector.ndim == 1:
		vector = vector[np.newaxis, :]

	if rotation_around is not None:
		vector = vector - rotation_around

	vector = vector.T

	rotation_matrix = np.array([
		[np.cos(theta), -np.sin(theta)],
		[np.sin(theta), np.cos(theta)]
	])

	output: np.ndarray = (rotation_matrix @ vector).T

	if rotation_around is not None:
		output = output + rotation_around

	return output.squeeze()

def rescale_and_rotate(comp_freq):
	"""rescale_and_rotate.

	Parameters
	----------
	comp_freq : pandas dataframe
		a pandas dataframe containing several metrics computed from activation frequencies

	Returns
	-------
		a pandas dataframe with f_ctrl and f_trt  that have been rescaled, rescaled and rotated

	"""
	rids = list(comp_freq.index)
	#rescale f_ctrl and f_treatment
	comp_freq = comp_freq.assign(
		f_ctrl_rescale = (comp_freq['f_ctrl'] - 0.5) * 2,
		f_treatment_rescale = (comp_freq['f_treatment'] - 0.5) * 2
	)
	#Rotation of f_ctrl and f_treatment
	rotated_f = pd.DataFrame(rotate(comp_freq.iloc[:,1:3],math.pi/4), columns=['f_ctrl_r','f_treatment_r'], index = comp_freq.index)
	comp_freq = pd.concat([comp_freq,rotated_f],axis=1)

	#Rotation of f_ctrl_rescale and f_treatment_rescale
	rotated_frescale = pd.DataFrame(rotate(comp_freq.iloc[:,3:5],math.pi/4), columns=['f_ctrl_rescale_r','f_treatment_rescale_r'], index = comp_freq.index)
	comp_freq = pd.concat([comp_freq,rotated_frescale],axis=1)

	return comp_freq

def findCircleCenter(A,B,C):
	"""findCircleCenter.

	Parameters
	----------
	A : list
		a list containing x and y coordinates of point A, a point of the circle O
	B : list
		a list containing x and y coordinates of point B, a point of the circle O
	C : list
		a list containing x and y coordinates of point C, a point of the circle O

	Returns
	-------
		a pandas dataframe containing :calculated properties of the circle O

	"""
	Ax = A[0] #x1
	Ay = A[1] #y1
	Bx = B[0] #x2
	By = B[1] #y2
	Cx = C[0] #x3
	Cy = C[1] #y3
	xAB = Ax - Bx #x12
	xAC = Ax - Cx #x13
	yAB = Ay - By #y12
	yAC = Ay - Cy #y13
	xCA = Cx - Ax #x31
	xBA = Bx - Ax #x21
	yCA = Cy - Ay #y31
	yBA = By - Ay #y21

	sxAC = (Ax ** 2) - (Cx ** 2) #sx13
	syAC = (Ay ** 2) - (Cy ** 2) #sy13
	sxBA = (Bx ** 2) - (Ax ** 2) #sx21
	syBA = (By ** 2) - (Ay ** 2) #sy21
	try:
		f = ((sxAC) * (xAB)
		+ (syAC) * (xAB)
		+ (sxBA) * (xAC)
		+ (syBA) * (xAC)) / \
		(2 * ((yCA) * (xAB) - (yBA) * (xAC)))

		g = ((sxAC) * (yAB)
		+ (syAC) * (yAB)
		+ (sxBA) * (yAC)
		+ (syBA) * (yAC)) / \
		(2 * ((xCA) * (yAB) - (xBA) * (yAC)))  
	except ZeroDivisionError:
		#ZeroDivision Error raised, return nan values
		return pd.DataFrame([[np.nan,np.nan,np.nan, np.nan]], columns = ['x','y','d','dist_to_OO'])
	c = -1*(Ax ** 2) - (Ay ** 2) - 2 * g * Ax - 2 * f *Ay
	sqr_of_r = -g * (-g) + -f * (-f) -c
	r  = math.sqrt(sqr_of_r)

	df = pd.DataFrame(list(zip([0,-g],[0,-f])))
	dist_to_OO = distance_matrix(df.values,df.values).max()
	return pd.DataFrame([[-g,-f,r*2, dist_to_OO]], columns = ['x','y','d','dist_to_OO'])

def compute_scores(comp_freq,crossing_point=1,crossing_point_1_2=1.2,b=1):
	"""compute_scores.

	Parameters
	----------
	comp_freq : pandas dataframe
		a pandas dataframe containing several metrics computed from activation frequencies
	crossing_point : int
		a crossing point factor for the circle center calculation
	crossing_point_1_2 : float
		another crossing point factor the circle center calculation
	b : int
		a parameter for the calculation of ellipses 
	Returns
	-------
		a pandas dataframe with several metrics computed: R2, center of circle, dist_to_OO, center of circle 1.2

	"""
	#instantiate theta
	theta = math.pi/4
	#create the data_id col with Reactions ids
	comp_freq = comp_freq.assign(
		data_id=comp_freq.index
		).assign(
			R2 = lambda x: (comp_freq['f_ctrl'] - comp_freq['f_treatment']) **2
		).assign(
			X = comp_freq['f_ctrl_rescale'] * np.cos(theta) + comp_freq['f_treatment_rescale'] * np.sin(theta)
		).assign(
			Y = comp_freq['f_treatment_rescale'] * np.cos(theta) - comp_freq['f_ctrl_rescale'] * np.cos(theta)
		)
	comp_freq = comp_freq.assign(
			ellipse = (comp_freq['X']/math.pi) ** 2 + (comp_freq['Y']/b) ** 2
		)
	circle_center_crossing_point = pd.DataFrame()
	for i in range(0,comp_freq.shape[0]):
		circle_center_crossing_point = pd.concat([circle_center_crossing_point,findCircleCenter(A=[-1*crossing_point,-1*crossing_point],
										B=[crossing_point,crossing_point],
										C=list(comp_freq.iloc[i,3:5]))])
	#rename cols and rows
	circle_center_crossing_point.index = comp_freq.index
	circle_center_crossing_point.columns = ['x','y','center_of_circle','dist_to_OO']
	comp_freq = pd.concat([comp_freq,circle_center_crossing_point.loc[:,['center_of_circle','dist_to_OO']]],axis=1)

	comp_freq = comp_freq.assign(
			center_of_circle_sqrt = np.sqrt(comp_freq['center_of_circle'])
		).assign(
			center_of_circle_log = np.log(comp_freq['center_of_circle'])
		)
	circle_center_crossing_point_1_2 = pd.DataFrame()
	for i in range(0,comp_freq.shape[0]):
		circle_center_crossing_point_1_2 = pd.concat([circle_center_crossing_point_1_2,findCircleCenter(A=[-1*crossing_point_1_2,-1*crossing_point_1_2],
										B=[crossing_point_1_2,crossing_point_1_2],
										C=list(comp_freq.iloc[i,3:5]))])
	#rename cols and rows
	circle_center_crossing_point_1_2.index = comp_freq.index
	circle_center_crossing_point_1_2.columns = ['x','y','center_of_circle_1_2','dist_to_OO_1_2']
	comp_freq = pd.concat([comp_freq,circle_center_crossing_point_1_2.loc[:,['center_of_circle_1_2','dist_to_OO_1_2']]],axis=1)
	comp_freq = comp_freq.assign(
			center_of_circle_1_2_sqrt = np.sqrt(comp_freq['center_of_circle_1_2'])
		).assign(
			center_of_circle_1_2_log = np.log(comp_freq['center_of_circle_1_2'])
		)
	return comp_freq