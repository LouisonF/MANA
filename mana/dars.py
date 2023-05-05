import math
import os
import pandas as pd
import numpy as np
import scipy.cluster.hierarchy as hc 

from scipy.spatial import distance_matrix

def calculate_frequencies(data, name):
	prob_vec = [sum(data.iloc[:, col])/data.shape[0] for col in range(data.shape[1])]
	return pd.Series(prob_vec, name=name)
def calculate_frequencies_for_dir(full_enum_path,rList,output_file= ""):
	csv_files = os.listdir(full_enum_path)
	prob_table = pd.DataFrame()
	prob_table = pd.concat([calculate_frequencies(pd.read_csv(full_enum_path+"/"+csv_file,index_col=0), csv_file.split('_')[0])
							for csv_file in csv_files], axis=1).transpose()
	prob_table.insert(0,"Barcode",prob_table.index)
	if "Barcode" not in rList:
		rList.insert(0,"Barcode")
	prob_table.columns = rList
	if len(output_file) > 0:
		prob_table.to_csv(output_file)
	return prob_table

def calculate_freq_ctrls(rListFile, all_cpds, time, pheno, working_path):
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
	"""
	reference: https://en.wikipedia.org/wiki/Rotation_matrix#In_two_dimensions
	:param vector: list of length 2 OR
				   list of list where inner list has size 2 OR
				   1D numpy array of length 2 OR
				   2D numpy array of size (number of points, 2)
	:param theta: rotation angle in radians (+ve value of anti-clockwise rotation)
	:param rotation_around: "vector" will be rotated around this point, 
					otherwise [0, 0] will be considered as rotation axis
	:return: rotated "vector" about "theta" degree around rotation
			 axis "rotation_around" numpy array
	"""
	vector = np.array(vector)

	if vector.ndim == 1:
		vector = vector[np.newaxis, :]

	if rotation_around is not None:
		vector = vector - rotation_around

	vector = vector.T

#     theta = np.radians(theta)

	rotation_matrix = np.array([
		[np.cos(theta), -np.sin(theta)],
		[np.sin(theta), np.cos(theta)]
	])

	output: np.ndarray = (rotation_matrix @ vector).T

	if rotation_around is not None:
		output = output + rotation_around

	return output.squeeze()

def rescale_and_rotate(comp_freq):
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

def extract_reactions_from_clusters(matrix,title,write_files=False,file_name='cluster',header=True):
	#Show dendro
	plt.figure(figsize=(20, 7))
	plt.title(title)
	colors = []
	#populate the colors_vector with black hexa code
	[colors.append("#000000") for x in range(0,matrix.shape[0]*matrix.shape[1])]
	# Create dendrogram
	linkage = hc.linkage(matrix, method='ward')
	dendro = hc.dendrogram(linkage,labels=list(matrix.index),link_color_func=lambda k: colors[k])

	plt.xlabel('Reaction ID')
	plt.ylabel('Euclidean distance')
	plt.show()
	#ask user for the number of clusters:
	nclusters = int(input("Enter the number of clusters to extract"))
	cutree = hc.cut_tree(linkage,n_clusters=nclusters).T
	labels = list(matrix.index)
	print(cutree)
	#for each cluster return reaction ids in a dict
	#init dict
	clusters_dict = {}
	for i in range(0,nclusters):
		#for each cluster, define
		print([labels[x] for x in np.where(cutree == i)[1]])
		clusters_dict["cluster"+str(i)] = [labels[x] for x in np.where(cutree == i)[1]]
	if write_files:
		j=1
		for key in clusters_dict.keys():
			with open(file_name+str(j)+".tab",'w') as w_hdler:
				if header:
					w_hdler.write("Reaction ID \n")
				for elem in clusters_dict[key]:
					w_hdler.write(elem.strip('\t')+'\n')
			j=j+1
	return cutree

def compute_scores(comp_freq,crossing_point=1,crossing_point_1_2=1.2,b=1):
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