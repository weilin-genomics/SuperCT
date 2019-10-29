#!/usr/bin/env python
# coding: utf-8
import sys,getopt

usage='''Usage: python pred_types.py -i <input_dge.csv> -o <cell_type_output.csv> -p -m <model.h5> -s <human or mouse> -t <id2type.csv>
	or: python pred_types.py --input_file=<input_dge.csv> --output_file=<cell_type_output.csv> --plot_hist --model=<model.h5> --species=<human or mouse> --types=<id2type.csv>
	'''
try:
	opts,args = getopt.getopt(sys.argv[1:], "hpi:o:m:s:t:", ["help",'plot','input_file=','output_file=','model=','species=','types='])
except getopt.GetoptError:
	print('Error\n',usage)
	sys.exit(2)
output_file='cell_type_output.csv'
input_file=''
model_file=''
plot_hist=False
species='human'
types=''
for opt,arg in opts:
	if opt in ("-h", "--help"):
		print(usage)
		sys.exit()
	elif opt in ("-i", "--input_file"):
		input_file = arg
	elif opt in ("-o", "--output_file"):
		output_file = arg
	elif opt in ("-m", "--model"):
		model_file = arg
	elif opt in ('-p','--plot_hist'):
		plot_hist=True
	elif opt in ('-s','--species'):
		species = arg
	elif opt in ('-t','--types'):
		types=arg
if len(input_file)==0:
	print('expression file must be supplied!')
	print(usage)
	sys.exit()
if len(model_file)==0:
	print('model file must be supplied!')
	print(usage)
	sys.exit()
#make sure these module installed
import numpy as np
import pandas as pd
from keras.models import load_model
import matplotlib.pyplot as plt

def load_genes(genes_file,species):
	if species=='human':
		genes_used = pd.read_csv(genes_file,usecols=[0],header = None)
	elif species=='mouse':
		genes_used = pd.read_csv(genes_file,usecols=[1],header = None)
	genes_used = genes_used.values.tolist()
	for i,g in enumerate(genes_used):
		genes_used[i]=g[0]
	return genes_used
def load_data(file_name,genes_used):
	data = pd.read_csv(file_name, header = 0, index_col = 0)
	data = data[~data.index.duplicated(keep='first')]
	genes_valid=np.intersect1d(data.index,genes_used)
	if len(genes_valid)==len(genes_used):
		data=data.loc[genes_used,:]
	elif len(genes_valid) < len(genes_used):
		genes_not_found=np.setdiff1d(genes_used,genes_valid)
		nonezero_df=pd.DataFrame(np.zeros((len(genes_not_found),data.shape[1])).astype('int64'))
		nonezero_df.columns=data.columns
		data=pd.concat([data,nonezero_df])
		data=data.loc[genes_used,:]
		print('Warning:',genes_not_found,'not found in your expression file,please check them.but the prediction will go on by setting these genes with zero count.')
	return data
def load_cell_types(types):
	cell_types={}
	all_types=pd.read_csv(types)
	for i,each in enumerate(all_types['celltype']):
		cell_types[i]=each
	return cell_types
def predict_type(model,data,cell_types):
	data = data > 0
	data = data.astype(np.uint8)
	result = model.predict(data.values.T)
	result_df = pd.DataFrame(result)
	result_df.columns = list(cell_types.values())
	cell_type_pred = [cell_types[i] for i in list(np.argmax(result, axis = 1))]
	result_df["pred_type"] = cell_type_pred
	return result_df

model=load_model(model_file)
genes_used=load_genes('inthomgenes.csv',species)
data=load_data(input_file,genes_used)
cell_types=load_cell_types(types)

result_df=predict_type(model,data,cell_types)
result_out=pd.DataFrame({'cell_id': data.columns,'pred_type':result_df['pred_type']})
result_out.to_csv(output_file,index=False)

if plot_hist:
	d={}
	res_list=list(result_out['pred_type'])
	for i in res_list:
		d.update({i:res_list.count(i)})
	plt.bar(d.keys(),d.values())
	plt.title('Frequency of predicted cell types')
	plt.xticks(rotation=90)
	plt.savefig('predicted_types.pdf',bbox_inches='tight')
