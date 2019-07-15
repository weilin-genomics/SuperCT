import sys,getopt
usage='Usage: python transform_10x_to_dge.py -m matrix.mtx -g genes.tsv -b barcodes.tsv -o transformed_dge.csv'
try:
	opts,args = getopt.getopt(sys.argv[1:], "hm:g:b:o:")
except getopt.GetoptError:
	print('Error\n',usage)
	sys.exit(2)

out_file='transformed_dge.csv'
for opt,arg in opts:
	if opt in ("-h"):
		print(usage)
		sys.exit()
	if opt in ('-m'):
		matrix_file=arg
	if opt in ('-g'):
		genes_file=arg
	if opt in ('-b'):
		barc_file=arg
	if opt in ('-o'):
		out_file=arg

from scipy.io import mmread
import numpy as np
import pandas as pd

data=mmread(matrix_file)
data=pd.DataFrame(data.todense())

genes=pd.read_csv(genes_file,header=None,usecols=[0],sep='\t')
genes_list=[each[0] for i,each in enumerate(genes.values.tolist())]
barc=pd.read_csv(barc_file,header=None,usecols=[0])
barc_list=[each[0] for i,each in enumerate(barc.values.tolist())]
data.index=genes_list
data.columns=barc_list
data.to_csv(out_file)

print("transforming expression profile...")
print("your data shape:",data.shape)
print('Done!data saved as ',out_file,'.')

