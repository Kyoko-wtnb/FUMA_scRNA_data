#!/usr/bin/python
### USAGE: python loom2txt.py <input loom file> <output prefix>
### Output:	Tab delimited plain text files
###			<prefix>_col_attrs.txt
###			<prefix>_row_attrs.txt
###			<prefix>_matrix.txt

import numpy as np
import h5py
import sys

def main():
	if len(sys.argv) < 3:
		print "ERROR: not enough arguments.\nUSAGE: python loom2txt.py <input loom file> <output prefix>"
		sys.exit()
	fin = sys.argv[1]
	fout = sys.argv[2]
	f = h5py.File(fin, 'r')

	### col_attrs
	g = "col_attrs"
	k = f[g].keys()
	cols = []
	for i in k:
		if len(cols)==0:
			cols = f[g][i][:]
		else:
			cols = np.c_[cols, f[g][i][:]]
	with open(fout+"_col_attrs.txt", 'w') as o:
		o.write("\t".join(k)+"\n")
	with open(fout+"_col_attrs.txt", 'a') as o:
		np.savetxt(o, cols, fmt="%s", delimiter="\t")

	### row_attrs
	g = "row_attrs"
	k = f[g].keys()
	rows = []
	for i in k:
		if len(rows)==0:
			rows = f[g][i][:]
		else:
			rows = np.c_[rows, f[g][i][:]]
	with open(fout+"_row_attrs.txt", 'w') as o:
		o.write("\t".join(k)+"\n")
	with open(fout+"_row_attrs.txt", 'a') as o:
		np.savetxt(o, rows, fmt="%s", delimiter="\t")

	### data
	with open(fout+"_matrix.txt", 'a') as o:
		np.savetxt(o, f["matrix"][:,:], fmt="%s", delimiter="\t")

if __name__=="__main__": main()
