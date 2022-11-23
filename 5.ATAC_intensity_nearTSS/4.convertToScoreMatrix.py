## Written by KyungTae Lee
## Last update : 2022-02-01
## Python script : "4.convertToScoreMatrix.py"
## Script to convert baitch-correcdted ATAC-seq read counts to file of different format that can be used as input for R genomation package

#!/usr/bin/env python

import sys
import pandas as pd

def getScoreD(br_matrix):
## batch corrected probe - sample matrix file
	scoreD= dict()
	rownumD= dict()
	colnumD= dict()
	with open(br_matrix, "r") as matrix_open:
#		matrixlineL= matrix_open.readlines()
		i_line= matrix_open.readline()
		header= i_line
		colL= header.strip().split("\t")
		sampleL= colL[1:]
	
		while i_line!= "\n" and i_line!="":
			i_line= matrix_open.readline()
			if i_line!= "\n" and i_line!= "":
		
				infoL= i_line.strip().split("\t")
				probeid= infoL[0]
				rownum, region_num= probeid.split("_")
				rownum= rownum.lstrip('"').rstrip('"')
				region_num= region_num.lstrip('"').rstrip('"')
	
				for i_index in xrange(len(sampleL)):
					i_sample= sampleL[i_index]
					i_sample= i_sample.lstrip('"').rstrip('"')
					i_count= infoL[i_index+1]
					if not scoreD.has_key(i_sample):
						scoreD[i_sample]= dict()
					if not scoreD[i_sample].has_key(rownum):
						scoreD[i_sample][rownum]= dict()
					if not scoreD[i_sample][rownum].has_key(region_num):
						scoreD[i_sample][rownum][region_num]= i_count
			else:
				pass

	return scoreD

def getRowColNum(original_sm):
	sm_df = pd.read_csv(original_sm, sep="\t")
	rownum= len(sm_df)
	colnum= len(sm_df.columns)-1
	return rownum, colnum

def writeScoreMatrix(scoreD, colnum, rownum,output_dir):
	sampleL= scoreD.keys()
	colnumL= range(1,colnum+1)
	regionL= map(lambda x: "V"+ str(x), colnumL)
	rownumL= map(lambda x: str(x), range(1,rownum+1))
	for i_sample in sampleL:
		outputname= output_dir+ "%s.batchCorrect.scoreMatrix.txt"%i_sample
		with open(outputname, "w") as outputopen:
			header= "generow\t"+ "\t".join(regionL)+ "\n"
			outputopen.write(header)
			for i_row in rownumL:
				countL= list()
				for i_region in regionL:
					try:
						i_count= scoreD[i_sample][i_row][i_region]
						countL.append(i_count)
					except KeyError:
						continue
				if countL==[]:
					continue
				outputline= i_row+ "\t"+ "\t".join(countL)+ "\n"
				outputopen.write(outputline)

def main():
	br_matrix= sys.argv[1]
	original_sm= sys.argv[2]
	## Any of the original sm file to get row and colnums
	output_dir= sys.argv[3]
	if not output_dir.endswith("/"): outptudir+="/"
	
	scoreD= getScoreD(br_matrix)
	rownum, colnum= getRowColNum(original_sm)
	writeScoreMatrix(scoreD, colnum, rownum, output_dir)

if __name__=="__main__":
	main()





				
			
