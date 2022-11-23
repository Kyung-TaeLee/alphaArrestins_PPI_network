## Written by KyungTae Lee
## Last update : 2022-04-17
## Python script : "2.convertToGeneSampleTable.py"
## Script to convert ATAC-seq matrix, which is formatted as rows - genes, columns genomic positions relative to TSS, to 2D matrix which is formatted as rows - gene_relativeCoordnaites and coluimns- sample names to correct ATAC-seq read counts from different batches of batch effects.

#!/usr/bin/env python

import sys
import os
import programExecuteModule as execute

def getCountD(score_mat):
## genomation score matrix output file
	countD= dict()
#	probeidD= dict()
	with open(score_mat, "r") as mat_open:
		mat_lineL= mat_open.readlines()
		mat_header= mat_lineL[0]
		coordL= mat_header.strip().split("\t")[1:]
		for i_line in mat_lineL[1:]:
			infoL= i_line.strip().split("\t")
			generow= infoL[0]
			for i_index in xrange(len(coordL)):
				i_coord= coordL[i_index]
				gene_count= infoL[1+i_index] 
				probeid= "%s_%s"%(generow, i_coord)
				countD[probeid]= gene_count
#				probeidL.append(probeid)
#				probeidD[probeid]=""
#	return countD, probeidD
	return countD

def getFileName(filepath):
	if filepath.find("/")== -1:
		return filepath
	else:
		filename= filepath.split("/")[-1]
		return filename

def getSampleCountD_multi(inputarg, out_q):
	sample_countD= dict()
	for i_arg in inputarg:
		filepath= i_arg
		filename= getFileName(filepath)
		samplename= filename.split(".scoreMatrix.txt")[0]
		i_countD= getCountD(filepath)
		sample_countD[samplename]= i_countD
	out_q.put(sample_countD)	

def writeProbeSampleMat(inputdir, outputname):
## inputdir contaiing score matrix outputs
	inputfileL= os.listdir(inputdir)
	inputfileL= filter(lambda x: x.endswith("scoreMatrix.txt"), inputfileL)
	inputfileL= map(lambda x: inputdir+ x, inputfileL)
	
	sample_countD= dict()
	for i_file in inputfileL:
		print i_file
		filename= getFileName(i_file)
		samplename= filename.split(".scoreMatrix.txt")[0]
		i_countD = getCountD(i_file)
		sample_countD[samplename]= i_countD
#		for i_probe in i_probeidD.keys():
#			if not probeidD.has_key(i_probe):
#				probeidD[i_probe]= ""
#	sample_countD= execute.runMultipleQueueJobs(getSampleCountD_multi, inputfileL, data_type= "dict", proc_num= 2)	
	with open(outputname, "w") as outputopen:
		sampleL= sample_countD.keys()
		header= "probe\t"+ "\t".join(sampleL)+ "\n"
		outputopen.write(header)
		probeidL= sample_countD[sampleL[0]].keys()
		
		for i_probe in probeidL:
			try:
				countL= map(lambda x: sample_countD[x][i_probe], sampleL)
				outputline= i_probe+ "\t"+ "\t".join(countL)+ "\n"
				outputopen.write(outputline)
			except KeyError:
				continue
	

def main():
	inputdir= sys.argv[1]
	if not inputdir.endswith("/"):
		inputdir+= "/"
	
	outputname= sys.argv[2]
	writeProbeSampleMat(inputdir, outputname)

if __name__=="__main__":
	main()
