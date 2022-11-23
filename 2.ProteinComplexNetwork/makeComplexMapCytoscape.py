#!/usr/bin/env python

## Written by KyungTae Lee
## Last update : 2022-03-10
## Python script : "makeComplexMapCytoscape.py"
## Script to 1. cluster COMPLEAT protein complexes and GO cellular components based on overlap coefficient (>= 0.5), 2. Format the clustered complex information to be used in Cytoscape

import numpy as np
import sys
import os
import re
import scipyModule as statmodule
import flybaseGeneAnnotation as flybase

def getFileLineL(file_name):
    file_open= open(file_name, "r")
    file_lineL= file_open.readlines()
    file_open.close()
    return file_lineL

def getCompleatSummaryD(compleat_summaryF, species):
	summary_lineL= getFileLineL(compleat_summaryF)
	bait_flag=0
	complex_flag=1
	score_flag=3
	pvalue_flag= 4
	name_flag=6
	if species=="Human":
		gene_flag=7
	elif species== "Drosophila":
		gene_flag=10
		
	bait_compleatD= dict()

	compleat_geneD= dict()
	compleat_nameD= dict()
	for i_line in summary_lineL[1:]:
		i_line= i_line.strip()
		if i_line=="":
			continue
		i_line= re.sub("\"","",i_line)
		infoL= i_line.split("\t")
		i_bait= infoL[bait_flag]	
		i_complex= infoL[complex_flag]
		i_score= infoL[score_flag]
		i_pvalue= infoL[pvalue_flag]
		i_name= infoL[name_flag]
		i_gene= infoL[gene_flag]
		i_geneL= map(lambda x: x.strip(), i_gene.split(","))
		i_geneL= filter(lambda x: x!= "", i_geneL)
		i_geneL= sorted(i_geneL)

		if not bait_compleatD.has_key(i_complex):
			bait_compleatD[i_complex]= dict()
		
		bait_compleatD[i_complex][i_bait]= [i_score, i_pvalue]

		if not compleat_geneD.has_key(i_complex):
			compleat_geneD[i_complex]= i_geneL
		elif compleat_geneD[i_complex]!= i_geneL:
			print "Gene components do not match for %s"%i_complex
			sys.exit(1)
		if not compleat_nameD.has_key(i_complex):
			compleat_nameD[i_complex]= i_name
	return bait_compleatD, compleat_geneD, compleat_nameD 

def getCompleatInterD(compleat_interF):
## Interction file of no redundancy (hide redundant) should be given
	inter_lineL= getFileLineL(compleat_interF)
	compleat_interD= dict()
	for i_line in inter_lineL[1:]:
		i_line= i_line.strip()
		if i_line=="":
			continue
		i_line= re.sub("\"","",i_line)
		infoL= i_line.split("\t")
		gene1= infoL[0]
		gene2= infoL[1]
		i_complex= infoL[4]
		print i_complex
		complexL= i_complex.split(",")
		complexL= map(lambda x: x.strip(), complexL)
		complexL= filter(lambda x: x!= "", complexL)
		geneL= sorted([gene1, gene2])
		if gene1.find("%")!= -1 or gene2.find("%")!= -1:
			print "%s or %s has % in their Ids"%(gene1, gene2)
			sys.exit(1)
		i_inter= "%".join(geneL)
		for i_complex in complexL:
			if not compleat_interD.has_key(i_complex):
				compleat_interD[i_complex]= list()
			if not i_inter in compleat_interD[i_complex]:
				compleat_interD[i_complex].append(i_inter)
	return compleat_interD

def getInterL(geneL):
	geneL= sorted(geneL)
	gene_interL= list()
	for i_index in xrange(len(geneL)-1):
		gene1= geneL[i_index]
		for n_index in xrange(i_index+1, len(geneL)):
			gene2= geneL[n_index]
			if gene1.find("%")!= -1 or gene2.find("%")!= -1:
				print "%s or %s has % in their Ids"%(gene1, gene2)
				sys.exit(1)

			gene_inter= gene1+ "%"+ gene2
			if not gene_inter in gene_interL:
				gene_interL.append(gene_inter)
	return gene_interL

def getGoTermName(go_term):
	check_wave=re.findall("~", go_term)
	if len(check_wave)>1 or len(check_wave)==0:
		print "%s is not correctly formatted"%go_term
		sys.exit(1)
	term_id, term_name= go_term.split("~")
	return term_id, term_name

def get_GOCC_infoD(go_ccF, species):
	gocc_lineL= getFileLineL(go_ccF)
		
	complexbait_pvalueD= dict()
	complex_geneD= dict()
	complex_interD= dict()
	go_termD= dict()

	for i_line in gocc_lineL[1:]:
		i_line= i_line.strip()
		if i_line=="":
			continue
		i_line= re.sub("\"","",i_line)
		infoL= i_line.split("\t")
		i_bait= infoL[0]
		i_term = infoL[2]
		term_id, term_name= getGoTermName(i_term)
		i_gene= infoL[6]
		geneL= i_gene.split(",")
		geneL= map(lambda x: x.strip(), geneL)
		geneL= filter(lambda x: x!= "", geneL)
		if len(geneL) <3:
			continue
		if species== "Drosophila":
			i_pvalue=infoL[14]
		elif species== "Human":
			i_pvalue=infoL[12]

		if not complexbait_pvalueD.has_key(term_id):
			complexbait_pvalueD[term_id]= dict()
		if not complexbait_pvalueD[term_id].has_key(i_bait):
			complexbait_pvalueD[term_id][i_bait]= i_pvalue
		
		if not complex_geneD.has_key(term_id):
			complex_geneD[term_id]= list()
		complex_geneD[term_id]+= filter(lambda x: not x in complex_geneD[term_id], geneL)
		if not go_termD.has_key(term_id):
			go_termD[term_id]= term_name		

	for i_term in complex_geneD.keys():
		i_geneL= complex_geneD[i_term]
		i_interL= getInterL(i_geneL)			
		complex_interD[i_term]= i_interL

	return complexbait_pvalueD, complex_geneD, complex_interD, go_termD

def getOverlapCoef(list1, list2):
## Calculate overlap coefficent between protein complexes
	overlapL= list(set(list1)&set(list2))
	overlap_count= len(overlapL)
	min_count= sorted([len(list1), len(list2)], reverse=False)[0]
	overlap_coef= float(overlap_count)/min_count
	return overlap_coef

def getComplexGeneL(complex_name, compleat_geneD, goterm_geneD):
	complexL= complex_name.split("-")
	complex_geneL= list()
	for i_complex in complexL:
		if i_complex.startswith("GO"):
			complex_geneL+= filter(lambda x: not x in complex_geneL,  goterm_geneD[i_complex])
		else:
			complex_geneL+= filter(lambda x: not x in complex_geneL, compleat_geneD[i_complex])
	return complex_geneL	

def getComplexInterD(complex_name, compleat_interD, gocc_interD, complex_prefix="", complex_prefix_flag=False):
	complexL= complex_name.split("-")
	complex_interD= dict()
	if complex_prefix_flag:
		gene_prefix= 'complex'+str(complex_prefix)+":"
	else:
		gene_prefix=""

	for i_index in xrange(len(complexL)):
		i_complex= complexL[i_index]
		if i_complex.startswith("GO"):
		
			complex_interL=  gocc_interD[i_complex]
			i_flag= "GOCC"
		else:
			try:
				complex_interL= compleat_interD[i_complex]
				i_flag= "Compleat"
			except KeyError:
				continue
		for i_inter in complex_interL:
			inter_geneL	= sorted(i_inter.split("%"))
			inter_geneL= map(lambda x: gene_prefix+x, inter_geneL)
			inter_genes= "%".join(inter_geneL)
			if not complex_interD.has_key(inter_genes):
				complex_interD[inter_genes]= list()
			if not i_flag in complex_interD[inter_genes]:
				complex_interD[inter_genes].append(i_flag)
	return complex_interD
	
def getAverageSaintScoreD(geneL, saintscoreD, baitL, targetD, complex_prefix="", complex_prefix_flag=True):
## saintscoreD should have bait as first key, gene id as second key and its saintscore as value
	bait_count= len(baitL)
	for original_geneid in geneL:
		if complex_prefix_flag:
			gene_id= "complex"+str(complex_prefix)+":"+ original_geneid
		else:
			gene_id= original_geneid
		gene_scoreL= []
		for i_bait in baitL:
			bait_scoreD= saintscoreD[i_bait]
			if not bait_scoreD.has_key(original_geneid):
				gene_score= 0.0
			else:
				gene_score= float(bait_scoreD[original_geneid])
			gene_scoreL.append(gene_score)
		average_score= float(sum(gene_scoreL))/bait_count	
		if not targetD.has_key(gene_id):
			targetD[gene_id]= average_score
		else:
			print "%s exist more than once"%gene_id
			sys.exit(1)
		
def getComplexFunctionL(complex_name, compleat_functionD, gocc_functionD):
	complexL= complex_name.split("-")
	complex_nameL= list()
	for i_complex in complexL:
		if i_complex.startswith("GO"):
			complex_name= gocc_functionD[i_complex]
		else:
			complex_name= compleat_functionD[i_complex]
		complex_nameL.append(complex_name)
	return complex_nameL		
	
def clusterSingleRound(candidate_complexL, compleat_geneD, goterm_geneD):
	candidate_complexL= sorted(candidate_complexL)
	overlap_coefD= dict()
	for i_index in xrange(len(candidate_complexL)-1):
		complex1= candidate_complexL[i_index]
		complex1_geneL= getComplexGeneL(complex1, compleat_geneD, goterm_geneD)
		for n_index in xrange(i_index+1, len(candidate_complexL)):
			complex2= candidate_complexL[n_index]
			complex2_geneL= getComplexGeneL(complex2, compleat_geneD, goterm_geneD)
			if complex1.find("%")!= -1 or complex2.find("%")!= -1:
				print "% should not be used as separator for %s and %s"%(complex1, complex2)
				sys.exit(1)
			merged_complex= complex1+ "%"+ complex2
			overlap_coef= getOverlapCoef(complex1_geneL, complex2_geneL)
			overlap_coefD[merged_complex]= overlap_coef
	overlap_coefL= sorted(overlap_coefD.items(), key=lambda t: t[1], reverse=True)
	highest_coef= overlap_coefL[0][1]
	if highest_coef <= 0.5:
		return candidate_complexL, "ClusteringDone"
	else:
		highest_merged_complex= overlap_coefL[0][0]
		merged_complexL= highest_merged_complex.split("%")
		non_clusteredL= filter(lambda x: not x in merged_complexL, candidate_complexL)
		merged_complex_name= "-".join(highest_merged_complex.split("%"))
		non_clusteredL.append(merged_complex_name)
		total_clusteredL= non_clusteredL
		return total_clusteredL, "Clustering"

def calculateIQM(valueL):
## Clcaulte InterQuartileMean based on spectral counts
	sort_valueL= sorted(valueL, reverse=False)
	value_len= len(sort_valueL)
	weight_flag= value_len%4
	quartile_index= float(value_len)/ 4
	if weight_flag==0:
		quartile_index= int(quartile_index)
		iq_valueL= sort_valueL[quartile_index:value_len-quartile_index]			
		iq_len= len(iq_valueL)
		iq_sum= sum(iq_valueL)
		iqm= float(iq_sum)/ iq_len
	elif weight_flag!=0 and value_len >3:
		
		weight= 1-(float(quartile_index) - int(quartile_index))
		quartile_index= int(quartile_index)
		iq_valueL= sort_valueL[quartile_index+1: value_len-quartile_index-1]
		iq_len= len(iq_valueL)
		iq_weight_len= weight*2+ iq_len
		weight_value= float(weight)*(sort_valueL[quartile_index]+sort_valueL[value_len-quartile_index])
		iq_valueL.append(weight_value)
		iqm= float(sum(iq_valueL))/ iq_weight_len
	elif weight_flag!=0 and value_len==3:
		iqm= (sort_valueL[1]+ 0.75*(sort_valueL[0]+ sort_valueL[2]))/2.5
	else:
		print sort_valueL
		print "Wrong condition to calculate IQM"
		sys.exit(1)
	return iqm

def getSaintScoreL(saintscore_summaryF):
	saintscore_lineL=getFileLineL(saintscore_summaryF)			
	bait_id= os.path.splitext(saintscore_summaryF.split("/")[-1])[0]
	bait_saintscoreD= dict()
	for i_line in saintscore_lineL:
		i_line= i_line.strip()
		if i_line=="":
			continue
		i_line= re.sub("\"","",i_line)
		prey_id= i_line.split("\t")[0]
		saintscore= float(i_line.split("\t")[1])
		bait_saintscoreD[prey_id]= saintscore
	return bait_id, bait_saintscoreD
	
def getRandomIqmL(scoreL, sample_size, samplingN=1000):
## generate random protein complexs of size 1000 to be used to derive empirical P-value
	np.random.seed(1111)	
	randomL= list()
	for i_index in xrange(samplingN):
		random_out= np.random.choice(scoreL, size= sample_size,replace=False)
		random_out= list(random_out)
		randomL.append(random_out)
	
	random_iqmL= list()
	for i_random in randomL:
		i_iqm= calculateIQM(i_random)
		random_iqmL.append(i_iqm)
	return random_iqmL

def getPvalue(target_value, background_valueL):
	extremeL= filter(lambda x: x>= target_value, background_valueL)
	extreme_count= len(extremeL)
	background_count= len(background_valueL)
	p_value= float(extreme_count)/background_count
	return p_value

def writeF(output_name, lineL, header):
	output_open= open(output_name, "w")
	output_open.write(header)
	for i_line in lineL:
		output_open.write(i_line)
	output_open.close()

def main():
	compleat_summaryF= sys.argv[1]
	compleat_interF= sys.argv[2]		
	## The one with "hide_redundancy"
	gocc_supplementaryF= sys.argv[3]
	saintscore_dir= sys.argv[4]
	if not saintscore_dir.endswith("/"):
		saintscore_dir+= "/"
	species= sys.argv[5]
	output_prefix= sys.argv[6]

	compleat_bait_scoreD, compleat_geneD, compleat_nameD=	getCompleatSummaryD(compleat_summaryF, species)
	compleat_interD= getCompleatInterD(compleat_interF) 	
				
	gocc_bait_scoreD, gocc_geneD, gocc_interD, gocc_termD= get_GOCC_infoD(gocc_supplementaryF, species)

	compleat_complexL= compleat_bait_scoreD.keys()
	gocc_complexL= gocc_bait_scoreD.keys() 				

	total_complexL= compleat_complexL + gocc_complexL

	print len(total_complexL)
	clusteredL, cluster_flag= clusterSingleRound(total_complexL, compleat_geneD, gocc_geneD)	
	print len(clusteredL)
	print "First clustering"
	cluster_index=1
	while cluster_flag != "ClusteringDone":
		clusteredL, cluster_flag= clusterSingleRound(clusteredL, compleat_geneD, gocc_geneD)
		cluster_index+=1
		print "Round"+ str(cluster_index)+ " Clustering : "+ str(len(clusteredL))
	print "ClusteringDone"
		
	
	clusteredD= dict()
	complex_sizeL= list()
	for i_complex in clusteredL:
		i_geneL= getComplexGeneL(i_complex, compleat_geneD, gocc_geneD)
		clusteredD[i_complex]= i_geneL
		i_size= len(i_geneL)
		if not i_size in complex_sizeL:
			complex_sizeL.append(i_size)

	bait_scoreFL= os.listdir(saintscore_dir)
	bait_scoreFL= filter(lambda x: os.path.splitext(x)[-1]== ".txt",bait_scoreFL)
	bait_scoreFL= map(lambda x: saintscore_dir+ x, bait_scoreFL)

	bait_random_iqmD= dict()
	bait_saintscoreD= dict()
	bait_idL= list()

	for i_baitF in bait_scoreFL:
		bait_id, bait_scoreD=getSaintScoreL(i_baitF)
		bait_scoreL= map(lambda x: x[1], bait_scoreD.items())
		bait_saintscoreD[bait_id]= bait_scoreD
		bait_idL.append(bait_id)
		bait_random_iqmD[bait_id]= dict()

		for i_size in complex_sizeL:
			i_size_random= getRandomIqmL(bait_scoreL, i_size, samplingN=1000)
			bait_random_iqmD[bait_id][i_size]= i_size_random

	complex_pvalueD= dict()
	raw_pvalueL= list()

	for i_complex in clusteredL:
		complex_pvalueD[i_complex]= dict()
		i_geneL= clusteredD[i_complex]
		i_size= len(i_geneL)
		for i_bait in bait_idL:
			bait_complex_scoreL= list()
			for i_gene in i_geneL:
				if not bait_saintscoreD[i_bait].has_key(i_gene):
					bait_complex_scoreL.append(0)
				else:
					i_score=bait_saintscoreD[i_bait][i_gene]
					bait_complex_scoreL.append(i_score)
			i_complex_iqm= calculateIQM(bait_complex_scoreL)
			i_random_iqmL= bait_random_iqmD[i_bait][i_size]
			i_pvalue=getPvalue(i_complex_iqm, i_random_iqmL)
			complex_pvalueD[i_complex][i_bait]= [i_pvalue]
			raw_pvalueL.append(i_pvalue)

	corrected_pvalueL= statmodule.fdrCorrection(raw_pvalueL, method="fdr_bh")	
	p_index=0	
	complex_prefix_=0

	complex_info_lineL= list()
	complex_inter_lineL= list()
	bait_complex_inter_lineL= list()	
	complex_saintscoreD=dict()

	total_complex_counts= 0
	for i_complex in clusteredL:
		significant_baitL= list()
		for i_bait in bait_idL:
			corrected_pvalue= corrected_pvalueL[p_index]
			raw_pvalue= float(complex_pvalueD[i_complex][i_bait][0])
			complex_pvalueD[i_complex][i_bait].append(corrected_pvalue)
			p_index+=1
#			if corrected_pvalue<0.05:
			if raw_pvalue <0.05: 
				significant_baitL.append(i_bait)

		if len(significant_baitL)==0:
			continue
		complex_interD= getComplexInterD(i_complex, compleat_interD, gocc_interD, complex_prefix= complex_prefix_, complex_prefix_flag=True)
		inter_counts= len(complex_interD.keys())
		if inter_counts <2:
			continue
		total_complex_counts+=1
		i_geneL= clusteredD[i_complex]	
		i_size= len(i_geneL)	

		i_complex_saintscoreD= getAverageSaintScoreD(i_geneL, bait_saintscoreD, significant_baitL,complex_saintscoreD,complex_prefix= complex_prefix_, complex_prefix_flag=True)
		
		output_line= i_complex+ "\t"
		complex_functionL= getComplexFunctionL(i_complex, compleat_nameD, gocc_termD)
		output_line += " & ".join(complex_functionL)+ "\t"
		output_line += ",".join(i_geneL)+ "\t"+ str(i_size)+ "\t"
		pvalueL= map(lambda x: complex_pvalueD[i_complex][x], bait_idL)
#		pvalueL= map(lambda x: str(x[0])+"/"+str(x[1]), pvalueL)
		pvalueL= map(lambda x: str(x[0]), pvalueL)
		output_line+= "\t".join(pvalueL)+ "\n"
		complex_info_lineL.append(output_line)
				
		complex_interD= getComplexInterD(i_complex, compleat_interD, gocc_interD, complex_prefix= complex_prefix_, complex_prefix_flag=True)		
		complex_prefix_ +=1
		
		if complex_interD=={}:
			continue
		complex_inter_list= complex_interD.keys()
		for i_index in xrange(len(complex_inter_list)):
			i_inter= complex_inter_list[i_index]
			if i_index==0:
				example_prey= i_inter.split("%")[0]
			i_flagL= complex_interD[i_inter]
			if "GOCC" in i_flagL and "Compleat" in i_flagL:
				edge_flag= "Compleat_GOCC"
			elif "GOCC" in i_flagL and not "Compleat" in i_flagL:
				edge_flag= "GOCC"
			elif not "GOCC" in i_flagL and "Compleat" in i_flagL:
				edge_flag= "Compleat" 		
			else:
				print i_flagL
				print "Wrong Flag information"
				sys.exit(1)
			cytoscape_interaction= "\t".join(i_inter.split("%"))
			interaction_info= cytoscape_interaction+ "\t"+ edge_flag+ "\t"+ i_complex+ "\t"+ " & ".join(complex_functionL)+  "\n"
			complex_inter_lineL.append(interaction_info)
			
		bait_complex_interL= map(lambda x: x+ "\t"+ example_prey+"\t"+ str(complex_pvalueD[i_complex][x][0])+ "\t"+ str(complex_pvalueD[i_complex][x][1])+ "\n", significant_baitL)

		bait_complex_inter_lineL+= bait_complex_interL
	
	print "Number of clustered complexes in %s : %d"%(species, total_complex_counts)  		
	complex_summaryT= output_prefix+ "_ClusteredComplex.Compleat_GOCC.txt"
#	bait_pvalue_headerL= map(lambda x: x+ " (Raw-pvalue/BH-adjusted-pvalue", bait_idL)
#	complex_summary_header= "ClusteredComplex\tFunctions\tGenes\tSize\t"+ "\t".join(bait_pvalue_headerL)+ "\n"
	complex_summary_header= "ClusteredComplex\tFunctions\tGenes\tSize\t"+ "\t".join(bait_idL)+ "\n"
	writeF (complex_summaryT, complex_info_lineL, complex_summary_header)
 		
	complex_interT= output_prefix+ "_Interactions_within_complex.txt"
	complex_inter_header= "Gene1\tGene2\tEvidence\tComplex\tFunctions\n"
	writeF(complex_interT, complex_inter_lineL, complex_inter_header)
		
	bait_complex_interT= output_prefix+ "_Bait_complex_interactions.txt"
	bait_complex_inter_header= "Bait\tComplex(Representative_component)\tRawPvalue\tBH-adjusted-pvalue\n"
	writeF(bait_complex_interT, bait_complex_inter_lineL, bait_complex_inter_header)

	complex_gene_scoreT= output_prefix+ "_AverageSaintScores.txt"
	if species=="Human":
		complex_gene_header= "GeneId (ccomplexN)\tOriginalId\tAvergeSAINTexpressScore\n"
		gene_score_lineL= map(lambda x: x+ "\t"+ x.split(":")[1]+"\t"+str(complex_saintscoreD[x])+ "\n", complex_saintscoreD.keys())
	elif species== "Drosophila":
		flybase_annoD, flybase_idL= flybase.getFlybaseDL(primary_id_= "gene_symbol")
		complex_gene_header= "GeneId (ccomplexN)\tOriginalId\tGeneSymbol\tAvergeSAINTexpressScore\n"
		gene_score_lineL= map(lambda x: x+ "\t"+ x.split(":")[1]+"\t"+ flybase.convertFlySymbol(x.split(":")[1], flybase_annoD, flybase_idL, primary_id= "gene_symbol") +"\t"+str(complex_saintscoreD[x])+ "\n", complex_saintscoreD.keys())
	writeF(complex_gene_scoreT, gene_score_lineL, complex_gene_header)

if __name__=="__main__":
	main()
