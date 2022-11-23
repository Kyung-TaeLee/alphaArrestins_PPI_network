## Written by KyungTae Lee
## Last update : 2021-09-10
## Python script : "getGeneL2fcByAcs.py"
## Script to obtain log2 fold changes of genes or ACRs in TXNIP knock down condition compared to control. Genes or peaks will be gruped based on the changes of gene expression or chromatin accessiblity. Output files will be used to plot CDF of changes of gene expression or chromatin accessiblity.

#!/usr/bin/env python

import sys
import statistics

def getFileLineL(inputfile):
	with open(inputfile, "r") as inputopen:
		inputlineL= inputopen.readlines()
		return inputlineL
	

def getEdgerResult(edger_result, feature="gene", symbol_flag=False):
## edgeR result file containing information of all genes or peaks (diffbind)
## For feature= gene, Should contain geneID br_siNeg1_1 br_siNeg_2 br_siTXNIP3_1 br_siTXNIP3_2
## For feature= peak, Should contain seqnames start end br_siNeg_2 br_siTXNIP_1 br_siTXNIP_2
	outputD= dict()

	filelineL= getFileLineL(edger_result)
	header= filelineL[0]
	colL= header.strip().split("\t")
	
	if feature=="gene":

		id_index= colL.index("geneId")
		symbol_index= colL.index("geneSymbol")
		fold_index= colL.index("logFC")
		fdr_index= colL.index("FDR")
		
	elif feature=="peak":
		chr_index= colL.index("chr")
		start_index= colL.index("start")
		end_index= colL.index("end")
		fold_index= colL.index("Fold")
		fdr_index= colL.index("FDR")
#		logmean_index= colL.index("Conc")
		ctrl_index= colL.index("Conc_siNeg")
		sitxnip_index= colL.index("Conc_siTXNIP")
	
	for i_line in filelineL[1:]:
		infoL= i_line.strip().split("\t")
		l2fc= float(infoL[fold_index])

		if feature=="gene":
	
			if not symbol_flag:
				featureid= infoL[id_index]
				symbol= infoL[symbol_index]
				fdr= infoL[fdr_index]
				value= [l2fc, fdr,symbol]
				if not outputD.has_key(featureid):
					outputD[featureid]= value
				else:
					print "%s exist more than once in %s"%(featureid, edger_result)
					sys.exit(1)

			else:
				featureid= infoL[symbol_index]
				value= [l2fc,fdr]

				if not outputD.has_key(featureid):
					outputD[featureid]= list()
				outputD[featureid].append(value)		
			
		elif feature=="peak":
			chr= infoL[chr_index]
			start= infoL[start_index]
			end= infoL[end_index]
			featureid= "%s_%s_%s"%(chr, start, end)
			fdr= float(infoL[fdr_index])
#			logmean= float(infoL[logmean_index])
			ctrl_conc= float(infoL[ctrl_index])
			sitxnip_conc= float(infoL[sitxnip_index])
			l2fc= float(infoL[fold_index])

			value= [ctrl_conc, sitxnip_conc, fdr,l2fc]
			
			if not outputD.has_key(featureid):
				outputD[featureid]= value
			else:
				print "%s exist more than once in %s"%(featureid, edger_result)
				sys.exit(1)

	return outputD

def getGenePeakAnno(chippeak_anno):
## ChIPPeakAnno result file containig colmuns of seqnames start end feature
	gene_peak_interD= dict()

	filelineL= getFileLineL(chippeak_anno)
	header= filelineL[0]
	colL= header.strip().split("\t")
	chr_index= colL.index("seqnames")
	start_index= colL.index("start")
	end_index= colL.index("end")
	feature_index= colL.index("feature")
	distance_index= colL.index("distance")
		
	for i_line in filelineL[1:]:
		infoL= i_line.strip().split("\t")
		chr= infoL[chr_index]
		start= infoL[start_index]
		end= infoL[end_index]
		peakid= "%s_%s_%s"%(chr, start, end)
		geneid= infoL[feature_index]
		distance= infoL[distance_index]
		if not gene_peak_interD.has_key(geneid):
			gene_peak_interD[geneid]= list()
		peakinfo= [peakid, distance]
		if not peakinfo in gene_peak_interD[geneid]:
			gene_peak_interD[geneid].append(peakinfo)
	return gene_peak_interD

def getGenePeakInter(chipseeker):
## ChipSeeker result
	gene_peak_interD= dict()
	with open(chipseeker, "r") as fileopen:
		filelineL= fileopen.readlines()
		for i_line in filelineL[1:]:
			infoL= i_line.strip().split("\t")
			peakid= infoL[0]
			geneid= infoL[1]
			if not gene_peak_interD.has_key(geneid):
				gene_peak_interD[geneid]= list()
			if not peakid in gene_peak_interD[geneid]:
				gene_peak_interD[geneid].append(peakid)
	return gene_peak_interD

def getEnhancerInter(hela_enhancer_inter):
## Hela Enhancer intearction file from enhancer atlas V2
	gene_enhancer_interD= dict()
	with open(hela_enhancer_inter, "r") as enhancer_open:
		enhancer_lineL= enhancer_open.readlines()
		for i_line in enhancer_lineL:
			infoL= i_line.strip().split("\t")
			chr= infoL[0]
			start= str(int(infoL[1])+1)
			end= infoL[2]
			peakid= "%s_%s_%s"%(chr, start, end)
			genesymbol= infoL[9]
			if not gene_enhancer_interD.has_key(genesymbol):
				gene_enhancer_interD[genesymbol]= list()
			if not peakid in gene_enhancer_interD[genesymbol]:
				gene_enhancer_interD[genesymbol].append(peakid)
	return gene_enhancer_interD	

def writeGeneL2fcDensity (gene_outputD, peak_outputD, gene_peak_interD, outputname, peak_measure= "max", group_based="peak", fdr_flag=True):

## Currently, peak_measure should be either "logmean", "median", "mean"
## group_based should be either "peak" or "gene". it will be the feature on which group of other features will be based
	
	def getGroupName(l2fc, fdr=None):
		l2fc= float(l2fc)
		if fdr is None:
		## Not considering FDR
			if l2fc >= 1:
				group="up"
			elif l2fc<=-1:
				group= "down"
			elif l2fc >-0.5 and l2fc < 0.5 :
				group="nochange"
			else:
				return None
		else:
		## Considering FDR
			fdr= float(fdr)
			if l2fc >= 1 and fdr<= 0.05 :
				group="up"
			elif l2fc<=-1 and fdr <=0.05:
				group= "down"
			elif l2fc >-0.5 and l2fc < 0.5 :
				group="nochange"
			else:
				return None
		
		return group

	def getMean(valueL):
		mean= float(sum(valueL))/ len(valueL)
		return mean

	with open(outputname, "w") as outputopen:
		if group_based=="gene":
			header="peak_l2fc\tgene_group\n"
		elif group_based=="peak":
			header= "gene_l2fc\tpeak_group\n"
		outputopen.write(header)
		for i_gene in gene_peak_interD.keys():
			try:
				gene_l2fc= gene_outputD[i_gene][0]
			except KeyError:
				continue
			gene_fdr= gene_outputD[i_gene][1]
#			peakannoL= gene_peak_interD[i_gene]
			peakidL= gene_peak_interD[i_gene]
#			if len(peakannoL)==1:
			if len(peakidL)==1:
#				peak= peakannoL[0][0]
				peak= peakidL[0]
				try:
					peak_l2fc= float(peak_outputD[peak][3])
				except KeyError:
					continue
			else:
				gene_peakL= list()
#				for i_peak in peakannoL:
				for i_peak in peakidL:

#					peak_id, peak_dist= i_peak
					peak_id= i_peak
					try:
						peakinfo= peak_outputD[peak_id]		
					except KeyError:
						continue
					ctrl_conc, sitxnip_conc, peak_fdr,peak_l2fc= peakinfo
					peak_info= [peak_id, ctrl_conc, sitxnip_conc, peak_l2fc]
			
					gene_peakL.append(peak_info)
			
				if gene_peakL==[]:
					continue
	
				if peak_measure=="max":
					gene_maxpeakL= map(lambda x: [x[0], max(x[1], x[2]), x[3]], gene_peakL)
#					gene_maxpeakL= map(lambda x: [x[0], getMean([x[1], x[2]]), x[3]], gene_peakL)
					gene_maxpeakL= sorted(gene_peakL, key=lambda t: float(t[1]), reverse=True)
				## If multiple peaks are assigned to single gene, peaks of the highest mean count will be chosen
#				elif peak_measure=="distance":
#					gene_peakL= sorted(gene_peakL, key=lambda t: abs(float(t[1])), reverse=False)
				## If multiple peaks are assigned to single gene, peaks that are closest to TSS will be chosen
					peak_l2fc= float(gene_maxpeakL[0][3])
				elif peak_measure=="median":
				## median log2 fold change
					peak_l2fcL= map(lambda x: float(x[3]), gene_peakL)
					peak_l2fc= statistics.median(peak_l2fcL) ## median peak L2fc
				elif peak_measure=="mean":
				## mean siTXNIP log2 concentration / mean siNegative log2 concentration
					sineg_meanconc= getMean(map(lambda x: x[1], gene_peakL))
					sitxnip_meanconc= getMean(map(lambda x: x[2], gene_peakL))
					peak_l2fc= sitxnip_meanconc - sineg_meanconc
				else:
					print "peak_measure option should be either logmean or median"
					sys.exit(1)
			if group_based=="gene":
				if not fdr_flag:
				## not considering gene FDR
					group= getGroupName(gene_l2fc)
				elif fdr_flag:
				## considering gene FDR
					group= getGroupName(gene_l2fc, fdr= gene_fdr)
		
				if group is None:
					continue
				outputline= str(peak_l2fc)+ "\t"+ group+ "\n"
 
			elif group_based=="peak":
				group= getGroupName(peak_l2fc)
				if group is None:
					continue
				outputline= str(gene_l2fc)+ "\t"+ group+ "\n"
			else:
				"group_based should be either \"gene\" or \"peak\""
				sys.exit(1)
			
						
			outputopen.write(outputline)

def writeGeneAcsL2fcScatter (gene_outputD, peak_outputD, gene_peak_interD, outputname, peak_measure="logmean"):
## Currently, peak_measure should be either "logmean" or "distance"
	with open(outputname, "w") as outputopen:
		header= "geneid\tgenesymbol\tgene_l2fc\tgene_fdr\tacs_l2fc\tacs_fdr\n"
		outputopen.write(header)
		for i_gene in gene_peak_interD.keys():
			try:
				gene_l2fc= gene_outputD[i_gene][0]
				gene_fdr= gene_outputD[i_gene][1]
				genesymbol= gene_outputD[i_gene][2]
			except KeyError:
				continue

#			peakannoL= gene_peak_interD[i_gene]
			peakidL= gene_peak_interD[i_gene]
#			if len(peakannoL)==1:
			if len(peakidL)==1:
#				peak= peakannoL[0][0]
				peak= peakidL[0]
			else:
				gene_peakL= list()
#				for i_peak in peakannoL:
				for i_peak in peakidL:

#					peak_id, peak_dist= i_peak
					peak_id= i_peak
					try:
						peakinfo= peak_outputD[peak_id]		
					except KeyError:
						continue
#					peak_logmean, peak_fdr,peak_l2fc= peakinfo
					if peak_measure=="logmean":
						peak_logmean= peakinfo[-3]
						peak_info= [peak_id, peak_logmean]
#					elif peak_measure=="distance":
#						peak_info= [peak_id, peak_dist]
					gene_peakL.append(peak_info)
			
				if gene_peakL==[]:
					continue
	
				if peak_measure=="logmean":
					gene_peakL= sorted(gene_peakL, key=lambda t: float(t[1]), reverse=True)
				## If multiple peaks are assigned to single gene, peaks of the highest mean count will be chosen
#				elif peak_measure=="distance":
#					gene_peakL= sorted(gene_peakL, key=lambda t: abs(float(t[1])), reverse=False)
				## If multiple peaks are assigned to single gene, peaks that are closest to TSS will be chosen
				peak= gene_peakL[0][0]

			try:
				peak_fdr= peak_outputD[peak][-2]
				peak_l2fc= peak_outputD[peak][-1]
			except KeyError:
				continue
			outputline= i_gene+ "\t"+ genesymbol+ "\t"+str(gene_l2fc)+ "\t"+ str(gene_fdr)+ "\t"+ str(peak_l2fc)+"\t"+ str(peak_fdr)+ "\n"
			outputopen.write(outputline)

def getMedianD(l2fcD):
	medianD= dict()
	for i_gene in l2fcD.keys():
		l2fcL= l2fcD[i_gene]
		median_l2fc= statistics.median(l2fcL)
		medianD[i_gene]= [median_l2fc]
	return medianD
		
def main():
	edger_result= sys.argv[1]
	## gene edger result
	diffbind_result= sys.argv[2]
	## peak edger result
#	inter_anno_type= sys.argv[3]
	## should be either tss or enhancer
	anno_file= sys.argv[3]
	## inter_anno_type == tss : chipPeakAnno file
	## inter_anno_type == enhancer : enhancer bedtools intersect output 
	## Currently, ChipSeeker interaction txt file that contains peak and gene columns should be given (2021-09-10)
	outputname= sys.argv[4]
	analysis= sys.argv[5]
	## should be either "density" or "scatter"
	peak_measure= sys.argv[6]
	## should be either "max", "median", "mean"
	group_based= sys.argv[7]
	## should be either "gene" or "peak"
	fdr_flag= sys.argv[8]
	if fdr_flag=="true":
		fdr_flag=True
	elif fdr_flag=="false":
		fdr_flag= False
	else:
		print "fdr_flag should be either true or false"
		sys.exit(1)
	
#	if inter_anno_type=="tss":
	genel2fcD =getEdgerResult(edger_result, feature="gene")
#		gene_peak_interD=getGenePeakAnno(anno_file)
	gene_peak_interD= getGenePeakInter(anno_file)
	
	if analysis=="density":
		peak_infoD =getEdgerResult(diffbind_result, feature="peak")
		writeGeneL2fcDensity (genel2fcD, peak_infoD, gene_peak_interD, outputname, peak_measure=peak_measure, group_based= group_based, fdr_flag=fdr_flag)
	elif analysis== "scatter":
#			peak_infoD =getEdgerResult(diffbind_result, feature="peak", peak_group= False)
		peak_infoD =getEdgerResult(diffbind_result, feature="peak")
		writeGeneAcsL2fcScatter(genel2fcD, peak_infoD, gene_peak_interD, outputname, peak_measure= peak_measure)
	else:
		print "analysis should be either density or scatter"


if __name__=="__main__":
	main()	

