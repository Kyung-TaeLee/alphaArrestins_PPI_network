## Written by KyungTae Lee
## Last update : 2022-05-25
## Python script : "getOrthoInterTable.py"
## Script to summarize interaction between alpha-arrestin and ortholog proteins between human and fly. Manually curated biological functions of orthologous proteins are also summarized by this script

#!/usr/bin/env python

import sys
import copy

def getPreyBaitSpecD(saintexp_output):
## SAINTexpress output should be given

	saintD= dict()
	with open(saintexp_output, "r") as saintopen:
		saintlineL= saintopen.readlines()
		saintheader= saintlineL[0]
		colL= saintheader.strip().split("\t")
		bait_index= colL.index("Bait")
		prey_index= colL.index("Prey")
		spec_index= colL.index("AvgSpec")
		
		for i_line in saintlineL[1:]:
			infoL= i_line.strip().split("\t")
			bait= infoL[bait_index]
			prey= infoL[prey_index]
			spec= infoL[spec_index]
			if not saintD.has_key(prey):
				saintD[prey]= dict()
			if not saintD[prey].has_key(bait):
				saintD[prey][bait]= spec
			else:
				print "information of interaction between %s and %s exist more than once in %s"%(bait,prey, i_line)
				sys.exit(1)
	return saintD

def getOrthoD(ortho_table, gotermD):
	bait_list= ["CG7047","CG2993","CG10086","ARRDC3","CG18748","CG1105","ARRDC2","CG14696","ARRDC4","ARRDC5","ARRDC1","CG18746","CG18745","CG4674","TXNIP","CG2641"]
	
	def checkListOverlap(list1,list2):
		unionL= list1+list2
		uniqueL= list(set(unionL))
		if len(unionL) != len(uniqueL):
			return True
		else:
			return False
	def getGroupName(geneL):
		for i_gene in geneL:
			if i_gene.startswith("CG"):
				continue
			else:
				return "ortho_"+ i_gene
	def getGoTerm(geneL):
		gotermL= list()
		for i_gene in geneL:
			try:
				goterm= gotermD[i_gene]
			except KeyError:
				continue
			if not goterm in gotermL:
				gotermL.append(goterm)
		if len(gotermL)!= 1:
			print "Multiple go terms", gotermL
			collapse_goterm= "other"
		else:		
			collapse_goterm= gotermL[0]
		return collapse_goterm

	with open(ortho_table, "r") as ortho_open:
		ortho_lineL= ortho_open.readlines()
		ortho_list= list()
		for i_line in ortho_lineL[1:]:
			geneL= i_line.strip().split("\t")
			bait_flag=False
			for i_gene in geneL:
				if i_gene in bait_list:
					bait_flag=True
				else:
					pass
			if bait_flag:
				continue
			else:
				ortho_list.append(geneL)

		orthoD= dict()		
		ortho_gotermD= dict()
		group=0
		while ortho_list!= []:
			query_geneL = ortho_list[0]
			other_geneL= ortho_list[1:]

			overlapflag=True
			while overlapflag:
				overlapflag=False
				nonoverlapL= list()
				for i_gene in other_geneL:
					listoverlap= checkListOverlap(query_geneL, i_gene)
					if listoverlap:
						query_geneL+= filter(lambda x: not x in query_geneL, i_gene)
						overlapflag=True
					else:
						nonoverlapL.append(i_gene)
				other_geneL= copy.deepcopy(nonoverlapL)

			ortho_group= getGroupName(query_geneL)
			goterm= getGoTerm(query_geneL)
			orthoD[ortho_group]= query_geneL
			ortho_gotermD[ortho_group]= goterm
			ortho_list= nonoverlapL
			group+=1
		print "Number of groups : %d"%group
		print orthoD	
	return orthoD, ortho_gotermD

def getTermD(collapse_goterm):
## collapsed GO term file containing name of prey and collapsed GO terms
	gotermD= dict()
	with open(collapse_goterm, "r") as fileopen:
		filelineL= fileopen.readlines()
		header= filelineL[0]
		colL= header.strip().split("\t")
		goterm_index= colL.index("Collapsed_GOterm")
		gene_index= colL.index("name")
		
		for i_line in filelineL[1:]:
			print i_line
			infoL= i_line.strip().split("\t")
			goterm= infoL[goterm_index]
			gene= infoL[gene_index]
			gotermD[gene]= goterm
	return gotermD

def getInterTable(saintD, orthoD, outputname):
		
	baitorthoD= dict()		
	baitL= list()			
	for i_ortho in orthoD.keys():
		for i_prey in orthoD[i_ortho]:
			i_saintD= saintD[i_prey]
			for i_bait in i_saintD.keys():
				i_spec= float(i_saintD[i_bait])
				if not baitorthoD.has_key(i_ortho):
					baitorthoD[i_ortho]= dict()
				if not baitorthoD[i_ortho].has_key(i_bait):
					baitorthoD[i_ortho][i_bait]= list()
				baitorthoD[i_ortho][i_bait].append(i_spec)

				if not i_bait in baitL:
					baitL.append(i_bait)
	with open(outputname, "w") as outputopen:
		orthoL= baitorthoD.keys()
		header= "ortho\t"+ "\t".join(baitL)+ "\n"
		outputopen.write(header)
		for i_ortho in orthoL:
			specs= list()
			for i_bait in baitL:
				try:
					specL= baitorthoD[i_ortho][i_bait]
					avgspec= str(float(sum(specL))/ len(specL))
				except KeyError:
					avgspec="0.0"
				specs.append(avgspec)
			outputline= i_ortho+ "\t"+ "\t".join(specs)+ "\n"
			outputopen.write(outputline)
					
def main():
	human_saint= sys.argv[1]
	fly_saint= sys.argv[2]	
	ortho_table= sys.argv[3]
	goterm_table= sys.argv[4]
	outputprefix= sys.argv[5]

	human_saintD= getPreyBaitSpecD(human_saint)
	fly_saintD= getPreyBaitSpecD(fly_saint)			
	mergedD= dict()
	mergedD.update(human_saintD)
	mergedD.update(fly_saintD)
	
	gotermD= getTermD(goterm_table)
	orthoD, ortho_gotermD= getOrthoD(ortho_table, gotermD)
	outputname= outputprefix+ ".txt"
	getInterTable(mergedD, orthoD, outputname)

	group_anno= outputprefix+ "group_anno.txt"
	with open(group_anno, "w") as anno_open:
		header= "group\tcollapsed_goterm\n"
		anno_open.write(header)
		for i_group in ortho_gotermD.keys():
			goterm= ortho_gotermD[i_group]
			outputline= i_group+ "\t"+ goterm+ "\n"
			anno_open.write(outputline)
	

if __name__=="__main__":
	main()
