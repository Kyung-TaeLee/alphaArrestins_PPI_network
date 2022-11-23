#!/usr/bin/env python

## Written by KyungTae Lee
## Last update : 2022-02-01
## Python script : "flybaseGeneAnnotation.py"
## Containing modules to convery Fly identifiers to primary identifier (Gene Id : fbgn--- or gene symbol : CG--).

## Example usage of the script
#flybase_map= dataset.getFlybaseMap()
#flybase_annoD, flybase_annoIdL= flybase.getFlybaseD(flybase_map, primary_id= "annotation_id")
## Above should be decalred at the beginning
#flybase.convertFlySymbol(theGeneId, theFlybaseAnnoD, theFlybaseAnnoIdL, primary_id= "annotaton_id")
## The parameters, primary_id, in the "flybase.getFlybaseD()" and "flybase.convertFlySymbol()" should be the same

import sys

def getFlybaseD(flybase_anno_table, primary_id= "annotation_id"):	
## Input1 ("flybase_anno_table")= Flybase annotation file of gene Id and symbol (output from the function, \
## "getFlybaseGeneAnnoFile()"
## Input2 ("primary_id") =  Should be either "annotation_id"(default),"gene_symbol"  or "gene_id"           
## Output1 ("flybase_annoD")= Dictionary that has gene identifiers other than the primary ID (Primary gene symbol - CG--- or primary gene ID - fbgn---) as the key and list of primary ID(s) as the value. 
## Output2 ("flybase_anno_idL")= List of primary gene ID (or symbols) 
## This function mainly focuses on gene symbol. If primary_id is set to "gene_id", Only the primary gene symbol
## can be converted to primary gene ID (Should be aware of it!)
	flybase_open= open(flybase_anno_table, 'r')
	flybase_lines= filter(lambda x: not x.startswith("#")  and x!= '' and  x!= '\n' and x!= '\r\n' and \
							x.split('\t')[1]=="Dmel", flybase_open.readlines())
	flybase_open.close()

	flybase_annoD= dict()
	flybase_anno_idL= []
	def addToDicString(gene_id, anno_id, targetD):
		if gene_id==anno_id:
			return None
		if not targetD.has_key(gene_id):
			targetD[gene_id]= [anno_id]
		elif not anno_id in targetD[gene_id]:
			targetD[gene_id].append(anno_id)
#			print "%s  origin is ambiguous"% gene_id
	def addToDicList(gene_idL, anno_id, targetD):
		for i_gene in gene_idL:
			if i_gene=='':
				continue
			addToDicString(i_gene, anno_id, targetD)

	for i_line in flybase_lines:
		line= i_line.strip('\n').strip('\r\n')
		identifiers= line.split('\t')
		try :
			gene_symbol= identifiers[0]
			prim_fbgn= identifiers[2]
			second_fbgns= identifiers[3].split(',')
			primary_symbol= identifiers[4]
			second_annos= identifiers[5].split(',')
		except IndexError:
			print i_line
			sys.exit(1)
		if primary_id== "annotation_id":
			anno_id= primary_symbol
			second_prim_id= prim_fbgn
		elif primary_id== "gene_id":
			anno_id= prim_fbgn
			second_prim_id= primary_symbol
		elif primary_id== "gene_symbol":
			anno_id= gene_symbol
			second_prim_id= primary_symbol
		else:
			print "Type of gene identifiers for drosophila shoud be either %s , %s or %s"\
					%("annotataion_id", "gene_id", "gene_symbol")
			sys.exit(1)
		if not anno_id in flybase_anno_idL:
			flybase_anno_idL.append(anno_id)
		else:
			print "%s  present more than once as primary annotation Id"%anno_id
			sys.exit(1)
		if primary_id== "annotation_id":
			addToDicString(gene_symbol, anno_id, flybase_annoD)
			addToDicList(second_fbgns, anno_id, flybase_annoD)
			addToDicList(second_annos, anno_id, flybase_annoD)
			addToDicString(second_prim_id, anno_id, flybase_annoD)
		elif primary_id== "gene_id" or primary_id== "gene_symbol":
			addToDicString(second_prim_id, anno_id, flybase_annoD)
			addToDicString(prim_fbgn, anno_id, flybase_annoD)		
		
	return flybase_annoD, flybase_anno_idL	

def getFlybaseD_ver2(flybase_anno_table, primary_id= "annotation_id", secondary_id="gene_symbol"):	
## Input1 ("flybase_anno_table")= Flybase annotation file of gene Id and symbol (output from the function,  "getFlybaseGeneAnnoFile()"
## "primary_id" and "secondary_id" should be either "annotation_id"(default),"gene_symbol"  or "gene_id". Primary Ids will be key and secondary ids will be values      
## Output ("flybase_annoD")= Dictionary that has the primary ID as key ans list of secondary Ids as values. 
	flybase_open= open(flybase_anno_table, 'r')
	flybase_lines= filter(lambda x: not x.startswith("#")  and x!= '' and  x!= '\n' and x!= '\r\n' and x.split('\t')[1]=="Dmel", flybase_open.readlines())
	flybase_open.close()

	flybase_annoD= dict()
	def addToDicString_2(gene_id, anno_id, targetD):
#		if gene_id==anno_id:
#			return None
		if not targetD.has_key(gene_id):
			targetD[gene_id]= anno_id
		else:
			print "%s exist more than once in %s"%(gene_id, flybase_anno_table)

	for i_line in flybase_lines:
		line= i_line.strip('\n').strip('\r\n')
		identifiers= line.split('\t')
		try :
			gene_symbol= identifiers[0]
			prim_fbgn= identifiers[2]
			second_fbgns= identifiers[3].split(',')
			prim_annoid= identifiers[4]
			second_annoid= identifiers[5].split(',')
		except IndexError:
			print i_line
			sys.exit(1)
		if primary_id== "annotation_id":
			prim_id= prim_annoid
		elif primary_id== "gene_id":
			prim_id= prim_fbgn
		elif primary_id== "gene_symbol":
			prim_id= gene_symbol
		else:
			print "primary_id should be either annotation_id gene_id gene_symbol"
			sys.exit(1)
		if secondary_id== "annotation_id":
			second_id= prim_annoid
		elif secondary_id== "gene_id":
			second_id= prim_fbgn
		elif secondary_id== "gene_symbol":
			second_id= gene_symbol
		else:
			print "secondary_id should be either annotation_id gene_id gene_symbol"
			sys.exit(1)
		addToDicString_2(prim_id, second_id, flybase_annoD)
	
	return flybase_annoD


def getDefaultFlybaseDL(flybase_map):
## flybase_map is gene annotation file downloaded from Flybase
	flybaseD, flybase_idL= getFlybaseD(flybase_map)
	return flybaseD, flybase_idL

def getFlybaseDL(primary_id_= "annotation_id"):
## flybase_map is gene annotation file downloaded from Flybase
	flybaseD, flybase_idL= getFlybaseD(flybase_map, primary_id= primary_id_)
	return flybaseD, flybase_idL	

def convertFlySymbol(fly_id, gene_name_symbolD, gene_orf_names, primary_id= "annotation_id"):
## Input1 ("fly_id")= fly_id (fbgn---, CG---, gene symbol, etc)
## Input2 ("gene_name_symbolD")= Dictionary generated from the function, "getFlybaseD()" as the first output
## Input3 ("gene_orf_names")= List of gene Id (or symbols) generated from the function, "getFlybaseD()" as the \
## second output
## Input4 ("primary_id")= Should be either "annotation_id", "gene_symbol" or "gene_id". It should be matchec with \
## the same primary_id set in the function, "getFlybaseD()" 
## Output ("fly_id")= Primary gene identifer (ID or symbol) of the query ID will be returned

	if fly_id.startswith("Dmel_"):
		fly_id= fly_id.split("Dmel_")[1]
	elif fly_id.startswith("Dmel/"):
		fly_id= fly_id.split("Dmel/")[1] 

	if primary_id== "annotation_id":
		exceptionD= {"lack":"CG4943","eIF2B-beta":'CG2677', "CG3884":"CG3884","FBpp0084874":"CG7883",\
						"FBpp0085012":"CG31022"}
	elif primary_id== "gene_id":
		exceptionD= {"lack":"FBgn0029006","eIF2B-beta":'FBgn0024996',"CG1135":"FBgn0263832"}
		gene_name_symbolD["CG1135"]= ["FBgn0263832"]	
		gene_name_symbolD["CG1453"]= ["FBgn0030268"]
	elif primary_id== "gene_symbol":
		exceptionD={"lack":"Smurf", "FBpp0084874":"eIF2Balpha", "eIF2B-beta": "eIF2Bbeta",\
					"CG3884":"CG3884", "FBpp0085012":"PH4alphaEFB"} 
#       if fly_id== "CG4674":
#           fly_id="Leash" 
#       elif fly_id== "CG7047":
#           "Vdup1"
	else:
		print primary_id
		print "primary_id option should be set to either \"annotation_id\", \"gene_id\" or \"gene_symbol\""
		sys.exit(1)
	if fly_id in gene_orf_names :
		return fly_id
	elif gene_name_symbolD.has_key(fly_id):
		fly_symbol= gene_name_symbolD[fly_id]
		if len(fly_symbol) >1:
			print "%s  has more than one annotation IDs"%fly_id
			return fly_id
		else:
			fly_symbol= fly_symbol[0]
			return fly_symbol
	elif exceptionD.has_key(fly_id):
		fly_symbol= exceptionD[fly_id]
		return fly_symbol
	else:
#		print "\"%s\" is not recognizable"%fly_id
		return fly_id	

