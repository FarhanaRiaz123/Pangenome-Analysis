import sys
import argparse
import re
import pandas as pd
from Bio import SeqIO

parser = argparse.ArgumentParser(
    description="This script parses orthoMCL output to form count table"+
    ", presence-absence matrix and prints core, accessory and unique genes"+
    " in fasta file.")
parser.add_argument("-g",dest="groupFile",default="group.txt",help="Group file created by orthoMCL"+
                    ". Usually located at my_orthomcl_dir/groups/groups.txt")
parser.add_argument("-f",dest="fastaFile",default="goodProteins.fasta",help="goodProteins file"+
                    ". Usually located at my_orthomcl_dir/blast_dir/goodProteins.fasta")
parser.add_argument("-n",dest="nameFile",default="names.txt",help="List of names of genomes used by orthoMCL"+
                    ". List of names of genomes/organisms used while running orthomcl")
parser.add_argument("-s",dest="singletonFile",default="singletons.txt",help="Singletons file"+
                    ". Created by running orthomclSingletons")

#if(len(sys.argv) == 1):
#    parser.print_help(sys.stderr)
#    sys.exit()

#parsing arguments
#args = parser.parse_args()
#groupFile = args.groupFile
#fastaFile = args.fastaFile
#nameFile = args.nameFile
#singletonFile = args.singletonFile

groupFile = "groups.txt"
fastaFile = "goodProteins.fasta"
nameFile = "names.txt"
singletonFile = "singletons.txt"


############ definitions #######3
def addGenesToDict(dictToReturn,groupList,refDict):
    """
    :param groupList:
    :param refDict:
    :return:
    """
    for i in groupList:
        dictToReturn[i] = refDict[i]
    return dictToReturn

#make list of names
nameList = []
with open(nameFile,"r") as f:
    for i in f.readlines():
        nameList.append(i.strip())

#making group dictionary and empty dataframe
groupDict = dict()
df = pd.DataFrame()

with open(groupFile,"r") as f:
    for i in f.readlines():
        #make dummy dict
        dummyDict = {x:0 for x in nameList}
        members = i.split(": ")[1].strip().split(" ") #all genes in group
        orthoGroup = i.split(": ")[0] #ortho group name
        groupDict[orthoGroup] = members #insert into dictionary
        #calculate count
        temp = [re.sub("\|\S+","",x) for x in members]
        tempDict = {x:temp.count(x) for x in temp}
        dummyDict.update(tempDict)
        df = df.append(dummyDict,ignore_index=True)
df_index = list(groupDict.keys())

#add singletons to matrix
#with open(singletonFile,"r") as f:
#    for i in f.readlines():
#        # make dummy dict
#        dummyDict = {x: 0 for x in nameList}
#        genome = i.strip().split("|")[0]
#        dummyDict[genome] = 1
#        df = df.append(dummyDict,ignore_index=True)
#        gene = i.strip().split("|")[1]
#        df_index.append(gene)

df.index = df_index

######################################################################
# Analysis of  ALL Genes (Not Only Core Genes)
#
Genes_FullNames_ALL = []   # ALL Genes (Not Only Core Genes)
for i in df.index:
    Genes_FullNames_ALL = Genes_FullNames_ALL + [str(x) for x in groupDict[i]]  
        #coreGenes_FullNames_temp = list(set(coreGenes_FullNames_temp))             # remove duplicates from FullNames Genes
        #####coreGenesDict = addGenesToDict(coreGenesDict,groupDict[i],seqDict)
Genes_FullNames = list(set(Genes_FullNames_ALL))             # remove duplicates from FullNames Genes


ProteinsByGenomeALLDict = dict()   # All Proteins including duplicates
for i in Genes_FullNames_ALL:       # All Proteins including duplicates
    temp_ProteinName = str(i).split("|")[1].split(":")[0]   
    ProteinsByGenomeALLDict[temp_ProteinName]=[]
  
for i in Genes_FullNames_ALL:       # ALL Proteins including duplicates
    temp_ProteinName = str(i).split("|")[1].split(":")[0]   
    temp_GenomeName = str(i).split("|")[0]  
    ProteinsByGenomeALLDict[temp_ProteinName].append(temp_GenomeName)

ProteinsByGenomeALL_CountDict=dict()  # ALL Proteins including duplicates 
ProteinsByGenomeDuplicatesALL_CountDict=dict()   # Only Duplicates
for k,v in ProteinsByGenomeALLDict.items():
    ProteinsByGenomeALL_CountDict[k] = [len(v),v]
    if len(v)>1:
        ProteinsByGenomeDuplicatesALL_CountDict[k] = [len(v),v]

temp_count = [len(v[1]) for v in list(ProteinsByGenomeALL_CountDict.values())]
ProteinsByGenomeALL_CountDict_Sum = {x:temp_count.count(x) for x in temp_count}


ProteinsByGenomeDict = dict()   # Proteins Without Duplicates
for i in Genes_FullNames:       # Proteins Without Duplicates
    temp_ProteinName = str(i).split("|")[1].split(":")[0]   
    ProteinsByGenomeDict[temp_ProteinName]=[]
  
for i in Genes_FullNames:       # Without Duplicates
    temp_ProteinName = str(i).split("|")[1].split(":")[0]   
    temp_GenomeName = str(i).split("|")[0]  
    ProteinsByGenomeDict[temp_ProteinName].append(temp_GenomeName)

ProteinsByGenome_CountDict=dict()  # Without Duplicates 
ProteinsByGenomeDuplicates_CountDict=dict()   # Only Duplicates
for k,v in ProteinsByGenomeDict.items():
    ProteinsByGenome_CountDict[k] = [len(v),v]
    if len(v)>1:
        ProteinsByGenomeDuplicates_CountDict[k] = [len(v),v]

with open("Proteins.txt",'w') as f:
    for key,value in ProteinsByGenome_CountDict.items():
        f.write("{0:20} : {1} : ".format(key,str(value[0])))
        for v in value[1]:
            f.write(" " + str(v))
        f.write("\n")

with open("DuplicateProteins.txt",'w') as f:
    for key,value in ProteinsByGenomeDuplicates_CountDict.items():
        f.write("{0:20} : {1} : ".format(key,str(value[0])))
        for v in value[1]:
            f.write(" " + str(v))
        f.write("\n")


#############################################################3
# SCO groups
scoGroupsDict = dict()
scoGroups_Names = []
for i in df.index:
    scofound = True
    for j in df.loc[i]:
        if(j!=1):
            scofound = False
            break
    if(scofound == True):
        scoGroups_Names.append(i)
        scoGroupsDict[i] = groupDict[i]

# write SCO groups with count
with open("scoGroups_Names.txt",'w') as f:
    for i in scoGroups_Names:
        f.writelines(i+"\n")

with open("scoGroups.txt",'w') as f:
    for key,value in scoGroupsDict.items():
        f.write(key+":")
        for v in value:
            f.write(" " + v)
        f.write("\n")

#extract specific Protein from SCO Groups, and call it SCO Group GeneName
scoGroups_GeneNameDict = dict()
for key,value  in scoGroupsDict.items():
    for v in value:
        tempgene = str(v).split("|")[1].split(":")[0]  #gene
        if tempgene.lower().find('pf',0,2) != -1:   # if 'pf' found
            scoGroups_GeneNameDict[key] = tempgene
            break

with open("scoGroups_GeneNames.txt",'w') as f:
    for key,value in scoGroups_GeneNameDict.items():
        f.write(key + ": " + value + "\n")
        ###f.write(value + "\n")

#read 3D7_AnnotatedCDS.fasta
seq3D7Dict = SeqIO.to_dict(SeqIO.parse("3D7_AnnotatedCDS.fasta",'fasta'))

scoGroups_Genes = dict()
for i in scoGroups_GeneNameDict.values():
    scoGroups_Genes[i] = seq3D7Dict[i] 

with open("scoGroups_Genes.fasta",'w') as f:
    SeqIO.write(scoGroups_Genes.values(),f,'fasta')   

##########################################################################################


#write count Matrix    
df.to_csv("count_proteins.txt",sep="\t") #count data

#write presence_absense Matrix
df[df>0] = 1
df.to_csv("presence_absence_proteins.txt",sep="\t")

#read goodproteins.fasta
seqDict = SeqIO.to_dict(SeqIO.parse(fastaFile,'fasta'))

coreGenesDict = dict()
coreGenesDict_Full = dict()
accessoryGenesDict = dict()
accessoryGenesDict_Full = dict()

uniqGenesDict = dict()

#print core and accessory genes
coreGenes_Names = []                # only Gene name  (wihtout duplicates)
coreGenes_FullNames = []            # includes Genome|Gene name  (wihtout duplicates)
coreGenes_FullNames_temp = []       # includes Genome|Gene name  (With duplicates)

accessoryGenes_Names = []           # only Gene name (without duplicates)
accessoryGenes_FullNames = []       # includues Genome|Gene name (without duplicates)
accessoryGenes_FullNames_temp = []  # includues Genome|Gene name (With duplicates)

#for i in df.index:
#    if(df.loc[i].sum() == len(nameList)):
        #coreGenesDict = addGenesToDict(coreGenesDict,groupDict[i],seqDict)
#    elif(df.loc[i].sum() < len(nameList) and df.loc[i].sum() > 1):
#        accessoryGenesDict = addGenesToDict(accessoryGenesDict,groupDict[i],seqDict)
for i in df.index:
    if(df.loc[i].sum() == len(nameList)):
        coreGenes_FullNames_temp = coreGenes_FullNames_temp + [str(x) for x in groupDict[i]]  
        #coreGenes_FullNames_temp = list(set(coreGenes_FullNames_temp))             # remove duplicates from FullNames Genes
        #####coreGenesDict = addGenesToDict(coreGenesDict,groupDict[i],seqDict)
    elif(df.loc[i].sum() < len(nameList) and df.loc[i].sum() > 1):
        accessoryGenes_FullNames_temp = accessoryGenes_FullNames_temp + [str(x) for x in groupDict[i]]
        #accessoryGenes_FullNames_temp = list(set(accessoryGenes_FullNames_temp))   # remove duplicates from FullNames Genes
        #####accessoryGenesDict = addGenesToDict(accessoryGenesDict,groupDict[i],seqDict)

coreGenes_FullNames = list(set(coreGenes_FullNames_temp))             # remove duplicates from FullNames Genes
accessoryGenes_FullNames = list(set(accessoryGenes_FullNames_temp))   # remove duplicates from FullNames Genes



# Fill Core Genes Dictionary with its Sequences
for i in coreGenes_FullNames:       
    temp_Name = str(i).split("|")[1].split(":")[0]   # extract Gene only from FullNames
    coreGenes_Names.append(temp_Name)
    coreGenesDict_Full[temp_Name] = seqDict[i]       




coreProteinsByGenomeDict = dict()
for i in coreGenes_FullNames:       
    temp_ProteinName = str(i).split("|")[1].split(":")[0]   
    coreProteinsByGenomeDict[temp_ProteinName]=[]
  
for i in coreGenes_FullNames:       
    temp_ProteinName = str(i).split("|")[1].split(":")[0]   
    temp_GenomeName = str(i).split("|")[0]  
    coreProteinsByGenomeDict[temp_ProteinName].append(temp_GenomeName)

coreProteinsByGenome_CountDict=dict()
coreProteinsByGenomeDuplicates_CountDict=dict()
for k,v in coreProteinsByGenomeDict.items():
    coreProteinsByGenome_CountDict[k] = [len(v),v]
    if len(v)>1:
        coreProteinsByGenomeDuplicates_CountDict[k] = [len(v),v]

with open("coreDuplicateProteins.txt",'w') as f:
    for key,value in coreProteinsByGenomeDuplicates_CountDict.items():
        #f.write(key+" \t: "+str(value[0])+" : ")
        f.write("{0:20} : {1} : ".format(key,str(value[0])))
        for v in value[1]:
            f.write(" " + str(v))
        f.write("\n")

##coreGenes_Names_Duplicates =[]
##coreGenes_Names_Duplicates = coreGenes_Names
##coreGenes_Names_DuplicatesDict={x:coreGenes_Names_Duplicates.count(x) for x in coreGenes_Names_Duplicates}


coreGenes_Names = list(set(coreGenes_Names))    # remove duplicates from Genes
for key,value in coreGenesDict_Full.items():    # remove duplicates from Genes Sequence Dictionary
    if key not in coreGenesDict.keys():
        coreGenesDict[key] = value


# Fill Accessoru Genes Dictionary with its Sequences
for i in accessoryGenes_FullNames:       
    temp_Name = str(i).split("|")[1].split(":")[0]   # extract Gene only from FullNames
    accessoryGenes_Names.append(temp_Name)
    accessoryGenesDict_Full[temp_Name] = seqDict[i]       

accessoryGenes_Names = list(set(accessoryGenes_Names))    # remove duplicates from Genes
for key,value in accessoryGenesDict_Full.items():    # remove duplicates from Genes Sequence Dictionary
    if key not in accessoryGenesDict.keys():
        accessoryGenesDict[key] = value



#get unique genes (singletons)
with open(singletonFile,'r') as f:
    for i in f.readlines():
        uniqGenesDict[i.strip()] = seqDict[i.strip()]


# Core Genes Names
###coreGenes_Names = []
#for i in list(coreGenes.keys()):
#    coreGenes_Names.append(str(i).split("|")[1].split(":")[0] + "\n")


coreGenes_Names_Count = dict()
coreGenes_Names_Count = {x:coreGenes_Names.count(x) for x in coreGenes_Names} 

with open("core_proteins_names.txt",'w') as f:
    for i in coreGenes_Names:
        f.writelines(i+"\n")


#writing files
with open("core_proteins.fasta",'w') as f:
    SeqIO.write(coreGenesDict.values(),f,'fasta')
with open('accessory_proteins.fasta','w') as f:
    SeqIO.write(accessoryGenesDict.values(),f,'fasta')
with open('uniq_proteins.fasta','w') as f:
    SeqIO.write(uniqGenesDict.values(),f,'fasta')





#writing files manually
#with open("core_genes_manual.fasta",'w') as f:
#    for key,value in coreGenesDict.items():
#        f.writelines(">" + key + "\n")
#        f.writelines(value.)
#    SeqIO.write(coreGenesDict.values(),f,'fasta')