from Bio.Blast import NCBIXML
import xml.etree.ElementTree as ET
import numpy as np

def getEVal(alScore,K,Lambda,effSearchSpace):
    #expected number of HSPs with score>=S
    return K*effSearchSpace*np.exp(-1.0*Lambda*alScore)

lengthVec=[]
with open('sequencesToBlast_20190305.fasta','rU') as f:
    for counter, eachLine in enumerate(f):
        if eachLine[0]!='>':
            lengthVec.append(len(eachLine.strip()))
    
alScore=[]
with open('blastresults6_outfmt5_20190305.xml','rU') as f:
    records=NCBIXML.parse(f)
    for eachRec in records:
        alScore.append(eachRec.alignments[0].hsps[0].score)

effHSPLen=[]
tree=ET.parse('blastresults5_outfmt5_20190305.xml')
root=tree.getroot()
a=root.findall('./BlastOutput_iterations/Iteration/Iteration_stat/Statistics/Statistics_hsp-len')
for eachEl in a:
    effHSPLen.append(int(eachEl.text))


lambdaPar=1.37406
k=0.710603

numberSeq=46882715
totalLength=174044644328

searchSpace=[]
for counter, eachEl in enumerate(alScore):
    searchSpace.append((totalLength-numberSeq*effHSPLen[counter])*(lengthVec[counter]-effHSPLen[counter]))    

eValVector=[]
for counter, eachEl in enumerate(alScore):
    eValVector.append(getEVal(eachEl,k,lambdaPar,searchSpace[counter]))



