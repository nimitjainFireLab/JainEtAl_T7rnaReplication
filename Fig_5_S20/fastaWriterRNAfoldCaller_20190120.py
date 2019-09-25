import cPickle
import sys
from moduleUsefulFunctions_20180215 import *
from random import randint
import subprocess
import forgi
import forgi.graph.bulge_graph as fgb

def returnDotBracket(sequence,sampleLabel):
    with open('outputFastaWriterRNAfoldCaller_20190120_'+sampleLabel+'.txt','w') as f:
        processCall=subprocess.Popen(['RNAfold','--noPS','--noGU'],stdout=f,stdin=subprocess.PIPE)
        processCall.communicate(input=sequence)


sortWithPCRDuplicates=0 #parameter to dictate whether sorting is happening with or without subtraction of PCR duplicates, 0 means that each barcode is only counted once for each insert sequence, 1 means that each barcode is counted as many times as it is observed even if the insert sequence is the same
k=20 #used to define identity of a reference sequence

with open(sys.argv[1],'rb') as f:
    testDict=cPickle.load(f)

sampleName=sys.argv[1]
sampleName=sampleName[sampleName.rfind('/')+1:sampleName.rfind('_testDict.pckl')]

seqDict={}
for eachEntry in testDict:
    seqDict[eachEntry]=testDict[eachEntry][sortWithPCRDuplicates]

sortedDictionary=dictSorter(seqDict)
print len(seqDict)
numberSeqProcessed=0
fastaString=''

with open(sampleName+'_toFold_20190120.fasta','w') as f:
    for counter, eachEl in enumerate(sortedDictionary):
        if len(eachEl)>=k:
            line1='>'+str(numberSeqProcessed)+'\n'
            line2=eachEl+'\n'
            f.write(line1)
            f.write(line2)
            fastaString+=line1
            fastaString+=line2            
        numberSeqProcessed+=1
        
returnDotBracket(fastaString,sampleName)
