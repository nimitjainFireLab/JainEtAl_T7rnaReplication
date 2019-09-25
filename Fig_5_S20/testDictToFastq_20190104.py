import sys
from moduleUsefulFunctions_20180215 import *

for i in range(1,17):
    with open(str(i)+'_trimmomaticP_thresholdTestDict_1_seqErrorCorrect_1_testDict.pckl','rb') as f:
        testDict=cPickle.load(f)

    seqDictCollapseBarcode={}

    for eachEntry in testDict:
        seqDictCollapseBarcode[eachEntry]=testDict[eachEntry][0]

    sortedDictionary=dictSorter(seqDictCollapseBarcode)
    with open(str(i)+'_ms2ToAlign.fastq','w') as f:
        for counter, eachEl in enumerate(sortedDictionary):
            seqName=str(i)+'_'+str(counter+1)
            f.write('@'+seqName+'\n')
            f.write(eachEl+'\n')
            f.write('+\n')
            f.write('I'*len(eachEl)+'\n')
    
