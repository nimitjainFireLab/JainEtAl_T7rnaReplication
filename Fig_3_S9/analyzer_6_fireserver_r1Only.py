import sys

from moduleUsefulFunctions_20180215 import *

###########################

#PARAMETER SECTION

###########################

kmerLengthToFind=11

###########################

#END OF PARAMETER SECTION

###########################

###########################

#DATA SECTION

###########################
#Library data entry
seqErrorCorrect=1 #1 if want to correct for sequencing error

adapterR1='AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'

adapterR1dict={}

for i in range(0,len(adapterR1)+1-kmerLengthToFind):
    adapterR1dict[adapterR1[i:i+kmerLengthToFind]]=i
                        
thresholdTestDict=1
barcodeL=6 #random hexamer
totalReadNum = 0
R1shorterThanBarcode=0
R1hasN=0
seqDict={}
seqDictCollapseBarcode={}
seqDictIncludingPCR={}
adapterCount=0

fileNameR1=sys.argv[1]
sampleName=fileNameR1[fileNameR1.rfind('/')+1:fileNameR1.rfind('_R1.fastq')]+'_thresholdTestDict_'+str(thresholdTestDict)+'_seqErrorCorrect_'+str(seqErrorCorrect)

alignmentMetrics={}
alignmentMetrics['thresholdTestDict']=thresholdTestDict
alignmentMetrics['fileNameR1']=fileNameR1
alignmentMetrics['adapterR1']=adapterR1
alignmentMetrics['seqErrorCorrect']=seqErrorCorrect
alignmentMetrics['barcodeL']=barcodeL
alignmentMetrics['kmerLengthToFind']=kmerLengthToFind

testDictPickleName=sampleName+'_testDict.pckl'
metricsName=sampleName+'_metrics.pckl'
outputDir=''
testOut1=outputDir+sampleName+'_out1.txt'
testOut2=outputDir+sampleName+'_out2.txt'
testOut3=outputDir+sampleName+'_out3.txt'
testOut4=outputDir+sampleName+'_out4.txt'
testOut5=outputDir+sampleName+'_out5.txt'

with open(fileNameR1,'rU') as f1, open(testOut1,'w') as w1, open(testOut2,'w') as w2, open(testOut3,'w') as w3, open(testOut4,'w') as w4, open(testOut5,'w') as w5:
    for counter, eachLine1 in enumerate(f1):
        if counter%10000==0:
            print counter
        if counter%4==1:
            totalReadNum+=1
            R1seq=eachLine1.strip()
            if len(R1seq)>barcodeL:                
                if checkNoN(R1seq):
                    #first check for adapter contamination--here assuming that the adapter is somewhere in the beginning of the R1 and R2 reads, should be less than a barcode away from start of the reads, assumption that R1 and R2 read lengths more than kmerLengthToFind otherwise this code would just be skipped
                    adapterFound=0
                    for i in range(0,len(R1seq)+1-kmerLengthToFind):
                        if R1seq[i:i+kmerLengthToFind] in adapterR1dict and i<=(adapterR1dict[R1seq[i:i+kmerLengthToFind]]+barcodeL):
                            adapterCount+=1
                            adapterFound=1
                            w3.write(R1seq+'\n')
                            w3.write('\n\n')
                            break
                    if adapterFound==0:                
                        dictIncrement(R1seq,seqDict)
                else:
                    R1hasN+=1
            else:
                R1shorterThanBarcode+=1


for eachEntry in seqDict:
    dictIncrement(eachEntry[:-barcodeL],seqDictCollapseBarcode)
    if eachEntry[:-barcodeL] in seqDictIncludingPCR:
        seqDictIncludingPCR[eachEntry[:-barcodeL]]+=seqDict[eachEntry]
    else:
        seqDictIncludingPCR[eachEntry[:-barcodeL]]=seqDict[eachEntry]

testDict={}
for eachEntry in seqDictCollapseBarcode:
    if seqDictCollapseBarcode[eachEntry]>=thresholdTestDict:
        testDict[eachEntry]=np.array([seqDictCollapseBarcode[eachEntry], seqDictIncludingPCR[eachEntry]])

alignmentMetrics['totalReadNum']=totalReadNum
alignmentMetrics['R1shorterThanBarcode']=R1shorterThanBarcode
alignmentMetrics['R1hasN']=R1hasN
alignmentMetrics['adapterCount']=adapterCount
alignmentMetrics['len(seqDict)']=len(seqDict)
alignmentMetrics['len(seqDictCollapseBarcode)']=len(seqDictCollapseBarcode)
alignmentMetrics['len(seqDictIncludingPCR)']=len(seqDictIncludingPCR)
alignmentMetrics['sum(seqDict.values())']=sum(seqDict.values())
alignmentMetrics['sum(seqDictCollapseBarcode.values())']=sum(seqDictCollapseBarcode.values())
alignmentMetrics['sum(seqDictIncludingPCR.values())']=sum(seqDictIncludingPCR.values())
alignmentMetrics['len(testDict)']=len(testDict)
alignmentMetrics['sum(testDict.values())']=np.sum(np.array(testDict.values()),axis=0)
alignmentMetrics['codeUsedToGenerateFiles']=sys.argv[0]

###########################

#END OF DATA SECTION

###########################


with open(testDictPickleName,'wb') as f:
    cPickle.dump(testDict,f,protocol=cPickle.HIGHEST_PROTOCOL)

with open(metricsName,'wb') as f:
    cPickle.dump(alignmentMetrics,f,protocol=cPickle.HIGHEST_PROTOCOL)




