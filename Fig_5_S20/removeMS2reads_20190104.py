import sys
from moduleUsefulFunctions_20180215 import *
import pysam

for i in range(1,17):
    fileName=str(i)+'.Aligned.sorted.bam'
    alignmentBamFile=pysam.Samfile(fileName,'rb') #pass coordinate sorted bam file---important to have indexing already done for the fetch function in pysam
    sampleName=str(i)
    with open(str(i)+'_trimmomaticP_thresholdTestDict_1_seqErrorCorrect_1_testDict.pckl','rb') as f:
        testDict=cPickle.load(f)
    nameSeqDict={}
    with open(str(i)+'_ms2ToAlign.fastq','rU') as f:
        for counter, eachLine in enumerate(f):
            if counter%4==0:
                seqName1=eachLine.strip()[1:]
            if counter%4==1:
                seqID=eachLine.strip()
                nameSeqDict[seqName1]=seqID

    
    mappingQualityThreshold=60
    numberEditsAllowed=4
    mCigarRequired=31
    pcrDuplicates=0

    totalReadsMS2=0
    totalReadsFilteredMS2=0
    listReadNamesToRemove={}
    outputFile=sampleName+'_alignmentStats_60.txt'
    with open(outputFile,'w') as f:
        f.write('Filtering param mapQ threshold = '+str(mappingQualityThreshold))
        f.write('\n')
        f.write('Filtering param number edits allowed = '+str(numberEditsAllowed))
        f.write('\n')
        f.write('Filtering param M in cigar minimum required = '+str(mCigarRequired))
        f.write('\n')
        f.write('PCR duplicates for counting = '+str(pcrDuplicates))
        f.write('\n')
        #NC_001417.2 is MS2 genome name, length of genome is 3569
        for counter, alignedread in enumerate(alignmentBamFile.fetch('NC_001417.2',0,3569)):
            seqName=alignedread.qname
            totalReadsMS2+=testDict[nameSeqDict[seqName]][pcrDuplicates]
            tagInfo=alignedread.tags
            numberEdits=0
            for eachTag in tagInfo:
                if eachTag[0]=='NM':
                    numberEdits=eachTag[1]
                    break #assume that NM tag only occurs once in list of tags--tested this assumption as well and was good
            mCigarNumber=0
            cigarList=alignedread.cigar
            for eachPart in cigarList:
                if eachPart[0]==0:
                    mCigarNumber+=eachPart[1] #note that there can be multiple tuples in the cigarList where the first coordinate is 0
            if alignedread.mapq>=mappingQualityThreshold and numberEdits<=numberEditsAllowed and mCigarNumber>=mCigarRequired:                   
                totalReadsFilteredMS2+=testDict[nameSeqDict[seqName]][pcrDuplicates]
                listReadNamesToRemove[seqName]=0
        f.write('Total reads MS2: '+str(totalReadsMS2))
        f.write('\n')
        f.write('Total reads MS2 after filtering: '+str(totalReadsFilteredMS2))
        f.write('\n')        
    alignmentBamFile.close()        
    newTestDict={}
    for eachEl in nameSeqDict:
        if eachEl not in listReadNamesToRemove:
            newTestDict[nameSeqDict[eachEl]]=testDict[nameSeqDict[eachEl]]
    with open(sampleName+'_afterMS2removal_testDict.pckl','wb') as f:
        cPickle.dump(newTestDict,f,protocol=cPickle.HIGHEST_PROTOCOL)

#NOTE THAT ALL THE ALIGNEDREAD property names have now been deprecated in the new version of pysam but the old ones still seem to work
