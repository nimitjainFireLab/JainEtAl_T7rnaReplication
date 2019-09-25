import sys

from moduleUsefulFunctions_20180215 import *

###########################

#FUNCTION SECTION

###########################

def needlemanWunschShortInsert(seq1,seq2,barcodeL,m=2,s=-1,d=-2,n=0): #MODIFIED FOR ALIGNMENTS BETWEEN R1 and R2 reads
    #m is for matches, s for mismatches, d for gaps, n is for ns in reference--sites where diversity introduced
    #ensure that R1 and R2 reads don't have Ns in them
    seq1=seq1.upper()
    seq2=seq2.upper()
    dpTable=np.zeros((len(seq1)+1,len(seq2)+1))  #assume R1 is seq1 and revComp(R2) is seq2
    pointerTable=[]
    seqEmitTable=[]
    dummyVector=[]
    for i in range(0,len(seq2)+1):
        dummyVector.append(0)
    for j in range(0,len(seq1)+1):
        pointerTable.append(list(dummyVector))
        seqEmitTable.append(list(dummyVector))
    #initialization
    for i in range(1,len(seq1)+1):
        dpTable[i,0]=0 #changed compared to aligner_20170227_ver2.py
        pointerTable[i][0]=(i-1,0)
        seqEmitTable[i][0]=(seq1[i-1],'_')
    for j in range(1,len(seq2)+1):
        dpTable[0,j]=0 
        pointerTable[0][j]=(0,j-1)
        seqEmitTable[0][j]=('_',seq2[j-1])

    #update dynamic programming table
    for i in range(1,len(seq1)+1):
        for j in range(1,len(seq2)+1):
            if seq1[i-1]==seq2[j-1]:
                subScore=m
            else:
                subScore=s
            diagScore=dpTable[i-1,j-1]+subScore
            topDownScore=dpTable[i-1,j]+d
            leftRightScore=dpTable[i,j-1]+d

            if diagScore>=topDownScore and diagScore>=leftRightScore:
                dpTable[i,j]=diagScore
                pointerTable[i][j]=(i-1,j-1)
                seqEmitTable[i][j]=(seq1[i-1],seq2[j-1])
            elif topDownScore>=diagScore and topDownScore>=leftRightScore:
                dpTable[i,j]=topDownScore
                pointerTable[i][j]=(i-1,j)
                seqEmitTable[i][j]=(seq1[i-1],'_')
            elif leftRightScore>=diagScore and leftRightScore>=topDownScore:
                dpTable[i,j]=leftRightScore
                pointerTable[i][j]=(i,j-1)
                seqEmitTable[i][j]=('_',seq2[j-1])
                
    maxScore=d*(max(len(seq1),len(seq2))+1)
    hexamerMismatch=0
    combString=''   
    
    for i in range(len(seq1),-1,-1): #checking last column
        if dpTable[i,len(seq2)]>maxScore:
            maxScore=dpTable[i,len(seq2)]
            endCoordinate=(i,len(seq2))
    currentCoordinate=endCoordinate
    while currentCoordinate[1]>=len(seq2)-(barcodeL-1):
        currentCoordinate=pointerTable[currentCoordinate[0]][currentCoordinate[1]]
    return (seq1[:currentCoordinate[0]+barcodeL], maxScore)

def needlemanWunsch(seq1,seq2,seqErrorCorrect,barcodeL,m=2,s=-1,d=-2,n=0): #MODIFIED FOR ALIGNMENTS BETWEEN R1 and R2 reads--penalizes alignments that are shorter than the read length
    #m is for matches, s for mismatches, d for gaps, n is for ns in reference--sites where diversity introduced
    #ensure that R1 and R2 reads don't have Ns in them
    seq1=seq1.upper()
    seq2=seq2.upper()
    dpTable=np.zeros((len(seq1)+1,len(seq2)+1))  #assume R1 is seq1 and revComp(R2) is seq2
    pointerTable=[]
    seqEmitTable=[]
    dummyVector=[]
    for i in range(0,len(seq2)+1):
        dummyVector.append(0)
    for j in range(0,len(seq1)+1):
        pointerTable.append(list(dummyVector))
        seqEmitTable.append(list(dummyVector))
    #initialization
    for i in range(1,len(seq1)+1):
        dpTable[i,0]=0 #changed from i*d
        pointerTable[i][0]=(i-1,0)
        seqEmitTable[i][0]=(seq1[i-1],'_')
    for j in range(1,len(seq2)+1):
        dpTable[0,j]=j*d #penalizes alignments shorter than read length
        pointerTable[0][j]=(0,j-1)
        seqEmitTable[0][j]=('_',seq2[j-1])

    #update dynamic programming table
    for i in range(1,len(seq1)+1):
        for j in range(1,len(seq2)+1):
            if seq1[i-1]==seq2[j-1]:
                subScore=m
            else:
                subScore=s
            diagScore=dpTable[i-1,j-1]+subScore
            topDownScore=dpTable[i-1,j]+d
            leftRightScore=dpTable[i,j-1]+d

            if diagScore>=topDownScore and diagScore>=leftRightScore:
                dpTable[i,j]=diagScore
                pointerTable[i][j]=(i-1,j-1)
                seqEmitTable[i][j]=(seq1[i-1],seq2[j-1])
            elif topDownScore>=diagScore and topDownScore>=leftRightScore:
                dpTable[i,j]=topDownScore
                pointerTable[i][j]=(i-1,j)
                seqEmitTable[i][j]=(seq1[i-1],'_')
            elif leftRightScore>=diagScore and leftRightScore>=topDownScore:
                dpTable[i,j]=leftRightScore
                pointerTable[i][j]=(i,j-1)
                seqEmitTable[i][j]=('_',seq2[j-1])
                
    maxScore=d*(max(len(seq1),len(seq2))+1)
    hexamerMismatch=0
    combString=''
    
    for i in range(len(seq2),-1,-1): #checking last row
        if dpTable[len(seq1),i]>maxScore:
            maxScore=dpTable[len(seq1),i]
            endCoordinate=(len(seq1),i)
    #maxScore is alignment score    
    if endCoordinate[1]==len(seq2): #max is in last column and last row--this means that random hexamer has a mismatch between R1 and R2
        hexamerMismatch=1
        if seqErrorCorrect==1:
            return (0,combString,hexamerMismatch, maxScore) #the first element of the tuple indicates that it is discordant
        elif seqErrorCorrect==0: #get the R1 read sequence after defining where the barcode aligns
            (alignedSeq, alignmentScore)=needlemanWunschShortInsert(seq1,seq2,barcodeL)
            return(1,alignedSeq,hexamerMismatch, alignmentScore)
    else: #max is in last row but not in last column--this means that the insert size is longer than the read length
        x=len(seq2)-1
        while x!=endCoordinate[1]-1:
            combString=seq2[x]+combString            
            x=x-1
        if seqErrorCorrect==0: #have enough info to return at this point
            return (1,seq1+combString,hexamerMismatch, maxScore)
        elif seqErrorCorrect==1:            
            currentCoordinate=endCoordinate
            while currentCoordinate[1]!=0 and currentCoordinate[0]!=0: #and currentCoordinate[0]!=0 checks for those niche cases where R2 is greater in length than R1 on the 5' end of the insert
                emitted=seqEmitTable[currentCoordinate[0]][currentCoordinate[1]]
                if emitted[0]!=emitted[1]:
                    return (0,'',hexamerMismatch, maxScore)
                combString=emitted[0]+combString
                currentCoordinate=pointerTable[currentCoordinate[0]][currentCoordinate[1]]
            combString=seq1[:currentCoordinate[0]]+combString
            return (1,combString,hexamerMismatch, maxScore)      

###########################

#END OF FUNCTION SECTION

###########################


###########################

#PARAMETER SECTION

###########################

overlapLengthMin=10
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
adapterR2='AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT'

adapterR1dict={}
adapterR2dict={}
for i in range(0,len(adapterR1)+1-kmerLengthToFind):
    adapterR1dict[adapterR1[i:i+kmerLengthToFind]]=i
                        
for i in range(0,len(adapterR2)+1-kmerLengthToFind):
    adapterR2dict[adapterR2[i:i+kmerLengthToFind]]=i

thresholdTestDict=1
barcodeL=6 #random hexamer
totalReadNum = 0
R2shorterThanBarcode=0
hexamerNotFound=0
hexamerAtBeginning=0
eitherR1orR2hasN=0
seqDict={}
seqDictCollapseBarcode={}
seqDictIncludingPCR={}

numberDiscordant=0
hexamerMismatchCounter=0
longInsertAboveLength=0
longInsertBelowLength=0
adapterCount=0

fileNameR1=sys.argv[1]
fileNameR2=sys.argv[2]
sampleName=fileNameR1[fileNameR1.rfind('/')+1:fileNameR1.rfind('_R1.fastq')]+'_thresholdTestDict_'+str(thresholdTestDict)+'_seqErrorCorrect_'+str(seqErrorCorrect)

alignmentMetrics={}
alignmentMetrics['thresholdTestDict']=thresholdTestDict
alignmentMetrics['fileNameR1']=fileNameR1
alignmentMetrics['fileNameR2']=fileNameR2
alignmentMetrics['adapterR1']=adapterR1
alignmentMetrics['adapterR2']=adapterR2
alignmentMetrics['seqErrorCorrect']=seqErrorCorrect
alignmentMetrics['barcodeL']=barcodeL
alignmentMetrics['overlapLengthMin']=overlapLengthMin
alignmentMetrics['kmerLengthToFind']=kmerLengthToFind

testDictPickleName=sampleName+'_testDict.pckl'
metricsName=sampleName+'_metrics.pckl'
outputDir=''
testOut1=outputDir+sampleName+'_out1.txt'
testOut2=outputDir+sampleName+'_out2.txt'
testOut3=outputDir+sampleName+'_out3.txt'
testOut4=outputDir+sampleName+'_out4.txt'
testOut5=outputDir+sampleName+'_out5.txt'

with open(fileNameR1,'rU') as f1, open(fileNameR2,'rU') as f2, open(testOut1,'w') as w1, open(testOut2,'w') as w2, open(testOut3,'w') as w3, open(testOut4,'w') as w4, open(testOut5,'w') as w5:
    for counter, eachLine2 in enumerate(f2):
        if counter%10000==0:
            print counter
        eachLine1=f1.next().strip()
        if counter%4==1:
            totalReadNum+=1
            R1seq=eachLine1
            R2seq=eachLine2.strip()
            if len(R2seq)>barcodeL:                
                if checkNoN(R1seq) and checkNoN(R2seq):
                    #first check for adapter contamination--here assuming that the adapter is somewhere in the beginning of the R1 and R2 reads, should be less than a barcode away from start of the reads, assumption that R1 and R2 read lengths more than kmerLengthToFind otherwise this code would just be skipped
                    adapterFound=0
                    for i in range(0,len(R1seq)+1-kmerLengthToFind):
                        if R1seq[i:i+kmerLengthToFind] in adapterR1dict and i<=(adapterR1dict[R1seq[i:i+kmerLengthToFind]]+barcodeL):
                            adapterCount+=1
                            adapterFound=1
                            w3.write(R1seq+'\n')
                            w3.write(revComp(R2seq)+'\n')
                            w3.write('\n\n')
                            break
                    if adapterFound==0:
                        for i in range(0,len(R2seq)+1-kmerLengthToFind):
                            if R2seq[i:i+kmerLengthToFind] in adapterR2dict and i<=(adapterR2dict[R2seq[i:i+kmerLengthToFind]]+barcodeL):
                                adapterCount+=1
                                adapterFound=1
                                w3.write(R1seq+'\n')
                                w3.write(revComp(R2seq)+'\n')
                                w3.write('\n\n')
                                break
                    if adapterFound==0:                
                        barcode = R2seq[:barcodeL]
                        pos = R1seq.rfind(revComp(barcode))                    
                        if pos>0:                       
                            #for sequencing error correction
                            if seqErrorCorrect==0:
                                dictIncrement(R1seq[:pos+barcodeL],seqDict)
                            else:
                                seqErrorIndex=min(pos+barcodeL,len(R2seq)) #assuming that R1 read is >= in length than R2 read, what if R1 read is shorter than R2 read? This code should still take care of that situation
                                if revComp(R2seq[:seqErrorIndex])==R1seq[pos+barcodeL-seqErrorIndex:pos+barcodeL]:
                                    dictIncrement(R1seq[:pos+barcodeL],seqDict)
                                else:
                                    numberDiscordant+=1                            
                        elif pos==0:
                            hexamerAtBeginning+=1 #expected to be 0 because the first check is for len(R2seq)>barcodeL as opposed to len(R2seq)>=barcodeL
                        elif pos==-1:  #could be because hexamers have mismatches between R1 and R2 reads or because insert is longer than R1 read length                        
                            hexamerNotFound+=1
                            alignerResult=needlemanWunsch(R1seq,revComp(R2seq),seqErrorCorrect,barcodeL)
                            alignedSeq=alignerResult[1]
                            if alignerResult[0]==0: #when seqErrorCorrect=0, this should never be triggered; when seqErrorCorrect=1, this should be triggered for hexamerMismatch cases and for a subset of long inserts when there are errors between R1 and R2
                                numberDiscordant+=1
                                w1.write(R1seq+'\n')
                                w1.write(revComp(R2seq)+'\n')
                                w1.write(alignedSeq+'\n')
                                w1.write('Hexamer mismatch= '+str(alignerResult[2])+'\n')
                                w1.write('Alignment Score= '+str(alignerResult[3])+'\n')
                                w1.write('\n\n')
                            else: #when seqErrorCorrect=0, this should always be triggered; when seqErrorCorrect=1, this should only be triggered for a subset of the long inserts                            
                                if alignerResult[2]==1: #when seqErrorCorrect=0, this can be triggered for hexamerMismatch cases, when seqErrorCorrect=1, this should never be triggered
                                    hexamerMismatchCounter+=1
                                    dictIncrement(alignedSeq,seqDict)
                                    w2.write(R1seq+'\n')
                                    w2.write(revComp(R2seq)+'\n')
                                    w2.write(alignedSeq+'\n')
                                    w2.write('Hexamer mismatch= '+str(alignerResult[2])+'\n')
                                    w2.write('Alignment Score= '+str(alignerResult[3])+'\n')
                                    w2.write('\n\n')                            
                                else: #long insert cases for both seqErrorCorrect=0 and seqErrorCorrect=1
                                    #remove adapter contamination first
                                    #assuming that length of R1seq and R2seq is greater than kmerLengthToFind as long insert case
                                    adapterFound1=0
                                    for i in range(0,len(R1seq)+1-kmerLengthToFind):
                                        if R1seq[i:i+kmerLengthToFind] in adapterR1dict:
                                            adapterCount+=1
                                            adapterFound1=1
                                            w3.write(R1seq+'\n')
                                            w3.write(revComp(R2seq)+'\n')
                                            w3.write(alignedSeq+'\n')
                                            w3.write('Hexamer mismatch= '+str(alignerResult[2])+'\n')
                                            w3.write('Alignment Score= '+str(alignerResult[3])+'\n')
                                            w3.write('\n\n')
                                            break
                                    if adapterFound1==0:
                                        for i in range(0,len(R2seq)+1-kmerLengthToFind):
                                            if R2seq[i:i+kmerLengthToFind] in adapterR2dict:
                                                adapterCount+=1
                                                adapterFound1=1
                                                w3.write(R1seq+'\n')
                                                w3.write(revComp(R2seq)+'\n')
                                                w3.write(alignedSeq+'\n')
                                                w3.write('Hexamer mismatch= '+str(alignerResult[2])+'\n')
                                                w3.write('Alignment Score= '+str(alignerResult[3])+'\n')
                                                w3.write('\n\n')
                                                break
                                    if adapterFound1==0:
                                        alignedLength=min(len(R2seq)-(len(alignedSeq)-len(R1seq)),len(R1seq))
                                        if alignedLength>=overlapLengthMin:                                                                 
                                            longInsertAboveLength+=1
                                            dictIncrement(alignedSeq,seqDict)
                                            w4.write(R1seq+'\n')
                                            w4.write(revComp(R2seq)+'\n')
                                            w4.write(alignedSeq+'\n')
                                            w4.write('Hexamer mismatch= '+str(alignerResult[2])+'\n')
                                            w4.write('Alignment Score= '+str(alignerResult[3])+'\n')
                                            w4.write('\n\n')
                                        else:
                                            longInsertBelowLength+=1
                                            w5.write(R1seq+'\n')
                                            w5.write(revComp(R2seq)+'\n')
                                            w5.write(alignedSeq+'\n')
                                            w5.write('Hexamer mismatch= '+str(alignerResult[2])+'\n')
                                            w5.write('Alignment Score= '+str(alignerResult[3])+'\n')
                                            w5.write('\n\n')
                else:
                    eitherR1orR2hasN+=1
            else:
                R2shorterThanBarcode+=1


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
alignmentMetrics['R2shorterThanBarcode']=R2shorterThanBarcode
alignmentMetrics['hexamerNotFound']=hexamerNotFound
alignmentMetrics['hexamerAtBeginning']=hexamerAtBeginning
alignmentMetrics['eitherR1orR2hasN']=eitherR1orR2hasN
alignmentMetrics['numberDiscordant']=numberDiscordant
alignmentMetrics['hexamerMismatchCounter']=hexamerMismatchCounter
alignmentMetrics['longInsertAboveLength']=longInsertAboveLength
alignmentMetrics['longInsertBelowLength']=longInsertBelowLength
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




