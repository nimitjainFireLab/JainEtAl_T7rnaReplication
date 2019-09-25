import sys
import itertools
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio.Phylo.TreeConstruction import _DistanceMatrix
from Bio import Phylo

from moduleUsefulFunctions_20180215 import *

def needlemanWunschAffine(seq1,seq2,m=2,s=-1,d=-2,n=1,e=-1): 
    #m is for matches, s for mismatches, d for gaps, n is for ns in reference--sites where diversity there, e is gap extend penalty
    seq1=seq1.upper()
    seq2=seq2.upper()
    dpTableV=np.zeros((len(seq1)+1,len(seq2)+1)) #assume seq1 is reference
    dpTableF=np.zeros((len(seq1)+1,len(seq2)+1))
    dpTableG=np.zeros((len(seq1)+1,len(seq2)+1))
    dpTableH=np.zeros((len(seq1)+1,len(seq2)+1))
    pointerTableV=[]
    pointerTableF=[]
    pointerTableG=[]
    pointerTableH=[]    
    seqEmitTableV=[]
    seqEmitTableF=[]
    seqEmitTableG=[]
    seqEmitTableH=[]
    
    dummyVector=[]
    for i in range(0,len(seq2)+1):
        dummyVector.append(0)
    for j in range(0,len(seq1)+1):
        pointerTableV.append(list(dummyVector))
        pointerTableF.append(list(dummyVector))
        pointerTableG.append(list(dummyVector))
        pointerTableH.append(list(dummyVector))        
        seqEmitTableV.append(list(dummyVector))
        seqEmitTableF.append(list(dummyVector))
        seqEmitTableG.append(list(dummyVector))
        seqEmitTableH.append(list(dummyVector))
         
    #initialization
    for i in range(1,len(seq1)+1):
        dpTableF[i][0]=2*d*max(len(seq1),len(seq2))
        pointerTableF[i][0]=(i-1,0,'V')
        seqEmitTableF[i][0]=(seq1[i-1],'_',i-1,'deletion')

        dpTableH[i][0]=2*d*max(len(seq1),len(seq2))
        pointerTableH[i][0]=(i-1,0,'V')
        seqEmitTableH[i][0]=(seq1[i-1],'_',i-1,'deletion')

        dpTableG[i][0]=d+(i-1)*e
        pointerTableG[i][0]=(i-1,0,'V')
        seqEmitTableG[i][0]=(seq1[i-1],'_',i-1,'deletion')     

        dpTableV[i][0]=d+(i-1)*e
        pointerTableV[i][0]=(i-1,0,'V')
        seqEmitTableV[i][0]=(seq1[i-1],'_',i-1,'deletion')
        
    for j in range(1,len(seq2)+1):
        dpTableF[0][j]=2*d*max(len(seq1),len(seq2))
        pointerTableF[0][j]=(0,j-1,'V')
        seqEmitTableF[0][j]=('_',seq2[j-1],'insertion',j-1)

        dpTableH[0][j]=d+(i-1)*e
        pointerTableH[0][j]=(0,j-1,'V')
        seqEmitTableH[0][j]=('_',seq2[j-1],'insertion',j-1)

        dpTableG[0][j]=2*d*max(len(seq1),len(seq2))
        pointerTableG[0][j]=(0,j-1,'V')
        seqEmitTableG[0][j]=('_',seq2[j-1],'insertion',j-1)

        dpTableV[0][j]=d+(i-1)*e
        pointerTableV[0][j]=(0,j-1,'V')
        seqEmitTableV[0][j]=('_',seq2[j-1],'insertion',j-1)  
        
    #update dynamic programming table
    
    for i in range(1,len(seq1)+1):
        for j in range(1,len(seq2)+1):
            #First take care of F matrix
            if seq1[i-1]=='N':
                subScore=n
            else:
                if seq1[i-1]==seq2[j-1]:
                    subScore=m
                else:
                    subScore=s
            fNewScore=dpTableV[i-1][j-1]+subScore
            dpTableF[i][j]=fNewScore
            pointerTableF[i][j]=(i-1,j-1,'V')  #need to know which matrix it came from
            seqEmitTableF[i][j]=(seq1[i-1],seq2[j-1],i-1,j-1)

            #Now address G matrix    
            gScore1=dpTableV[i-1][j]+d
            gScore2=dpTableG[i-1][j]+e
            if gScore1>=gScore2:
                dpTableG[i][j]=gScore1
                pointerTableG[i][j]=(i-1,j,'V')
                seqEmitTableG[i][j]=(seq1[i-1],'_',i-1,'deletion')
            elif gScore2>gScore1:
                dpTableG[i][j]=gScore2
                pointerTableG[i][j]=(i-1,j,'G')
                seqEmitTableG[i][j]=(seq1[i-1],'_',i-1,'deletion')

            #Now address H matrix
            hScore1=dpTableV[i][j-1]+d
            hScore2=dpTableH[i][j-1]+e
            if hScore1>=hScore2:
                dpTableH[i][j]=hScore1
                pointerTableH[i][j]=(i,j-1,'V')
                seqEmitTableH[i][j]=('_',seq2[j-1],'insertion',j-1)
            elif hScore2>hScore1:
                dpTableH[i][j]=hScore2
                pointerTableH[i][j]=(i,j-1,'H')
                seqEmitTableH[i][j]=('_',seq2[j-1],'insertion',j-1)

            #Finally the V matrix
            if dpTableF[i][j]>=dpTableG[i][j] and dpTableF[i][j]>=dpTableH[i][j]:
                dpTableV[i][j]=dpTableF[i][j]
                pointerTableV[i][j]=pointerTableF[i][j]
                seqEmitTableV[i][j]=seqEmitTableF[i][j]
            elif dpTableG[i][j]>=dpTableF[i][j] and dpTableG[i][j]>=dpTableH[i][j]:
                dpTableV[i][j]=dpTableG[i][j]
                pointerTableV[i][j]=pointerTableG[i][j]
                seqEmitTableV[i][j]=seqEmitTableG[i][j]
            else:
                dpTableV[i][j]=dpTableH[i][j]
                pointerTableV[i][j]=pointerTableH[i][j]
                seqEmitTableV[i][j]=seqEmitTableH[i][j]
                
    #construct alignment--first find max value in dpTableV
    maxIndexTuple=(len(seq1),len(seq2))               # always will be last row and column for global alignment
    #dpTableV[maxIndexTuple[0]][maxIndexTuple[1]] should be alignment score
    #alignmentList=[]
    alignmentInfo=[]
    currentCoordinate=maxIndexTuple
    currentPointer=pointerTableV[currentCoordinate[0]][currentCoordinate[1]]
    alignmentInfo.append(alignmentInfoEvent('End',currentCoordinate[0],'')) 
    if currentCoordinate[1]==len(seq2):
        pass
    else:
        alignmentInfo.append(alignmentInfoEvent('Insertion',currentCoordinate[0],seq2[currentCoordinate[1]:]))
    relevantSeqEmit=seqEmitTableV[currentCoordinate[0]][currentCoordinate[1]]
    if relevantSeqEmit[3]=='deletion':
        alignmentInfo.append(alignmentInfoEvent('Deletion',currentCoordinate[0],''))
    elif relevantSeqEmit[2]=='insertion':
        alignmentInfo.append(alignmentInfoEvent('Insertion',currentCoordinate[0],relevantSeqEmit[1]))
    else:
        if relevantSeqEmit[0]=='N':
            alignmentInfo.append(alignmentInfoEvent('Diversity',currentCoordinate[0],relevantSeqEmit[1]))
        elif relevantSeqEmit[0]==relevantSeqEmit[1]:
            alignmentInfo.append(alignmentInfoEvent('Match',currentCoordinate[0],relevantSeqEmit[1]))
        else:
            alignmentInfo.append(alignmentInfoEvent('Mismatch',currentCoordinate[0],relevantSeqEmit[1]))   
    while currentPointer!=0:
        currentCoordinate=(currentPointer[0],currentPointer[1])
        if currentPointer[2]=='V':
            relevantSeqEmit=seqEmitTableV[currentCoordinate[0]][currentCoordinate[1]]
            currentPointer=pointerTableV[currentCoordinate[0]][currentCoordinate[1]]
        elif currentPointer[2]=='G':
            relevantSeqEmit=seqEmitTableG[currentCoordinate[0]][currentCoordinate[1]]
            currentPointer=pointerTableG[currentCoordinate[0]][currentCoordinate[1]]
        elif currentPointer[2]=='H':
            relevantSeqEmit=seqEmitTableH[currentCoordinate[0]][currentCoordinate[1]]
            currentPointer=pointerTableH[currentCoordinate[0]][currentCoordinate[1]]
        if relevantSeqEmit!=0:
            if relevantSeqEmit[3]=='deletion':
                alignmentInfo.append(alignmentInfoEvent('Deletion',currentCoordinate[0],''))
            elif relevantSeqEmit[2]=='insertion':
                alignmentInfo.append(alignmentInfoEvent('Insertion',currentCoordinate[0],relevantSeqEmit[1]))
            else:
                if relevantSeqEmit[0]=='N':
                    alignmentInfo.append(alignmentInfoEvent('Diversity',currentCoordinate[0],relevantSeqEmit[1]))
                elif relevantSeqEmit[0]==relevantSeqEmit[1]:
                    alignmentInfo.append(alignmentInfoEvent('Match',currentCoordinate[0],relevantSeqEmit[1]))
                else:
                    alignmentInfo.append(alignmentInfoEvent('Mismatch',currentCoordinate[0],relevantSeqEmit[1]))
        else:
            if currentCoordinate[1]==0:
                pass
            else:
                alignmentInfo.append(alignmentInfoEvent('Insertion',currentCoordinate[0],seq2[:currentCoordinate[1]]))
            alignmentInfo.append(alignmentInfoEvent('Start',currentCoordinate[0],'')) 

    #alignmentList=alignmentList[::-1]
    alignmentInfo=alignmentInfo[::-1]

    #collapse insertions and deletions that occur together into one sequence    
    
    newAlignmentInfo=[]
    insertionFound=0
    insertCoord=0
    insertSeq=''    
    for eachEl in alignmentInfo:
        if eachEl.typeOfEvent=='Insertion':            
            insertionFound=1
            insertCoord=eachEl.position
            insertSeq=insertSeq+eachEl.notes
        else:
            if insertionFound==1:
                newAlignmentInfo.append(alignmentInfoEvent('Insertion',insertCoord,insertSeq))
                insertCoord=0
                insertSeq=''
                insertionFound=0
            newAlignmentInfo.append(eachEl)

    newAlignmentInfoDel=[]
    deletionFound=0
    deletionCoord=0
    for eachEl in newAlignmentInfo[::-1]:
        if eachEl.typeOfEvent=='Deletion':
            deletionFound=1
            deletionCoord=eachEl.position
        else:
            if deletionFound==1:
                newAlignmentInfoDel.append(alignmentInfoEvent('Deletion',deletionCoord,''))
                deletionCoord=0
                deletionFound=0
            newAlignmentInfoDel.append(eachEl)
    newAlignmentInfoDel=newAlignmentInfoDel[::-1]
    return [dpTableV[maxIndexTuple[0]][maxIndexTuple[1]],newAlignmentInfoDel]
    #printAlignmentInfoList(newAlignmentInfo)
    #print alignmentList
    #print ''

def returnNumberMutations(alignInfoList):
    numberMutations=0
    for eachEl in alignInfoList:
        if eachEl.typeOfEvent=='Start':
            startPosition=eachEl.position
        elif eachEl.typeOfEvent=='End':
            endPosition=eachEl.position
    for eachEl in alignInfoList:
        if eachEl.typeOfEvent=='Insertion' or eachEl.typeOfEvent=='Mismatch' or eachEl.typeOfEvent=='Deletion':
            numberMutations+=1
    return numberMutations

def check_symmetric(a, tol=1e-8):
    return np.allclose(a, a.T, atol=tol)

names1=[]
names2=[]
for i in range(1,13):
    names1.append('tf'+str(i))
    names2.append(str(i+8))

names=names1+names2

namesSeq=[]
actualSeq=[]
abundanceSeq=[]
threshold=5 #choose sequences greater than threshold percent abundance
for counter, eachEl in enumerate(names):
    with open(eachEl+'_trimmomaticP_thresholdTestDict_1_seqErrorCorrect_1_numberStoreDictVer2.pckl','rb') as f:
        numberStoreDict=cPickle.load(f)
    newRefDict={}
    for eachEl1 in numberStoreDict:
        if eachEl1!='Clashed' and eachEl1!='Unaligned':
            newRefDict[eachEl1]=numberStoreDict[eachEl1]
    sortedList=sorted(newRefDict,key=lambda x: np.sum(newRefDict[x]), reverse=True)
    for counter1, eachRef in enumerate(sortedList):
        if np.sum(newRefDict[eachRef])>=threshold:
            namesSeq.append(str(counter+1)+'.'+str(counter1+1))
            actualSeq.append(eachRef)
            abundanceSeq.append(newRefDict[eachRef])

namesSeq.append('Y')
actualSeq.append('GGAAAATTTCAAGATCAGGGCTTGAAATTTTCCAAAATTTCAAGCCCTGATCTTGAAATTTTGG')
distanceMatrix=np.zeros((len(actualSeq),len(actualSeq)))
print len(actualSeq)
for counter, seq1 in enumerate(actualSeq):    
    print counter
    for counter1, seq2 in enumerate(actualSeq):
        distanceMatrix[counter][counter1]=min([returnNumberMutations(needlemanWunschAffine(seq1,seq2)[1]),returnNumberMutations(needlemanWunschAffine(seq1,revComp(seq2))[1]),returnNumberMutations(needlemanWunschAffine(revComp(seq1),seq2)[1]),returnNumberMutations(needlemanWunschAffine(revComp(seq1),revComp(seq2))[1]),returnNumberMutations(needlemanWunschAffine(seq2,seq1)[1]),returnNumberMutations(needlemanWunschAffine(seq2,revComp(seq1))[1]),returnNumberMutations(needlemanWunschAffine(revComp(seq2),seq1)[1]),returnNumberMutations(needlemanWunschAffine(revComp(seq2),revComp(seq1))[1])])


print 'Distance Matrix is symmetric: '+str(check_symmetric(distanceMatrix))

with open('FinalDistanceMatrix_Fitch_20180315.txt','w') as f:
    f.write('{:10s}'.format(str(len(actualSeq))))
    f.write('\n')
    for counter, eachRow in enumerate(distanceMatrix):
        f.write('{:10s}'.format(namesSeq[counter]))
        for eachEl in eachRow:
            f.write('\t')
            f.write('{:4s}'.format(str(eachEl)))
        f.write('\n')
        
    
    

    
    





