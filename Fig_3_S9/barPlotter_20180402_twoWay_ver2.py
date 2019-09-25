from moduleUsefulFunctions_20180215 import *
import itertools
plt.close()

def returnStartPosition(alignInfoList):
    for eachEl in alignInfoList:
        if eachEl.typeOfEvent=='Start':
            return eachEl.position

def returnEndPosition(alignInfoList):
    for eachEl in alignInfoList:
        if eachEl.typeOfEvent=='End':
            return eachEl.position

def isFullLength(alignmentObject):
    refSeq=alignmentObject.refUsed
    startPos=returnStartPosition(alignmentObject.alignmentInfo[1])
    endPos=returnEndPosition(alignmentObject.alignmentInfo[1])
    if startPos==0 and endPos==len(refSeq):
        return True
    else:
        return False

def returnNumberMutations(alignInfoList):
    numberMutations=0
    for eachEl in alignInfoList:
        if eachEl.typeOfEvent=='Start':
            startPosition=eachEl.position
        elif eachEl.typeOfEvent=='End':
            endPosition=eachEl.position
    for eachEl in alignInfoList:
        if eachEl.typeOfEvent=='Insertion':
            if eachEl.position==startPosition or eachEl.position==endPosition: #only count insertions in the middle of the sequence as mutations
                pass
            else:
                numberMutations+=1 #as insertions are not being collapsed together (collapseInsertions=False in the smithWaterman calls below), each base inserted is counted separately
        elif eachEl.typeOfEvent=='Mismatch' or eachEl.typeOfEvent=='Deletion':
            numberMutations+=1
    return numberMutations

def diversityStringReturn(alignmentInfoList):
    stringDiversity=''
    for eachEl in alignmentInfoList:
        if eachEl.typeOfEvent=='Diversity':
            stringDiversity+=eachEl.notes
    return stringDiversity

def createEmptyDict(lenComb):
    dictToReturn={}
    allPerms=itertools.product(nucLetters.keys(),repeat=lenComb)
    for eachEl in allPerms:
        stringToAdd=''
        for eachEl1 in eachEl:
            stringToAdd+=eachEl1
        dictToReturn[stringToAdd]=0.0
    return dictToReturn
    
nucLetters={'A':0,'C':0,'G':0,'T':0}
pcrDuplicates=0
numberMutationLimit=1
thresholdBarToReport=1.0 #in percentage
toPlotComb=['AT','TA','CG','GC']
#For Watson Crick combinations--colors are in order red, blue, yellow, green
colors=["#e74c3c","#3498db",(211.0/255.0,161.0/255.0,0.0/255.0),"#2ecc71"]
grayColor="#95a5a6" 
#palette chosen from https://seaborn.pydata.org/tutorial/color_palettes.html, slightly modified
lenComb=len(toPlotComb[0])

newObjDictNames=[
    'AF_SOL_826_257_1_trimmomaticP_thresholdTestDict_1_seqErrorCorrect_1_newObjDict.pckl',
    'AF_SOL_826_257_9_trimmomaticP_thresholdTestDict_1_seqErrorCorrect_1_newObjDict.pckl',
    'AF_SOL_826_257_2_trimmomaticP_thresholdTestDict_1_seqErrorCorrect_1_newObjDict.pckl',
    'AF_SOL_826_257_10_trimmomaticP_thresholdTestDict_1_seqErrorCorrect_1_newObjDict.pckl',
    'AF_SOL_826_258_1_trimmomaticP_thresholdTestDict_1_seqErrorCorrect_1_newObjDict.pckl',
    'AF_SOL_826_258_7_trimmomaticP_thresholdTestDict_1_seqErrorCorrect_1_newObjDict.pckl',
    'AF_SOL_826_258_2_trimmomaticP_thresholdTestDict_1_seqErrorCorrect_1_newObjDict.pckl',
    'AF_SOL_826_258_8_trimmomaticP_thresholdTestDict_1_seqErrorCorrect_1_newObjDict.pckl',
    'AF_SOL_722_7_trimmomaticP_thresholdTestDict_1_seqErrorCorrect_1_newObjDict.pckl',
    'AF_SOL_722_11_trimmomaticP_thresholdTestDict_1_seqErrorCorrect_1_newObjDict.pckl',
    'AF_SOL_722_13_trimmomaticP_thresholdTestDict_1_seqErrorCorrect_1_newObjDict.pckl',
    'AF_SOL_722_17_trimmomaticP_thresholdTestDict_1_seqErrorCorrect_1_newObjDict.pckl',
    'AF_SOL_722_4_trimmomaticP_thresholdTestDict_1_seqErrorCorrect_1_newObjDict.pckl',
    'AF_SOL_713_Sample5_trimmomatic.fast_thresholdTestDict_1_seqErrorCorrect_1_newObjDict.pckl'
]

twoWayDiversityCoord=[[1,4],[1,4],[1,4],[1,4],[1,4],[1,4],[1,4],[1,4],[1,4],[1,4],[1,4],[1,4],[0,4],[0,4]]
numberDiversityPositions=[6,6,6,6,6,6,6,6,6,6,6,6,5,5]

f, axes=plt.subplots(nrows=1, ncols=1)

for counter, eachNewObjDictName in enumerate(newObjDictNames):
    print 'Sample: '+str(counter+1)
    with open(eachNewObjDictName,'rb') as f:
        newObjDict=cPickle.load(f)
    dataDict=createEmptyDict(lenComb)
    relevantCoord=twoWayDiversityCoord[counter]
    relevantDiversityPositions=numberDiversityPositions[counter]
    for eachRef in newObjDict:
        if eachRef[0]=='G': #use 2-way combinations as defined on G strand but sum for both G and C strands
            takeRevComp=0
        elif eachRef[0]=='C':
            takeRevComp=1
        for newObj in newObjDict[eachRef]:
            if isFullLength(newObj):
                divString=diversityStringReturn(newObj.alignmentInfo[1])
                if len(divString)==relevantDiversityPositions and returnNumberMutations(newObj.alignmentInfo[1])<numberMutationLimit:
                    twoWayComb=''
                    for eachEl in relevantCoord:
                        twoWayComb+=divString[eachEl]
                    if takeRevComp==0:
                        dataDict[twoWayComb]+=newObj.count[pcrDuplicates]
                    elif takeRevComp==1:
                        dataDict[revComp(twoWayComb)]+=newObj.count[pcrDuplicates]
    totalReads=float(sum(dataDict.values()))
    for eachEl in dataDict:
        dataDict[eachEl]=(float(dataDict[eachEl])/totalReads)*100
    bottomVal=0.0
    for i, eachEl in enumerate(toPlotComb):
        axes.bar(counter+1,dataDict[eachEl],bottom=bottomVal,color=colors[i], align='center')
        bottomVal+=dataDict[eachEl]
    sortedDataDict=dictSorter(dataDict)
    for eachEl in sortedDataDict:
        if eachEl not in toPlotComb:
            if dataDict[eachEl]>thresholdBarToReport:
                axes.bar(counter+1,dataDict[eachEl],bottom=bottomVal,color=grayColor,align='center')
                bottomVal+=dataDict[eachEl]
    axes.bar(counter+1,100.0-bottomVal,bottom=bottomVal,color='white',align='center')

axes.set_ylim([0,105])
axes.set_xlim([0.5,len(newObjDictNames)+0.5])
plt.savefig('barPlot_20180402_twoWay_ver2.pdf',dpi=300)
                       
                       
                       
                       
        

    
    

                
            
            
