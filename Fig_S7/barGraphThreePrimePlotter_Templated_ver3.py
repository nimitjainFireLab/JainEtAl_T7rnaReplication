import sys
from moduleUsefulFunctions_20180215 import *
plt.close()

def returnEndPosition(alignInfoList):
    for eachEl in alignInfoList:
        if eachEl.typeOfEvent=='End':
            return eachEl.position

def ThreePrimeClipSequence(alignInfoList):
    endPosition=returnEndPosition(alignInfoList)
    for eachEl in alignInfoList:
        if eachEl.typeOfEvent=='Insertion' and eachEl.position==endPosition:
            return eachEl.notes
    return ''

colors=["#9b59b6", "#3498db", "#95a5a6", "#e74c3c", (211.0/255.0,161.0/255.0,0.0/255.0), "#2ecc71"] #palette chosen from https://seaborn.pydata.org/tutorial/color_palettes.html, slightly modified
pcrDuplicates=1 #this chosen because otherwise lose lot of reads for RNA oligo
numberDataToPlot=6

reactionsToProcess=['18','5','19','9']

referenceY2G='GGAAAATTTCAAGATCAGGGCTTGAAATTTTACAAAATTTCAAGCCCTGATCTTGAAATTTTGG'
#refSeqExtraBase=[referenceY2G,referenceY2G,revComp(referenceY2G),revComp(referenceY2G)]

dataBarPlot=[]
for i in range(0,numberDataToPlot): # there are 6 things to plot in stacked bar chart in order, end 0 with -, end 0 with A, end 0 with C, end 0 with G, end 0 with T, end 0 with more than 1 extra base
    dataBarPlot.append([])

gStrandNumbers=[]
cStrandNumbers=[]
for counter, eachEl in enumerate(reactionsToProcess):
    with open(eachEl+'_trimmomaticP_thresholdTestDict_1_seqErrorCorrect_1_numberStoreDict_alignmentTemplatedThreePrimeBases.pckl','rb') as f:
        numberStoreDict=cPickle.load(f)
    gStrandNumbers.append(numberStoreDict[referenceY2G][0])
    cStrandNumbers.append(numberStoreDict[referenceY2G][1])
    with open(eachEl+'_trimmomaticP_thresholdTestDict_1_seqErrorCorrect_1_newObjDict.pckl','rb') as f: #should only contain two sequences
        newObjDict=cPickle.load(f)
    for eachSeq in [referenceY2G, revComp(referenceY2G)]:
        dataForSeq=np.zeros(numberDataToPlot)
        if eachSeq in newObjDict:
            total=0.0
            refLength=len(eachSeq)
            for newObj in newObjDict[eachSeq]:
                endPos=returnEndPosition(newObj.alignmentInfo[1])
                if endPos==refLength:
                    total+=newObj.count[pcrDuplicates]
            print total
            for newObj in newObjDict[eachSeq]:       
                endPos=returnEndPosition(newObj.alignmentInfo[1])
                if endPos==refLength:                        
                    insertionSeq=ThreePrimeClipSequence(newObj.alignmentInfo[1])
                    if insertionSeq=='':
                        dataForSeq[0]+=(newObj.count[pcrDuplicates]/total)*100
                    elif insertionSeq=='A':
                        dataForSeq[1]+=(newObj.count[pcrDuplicates]/total)*100
                    elif insertionSeq=='C':
                        dataForSeq[2]+=(newObj.count[pcrDuplicates]/total)*100
                    elif insertionSeq=='G':
                        dataForSeq[3]+=(newObj.count[pcrDuplicates]/total)*100
                    elif insertionSeq=='T':
                        dataForSeq[4]+=(newObj.count[pcrDuplicates]/total)*100
                    else:
                        dataForSeq[5]+=(newObj.count[pcrDuplicates]/total)*100           
        for i in range(0,numberDataToPlot):
            dataBarPlot[i].append(dataForSeq[i])

f, axes=plt.subplots(nrows=1,ncols=2)

axes[0].bar(range(1,1+len(reactionsToProcess)),gStrandNumbers,color='black',align='center')
axes[0].bar(range(1,1+len(reactionsToProcess)),cStrandNumbers,bottom=gStrandNumbers,color='gray',align='center')
axes[0].set_ylim([0,105])
axes[0].set_xlim([0.5,len(reactionsToProcess)+0.5])

bottomVal=np.zeros(len(reactionsToProcess)*2)        
for counter, eachArray in enumerate(dataBarPlot):
    #axes[1].bar(range(1,1+len(reactionsToProcess)),eachArray,bottom=bottomVal,color=np.array(colors[counter])/255.0,align='center')
    axes[1].bar(range(1,1+len(reactionsToProcess)*2),eachArray,bottom=bottomVal,color=colors[counter],align='center')
    bottomVal+=np.array(eachArray)
axes[1].set_ylim([0,105])
axes[1].set_xlim([0.5,len(reactionsToProcess)*2+0.5])
#plt.show()
plt.savefig('barGraphThreePrimePlotter_Templated_ver3_barPlot_20180323.pdf',dpi=300)
    
            
            
            
        
        
                    

            
                

    
