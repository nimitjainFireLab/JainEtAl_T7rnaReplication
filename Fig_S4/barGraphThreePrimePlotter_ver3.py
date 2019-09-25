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

colors=["#95a5a6","#34495e",(0.73094965079251462, 0.83947713375091548, 0.92132257293252384),(0.32628989885835086, 0.61862362903707169, 0.80279893524506507),(0.73714726672453035, 0.89551711503197162, 0.7108343117377337), (0.33882353901863099, 0.71172627631355734, 0.40584391180206747)]
names1=[]
names2=[]
for i in range(1,13):
    names1.append('tf'+str(i))
    names2.append(str(i+8))
names=names1+names2

mappingDict={}
seqStoreList=[]
for counter, eachEl in enumerate(names):
    mappingDict[eachEl]=counter+1
    seqStoreList.append([])

with open('sequencesToBlast_20180315.fasta','rU') as f:
    for counter, eachLine in enumerate(f):
        lineInfo=eachLine.strip()
        if counter%2==0:
            seqName=lineInfo[1:].replace('-','.')
            reactionNumber=int(seqName[:seqName.find('.')])            
        elif counter%2==1:
            sequence=lineInfo
            seqStoreList[reactionNumber-1].append((seqName,sequence))

pcrDuplicates=0
numberDataToPlot=2
additionSeqName={0:'', 1:'`'}

reactionsToProcess=names

namesSequencesX=[]
dataBarPlot=[]
for i in range(0,numberDataToPlot): # there are 2 things to plot in stacked bar chart in order, ends at -1, -2 and last but without insertion; ends at -1, -2 and last but with insertion
    dataBarPlot.append([])
    
for eachEl in reactionsToProcess:
    with open(eachEl+'_trimmomaticP_thresholdTestDict_1_seqErrorCorrect_1_newObjDict.pckl','rb') as f:
        newObjDict=cPickle.load(f)
    for eachRef in seqStoreList[mappingDict[eachEl]-1]:
        refSequenceName=eachRef[0]
        refSequence=eachRef[1]        
        for counter, eachSeq in enumerate([refSequence, revComp(refSequence)]):
            dataForSeq=np.zeros(numberDataToPlot)
            namesSequencesX.append(refSequenceName+additionSeqName[counter])
            print refSequenceName+additionSeqName[counter]            
            if eachSeq in newObjDict:
                total=0.0
                refLength=len(eachSeq)
                #print 'Reference length: '+str(refLength)
                for newObj in newObjDict[eachSeq]:
                    total+=newObj.count[pcrDuplicates]
                '''
                Print array of end abundances
                
                endAbundance=np.zeros(refLength)
                for newObj in newObjDict[eachSeq]:
                    endPos=returnEndPosition(newObj.alignmentInfo[1])
                    endAbundance[endPos-1]+=(newObj.count[pcrDuplicates]/total)*100                
                for i, eachVal in enumerate(endAbundance):
                    print str(i)+': '+str(eachVal)
                print ''
                print ''
                '''
                for newObj in newObjDict[eachSeq]:
                    endPos=returnEndPosition(newObj.alignmentInfo[1])
                    insertionSeq=ThreePrimeClipSequence(newObj.alignmentInfo[1])
                    if endPos>=refLength-2:
                        if insertionSeq=='': #no insertion
                            dataForSeq[0]+=(newObj.count[pcrDuplicates]/total)*100                            
                        else:
                            dataForSeq[1]+=(newObj.count[pcrDuplicates]/total)*100
            for i in range(0,numberDataToPlot):
                dataBarPlot[i].append(dataForSeq[i])

bottomVal=np.zeros(len(namesSequencesX))        
for counter, eachArray in enumerate(dataBarPlot):
    #plt.bar(range(1,1+len(namesSequencesX)),eachArray,bottom=bottomVal,color=np.array(colors[counter])/255.0,align='center')
    plt.bar(range(1,1+len(namesSequencesX)),eachArray,bottom=bottomVal,color=colors[counter],align='center')
    bottomVal+=np.array(eachArray)
plt.gca().set_ylim([0,105])
plt.gca().set_xlim([0.5,len(namesSequencesX)+0.5])
plt.gca().set_xticks(range(1,1+len(namesSequencesX)))
plt.gca().set_xticklabels(namesSequencesX)
plt.savefig('barGraphThreePrimePlotter_ver3_barPlot.pdf',dpi=300)
    
            
            
            
        
        
                    

            
                

    
