from moduleUsefulFunctions_20180215 import *
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

def FivePrimeClipSequence(alignInfoList):
    startPosition=returnStartPosition(alignInfoList)
    for eachEl in alignInfoList:
        if eachEl.typeOfEvent=='Insertion' and eachEl.position==startPosition:
            return eachEl.notes
    return ''

def ThreePrimeClipSequence(alignInfoList):
    endPosition=returnEndPosition(alignInfoList)
    for eachEl in alignInfoList:
        if eachEl.typeOfEvent=='Insertion' and eachEl.position==endPosition:
            return eachEl.notes
    return ''

def returnFivePrimeDict(newObjDictName,numberMutationLimit,pcrDuplicates):
    with open(newObjDictName,'rb') as f:
        newObjDict=cPickle.load(f)
    fivePrimeInsertions={}
    for eachRef in newObjDict:
        fivePrimeInsertions[eachRef]={}
        for newObj in newObjDict[eachRef]:
            if isFullLength(newObj) and returnNumberMutations(newObj.alignmentInfo[1])<numberMutationLimit:
                if ThreePrimeClipSequence(newObj.alignmentInfo[1])=='A': #restricting to 3' end being A
                    fivePrimeSeq=FivePrimeClipSequence(newObj.alignmentInfo[1])
                    if fivePrimeSeq in fivePrimeInsertions[eachRef]:
                        fivePrimeInsertions[eachRef][fivePrimeSeq]+=newObj.count[pcrDuplicates]
                    else:
                        fivePrimeInsertions[eachRef][fivePrimeSeq]=newObj.count[pcrDuplicates]    
    return fivePrimeInsertions

def returnCumList(dict1):
    toReturnCounts=np.zeros(len(listEl))
    for counter, eachEl in enumerate(listEl[:-1]):
        if eachEl in dict1:
            toReturnCounts[counter]+=float(dict1[eachEl])
    toReturnCounts=toReturnCounts/sum(dict1.values())
    toReturnCounts[-1]=1.0-np.sum(toReturnCounts)
    return toReturnCounts

pcrDuplicates=1
numberMutationLimit=1
listEl=['','A','C','G','T','AA','AC','AG','AT','CA','CC','CG','CT','GA','GC','GG','GT','TA','TC','TG','TT','Other']
colors1=[(211.0/255.0,161.0/255.0,0.0/255.0), "#3498db", "#95a5a6"]
colors2=["#e74c3c", "#9b59b6", "#2ecc71"] #palette chosen from https://seaborn.pydata.org/tutorial/color_palettes.html, slightly modified

newObjDictNames=[
    '18_trimmomaticP_thresholdTestDict_1_seqErrorCorrect_1_newObjDict.pckl',
    '5_trimmomaticP_thresholdTestDict_1_seqErrorCorrect_1_newObjDict.pckl',
    '19_trimmomaticP_thresholdTestDict_1_seqErrorCorrect_1_newObjDict.pckl',
    '9_trimmomaticP_thresholdTestDict_1_seqErrorCorrect_1_newObjDict.pckl'
]
oligo223=returnFivePrimeDict(newObjDictNames[0],numberMutationLimit,pcrDuplicates)
products223=returnFivePrimeDict(newObjDictNames[1],numberMutationLimit,pcrDuplicates)
oligo224=returnFivePrimeDict(newObjDictNames[2],numberMutationLimit,pcrDuplicates)
products224=returnFivePrimeDict(newObjDictNames[3],numberMutationLimit,pcrDuplicates)
keysList=products223.keys()
for eachEl in keysList:
    if eachEl[0]=='G':
        list1=returnCumList(oligo223[eachEl])
        list2=returnCumList(products223[eachEl])
        list3=returnCumList(products224[eachEl])
    else:
        list4=returnCumList(oligo224[eachEl])
        list5=returnCumList(products223[eachEl])
        list6=returnCumList(products224[eachEl])

f, ax=plt.subplots(nrows=2,ncols=1)
f.set_size_inches((4.25,2.5))

ax[0].plot(np.arange(1,len(listEl)+1),list1,color=colors1[0],label='223')
ax[0].plot(np.arange(1,len(listEl)+1),list2,color=colors1[1],label='Gprod223')
ax[0].plot(np.arange(1,len(listEl)+1),list3,color=colors1[2],label='Gprod224')

ax[1].plot(np.arange(1,len(listEl)+1),list4,color=colors2[0],label='224')
ax[1].plot(np.arange(1,len(listEl)+1),list5,color=colors2[1],label='Cprod223')
ax[1].plot(np.arange(1,len(listEl)+1),list6,color=colors2[2],label='Cprod224')

ax[0].set_xlim([0,len(listEl)+1])
ax[1].set_xlim([0,len(listEl)+1])

ax[0].set_xticks(np.arange(1,len(listEl)+1))
ax[1].set_xticks(np.arange(1,len(listEl)+1))
listElReplaced=[]
for eachEl in listEl:
    listElReplaced.append(eachEl.replace('T','U'))

ax[0].set_xticklabels(listElReplaced)
ax[1].set_xticklabels(listElReplaced)

ax[0].set_ylim([0,0.6])
ax[1].set_ylim([0,0.6])
ax[0].legend()
ax[1].legend()

plt.savefig('output_20190409_initiationNuc_v16Code.pdf',dpi=300)


