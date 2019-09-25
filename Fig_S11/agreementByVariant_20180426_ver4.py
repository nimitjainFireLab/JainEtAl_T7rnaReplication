from moduleUsefulFunctions_20180215 import *
import scikits.bootstrap as boot

plt.close()

def updateDictData(dictionaryToUpdate,sequence,counts):
    if sequence in dictionaryToUpdate:
        dictionaryToUpdate[sequence]+=counts
    else:
        dictionaryToUpdate[sequence]=counts

def listToDictConverter(listOfDimers):
    dictToReturn={}
    for eachEl in listOfDimers:
        if eachEl in dictToReturn:
            dictToReturn[eachEl]+=1
        else:
            dictToReturn[eachEl]=1
    return dictToReturn

def bootStrapDictReturn(indexList,dimerList):
    dictToReturn={}
    for eachEl in indexList:
        if dimerList[eachEl] in dictToReturn:
            dictToReturn[dimerList[eachEl]]+=1
        else:
            dictToReturn[dimerList[eachEl]]=1
    return dictToReturn

def getHomoDict(dictData):
    dictToReturn={}
    for eachEl in dictData:
        if eachEl[0]==eachEl[1]:
            dictToReturn[eachEl]=dictData[eachEl]
    return dictToReturn

def getHeteroDict(dictData):
    dictToReturn={}
    for eachEl in dictData:
        if eachEl[0]!=eachEl[1]:
            dictToReturn[eachEl]=dictData[eachEl]
    return dictToReturn

def getFirstHalf(dictData):
    dictToReturn={}
    for eachEl in dictData:
        if eachEl[0] in dictToReturn:
            dictToReturn[eachEl[0]]+=dictData[eachEl]
        else:
            dictToReturn[eachEl[0]]=dictData[eachEl]
    return dictToReturn

def getSecondHalf(dictData):
    dictToReturn={}
    for eachEl in dictData:
        if eachEl[1] in dictToReturn:
            dictToReturn[eachEl[1]]+=dictData[eachEl]
        else:
            dictToReturn[eachEl[1]]=dictData[eachEl]
    return dictToReturn

def getTotal(refName,dictData):
    total=0
    for eachEl in dictData:
        if eachEl[0][0]+eachEl[1][0]==refName:
            total+=dictData[eachEl]
    return total

def getCorrection(refName,firstHalfDict,secondHalfDict,totalsDict):
    #correction required because the total of hetero dicts does not include homodimers produced because of intermolecular template switching
    correction=1.0
    for eachEl in firstHalfDict:
        if eachEl in secondHalfDict:
            if eachEl[0]==refName[0]:
                correction=correction-(float(firstHalfDict[eachEl])/float(totalsDict[refName]))*(float(secondHalfDict[eachEl])/float(totalsDict[refName]))
    return correction

def getXRange(dataVector,expPower):
    xToUse=0.0
    xVector=[]
    for counter, eachEl in enumerate(dataVector):
        if counter==0:
            pass
        else:
            xToUse+=np.power(eachEl,expPower)
        xVector.append(xToUse)
    return xVector          

def numberDiffHalves(seq1,seq2):
    totalDiff=0
    for counter, eachEl in enumerate(seq1):
        if eachEl!=seq2[counter]:
            totalDiff+=1
    return totalDiff
            
with open(sys.argv[1],'rb') as f:
    dimerCounts=cPickle.load(f) #dimerDump_ver2.pckl file

def bootstrapFunction(data):
    #data is an indices list
    dictSimData=bootStrapDictReturn(data,listDimerNames)
    homoSimData=getHomoDict(dictSimData)
    heteroSimData=getHeteroDict(dictSimData)
    firstHalfWholeSimData=getFirstHalf(dictSimData)
    secondHalfWholeSimData=getSecondHalf(dictSimData)
    firstHalfHeteroSimData=getFirstHalf(heteroSimData)
    secondHalfHeteroSimData=getSecondHalf(heteroSimData)
    totalsWholeSim={}
    totalsHeteroSim={}
    totalsHeteroSubtractSim={}
    listToReturn=[]
    for eachEl in refsToUse:
        totalsWholeSim[eachEl]=getTotal(eachEl,dictSimData)
        totalsHeteroSim[eachEl]=getTotal(eachEl,heteroSimData)
        totalsHeteroSubtractSim[eachEl]=getCorrection(eachEl,firstHalfHeteroSimData,secondHalfHeteroSimData,totalsHeteroSim)

    for counter, eachEl in enumerate(dictSorter(homoOriginalExpDict)):
        refName=eachEl[0][0]+eachEl[1][0]
        if eachEl[0] in firstHalfWholeSimData and eachEl[1] in secondHalfWholeSimData:
            listToReturn.append((float(firstHalfWholeSimData[eachEl[0]])*float(secondHalfWholeSimData[eachEl[1]]))/float(totalsWholeSim[refName]))
        else:
            listToReturn.append(0.0)
    for counter, eachEl in enumerate(dictSorter(heteroOriginalExpDict)):
        refName=eachEl[0][0]+eachEl[1][0]
        if eachEl[0] in firstHalfHeteroSimData and eachEl[1] in secondHalfHeteroSimData:
            listToReturn.append(((float(firstHalfHeteroSimData[eachEl[0]])*float(secondHalfHeteroSimData[eachEl[1]]))/float(totalsHeteroSim[refName]))/float(totalsHeteroSubtractSim[refName]))
        else:
            listToReturn.append(0.0)
    return np.array(listToReturn)

sampleName=sys.argv[1]
sampleName=sampleName[:sampleName.rfind('_dimerDump_ver2.pckl')]

#first generate sample to bootstrap on, also restrict to GG and CC dimers
refsToUse=['GG','CC']
pcrDuplicates=0
listDimerNames=[]
for eachComb in dimerCounts:
    refType=eachComb[0][0]+eachComb[1][0]
    if refType in refsToUse:
        for i in range(0,dimerCounts[eachComb][pcrDuplicates]):
            listDimerNames.append(eachComb)

dictOriginalData=listToDictConverter(listDimerNames)
homoOriginalData=getHomoDict(dictOriginalData)
heteroOriginalData=getHeteroDict(dictOriginalData)
firstHalfWholeOriginalData=getFirstHalf(dictOriginalData)
secondHalfWholeOriginalData=getSecondHalf(dictOriginalData)
firstHalfHeteroOriginalData=getFirstHalf(heteroOriginalData)
secondHalfHeteroOriginalData=getSecondHalf(heteroOriginalData)
totalsWholeOriginal={}
totalsHeteroOriginal={}
totalsHeteroSubtractOriginal={}
for eachEl in refsToUse:
    totalsWholeOriginal[eachEl]=getTotal(eachEl,dictOriginalData)
    totalsHeteroOriginal[eachEl]=getTotal(eachEl,heteroOriginalData)
    totalsHeteroSubtractOriginal[eachEl]=getCorrection(eachEl,firstHalfHeteroOriginalData,secondHalfHeteroOriginalData,totalsHeteroOriginal)

f,ax=plt.subplots(nrows=1,ncols=2,figsize=(24,12))
#plot data and expected counts
obsHomoCounts=[]
expHomoCounts=[]
homoOriginalExpDict={}
for eachEl in homoOriginalData:
    refName=eachEl[0][0]+eachEl[1][0]
    homoOriginalExpDict[eachEl]=(float(firstHalfWholeOriginalData[eachEl[0]])*float(secondHalfWholeOriginalData[eachEl[1]]))/float(totalsWholeOriginal[refName])

for eachEl in dictSorter(homoOriginalExpDict):
    obsHomoCounts.append(homoOriginalData[eachEl])
    expHomoCounts.append(homoOriginalExpDict[eachEl])

obsHeteroCounts=[]
expHeteroCounts=[]
heteroOriginalExpDict={}

for eachEl in heteroOriginalData:
    refName=eachEl[0][0]+eachEl[1][0]
    heteroOriginalExpDict[eachEl]=((float(firstHalfHeteroOriginalData[eachEl[0]])*float(secondHalfHeteroOriginalData[eachEl[1]]))/float(totalsHeteroOriginal[refName]))/float(totalsHeteroSubtractOriginal[refName])

for eachEl in dictSorter(heteroOriginalExpDict):
    obsHeteroCounts.append(heteroOriginalData[eachEl])
    expHeteroCounts.append(heteroOriginalExpDict[eachEl])

#generate bootstrap data
confidenceIntervals=boot.ci(np.arange(len(listDimerNames)), statfunction=bootstrapFunction, alpha=0.05, n_samples=100000, method='bca', output='lowhigh')
#first row of confidenceIntervals is lowPercentile
#second row of confidenceIntervals is highPercentile

lowPercentile=confidenceIntervals[0,:]
highPercentile=confidenceIntervals[1,:]
obsHomoCounts=np.array(obsHomoCounts)
expHomoCounts=np.array(expHomoCounts)
obsHeteroCounts=np.array(obsHeteroCounts)
expHeteroCounts=np.array(expHeteroCounts)

#re-sort based on high percentile values
lowPercentileHomo=lowPercentile[:len(homoOriginalExpDict)]
lowPercentileHetero=lowPercentile[len(homoOriginalExpDict):]

highPercentileHomo=highPercentile[:len(homoOriginalExpDict)]
highPercentileHetero=highPercentile[len(homoOriginalExpDict):]

#homo sort
newIndexArray=np.argsort(highPercentileHomo)
newIndexArray=newIndexArray[::-1]
lowPercentileHomo=lowPercentileHomo[newIndexArray]
highPercentileHomo=highPercentileHomo[newIndexArray]
obsHomoCounts=obsHomoCounts[newIndexArray]
expHomoCounts=expHomoCounts[newIndexArray]

#hetero sort
newIndexArray=np.argsort(highPercentileHetero)
newIndexArray=newIndexArray[::-1]
lowPercentileHetero=lowPercentileHetero[newIndexArray]
highPercentileHetero=highPercentileHetero[newIndexArray]
obsHeteroCounts=obsHeteroCounts[newIndexArray]
expHeteroCounts=expHeteroCounts[newIndexArray]

firstXRange=getXRange(highPercentileHomo,0.25)
secondXRange=getXRange(highPercentileHetero,0.7)
ax[0].plot(firstXRange,obsHomoCounts,'bo',mec='b')
ax[0].fill_between(firstXRange,lowPercentileHomo,highPercentileHomo,color=(230.0/255.0,230.0/255.0,0.0/255.0))
ax[0].set_yscale('symlog',linthreshy=1.0)
ax[0].set_xlim([-1*(firstXRange[2]-firstXRange[1]),None])
ax[0].set_ylim([0,10**5])
ax[0].set_xlabel('Homodimer species')
ax[0].set_ylabel('Counts')
ax[1].plot(secondXRange,obsHeteroCounts,'bo',mec='b')
ax[1].fill_between(secondXRange,lowPercentileHetero,highPercentileHetero,color=(230.0/255.0,230.0/255.0,0.0/255.0))
ax[1].set_yscale('symlog',linthreshy=1.0)
ax[1].set_xlim([-1*(secondXRange[2]-secondXRange[1]),None])
ax[1].set_ylim([0,10**5])
ax[1].set_xlabel('Heterodimer species')
ax[1].set_ylabel('Counts')
#plt.show()
plt.savefig(sampleName+'_agreementByVariant_20180426_ver4.pdf',dpi=300)

    
    
