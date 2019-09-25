from moduleUsefulFunctions_20180215 import *
plt.close()

def updateDictData(dictionaryToUpdate,sequence,counts):
    if sequence in dictionaryToUpdate:
        dictionaryToUpdate[sequence]+=counts
    else:
        dictionaryToUpdate[sequence]=counts

def elementsSame(list1,list2):
    for eachEl in list1:
        if eachEl not in list2:
            return False
    return True


with open(sys.argv[1],'rb') as f:
    dimerCounts=cPickle.load(f) #dimerDump_ver2.pckl file

sampleName=sys.argv[1]
sampleName=sampleName[:sampleName.rfind('_dimerDump_ver2.pckl')]

pcrDuplicates=0
firstHalf={}
secondHalf={}
allHalves={}

for eachComb in dimerCounts:
    firstHalfSeq=eachComb[0]
    secondHalfSeq=eachComb[1]
    updateDictData(firstHalf,firstHalfSeq,dimerCounts[eachComb][pcrDuplicates])
    updateDictData(secondHalf,secondHalfSeq,dimerCounts[eachComb][pcrDuplicates])
    updateDictData(allHalves,firstHalfSeq,dimerCounts[eachComb][pcrDuplicates])
    updateDictData(allHalves,secondHalfSeq,dimerCounts[eachComb][pcrDuplicates])

variantsOnHeatMap=10

sortedFirst=dictSorter(allHalves)
listHeatmap=sortedFirst[:variantsOnHeatMap]
for i in range(0,variantsOnHeatMap):
    print 'i: '+str(i)+', '+listHeatmap[i]
bigData=[]
for eachEl in listHeatmap:
    arrayData=[]
    for eachEl1 in listHeatmap:
        if (eachEl,eachEl1) in dimerCounts:
            arrayData.append((float(dimerCounts[(eachEl,eachEl1)][pcrDuplicates])/float(firstHalf[eachEl])))
        else:
            arrayData.append(0.0)
    bigData.append(arrayData)

bigData=np.array(bigData)

#cmap=plt.get_cmap('Greens')        
cmap=plt.get_cmap('YlGnBu')        
plt.matshow(bigData, cmap=cmap, vmin=0.0, vmax=1.0) 
plt.grid(which="minor",ls="-",lw=1.5)
cbar=plt.colorbar()
cbar.solids.set_edgecolor("face")
plt.savefig(sampleName+'_agreementHeatmap_20180413_ver2.pdf',dpi=300)
#plt.show()
    


    

