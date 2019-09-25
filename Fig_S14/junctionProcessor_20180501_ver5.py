from moduleUsefulFunctions_20180215 import *
plt.close()

class dimerSixBaseObject:
    def __init__(self):
        self.readSeqDict={} #keys are tuples of the form ('Mutations in 1st half','Dimer junction','Mutations in 2nd half','Three Prime Base addition')

def junctionAgreement(junctionList,threePrimeList):
    total=0.0
    for counter, eachEl in enumerate(threePrimeList):
        if len(junctionList[counter])>=len(eachEl):
            if junctionList[counter][:len(eachEl)]==eachEl:
                total+=1.0
    return total

def mutationAgreement(firstHalfList,secondHalfList):
    total=0.0
    for counter, eachEl in enumerate(firstHalfList):
        if eachEl==secondHalfList[counter] and eachEl!='': #don't count cases where there is no mutation
            total+=1.0
    return total    

with open(sys.argv[1],'rb') as f:
    dimerPairDict=cPickle.load(f) #_junctions_20180430_ver2.pckl 

numberPermutations=100000

numberThreePrimeBaseSeq=3

numberReadsObsAgreeMutations=0.0
numberReadsObsAgreeJunctions=0.0

numberReadsSimAgreeMutations=np.zeros(numberPermutations)
numberReadsSimAgreeJunctions=np.zeros(numberPermutations)
totalReads=0.0
for counter1, eachComb in enumerate(dimerPairDict):
    #lists 1-4 are for storing all info from the tuples in readSeqDict
    list1=[]
    list2=[]
    list3=[]
    list4=[]
    for eachTup in dimerPairDict[eachComb].readSeqDict:
        if len(eachTup[3])==numberThreePrimeBaseSeq:
            for i in range(0,dimerPairDict[eachComb].readSeqDict[eachTup]):
                list1.append(eachTup[0])
                list2.append(eachTup[1])
                list3.append(eachTup[2])
                list4.append(eachTup[3])
    numberReadsObsAgreeMutations+=mutationAgreement(list1,list3)
    numberReadsObsAgreeJunctions+=junctionAgreement(list2,list4)
    totalReads+=len(list1)
    for i in range(0,numberPermutations):
        simList1=[]
        simList2=[]
        simList3=[]
        simList4=[]
        permutedIndices=np.random.permutation(len(list1))
        for counter, j in enumerate(permutedIndices):
            simList1.append(list1[j])
            simList2.append(list2[j])
            simList3.append(list3[counter])
            simList4.append(list4[counter])
        numberReadsSimAgreeMutations[i]+=mutationAgreement(simList1,simList3)
        numberReadsSimAgreeJunctions[i]+=junctionAgreement(simList2,simList4)

sampleName=sys.argv[1]
sampleName=sampleName[:sampleName.find('_junctions')]

with open(sampleName+'_junctionProcessor_20180501_ver5.pckl','wb') as f:
    cPickle.dump([totalReads, numberReadsObsAgreeMutations, numberReadsObsAgreeJunctions, numberReadsSimAgreeMutations, numberReadsSimAgreeJunctions, numberPermutations, numberThreePrimeBaseSeq],f,protocol=cPickle.HIGHEST_PROTOCOL)

          


