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

with open('AF_SOL_850_4_junctions_20180430_ver2.pckl','rb') as f:
    dimerPairDict=cPickle.load(f) #_junctions_20180430_ver2.pckl 

strandToPullOut='C'
diversityString='CGGCCG'

searchString=strandToPullOut+diversityString
counter=0
for eachTup in dictSorter(dimerPairDict[(searchString,searchString)].readSeqDict):
    if counter==10:
        break
    if eachTup[0]=='' and eachTup[2]=='': #get sequences without internal mutations
        print eachTup
        print dimerPairDict[(searchString,searchString)].readSeqDict[eachTup]
        print ''
        counter+=1

