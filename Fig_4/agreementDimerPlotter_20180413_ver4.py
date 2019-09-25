from moduleUsefulFunctions_20180215 import *
plt.close()

xRep1='AF_SOL_850_3_dimerDump_ver2.pckl'
xRep2='AF_SOL_850_4_dimerDump_ver2.pckl'
y2Rep1='AF_SOL_850_1_dimerDump_ver2.pckl'
y2Rep2='AF_SOL_850_2_dimerDump_ver2.pckl'

sampleNames=[xRep1,xRep2,y2Rep1,y2Rep2]

pcrDuplicates=0

resultsMatrix=np.zeros((len(sampleNames),2)) #1st column for agreeing dimers, 2nd column for all dimers

for counter, eachEl in enumerate(sampleNames):
    with open(eachEl, 'rb') as f:
        dimerCounts=cPickle.load(f)
    for eachComb in dimerCounts:
        resultsMatrix[counter,1]+=dimerCounts[eachComb][pcrDuplicates]
        if eachComb[0]==eachComb[1]:
            resultsMatrix[counter,0]+=dimerCounts[eachComb][pcrDuplicates]

print resultsMatrix

for counter, eachEl in enumerate(sampleNames):
    print eachEl+':\t'+str((resultsMatrix[counter,0]/resultsMatrix[counter,1])*100)
