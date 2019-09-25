from moduleUsefulFunctions_20180215 import *

refList=['GGAAAATTTCAAACATTATGTTGTAATTTGTTTGAAAATTTCAAACAAATTACAACATAATGTTTGAAATTTTGGGGGGAAAAT','CCCAATATCATCAATTGCTGACGAAGATGATATTGATAATATCATCTTCGTCAGCAATTGATGATATT','GGAAAATCAATGACTGGTCAATCTCATTGATTTTTGAAATCAATGAGATTGACCAGTCATTGATTTT'] # in order from 13, 16 to 18

outputNames=['13','16','18']
inputNames=['1','9','11']

stripPlotList=[]
for counter, eachName in enumerate(outputNames):
    with open(eachName+'_numberStoreDictVer2.pckl','rb') as f:
        numberStoreDict=cPickle.load(f)
    print refList[counter]+'\t'+str(np.sum(numberStoreDict[refList[counter]]))
    inputList=[]
    for counter1, eachName1 in enumerate(inputNames):
        with open(eachName1+'_numberStoreDictVer2.pckl','rb') as f1:
            numberStoreDict1=cPickle.load(f1)
        inputList.append(np.sum(numberStoreDict1[refList[counter]]))
    print inputList
    normalizedInputList=np.array(inputList)/sum(inputList)
    print normalizedInputList
    stripPlotList.append(np.array(normalizedInputList).reshape((1,3)))

f, ax = plt.subplots(nrows=3,ncols=1)
f.set_size_inches((3,4))

cmap=plt.get_cmap('Greens')

for counter, eachEl in enumerate(stripPlotList):
    ax[counter].matshow(eachEl,cmap=cmap,vmin=0.0,vmax=1.0)
    ax[counter].set_xticks([x-0.5 for x in range(1,len(inputNames)+1)],minor=True )
    ax[counter].grid(which="minor",ls="-",lw=0.5)

plt.savefig('mainOut_plotAlignmentQuantification_v1_20190503.pdf',dpi=300)

f, ax=plt.subplots(nrows=1,ncols=1)
f.set_size_inches((3,4))

imageVar=ax.matshow(np.array([0.1,0.4,0.5]).reshape((1,3)),cmap=cmap,vmin=0.0,vmax=1.0) #this plot is just for getting colorbar, dummy data plotted
cbar=f.colorbar(imageVar)
cbar.solids.set_edgecolor("face")
plt.savefig('colorbarOut_plotAlignmentQuantification_v1_20190503.pdf',dpi=300)
        
        
