from moduleUsefulFunctions_20180215 import *
plt.close()

with open(sys.argv[1],'rb') as f:
    data=cPickle.load(f)

numberReadsObsAgreeMutations=data[1]
numberReadsObsAgreeJunctions=data[2]
numberReadsSimAgreeMutations=data[3]
numberReadsSimAgreeJunctions=data[4]
numberPermutations=data[5]

f,ax=plt.subplots(nrows=1,ncols=1)

colors=["#9b59b6",(211.0/255.0,161.0/255.0,0.0/255.0)] #purple and yellow

ax.bar([1,4],[numberReadsObsAgreeJunctions,numberReadsObsAgreeMutations],color=colors[0],align='center')

ax.bar([2,5],[np.mean(numberReadsSimAgreeJunctions),np.mean(numberReadsSimAgreeMutations)],color=colors[1],align='center')
print numberPermutations
print np.percentile(numberReadsSimAgreeJunctions,q=2.5)
print np.percentile(numberReadsSimAgreeJunctions,q=97.5)
ax.set_xlim([0,6])
ax.set_ylabel('Read counts')
ax.set_ylim([0,35000])
#plt.show()
sampleName=sys.argv[1]
sampleName=sampleName[:sampleName.rfind('_junctionProcessor')]
plt.savefig(sampleName+'_junctionBarPlot_20180503_ver1.pdf',dpi=300)


