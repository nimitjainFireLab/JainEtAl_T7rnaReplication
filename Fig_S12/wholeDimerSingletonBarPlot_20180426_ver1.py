from moduleUsefulFunctions_20180215 import *
plt.close()

rep1=sys.argv[1]

with open(rep1,'rb') as f1:
    rep1Data=cPickle.load(f1)

agree=rep1Data[1]/rep1Data[4]

disagree=rep1Data[2]/rep1Data[4]

model=rep1Data[3]/rep1Data[4]
print 'Total reads'
print rep1Data[4]

#convert to percentage
agree=agree*100
disagree=disagree*100
model=model*100

f,ax=plt.subplots(nrows=1,ncols=1)

colors=["#e74c3c", "#3498db", "#2ecc71"] #red, blue and green in order
for i in range(0,3):
    if i==0:
        toPlot=agree
    elif i==1:
        toPlot=disagree
    elif i==2:
        toPlot=model
    halfWidth=len(toPlot)/2
    ax.bar(i+1,sum(toPlot[:halfWidth]),color=colors[i],align='center')
    ax.bar(i+5,sum(toPlot[halfWidth:]),color=colors[i],align='center')

ax.set_xlim([0,8])
ax.set_ylabel('Incidence [%]')
ax.set_ylim([0,35])
#plt.show()
sampleName=sys.argv[1]
sampleName=sampleName[:sampleName.rfind('.p')]
plt.savefig(sampleName+'_wholeDimerSingletonBarPlot_20180426_ver1.pdf',dpi=300)
