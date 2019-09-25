from moduleUsefulFunctions_20180215 import *
from Bio.Blast import NCBIXML
import xml.etree.ElementTree as ET
import sys

def getPVal(Evalue):
    #probability of finding at least 1 HSP with score>=S
    return -1.0*np.log10(1-np.exp(-1.0*Evalue))

def fixWithJitter(x,y,threshold):
    toFix=0
    for eachTup in alreadyPlottedData:
        a=np.array([eachTup[0]-x,eachTup[1]-y])
        if np.linalg.norm(a,ord=2)<threshold:
            toFix=1
            break
    if toFix==1:
        angleRand=np.random.uniform(low=0.0,high=2.0*np.pi)
        newX=threshold*np.cos(angleRand)+x
        newY=threshold*np.sin(angleRand)+y      
        alreadyPlottedData[(newX,newY)]=0
        return [newX,newY]
    else:
        alreadyPlottedData[(x,y)]=0
        return [x,y]

seqNames=[]
seqActual=[]
seqNamesIndex={}
counter=0
with open('allSeq_toBlast_blastRunner_v1_20190221.fasta','rU') as f:
    for eachLine in f:
        if eachLine[0]=='>':
            seqNames.append(eachLine.strip()[1:])
            seqNamesIndex[eachLine.strip()[1:]]=counter
            counter+=1
        else:
            seqActual.append(eachLine.strip())

alScore_refseq=[]
alScore_mySeeds=[]
eval_refseq=[]
eval_mySeeds=[]

with open('blastresults_allSeq_refseq_outfmt5_20190221.xml','rU') as f:
    records=NCBIXML.parse(f)
    for eachRec in records:
        alScore_refseq.append(eachRec.alignments[0].hsps[0].score)
        eval_refseq.append(eachRec.alignments[0].hsps[0].expect)

with open('blastresults_allSeq_mySeeds_outfmt5_20190221.xml','rU') as f:
    records=NCBIXML.parse(f)
    for eachRec in records:
        alScore_mySeeds.append(eachRec.alignments[0].hsps[0].score)
        eval_mySeeds.append(eachRec.alignments[0].hsps[0].expect)

#samplesToPlot='1_3_5_7'
samplesToPlot='2_4_6_8'

listSamples=samplesToPlot.split('_')

alScore_refseq_sample=[]
alScore_mySeeds_sample=[]
pval_refseq_sample=[]
pval_mySeeds_sample=[]
colorArray=[]
colorDict={'1':tuple(np.array([255.0,102.0,0.0])/255.0),'5':tuple(np.array([255.0,102.0,0.0])/255.0),'3':tuple(np.array([0.0,0.0,255.0])/255.0),'7':tuple(np.array([0.0,0.0,255.0])/255.0),'2':tuple(np.array([0.0,142.0,0.0])/255.0),'6':tuple(np.array([0.0,142.0,0.0])/255.0),'4':tuple(np.array([186.0,157.0,0.0])/255.0),'8':tuple(np.array([186.0,157.0,0.0])/255.0)}
for eachSample in listSamples:
    with open(eachSample+'_afterMS2removal_finalRefCollapsing_v1_20190221_allLists.pckl','rb') as f:
        [xToPlot1,yToPlot1,listAllRefs1,listAllClusters1,listdiGRef1]=cPickle.load(f)
    for counter, eachSeq in enumerate(listdiGRef1):
        alScore_refseq_sample.append(alScore_refseq[seqNamesIndex[eachSample+'.'+str(counter+1)]])
        alScore_mySeeds_sample.append(alScore_mySeeds[seqNamesIndex[eachSample+'.'+str(counter+1)]])
        pval_refseq_sample.append(getPVal(eval_refseq[seqNamesIndex[eachSample+'.'+str(counter+1)]]))
        pval_mySeeds_sample.append(getPVal(eval_mySeeds[seqNamesIndex[eachSample+'.'+str(counter+1)]]))
        colorArray.append(colorDict[eachSample])

keySampleDetails={'2_4_6_8':'Positive','1_3_5_7': 'Negative'}

#minLimit=min([min(alScore_refseq_sample),min(alScore_mySeeds_sample)])-2
#maxLimit=max([max(alScore_refseq_sample),max(alScore_mySeeds_sample)])+2
minLimit=15 #same limits across all plots
maxLimit=45 #same limits across all plots
f,ax=plt.subplots(nrows=1,ncols=2)
f.set_size_inches((4,2))
alreadyPlottedData={}
ax[0].plot(np.arange(minLimit,maxLimit+1),np.arange(minLimit,maxLimit+1),'k--',linewidth=0.25)
for i in range(len(alScore_refseq_sample)):
    toPlotPoints=fixWithJitter(alScore_refseq_sample[i],alScore_mySeeds_sample[i],0.35)
    ax[0].plot(toPlotPoints[0],toPlotPoints[1],markerfacecolor=colorArray[i],marker='o',markersize=3,markeredgecolor='none',linestyle='')

ax[0].set_xlim(minLimit,maxLimit)
ax[0].set_ylim(minLimit,maxLimit)
ax[0].set_xlabel('Alignment score (Best hit to Other Species)')
ax[0].set_ylabel('Alignment score (Best hit to Input Seeds)')
ax[0].set_title('Sample '+samplesToPlot+': '+keySampleDetails[samplesToPlot]+', numberSpecies='+str(len(alScore_refseq_sample)))

#minLimit=min([min(pval_refseq_sample),min(pval_mySeeds_sample)])-1
#maxLimit=max([max(pval_refseq_sample),max(pval_mySeeds_sample)])+1
minLimit=-1 #same limits across all plots
maxLimit=14 #same limits across all plots
alreadyPlottedData={}
ax[1].plot(np.arange(minLimit,maxLimit+1),np.arange(minLimit,maxLimit+1),'k--',linewidth=0.25)
for i in range(len(pval_refseq_sample)):
    toPlotPoints=fixWithJitter(pval_refseq_sample[i],pval_mySeeds_sample[i],0.35)
    ax[1].plot(toPlotPoints[0],toPlotPoints[1],markerfacecolor=colorArray[i],marker='o',markersize=3,markeredgecolor='none',linestyle='')

ax[1].set_xlim(minLimit,maxLimit)
ax[1].set_ylim(minLimit,maxLimit)
ax[1].set_xlabel('p-val (Best hit to Other Species)')
ax[1].set_ylabel('p-val (Best hit to Input Seeds)')
#ax[1].set_title('Sample '+samplesToPlot+': '+keySampleDetails[samplesToPlot]+', numberSpecies='+str(len(alScore_refseq_sample)))



#plt.show()
plt.savefig(samplesToPlot+'_compileAndPlot_v3_20190415.pdf',dpi=300)


