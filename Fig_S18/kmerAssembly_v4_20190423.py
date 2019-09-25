from moduleUsefulFunctions_20180215 import *
import networkx as nx
import pygraphviz as pgv
import matplotlib.pyplot as plt
import numpy as np
from networkx.drawing.nx_agraph import to_agraph
from scipy.interpolate import interp1d
import cPickle

sortWithPCRDuplicates=0

def returnTestDict(sampleName):
    with open(sampleName+'_trimmomaticP_thresholdTestDict_1_seqErrorCorrect_1_testDict.pckl','rb') as f:
        testDict=cPickle.load(f)
    seqDict={}
    for eachEntry in testDict:
        seqDict[eachEntry]=testDict[eachEntry][sortWithPCRDuplicates]
    sortedDictionary=dictSorter(seqDict)
    return [seqDict,sortedDictionary]

inputSample='8'
outputSample='15'
[inputTest,inputSorted]=returnTestDict(inputSample)
[outputTest,outputSorted]=returnTestDict(outputSample)

mergedTestDict={}
for eachSeq in inputTest:
    mergedTestDict[eachSeq]=inputTest[eachSeq]

for eachSeq in outputTest:
    if eachSeq in mergedTestDict:
        mergedTestDict[eachSeq]+=outputTest[eachSeq]
    else:
        mergedTestDict[eachSeq]=outputTest[eachSeq]

seqDict={}

with open('dummyR1_supercluster_v5_20190402.fastq','rU') as f:
    for counter, eachLine in enumerate(f):
        if counter%4==1:
            dictIncrement(eachLine.strip(),seqDict)

def returnRep(seq):
    seqRep=0
    revCompRep=0
    maxSeqRep=0
    maxRevCompRep=0
    for eachEl in mergedTestDict: #performs better than using seqDict here
        if seq in eachEl:
            seqRep+=mergedTestDict[eachEl]
            if mergedTestDict[eachEl]>maxSeqRep:
                maxSeqRep=mergedTestDict[eachEl]            
        if revComp(seq) in eachEl:
            revCompRep+=mergedTestDict[eachEl]
            if mergedTestDict[eachEl]>maxRevCompRep:
                maxRevCompRep=mergedTestDict[eachEl]
    return [seqRep,revCompRep,maxSeqRep,maxRevCompRep]        


def returnExtSequence(graphG,startingNode,startingSeq,direction,growthThresh):    
    currentNode=startingNode
    currentSeq=startingSeq
    addedSeq=''
    while 1:
        weightList=[]
        repList=[]
        repSeqList=[]
        repRevCompList=[]
        nodeList=[]
        seqList=[]
        print 'Grow to the '+direction+' now...'
        if direction=='left':
            iteratorList=G.predecessors(currentNode)
        elif direction=='right':
            iteratorList=G.successors(currentNode)
        for eachNode in iteratorList:
            if direction=='left':
                tempSeq=eachNode[0]+currentSeq
                weightList.append(G[eachNode][currentNode]['weight'])
            elif direction=='right':
                tempSeq=currentSeq+eachNode[-1]
                weightList.append(G[currentNode][eachNode]['weight'])
            resultRep=returnRep(tempSeq)
            repList.append(min(resultRep[2:]))
            repSeqList.append(resultRep[0])
            repRevCompList.append(resultRep[1])
            nodeList.append(eachNode)
            seqList.append(tempSeq)
        if len(weightList)==0:
            bestVal=-1
        else:
            tempWeights=[]
            tempRepSeqList=[]
            tempRepRevCompList=[]
            tempNodeList=[]
            tempSeqList=[]
            tempRepList=[]
            indexSorted=np.argsort(np.array(repList))
            bestRepVal=repList[indexSorted[-1]]
            for counter, eachEl in enumerate(repList):
                if eachEl==bestRepVal: 
                    tempWeights.append(weightList[counter])
                    tempRepSeqList.append(repSeqList[counter])
                    tempRepRevCompList.append(repRevCompList[counter])
                    tempNodeList.append(nodeList[counter])
                    tempSeqList.append(seqList[counter])
                    tempRepList.append(min([repSeqList[counter],repRevCompList[counter]]))
            indexSorted1=np.argsort(np.array(tempRepList))
            bestVal=tempRepList[indexSorted1[-1]]                               
        if bestVal<growthThresh:
            return addedSeq
        else:
            tempWeights2=[]
            tempRepSeqList2=[]
            tempRepRevCompList2=[]
            tempNodeList2=[]
            tempSeqList2=[]
            for counter, eachEl in enumerate(tempRepList):
                if eachEl==bestVal:
                    tempWeights2.append(tempWeights[counter])
                    tempRepSeqList2.append(tempRepSeqList[counter])
                    tempRepRevCompList2.append(tempRepRevCompList[counter])
                    tempNodeList2.append(tempNodeList[counter])
                    tempSeqList2.append(tempSeqList[counter])
            indexSorted2=np.argsort(np.array(tempWeights2))
            bestIndex=indexSorted2[-1]
            print tempSeqList2[bestIndex]+'\t'+str(tempRepSeqList2[bestIndex])+'\t'+str(tempRepRevCompList2[bestIndex])
            currentSeq=tempSeqList2[bestIndex]
            currentNode=tempNodeList2[bestIndex]
            if direction=='left':
                addedSeq=currentNode[0]+addedSeq
            elif direction=='right':
                addedSeq+=currentNode[-1]
                
'''
    if nx.algorithms.dag.is_directed_acyclic_graph(G):
        break
'''

def returnGraph(dictSequences,kmer):
    k=kmer
    seqDict=dictSequences
    G=nx.DiGraph() #make an empty directed graph, may not be ACYCLIC!
    for counter, eachSeq in enumerate(seqDict.keys()):
        if len(eachSeq)>=k:
            for i in range(0,len(eachSeq)+1-k,1):
                oneNode=eachSeq[i:i+k-1]
                secondNode=eachSeq[i+1:i+k]
                if G.has_edge(oneNode,secondNode):
                    G[oneNode][secondNode]['weight']+=seqDict[eachSeq]
                else:
                    G.add_edge(oneNode,secondNode,weight=seqDict[eachSeq])
    return G


k=20
G=returnGraph(seqDict,k)
print k
maxWeight=0
for eachEl in G.edges():
    if G[eachEl[0]][eachEl[1]]['weight']>maxWeight:
        maxWeight=G[eachEl[0]][eachEl[1]]['weight']
        maxEdge=(eachEl[0],eachEl[1])
        maxSeq=eachEl[0]+eachEl[1][-1]

print maxEdge
print maxWeight
print maxSeq

growthThreshold=float(maxWeight)/10.0 #counts

leftAdded=returnExtSequence(G,maxEdge[0],maxSeq,'left',growthThreshold)
rightAdded=returnExtSequence(G,maxEdge[1],leftAdded+maxSeq,'right',growthThreshold)

currentContig=leftAdded+maxSeq+rightAdded
#note that the currentContig can still be truncated because the graph is made using seqDict as opposed to mergedTestDict--extend the contig further using mergedTestDict now
newSeqDict={}
moreToLeft=[]
for eachEl in mergedTestDict: 
    if currentContig in eachEl or revComp(currentContig) in eachEl:
        newSeqDict[eachEl]=mergedTestDict[eachEl]
    if currentContig in eachEl:
        if eachEl.find(currentContig)>=(len(eachEl)-len(currentContig)-eachEl.find(currentContig)):
            moreToLeft.append(True)
        else:
            moreToLeft.append(False)

'''
for eachEl in dictSorter(newSeqDict):
    print eachEl+'\t'+str(newSeqDict[eachEl])
'''
k=len(currentContig)
G=returnGraph(newSeqDict,k)
print 'Is DAG: '+str(nx.algorithms.dag.is_directed_acyclic_graph(G)) #FOR CODE BELOW TO WORK, G DOESN'T HAVE TO BE ACYCLIC

maxEdge=(currentContig[:-1],currentContig[1:])
maxSeq=currentContig
newGrowthThreshold=5

#Decide which direction to grow in first based on location of currentContig in reads
if 2*sum(moreToLeft)>len(moreToLeft): #more bases are on left, extend left first
    leftAdded=returnExtSequence(G,maxEdge[0],maxSeq,'left',newGrowthThreshold)
    rightAdded=returnExtSequence(G,maxEdge[1],leftAdded+maxSeq,'right',newGrowthThreshold)    
else:
    rightAdded=returnExtSequence(G,maxEdge[1],maxSeq,'right',newGrowthThreshold)    
    leftAdded=returnExtSequence(G,maxEdge[0],maxSeq+rightAdded,'left',newGrowthThreshold)
    



    

    
