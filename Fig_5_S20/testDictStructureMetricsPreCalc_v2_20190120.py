import cPickle
import sys
from moduleUsefulFunctions_20180215 import *
from random import randint
import subprocess
import forgi
import forgi.graph.bulge_graph as fgb

def returnElNotation(dotBracket):
    bg=fgb.BulgeGraph(dotbracket_str=dotBracket)
    return bg.to_element_string()
    
def structureScoring(dotBracket,elementRep,s=1,h=0,i=-4,f=-2,t=-2,m=0,fThresh=2,tThresh=2):
    #s is score for stems, h is score for hairpin loop, m for multiloop segment, f for five prime unpaired, t for three prime unpaired, i for bulged out nucleotides
    bg=fgb.BulgeGraph(dotbracket_str=dotBracket)
    numberStems=0
    for eachEl in bg.stem_iterator():
        numberStems+=1
    if numberStems==0:
            return [i*len(dotBracket),False] #the second False passed here is to indicate that one hairpin not found        
    if numberStems>1: # could either be due to bulges or multiple stems
        if dotBracket.rfind('(')<dotBracket.find(')'): #all opening bases of hairpin are before
            pass
        else:
            return [i*len(dotBracket),False] #the second False passed here is to indicate that one hairpin not found
    scoreForStructure=0
    dictElements={}
    for eachEl in elementRep:
        dictIncrement(eachEl,dictElements)
    for eachEl in dictElements:
        if eachEl=='s':
            scoreForStructure+=s*dictElements[eachEl]
        elif eachEl=='h':
            scoreForStructure+=h*dictElements[eachEl]
        elif eachEl=='m':
            scoreForStructure+=m*dictElements[eachEl]
        elif eachEl=='i':
            scoreForStructure+=i*dictElements[eachEl]
        elif eachEl=='f':
            if dictElements[eachEl]>fThresh: #don't add to score otherwise
                scoreForStructure+=f*(dictElements[eachEl]-fThresh)
        elif eachEl=='t':
            if dictElements[eachEl]>tThresh: #don't add to score otherwise
                scoreForStructure+=t*(dictElements[eachEl]-tThresh)
    return [scoreForStructure,True] #the second True passed here is to indicate that one hairpin found

def dinucleotideRequirement(elementRep,sequence,diNucSeq):
    #assume here that sequence already checked for one hairpin forming
    fiveFree=''
    threeFree=''
    for counter, eachEl in enumerate(elementRep):
        if eachEl=='f':
            fiveFree+=sequence[counter]
        elif eachEl=='t':
            threeFree+=sequence[counter]
    diNucReq=(diNucSeq in fiveFree) and (diNucSeq in threeFree)
    return diNucReq

def ATHomopolymer(sequence):
    if 'C' in sequence or 'G' in sequence:
        return False
    else:
        return True

sortWithPCRDuplicates=0 #parameter to dictate whether sorting is happening with or without subtraction of PCR duplicates, 0 means that each barcode is only counted once for each insert sequence, 1 means that each barcode is counted as many times as it is observed even if the insert sequence is the same
k=20 #used to define identity of a reference sequence

with open(sys.argv[1],'rb') as f:
    testDict=cPickle.load(f)

sampleName=sys.argv[1]
sampleName=sampleName[sampleName.rfind('/')+1:sampleName.rfind('_testDict.pckl')]

seqDict={}
for eachEntry in testDict:
    seqDict[eachEntry]=testDict[eachEntry][sortWithPCRDuplicates]

sortedDictionary=dictSorter(seqDict)
print len(seqDict)
numberSeqProcessed=0
seqCount=[]
oneHairpin=[]
seqLength=[]
diGReq=[]
diCReq=[]
isATstretch=[]
dotBracketVector=[]
hairpinScoreVector=[]

dictionaryDotBracket={}
with open('outputFastaWriterRNAfoldCaller_20190120_'+sampleName+'.txt','rU') as f:
    for counter, eachLine in enumerate(f):
        if counter%3==0:
            indexToProbe=int(eachLine[1:])
        if counter%3==2:
            dotBracket=eachLine.split()[0]
            dictionaryDotBracket[indexToProbe]=dotBracket

for counter, eachEl in enumerate(sortedDictionary):
    numberSeqProcessed+=1
    if numberSeqProcessed%10000==0:
        print 'Processed '+str(numberSeqProcessed)+' of '+str(len(seqDict))
    seqCount.append(seqDict[eachEl])
    seqLength.append(len(eachEl))
    isATstretch.append(ATHomopolymer(eachEl))
    if len(eachEl)>=k:
        dotBracketNotation=dictionaryDotBracket[counter]
        dotBracketVector.append(dotBracketNotation)

        elRepNotation=returnElNotation(dotBracketNotation)
        scoreResult=structureScoring(dotBracketNotation,elRepNotation)
        hairpinScoreVector.append(float(scoreResult[0])/float(len(eachEl)))
        oneHairpin.append(scoreResult[1])

        if scoreResult[1]==True:
            diGReq.append(dinucleotideRequirement(elRepNotation,eachEl,'GG'))
            diCReq.append(dinucleotideRequirement(elRepNotation,eachEl,'CC'))
        else:
            diGReq.append(False)
            diCReq.append(False)           
    else:
        dotBracketVector.append('')
        hairpinScoreVector.append(-5.0)
        oneHairpin.append(False)
        diGReq.append(False)
        diCReq.append(False)
        
combinedList=[seqCount,oneHairpin,seqLength,diGReq,diCReq,dotBracketVector,isATstretch,hairpinScoreVector]    

with open(sampleName+'_combinedList_structureMetricsPreCalc_v2.pckl','wb') as f:
    cPickle.dump(combinedList,f,protocol=cPickle.HIGHEST_PROTOCOL)


