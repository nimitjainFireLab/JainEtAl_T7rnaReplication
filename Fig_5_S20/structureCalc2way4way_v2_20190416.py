from moduleUsefulFunctions_20180215 import *
from random import randint
import subprocess
import forgi
import forgi.graph.bulge_graph as fgb
import time
from Bio.Blast import NCBIXML

def returnDotBracket(fileName):
    with open(fileName,'rU') as f:
        for counter, eachLine in enumerate(f):
            if counter==1:
                dotBracket=eachLine.split()[0]
    return dotBracket

def returnElementNotation(sequence):
    with open('dummy_structureCalc2way4way_v2_20190416.txt','w') as f:
        processCall=subprocess.Popen(['RNAfold','--noPS','--noGU'],stdout=f,stdin=subprocess.PIPE)
        processCall.communicate(input=sequence)
    dotBracket=returnDotBracket('dummy_structureCalc2way4way_v2_20190416.txt')
    bg=fgb.BulgeGraph(dotbracket_str=dotBracket)
    return bg.to_element_string()

def extendTwoWayRepeat(elString,sequence):
    listElString=list(elString)
    #assume that a stem will be found as I curate replicating RNAs based on this
    curr_s5prime=elString.find('s')
    curr_s3prime=elString.rfind('s')
    curr_h5prime=elString.find('h')-1
    curr_h3prime=elString.rfind('h')+1
    possibleIndices=range(0,len(sequence))
    while (curr_s5prime-1) in possibleIndices and (curr_s3prime+1) in possibleIndices:
        if sequence[curr_s5prime-1]==revComp(sequence[curr_s3prime+1]): #extend
            curr_s5prime=curr_s5prime-1
            curr_s3prime-curr_s3prime+1
            listElString[curr_s5prime]='s'
            listElString[curr_s3prime]='s'
        else:
            break
    if curr_h5prime>0 and curr_h3prime>0: #takes care of case that no hairpin loop there by any chance
        while (curr_h3prime-1)>(curr_h5prime+1):
            if sequence[curr_h5prime+1]==revComp(sequence[curr_h3prime-1]): #extend
                curr_h5prime=curr_h5prime+1
                curr_h3prime=curr_h3prime-1
                listElString[curr_h5prime]='s'
                listElString[curr_h3prime]='s'
            else:
                break
    return [''.join(listElString),curr_s5prime,curr_s3prime,curr_h5prime,curr_h3prime] #take care here--the resultant element string may have lost all hairpin loop nucleotides!

def extendFourWayRepeat(i11,i12,i21,i22,i31,i32,i41,i42,sequence):
    possibleIndices=range(0,len(sequence))
    while (i11-1) in possibleIndices and (i42+1) in possibleIndices and (i22+1)<(i31-1):
        if sequence[i11-1]==revComp(sequence[i42+1]) and sequence[i22+1]==revComp(sequence[i11-1]) and sequence[i31-1]==revComp(sequence[i42+1]):
            i11=i11-1
            i22=i22+1
            i31=i31-1
            i42=i42+1
        else:
            break
    while (i12+1)<(i21-1) and (i32+1)<(i41-1):
        if sequence[i12+1]==revComp(sequence[i21-1]) and sequence[i32+1]==revComp(sequence[i41-1]) and sequence[i21-1]==revComp(sequence[i32+1]):
            i12=i12+1
            i21=i21-1
            i32=i32+1
            i41=i41-1
        else:
            break
    return [i11,i12,i21,i22,i31,i32,i41,i42]

def writeNewClustal(str1,str2,consStr):
    with open('dummyClustal_structureCalc2way4way_v2_20190416.aln','w') as wf:
        #note I am assuming that str1 and str2 will be below 60 characters, otherwise code below may break
        wf.write('CLUSTALW\n')
        wf.write('\n')
        wf.write('1\t'+str1+'\n')
        wf.write('2\t'+str2+'\n')
        wf.write(' \t'+consStr)
        wf.write('\n')    

def multipleStemsNotThere(dotBracket):
    bg=fgb.BulgeGraph(dotbracket_str=dotBracket)
    elementRep=bg.to_element_string()
    numberStems=0
    for eachEl in bg.stem_iterator():
        numberStems+=1
    if numberStems>1: # could either be due to bulges or multiple stems
        if dotBracket.rfind('(')<dotBracket.find(')'): #all opening bases of hairpin are before
            return True
        else:
            return False
    return True

def stemThere(dotBracket):
    bg=fgb.BulgeGraph(dotbracket_str=dotBracket)
    elementRep=bg.to_element_string()
    numberStems=0
    for eachEl in bg.stem_iterator():
        numberStems+=1
    if numberStems==0:
        return False
    else:
        return True


def multiAlWrite(elString,sequence,s5prime,s3prime,h5prime,h3prime):
    seq1=sequence[s5prime:h5prime+1]
    seq2=revComp(sequence[h3prime:s3prime+1])
    el1=elString[s5prime:h5prime+1]
    el2=elString[h3prime:s3prime+1]
    el2=el2[::-1]
    counter1=0
    counter2=0
    seq1ToWrite=''
    seq2ToWrite=''
    consString=''
    mappingDict1={}
    mappingDict2={}
    numberWritten=0
    #assume that there are no bulges at the ends of the stem loop
    while counter1!=len(seq1) and counter2!=len(seq2):
        mappingDict1[numberWritten]=counter1+s5prime
        mappingDict2[numberWritten]=s3prime-counter2
        if el1[counter1]=='s' and el2[counter2]=='s':
            seq1ToWrite+=seq1[counter1]
            seq2ToWrite+=seq2[counter2]
            consString+='*'
            counter1+=1
            counter2+=1
        elif el1[counter1]=='i' and el2[counter2]=='i':
            seq1ToWrite+=seq1[counter1]
            seq2ToWrite+=seq2[counter2]            
            consString+=' '
            counter1+=1
            counter2+=1            
        elif el1[counter1]=='s' and el2[counter2]=='i':
            seq1ToWrite+='-'
            seq2ToWrite+=seq2[counter2]
            consString+=' '
            counter2+=1
        elif el1[counter1]=='i' and el2[counter2]=='s':
            seq1ToWrite+=seq1[counter1]
            seq2ToWrite+='-'
            consString+=' '
            counter1+=1
        numberWritten+=1
    writeNewClustal(seq1ToWrite,seq2ToWrite,consString)
    return [seq1ToWrite,seq2ToWrite,mappingDict1,mappingDict2]

def isGU(index1,index2,sequence):
    if sequence[index1]=='G' and sequence[index2]=='T':
        return True
    elif sequence[index1]=='T' and sequence[index2]=='G':
        return True
    else:
        return False  
    
def removeGUWobble(seq1,seq2):
    with open('dummyAl_structureCalc2way4way_v2_20190416.txt','w') as f:
        processCall=subprocess.Popen('RNAalifold --noPS --noGU --temp=20 dummyClustal_structureCalc2way4way_v2_20190416.aln',shell=True, stdout=f,stdin=subprocess.PIPE) #temperature had to be lowered from default of 37C as otherwise some AT rich structures were not folding--note that this does change the folding for some RNAs as well
    time.sleep(1) #important otherwise dummyAl_structureCalc2way4way_v1_20190222.txt not created before being accessed
    dotBracket=returnDotBracket('dummyAl_structureCalc2way4way_v2_20190416.txt')
    bg=fgb.BulgeGraph(dotbracket_str=dotBracket)
    listDotBracket=list(dotBracket)
    for eachStem in bg.stem_iterator():
        for eachBp in bg.stem_bp_iterator(eachStem): #returns 1-index based numbering
            i1=eachBp[0]-1
            i2=eachBp[1]-1
            if isGU(i1,i2,seq1) and isGU(i1,i2,seq2):
                listDotBracket[i1]='.'
                listDotBracket[i2]='.'
    return ''.join(listDotBracket)  

def returnFourWayRepeatIndices(dotBracket,seq1,seq2,dict1,dict2,rnaSeq):
    #return a list giving four tuples--defining the four way repeat
    if multipleStemsNotThere(dotBracket) and stemThere(dotBracket):                
        listDotBracket=list(dotBracket)
        mappingDict={}
        counter1=0
        for counter, eachEl in enumerate(listDotBracket):
            if eachEl=='.':
                if seq1[counter]=='-' or seq2[counter]=='-':
                    listDotBracket[counter]=''
                    continue
            mappingDict[counter1]=counter
            counter1+=1
        bg=fgb.BulgeGraph(dotbracket_str=''.join(listDotBracket))
        maxLength=0
        for eachStem in bg.stem_iterator():
            if bg.stem_length(eachStem)>maxLength:
                maxLength=bg.stem_length(eachStem)
                bestStem=eachStem
        startNucPos=[]
        endNucPos=[]
        for eachBp in bg.stem_bp_iterator(bestStem):
            startNucPos.append(mappingDict[eachBp[0]-1])
            endNucPos.append(mappingDict[eachBp[1]-1])        
        i11=dict1[min(startNucPos)]
        i12=dict1[max(startNucPos)]
        i21=dict1[min(endNucPos)]
        i22=dict1[max(endNucPos)]
        i31=dict2[max(endNucPos)]
        i32=dict2[min(endNucPos)]
        i41=dict2[max(startNucPos)]
        i42=dict2[min(startNucPos)]
    
        newIndices=extendFourWayRepeat(i11,i12,i21,i22,i31,i32,i41,i42,rnaSeq)
        i11=newIndices[0]
        i12=newIndices[1]
        i21=newIndices[2]
        i22=newIndices[3]
        i31=newIndices[4]
        i32=newIndices[5]
        i41=newIndices[6]
        i42=newIndices[7]
        #the four way repeats are rnaSeq[i11:i12+1], rnaSeq[i21:i22+1], rnaSeq[i31:i32+1] and rnaSeq[i41:i42+1]
        return [i11,i12,i21,i22,i31,i32,i41,i42]           
    else:
        return []

def printTwoWayRepeat(coordVector,sequence):
    print sequence[coordVector[0]:coordVector[2]+1]
    print sequence[coordVector[3]:coordVector[1]+1]

def printFourWayRepeat(coordVector,sequence):
    if coordVector==[]:
        print 'Nothing'
    else:
        print sequence[coordVector[0]:coordVector[1]+1]
        print sequence[coordVector[2]:coordVector[3]+1]
        print sequence[coordVector[4]:coordVector[5]+1]
        print sequence[coordVector[6]:coordVector[7]+1]
        
def plotRNA(rnaName,rnaSeq,twoWayCoord,fourWayCoord,rnaNumber,blastCoordList,axH,totalNumberRNA,spacePerRNA,xLeftAlign):
    currentX=xLeftAlign
    currentY=(totalNumberRNA-rnaNumber)*spacePerRNA+(spacePerRNA/10)
    twoWayY=currentY+4*(spacePerRNA/10)
    fourWayY=currentY+2*(spacePerRNA/10)
    axH.text(1,currentY,rnaName,family='monospace',fontsize=6,weight='bold')
    axH.text(6,currentY,str(len(blastCoordList)),family='monospace',fontsize=10,weight='bold',color='red')
    for counter, eachChar in enumerate(rnaSeq):
        if counter in blastCoordList:
            axH.text(currentX,currentY,eachChar,family='monospace',fontsize=6,weight='bold',color='red')
        else:
            axH.text(currentX,currentY,eachChar,family='monospace',fontsize=6,weight='bold')
        currentX+=1
    axH.arrow(twoWayCoord[0]+xLeftAlign,twoWayY,twoWayCoord[2]-twoWayCoord[0],0,length_includes_head=True,head_width=0.05, head_length=0.5)
    axH.arrow(twoWayCoord[1]+xLeftAlign,twoWayY,-1.0*(twoWayCoord[1]-twoWayCoord[3]),0,length_includes_head=True,head_width=0.05, head_length=0.5)
    if fourWayCoord!=[]:
        axH.arrow(fourWayCoord[0]+xLeftAlign,fourWayY,fourWayCoord[1]-fourWayCoord[0],0,length_includes_head=True,head_width=0.05,head_length=0.5)
        axH.arrow(fourWayCoord[3]+xLeftAlign,fourWayY,-1.0*(fourWayCoord[3]-fourWayCoord[2]),0,length_includes_head=True,head_width=0.05,head_length=0.5)
        axH.arrow(fourWayCoord[4]+xLeftAlign,fourWayY,fourWayCoord[5]-fourWayCoord[4],0,length_includes_head=True,head_width=0.05,head_length=0.5)
        axH.arrow(fourWayCoord[7]+xLeftAlign,fourWayY,-1.0*(fourWayCoord[7]-fourWayCoord[6]),0,length_includes_head=True,head_width=0.05,head_length=0.5)
    
        
keySampleDetails={'1':'No seed','2': 'Seed','3': 'Seed+DNase','4': 'Seed+NaOH','5':'No seed','6': 'Seed','7': 'Seed+DNase','8': 'Seed+NaOH', '9':'No seed','10': 'Seed','11': 'Seed+DNase','12': 'Seed+NaOH','13':'No seed','14': 'Seed','15': 'Seed+DNase','16': 'Seed+NaOH'}

namesSeq=[]
actualSeq=[]
twoWayRepeatCoord=[]
fourWayRepeatCoord=[]
blastCoord=[]
counter1=0

with open('blastresults_allSeq_mySeeds_outfmt5_20190221.xml','rU') as f:
    records=NCBIXML.parse(f)
    for eachRec in records:
        listToAdd=[]
        startPos=eachRec.alignments[0].hsps[0].query_start  #1-index
        startPos=startPos-1
        matchString=eachRec.alignments[0].hsps[0].match #match and query have same length
        queryString=eachRec.alignments[0].hsps[0].query
        counter3=0
        for counter, eachEl in enumerate(matchString):
            if queryString[counter]=='-': #deletion in query
                pass
            else:
                if eachEl=='|': #only matching bases
                    listToAdd.append(startPos+counter3)
                counter3+=1 #increment counter3 for insertion, mismatch or matching base in query
        blastCoord.append(listToAdd)

with open('allSeq_toBlast_blastRunner_v1_20190221.fasta','rU') as f:
    seqName=''
    sequence=''
    for counter, eachLine in enumerate(f):
        if counter%2==0:
            seqName=eachLine[eachLine.find('>')+1:].strip()
        if counter%2==1:
            sequence=eachLine.strip()
            namesSeq.append(seqName)
            actualSeq.append(sequence)            
            twoWayResults=extendTwoWayRepeat(returnElementNotation(sequence),sequence)
            twoWayRepeatCoord.append(twoWayResults[1:])
            resultantStrings=multiAlWrite(twoWayResults[0],sequence,twoWayResults[1],twoWayResults[2],twoWayResults[3],twoWayResults[4])            
            newDotBracket=removeGUWobble(resultantStrings[0],resultantStrings[1])
            fourWayRepeatCoord.append(returnFourWayRepeatIndices(newDotBracket,resultantStrings[0],resultantStrings[1],resultantStrings[2],resultantStrings[3],sequence))
            #printFourWayRepeat(fourWayRepeatCoord[counter1],sequence)
            counter1+=1
            

for counter, eachSeq in enumerate(namesSeq):
    print eachSeq
    print actualSeq[counter]
    printTwoWayRepeat(twoWayRepeatCoord[counter],actualSeq[counter])
    printFourWayRepeat(fourWayRepeatCoord[counter],actualSeq[counter])
    print ''
    print ''


for eachSample in [1,2,3,4,5,6,7,8,9,10,11,12,13,16]:
    sampleToPlot=str(eachSample)
    repRNAToPlot={}
    #extract all the rep RNAs for a sample
    for counter, eachEl in enumerate(namesSeq):
        if eachEl[:eachEl.find('.')]==sampleToPlot:
            index=int(eachEl[eachEl.find('.')+1:])
            repRNAToPlot[index]=counter

    actualSeq_sample=[]
    namesSeq_sample=[]
    twoWayRepeatCoord_sample=[]
    fourWayRepeatCoord_sample=[]
    blastCoord_sample=[]
    for i in range(1,max(repRNAToPlot.keys())+1):
        indexToPlot=repRNAToPlot[i]
        actualSeq_sample.append(actualSeq[indexToPlot])
        namesSeq_sample.append(namesSeq[indexToPlot])
        twoWayRepeatCoord_sample.append(twoWayRepeatCoord[indexToPlot])
        fourWayRepeatCoord_sample.append(fourWayRepeatCoord[indexToPlot])
        blastCoord_sample.append(blastCoord[indexToPlot])
     
    maxLengthRNA=0
    for eachEl in actualSeq_sample:
        if len(eachEl)>maxLengthRNA:
                maxLengthRNA=len(eachEl)

    f,ax=plt.subplots(nrows=1,ncols=1)
    
    plt.xticks([])
    plt.yticks([])
    rnaNumberToPlot=len(actualSeq_sample)
    f.set_size_inches((8,rnaNumberToPlot))
    spacePerRNA=1.0
    xLeftAlign=11
    #useful link for formatting text: https://matplotlib.org/users/text_props.html

    for counter, eachSeq in enumerate(actualSeq_sample):
        plotRNA(namesSeq_sample[counter],eachSeq,twoWayRepeatCoord_sample[counter],fourWayRepeatCoord_sample[counter],counter+1,blastCoord_sample[counter],ax,rnaNumberToPlot,spacePerRNA,xLeftAlign)
    ax.set_title('Sample '+sampleToPlot+': '+keySampleDetails[sampleToPlot])
            
    ax.axis([0,2*xLeftAlign+maxLengthRNA,0,rnaNumberToPlot*spacePerRNA])
    plt.savefig(sampleToPlot+'_2way4wayBlastHits_structureCalc2way4way_v2_20190416.pdf',dpi=300)

