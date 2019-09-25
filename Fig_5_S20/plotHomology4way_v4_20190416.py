from moduleUsefulFunctions_20180215 import *
from random import randint
import subprocess
import forgi
import forgi.graph.bulge_graph as fgb
import time
from Bio.Blast import NCBIXML

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

def buildLocationFlanks(fourWayCoordinates,LengthSeq):
    list1=[]
    list1.append([0,fourWayCoordinates[0]])
    list1.append([fourWayCoordinates[0],fourWayCoordinates[1]+1])
    list1.append([fourWayCoordinates[1]+1,fourWayCoordinates[2]])
    list1.append([fourWayCoordinates[2],fourWayCoordinates[3]+1])
    list1.append([fourWayCoordinates[3]+1,fourWayCoordinates[4]])
    list1.append([fourWayCoordinates[4],fourWayCoordinates[5]+1])
    list1.append([fourWayCoordinates[5]+1,fourWayCoordinates[6]])
    list1.append([fourWayCoordinates[6],fourWayCoordinates[7]+1])
    list1.append([fourWayCoordinates[7]+1,LengthSeq])
    #note that some of the 2 coordinates in a list may be same--for example, fourWayCoordinates[3]+1 and fourWayCoordinates[4]--in this case a terminus will be mapped to the (fourWayCoordinates[4],fourWayCoordinates[5]+1) interval
    return list1

def returnCoordinates(terminusCoord,flankingLocs):
    #determine which flanks are there for terminus
    for counter, eachList in enumerate(flankingLocs):
        if terminusCoord>=eachList[0] and terminusCoord<eachList[1]:
            bestFlank=counter
            break
    coordinateToReturn=0.0
    for counter, eachEl in enumerate(locationLengths):
        if counter==bestFlank:
            break
        coordinateToReturn+=eachEl

    for counter1, i in enumerate(np.arange(flankingLocs[bestFlank][0],flankingLocs[bestFlank][1])):
        if i==terminusCoord:
            bestCounter1=counter1
            break

    getFullList=np.linspace(0,locationLengths[bestFlank],2+flankingLocs[bestFlank][1]-flankingLocs[bestFlank][0])
    coordinateToReturn+=getFullList[1+bestCounter1]
    return coordinateToReturn

with open('structureCalc2way4wayStoreData_v2_20190416_allLists.pckl','rb') as f:
    [namesSeq,actualSeq,twoWayRepeatCoord,fourWayRepeatCoord,blastCoord]=cPickle.load(f)

#Do all blastCoord go from 5' to 3' end?
for counter, eachEl in enumerate(blastCoord):
    if eachEl[0]<eachEl[-1]:
        pass
    else:
        print 'What??'



#9 locations: 1) before 1st 4-way repeat, 2) 1st 4-way repeat, 3) between 1st and 2nd 4-way repeat, 4) 2nd 4-way repeat, 5) between 2nd and 3rd 4-way repeat, 6) 3rd 4-way repeat, 7) between 3rd and 4th 4-way repeat, 8) 4th 4-way repeat, 9) after 4th 4-way repeat
locationLengths=[10,20,20,20,10,20,20,20,10] # in order--note the symmetry of this vector from left to right is important for the termini flipping in the code below

#samplesToUse=['1','3','5','7']
samplesToUse=['2','4','6','8']
rnasExcluded=0
terminus1List=[]
terminus2List=[]
thresholdLengthMatch=26
binWidth=5.0
for counter, eachName in enumerate(namesSeq):
    if eachName[:eachName.find('.')] in samplesToUse:
        if fourWayRepeatCoord[counter]==[]:
            rnasExcluded+=1
        else:
            if len(blastCoord[counter])>=thresholdLengthMatch:
                locationFlanks=buildLocationFlanks(fourWayRepeatCoord[counter],len(actualSeq[counter]))
                terminus1=returnCoordinates(blastCoord[counter][0],locationFlanks)
                terminus2=returnCoordinates(blastCoord[counter][-1],locationFlanks)
                #flip termini if needed, the next four lines can be commented out if flipping is not desired--flipping is only valid because of the left-right symmetry of locationLengths
                if (sum(locationLengths)-terminus2)<terminus1:
                    terminusTemp=terminus1
                    terminus1=sum(locationLengths)-terminus2
                    terminus2=sum(locationLengths)-terminusTemp
                terminus1List.append(terminus1)
                terminus2List.append(terminus2)

print len(terminus1List)
print 'RNAs without determined 4-way repeats: '+str(rnasExcluded)
f,ax=plt.subplots(nrows=1,ncols=1)
f.set_size_inches((3,2))

bins=np.arange(binWidth/2.0,sum(locationLengths),binWidth)
lowerBins=bins-binWidth/2.0
upperBins=bins+binWidth/2.0
density1=np.zeros(len(bins))
density2=np.zeros(len(bins))
for eachEl in terminus1List:
    for counter, eachEl1 in enumerate(bins):
        if eachEl>=lowerBins[counter] and eachEl<upperBins[counter]:
            density1[counter]+=1
            break

for eachEl in terminus2List:
    for counter, eachEl1 in enumerate(bins):
        if eachEl>=lowerBins[counter] and eachEl<upperBins[counter]:
            density2[counter]+=1
            break


density1=density1/np.sum(density1)
density2=density2/np.sum(density2)

max1=np.argsort(density1)[-1]
max2=np.argsort(density2)[-1]

upperLimPlot=0.8
for eachEl in np.cumsum(locationLengths):
    plt.plot([eachEl]*len(np.arange(0,upperLimPlot,0.01)),np.arange(0,upperLimPlot,0.01),'k--')

plt.plot([bins[max1]]*len(np.arange(0,density1[max1],0.01)),np.arange(0,density1[max1],0.01),'k--')
plt.plot([bins[max2]]*len(np.arange(0,density2[max2],0.01)),np.arange(0,density2[max2],0.01),'k--')
#plt.plot(bins,density1,'o-',markersize=5,color=(1.0,0.0,0.0))
#plt.plot(bins,density2,'o-',markersize=5,color=tuple(np.array([255.0,75.0,139.0])/255.0))
plt.bar(bins,density1,width=binWidth,align='center',color=(1.0,1.0,1.0))
plt.bar(bins,density2,width=binWidth,align='center',color=(0.5,0.5,0.5))

plt.xlim(0,sum(locationLengths))
plt.ylim(0,upperLimPlot)
plt.savefig('plot1dHomology4way_samples2468_v4_20190416.pdf',dpi=300)
        
   
