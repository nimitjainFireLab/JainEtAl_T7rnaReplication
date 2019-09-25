from moduleUsefulFunctions_20180215 import *
import subprocess

dictFirstTwoChars={}
#taxIDDict={}
giList=[]
#Note that under coli taxID--562, lots of plasmid sequences are also being included--good idea to exclude taxID 562 also because of matches to lambda and the pPD122.03 plasmid I used
#Note that lambda's taxID is 10710--and it is present in the ref_seq genomic database
#C. elegans taxID is 6239
#C. brenneri taxID is 135651
#C. remanei taxID is 31234
#S. cerevisiae taxID is 559292
# There are some other coli associated taxIDs--741093, 469008, 585396, 585395, 595496, 155864, 1048689, 536056, 511145, 316385, 316407, 1133853, 1134782, 444450, 866768, 544404, 1274814, 386585
# Some phages that are similar to lambda are still in database--decided not to exclude them at the current time
organismsOfInterest=['10710','6239','135651','31234','559292','562','741093', '469008', '585396', '585395', '595496', '155864', '1048689', '536056', '511145', '316385', '316407', '1133853', '1134782', '444450', '866768', '544404', '1274814','386585']
with open('allRefSeqTaxData.txt','rU') as f:
    for counter, eachLine in enumerate(f):
        a=eachLine.strip().split('\\t')
        dictIncrement(a[0][:2],dictFirstTwoChars)
        if a[0][:2]=='NC':
            #dictIncrement(a[-3],taxIDDict)
            if a[-3] not in organismsOfInterest:
                giList.append(a[1])

with open('giListInclude_refseqgenomic.txt','w') as wf:
    for eachEl in giList:
        wf.write(eachEl+'\n')




    
        
