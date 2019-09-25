baseDir ='DemultiplexedFiles/'
baseDir1=''
readFile1 = baseDir1 + 'Undetermined_S0_L001_R1_001.fastq'
readFile2 = baseDir1 + 'Undetermined_S0_L001_R2_001.fastq'
indexFile1 = baseDir1 + 'Undetermined_S0_L001_I1_001.fastq'
indexFile2 = baseDir1 + 'Undetermined_S0_L001_I2_001.fastq'

D701='ATTACTCG'
D702='TCCGGAGA'
D703='CGCTCATT'
D704='GAGATTCC'
D705='ATTCAGAA'
D706='GAATTCGT'

D501='TATAGCCT'
D502='ATAGAGGC'
D503='CCTATCCT'
D504='GGCTCTGA'
D505='AGGCGAAG'
D506='TAATCTTA'
D507='CAGGACGT'
D508='GTACTGAC'

indexSeq = {(D706,D501):[open(baseDir+'1_R1.fastq','w'),open(baseDir+'1_R2.fastq','w')],
            (D705,D502):[open(baseDir+'2_R1.fastq','w'),open(baseDir+'2_R2.fastq','w')],
            (D704,D503):[open(baseDir+'3_R1.fastq','w'),open(baseDir+'3_R2.fastq','w')],
            (D703,D504):[open(baseDir+'4_R1.fastq','w'),open(baseDir+'4_R2.fastq','w')],            
            (D702,D505):[open(baseDir+'5_R1.fastq','w'),open(baseDir+'5_R2.fastq','w')],
            (D701,D506):[open(baseDir+'6_R1.fastq','w'),open(baseDir+'6_R2.fastq','w')],
            (D706,D507):[open(baseDir+'7_R1.fastq','w'),open(baseDir+'7_R2.fastq','w')],
            (D705,D501):[open(baseDir+'8_R1.fastq','w'),open(baseDir+'8_R2.fastq','w')],
            (D704,D502):[open(baseDir+'9_R1.fastq','w'),open(baseDir+'9_R2.fastq','w')],
            (D703,D503):[open(baseDir+'10_R1.fastq','w'),open(baseDir+'10_R2.fastq','w')],
            (D702,D504):[open(baseDir+'11_R1.fastq','w'),open(baseDir+'11_R2.fastq','w')],
            (D701,D505):[open(baseDir+'12_R1.fastq','w'),open(baseDir+'12_R2.fastq','w')],
            (D706,D506):[open(baseDir+'13_R1.fastq','w'),open(baseDir+'13_R2.fastq','w')],
            (D705,D507):[open(baseDir+'14_R1.fastq','w'),open(baseDir+'14_R2.fastq','w')],
            (D704,D501):[open(baseDir+'15_R1.fastq','w'),open(baseDir+'15_R2.fastq','w')],
            (D703,D502):[open(baseDir+'16_R1.fastq','w'),open(baseDir+'16_R2.fastq','w')]}

with open(readFile1,'rU') as R1File, open(indexFile1,'rU') as I1File, open(readFile2,'rU') as R2File, open(indexFile2,'rU') as I2File:
    counter=0
    R1StringWrite=''
    R2StringWrite=''
    i1Index=''
    i2Index=''
    for eachIndex1 in I1File:                
        eachRead1 = R1File.next()
        eachRead2 = R2File.next()
        eachIndex2 = I2File.next()
        if counter==4: 
            if (i1Index[0:8],i2Index[0:8]) in indexSeq:
                indexSeq[(i1Index[0:8],i2Index[0:8])][0].write(R1StringWrite)
                indexSeq[(i1Index[0:8],i2Index[0:8])][1].write(R2StringWrite)
            counter=0
            R1StringWrite=''
            R2StringWrite=''
            i1Index=''
            i2Index=''
        if counter==1:
            i1Index=eachIndex1.strip()
            i2Index=eachIndex2.strip()
        R1StringWrite+=eachRead1
        R2StringWrite+=eachRead2
        counter+=1            
    #code to take care of last line of fastq files
    if (i1Index[0:8],i2Index[0:8]) in indexSeq:
        indexSeq[(i1Index[0:8],i2Index[0:8])][0].write(R1StringWrite)
        indexSeq[(i1Index[0:8],i2Index[0:8])][1].write(R2StringWrite)

for eachKey in indexSeq:
    for eachFile in indexSeq[eachKey]:
        eachFile.close()
    
    
    
        
    
    
    
