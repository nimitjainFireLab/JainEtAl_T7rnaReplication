import sys
import subprocess

dir1='AF_SOL_826/'
dir2='AF_SOL_850/'
dir3='AF_SOL_722/'
dir4='AF_SOL_713/'

names1=[]
for eachEl in ['257_1','257_2','257_9','257_10','258_1','258_2','258_7','258_8']:
    names1.append(dir1+eachEl)
names2=[]
for eachEl in ['3','4','1','2']:
    names2.append(dir2+eachEl)
names3=[]
for eachEl in ['4','5','7','11','13','17']:
    names3.append(dir3+eachEl)
names4=[]
for eachEl in ['5','6']:
    names4.append(dir4+'Sample'+eachEl)

names=names1+names2+names3

for counter, eachEl in enumerate(names):
    fileName1=eachEl+'_R1.fastq'
    fileName2=eachEl+'_R2.fastq'
    sampleName=eachEl[eachEl.rfind('/')+1:]
    if counter<8:
        addendum='AF_SOL_826_'
    elif counter<12:
        addendum='AF_SOL_850_'
    else:
        addendum='AF_SOL_722_'
    sampleName=addendum+sampleName
    trimmomatic1=sampleName+'_trimmomaticP_R1.fastq'
    trimmomatic2=sampleName+'_trimmomaticP_R2.fastq'
    trimmomaticu1=sampleName+'_trimmomaticUP_R1.fastq'
    trimmomaticu2=sampleName+'_trimmomaticUP_R2.fastq'
    print ''
    print ''
    subprocess.check_call('java -jar trimmomatic-0.36.jar PE -threads 5 '+fileName1+' '+fileName2+' '+trimmomatic1+' '+trimmomaticu1+' '+trimmomatic2+' '+trimmomaticu2+' ILLUMINACLIP:TruSeqAdaptersPETrimmomatic.fa:2:30:10:5:true',shell=True)

for eachEl in names4:
    fileName=eachEl+'.fastq'
    sampleName='AF_SOL_713_'+eachEl[eachEl.rfind('/')+1:]
    trimmomaticF=sampleName+'_trimmomatic.fastq'
    print ''
    print ''
    subprocess.check_call('java -jar trimmomatic-0.36.jar SE -threads 5 '+fileName+' '+trimmomaticF+' ILLUMINACLIP:TruSeqAdaptersSETrimmomatic.fa:2:30:10',shell=True)
    
    
