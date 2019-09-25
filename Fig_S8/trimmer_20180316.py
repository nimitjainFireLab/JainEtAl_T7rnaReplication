import sys
import subprocess

dir1='AF_SOL_771/'

names=['5','9','18','19']

for eachEl in names:
    fileName1=dir1+eachEl+'_R1.fastq'
    fileName2=dir1+eachEl+'_R2.fastq'
    sampleName=eachEl[eachEl.rfind('/')+1:]
    trimmomatic1=sampleName+'_trimmomaticP_R1.fastq'
    trimmomatic2=sampleName+'_trimmomaticP_R2.fastq'
    trimmomaticu1=sampleName+'_trimmomaticUP_R1.fastq'
    trimmomaticu2=sampleName+'_trimmomaticUP_R2.fastq'

    subprocess.check_call('java -jar trimmomatic-0.36.jar PE -threads 5 '+fileName1+' '+fileName2+' '+trimmomatic1+' '+trimmomaticu1+' '+trimmomatic2+' '+trimmomaticu2+' ILLUMINACLIP:TruSeqAdaptersPETrimmomatic.fa:2:30:10:5:true',shell=True)
