from moduleUsefulFunctions_20180215 import *
import subprocess

with open('allSeq_toBlast_blastRunner_v1_20190221.fasta','w') as wf:
    for i in [1,2,3,4,5,6,7,8,9,10,11,12,13,16]:
        with open(str(i)+'_afterMS2removal_finalRefCollapsing_v1_20190221_allLists.pckl','rb') as f:
            [xToPlot1,yToPlot1,listAllRefs1,listAllClusters1,listdiGRef1]=cPickle.load(f)
        for counter, eachSeq in enumerate(listdiGRef1):
            wf.write('>'+str(i)+'.'+str(counter+1)+'\n')
            wf.write(eachSeq+'\n')

subprocess.check_call('ncbi-blast-2.7.1+/bin/blastn -task blastn-short -out blastresults_allSeq_refseq_20190221.txt -num_threads 25 -db allSubsets_combined_refseqGenomic'+' -word_size 7 -query allSeq_toBlast_blastRunner_v1_20190221.fasta -outfmt 11 -evalue 100000 -num_alignments 10 -dust no -soft_masking false -show_gis -max_hsps 1 -perc_identity 95 -ungapped',shell=True)
subprocess.check_call('ncbi-blast-2.7.1+/bin/blast_formatter -archive blastresults_allSeq_refseq_20190221.txt -outfmt 7 -out blastresults_allSeq_refseq_outfmt7_20190221.txt',shell=True)
subprocess.check_call('ncbi-blast-2.7.1+/bin/blast_formatter -archive blastresults_allSeq_refseq_20190221.txt -outfmt 5 -out blastresults_allSeq_refseq_outfmt5_20190221.xml',shell=True)        
    
subprocess.check_call('ncbi-blast-2.7.1+/bin/blastn -task blastn-short -out blastresults_allSeq_mySeeds_20190221.txt -num_threads 25 -db mySeeds'+' -word_size 7 -query allSeq_toBlast_blastRunner_v1_20190221.fasta -outfmt 11 -evalue 100000 -num_alignments 10 -dust no -soft_masking false -show_gis -max_hsps 1 -perc_identity 95 -ungapped',shell=True)
subprocess.check_call('ncbi-blast-2.7.1+/bin/blast_formatter -archive blastresults_allSeq_mySeeds_20190221.txt -outfmt 7 -out blastresults_allSeq_mySeeds_outfmt7_20190221.txt',shell=True)
subprocess.check_call('ncbi-blast-2.7.1+/bin/blast_formatter -archive blastresults_allSeq_mySeeds_20190221.txt -outfmt 5 -out blastresults_allSeq_mySeeds_outfmt5_20190221.xml',shell=True) 
