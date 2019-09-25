import subprocess

subprocess.check_call('ncbi-blast-2.7.1+/bin/blastn -task blastn-short -out blastresults5_20190305.txt -num_threads 25 -db t7rp1_combined_nt'+' -word_size 7 -query sequencesToBlast_20190305.fasta -outfmt 11 -evalue 100000 -num_alignments 20 -dust no -soft_masking false -show_gis -max_hsps 3',shell=True)
subprocess.check_call('ncbi-blast-2.7.1+/bin/blast_formatter -archive blastresults5_20190305.txt -outfmt 7 -out blastresults5_outfmt7_20190305.txt',shell=True)
subprocess.check_call('ncbi-blast-2.7.1+/bin/blast_formatter -archive blastresults5_20190305.txt -outfmt 5 -out blastresults5_outfmt5_20190305.xml',shell=True) 
subprocess.check_call('ncbi-blast-2.7.1+/bin/blast_formatter -archive blastresults5_20190305.txt -outfmt 0 -out blastresults5_outfmt0_20190305.txt',shell=True) 
