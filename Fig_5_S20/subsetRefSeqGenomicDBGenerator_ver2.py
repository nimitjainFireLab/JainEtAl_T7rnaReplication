import subprocess

dbList='refseq_genomic.00 refseq_genomic.01 refseq_genomic.02 refseq_genomic.03 refseq_genomic.04 refseq_genomic.05 refseq_genomic.06 refseq_genomic.07 refseq_genomic.08 refseq_genomic.09 refseq_genomic.10 refseq_genomic.11 refseq_genomic.12 refseq_genomic.13 refseq_genomic.14 refseq_genomic.15 refseq_genomic.16 refseq_genomic.17 refseq_genomic.18 refseq_genomic.19 refseq_genomic.20 refseq_genomic.21 refseq_genomic.22 refseq_genomic.23 refseq_genomic.24 refseq_genomic.25 refseq_genomic.26 refseq_genomic.27 refseq_genomic.28 refseq_genomic.29 refseq_genomic.30 refseq_genomic.31 refseq_genomic.32 refseq_genomic.33 refseq_genomic.34 refseq_genomic.35 refseq_genomic.36 refseq_genomic.37 refseq_genomic.38 refseq_genomic.39 refseq_genomic.40 refseq_genomic.41 refseq_genomic.42 refseq_genomic.43 refseq_genomic.44 refseq_genomic.45 refseq_genomic.46 refseq_genomic.47 refseq_genomic.48 refseq_genomic.49 refseq_genomic.50 refseq_genomic.51 refseq_genomic.52 refseq_genomic.53 refseq_genomic.54 refseq_genomic.55 refseq_genomic.56 refseq_genomic.57 refseq_genomic.58 refseq_genomic.59 refseq_genomic.60 refseq_genomic.61 refseq_genomic.62 refseq_genomic.63 refseq_genomic.64 refseq_genomic.65 refseq_genomic.66 refseq_genomic.67 refseq_genomic.68 refseq_genomic.69 refseq_genomic.70 refseq_genomic.71 refseq_genomic.72 refseq_genomic.73 refseq_genomic.74 refseq_genomic.75 refseq_genomic.76 refseq_genomic.77 refseq_genomic.78 refseq_genomic.79 refseq_genomic.80 refseq_genomic.81 refseq_genomic.82 refseq_genomic.83 refseq_genomic.84 refseq_genomic.85 refseq_genomic.86 refseq_genomic.87 refseq_genomic.88 refseq_genomic.89 refseq_genomic.90 refseq_genomic.91 refseq_genomic.92 refseq_genomic.93 refseq_genomic.94 refseq_genomic.95 refseq_genomic.96 refseq_genomic.97 refseq_genomic.98 refseq_genomic.99 refseq_genomic.100 refseq_genomic.101 refseq_genomic.102 refseq_genomic.103 refseq_genomic.104 refseq_genomic.105 refseq_genomic.106 refseq_genomic.107 refseq_genomic.108 refseq_genomic.109 refseq_genomic.110 refseq_genomic.111 refseq_genomic.112 refseq_genomic.113 refseq_genomic.114 refseq_genomic.115 refseq_genomic.116 refseq_genomic.117 refseq_genomic.118 refseq_genomic.119 refseq_genomic.120 refseq_genomic.121 refseq_genomic.122 refseq_genomic.123 refseq_genomic.124 refseq_genomic.125 refseq_genomic.126 refseq_genomic.127 refseq_genomic.128 refseq_genomic.129 refseq_genomic.130 refseq_genomic.131 refseq_genomic.132 refseq_genomic.133 refseq_genomic.134 refseq_genomic.135 refseq_genomic.136 refseq_genomic.137 refseq_genomic.138 refseq_genomic.139 refseq_genomic.140 refseq_genomic.141 refseq_genomic.142 refseq_genomic.143 refseq_genomic.144 refseq_genomic.145 refseq_genomic.146 refseq_genomic.147 refseq_genomic.148 refseq_genomic.149 refseq_genomic.150 refseq_genomic.151 refseq_genomic.152 refseq_genomic.153 refseq_genomic.154 refseq_genomic.155 refseq_genomic.156 refseq_genomic.157 refseq_genomic.158 refseq_genomic.159 refseq_genomic.160 refseq_genomic.161 refseq_genomic.162 refseq_genomic.163 refseq_genomic.164 refseq_genomic.165 refseq_genomic.166 refseq_genomic.167 refseq_genomic.168 refseq_genomic.169 refseq_genomic.170 refseq_genomic.171 refseq_genomic.172 refseq_genomic.173 refseq_genomic.174 refseq_genomic.175 refseq_genomic.176 refseq_genomic.177 refseq_genomic.178 refseq_genomic.179 refseq_genomic.180 refseq_genomic.181 refseq_genomic.182 refseq_genomic.183 refseq_genomic.184 refseq_genomic.185 refseq_genomic.186 refseq_genomic.187 refseq_genomic.188 refseq_genomic.189 refseq_genomic.190 refseq_genomic.191 refseq_genomic.192 refseq_genomic.193 refseq_genomic.194 refseq_genomic.195 refseq_genomic.196 refseq_genomic.197 refseq_genomic.198 refseq_genomic.199 refseq_genomic.200 refseq_genomic.201 refseq_genomic.202 refseq_genomic.203 refseq_genomic.204 refseq_genomic.205 refseq_genomic.206 refseq_genomic.207 refseq_genomic.208 refseq_genomic.209 refseq_genomic.210 refseq_genomic.211 refseq_genomic.212 refseq_genomic.213 refseq_genomic.214 refseq_genomic.215 refseq_genomic.216 refseq_genomic.217 refseq_genomic.218 refseq_genomic.219 refseq_genomic.220 refseq_genomic.221 refseq_genomic.222 refseq_genomic.223 refseq_genomic.224 refseq_genomic.225 refseq_genomic.226 refseq_genomic.227 refseq_genomic.228 refseq_genomic.229 refseq_genomic.230 refseq_genomic.231 refseq_genomic.232 refseq_genomic.233 refseq_genomic.234 refseq_genomic.235 refseq_genomic.236 refseq_genomic.237 refseq_genomic.238 refseq_genomic.239 refseq_genomic.240 refseq_genomic.241 refseq_genomic.242 refseq_genomic.243 refseq_genomic.244 refseq_genomic.245 refseq_genomic.246 refseq_genomic.247 refseq_genomic.248 refseq_genomic.249 refseq_genomic.250 refseq_genomic.251 refseq_genomic.252 refseq_genomic.253 refseq_genomic.254 refseq_genomic.255 refseq_genomic.256 refseq_genomic.257 refseq_genomic.258 refseq_genomic.259 refseq_genomic.260 refseq_genomic.261 refseq_genomic.262 refseq_genomic.263 refseq_genomic.264 refseq_genomic.265 refseq_genomic.266 refseq_genomic.267 refseq_genomic.268 refseq_genomic.269 refseq_genomic.270 refseq_genomic.271 refseq_genomic.272 refseq_genomic.273 refseq_genomic.274 refseq_genomic.275 refseq_genomic.276 refseq_genomic.277 refseq_genomic.278 refseq_genomic.279 refseq_genomic.280 refseq_genomic.281 refseq_genomic.282 refseq_genomic.283 refseq_genomic.284 refseq_genomic.285 refseq_genomic.286 refseq_genomic.287 refseq_genomic.288 refseq_genomic.289 refseq_genomic.290 refseq_genomic.291 refseq_genomic.292 refseq_genomic.293 refseq_genomic.294 refseq_genomic.295 refseq_genomic.296 refseq_genomic.297 refseq_genomic.298 refseq_genomic.299 refseq_genomic.300 refseq_genomic.301 refseq_genomic.302 refseq_genomic.303 refseq_genomic.304 refseq_genomic.305 refseq_genomic.306 refseq_genomic.307 refseq_genomic.308 refseq_genomic.309 refseq_genomic.310'

dbList=dbList.split()
#print dbList

pathToAdd=''
#print newDBList

giDictOfInterest={}
with open('giListInclude_refseqgenomic.txt','rU') as f:
    for eachLine in f:
        giDictOfInterest[eachLine.strip()]=0

newVDBList=''

for counter, eachEl in enumerate(dbList):
    print counter
    subprocess.check_call('ncbi-blast-2.7.1+/bin/blastdbcmd -db BlastDB_20180312/'+eachEl+' -entry all -outfmt %g -out dummy.txt',shell=True)
    processDB=0
    with open('dummy.txt','rU') as f:
        for eachLine in f:
            if eachLine.strip() in giDictOfInterest:
                processDB=1
                break
    if processDB==1:
        subprocess.check_call('ncbi-blast-2.7.1+/bin/blastdb_aliastool -db BlastDB_20180312/'+eachEl+' -gilist giListInclude_refseqgenomic.txt -dbtype nucl -out subset_'+eachEl+' -title subset_'+eachEl,shell=True)   
        newVDBList+=' subset_'+eachEl

#newVDBList has a space at the beginning
newVDBList='"'+newVDBList[1:]+'"'
subprocess.check_call('ncbi-blast-2.7.1+/bin/blastdb_aliastool -dblist '+newVDBList+' -dbtype nucl -out allSubsets_combined_refseqGenomic -title allSubsets_combined_refseqGenomic',shell=True)


