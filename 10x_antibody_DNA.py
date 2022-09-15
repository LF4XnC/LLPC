# -*- coding: utf-8 -*-


path='file path'
file='file name'
out='file out'
ref='barcode list'


##import barcode list
import gzip
fh=open(ref,'r')
fh.readline()
sh={}
shture={}
for line in fh:
    line1=line.strip()
    line1=line1.split(',')
    sh[line1[4]]=line1[1]
    shture[line1[4]]=True
fh.close()

##cell primary filter    
count={}
true={}
fh= gzip.open(path+file+'_S1_L001_R1_001.fastq.gz', "rt")
fh2= gzip.open(path+file+'_S1_L001_R2_001.fastq.gz', "rt")
for line1 in fh:
    line2=fh.readline()
    line3=fh.readline()
    line4=fh.readline()       
    r2line1=fh2.readline()
    r2line2=fh2.readline()
    r2line3=fh2.readline()
    r2line4=fh2.readline()
    CB=line2[0:16]
    count[CB]=count.get(CB,0)+1
fh.close()
fh2.close()
for i in count:
    if count[i]>10:
        true[i]=True    
    

##count barcode for each cell barcode
fh= gzip.open(path+file+'_S1_L001_R1_001.fastq.gz', "rt")
fh2= gzip.open(path+file+'_S1_L001_R2_001.fastq.gz', "rt")

reads={}
umi={}
for line1 in fh:
    line2=fh.readline()
    line3=fh.readline()
    line4=fh.readline()       
    r2line1=fh2.readline()
    r2line2=fh2.readline()
    r2line3=fh2.readline()
    r2line4=fh2.readline()
    CB=line2[0:16]
    if true.get(CB,False):
        UMI=line2[16:26]
        if r2line2[15:28]=='CCCATATAAGAAA':
            feature=r2line2[0:15]
            if shture.get(feature,False):
                reads[CB]=reads.get(CB,{})
                reads[CB][feature]=reads[CB].get(feature,0)+1
                umi[CB]=umi.get(CB,{})
                umi[CB][feature]=umi[CB].get(feature,{})
                umi[CB][feature][UMI]=umi[CB][feature].get(UMI,0)+1

fhout=open(out+'reads.txt','w')
fhumi=open(out+'umi.txt','w')
title=''
for key in sh:
    title+='\t'+sh[key]
fhout.writelines(title+'\n')
fhumi.writelines(title+'\n')

for key in reads:
    write=key
    write2=key
    for key1 in sh:
        write+='\t'+str(reads[key].get(key1,0))
        write2+='\t'+str(len(umi[key].get(key1,{})))+','+str(umi[key].get(key1,{}))
    fhout.writelines(write+'\n')
    fhumi.writelines(write2+'\n')
fhout.close()
fhumi.close()
        
                
                
            



