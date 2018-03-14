#NOTE the first file (DATA_alignments.csv) should not have a header.

RemoveDups=function(bowtie_file_csv,fasta_file_csv){
 aligns=read.csv(paste(bowtie_file_csv),header=FALSE)
 misses=which(aligns$V1!=aligns$V3)
 
 lst3=vector();removes=vector();
 for (i in 1:length(misses)){  #first make sure your i is not in lst3
  lst1=vector()
  lst2=vector()
  if (i%in%lst3==FALSE){
   
   lst1=c(lst1,i)
   
   lst2=c(lst2,aligns[misses[i],]$V1,aligns[misses[i],]$V3)#put your aligns[misses[i],]$V1 and your aligns[misses[i],]$V3 in lst2
   
   ##***
   outz=1
   while(outz>0){
    outz=0;
    c=which(aligns[misses,]$V3%in%lst2)#see what aligns with your lst2
    d=which(aligns[misses,]$V1%in%lst2)
    
    #what are in those sets?
    
    if (length(c)>0){lst1=c(lst1,c)}#add these indices to list1 if they are new not already in lst1
    if (length(d)>0){lst1=c(lst1,d)}
    
    #now lst1 should have all the c and ds.
    lst1=unique(lst1)
    
    
    for (j in 1:length(lst1)){
     if((aligns[misses[lst1[j]],]$V1%in%lst2)==FALSE){lst2=c(lst2,aligns[misses[lst1[j]],]$V1)}#if all the pairs are already in lst2, list of loci
     outz=sum(((aligns[misses[lst1[j]],]$V1%in%lst2)==FALSE),outz)
     #check to see if that new lst2 is anywhere else
     if((aligns[misses[lst1[j]],]$V3%in%lst2)==FALSE){lst2=c(lst2,aligns[misses[lst1[j]],]$V3)} 
     outz=sum(((aligns[misses[lst1[j]],]$V1%in%lst2)==FALSE),outz)
    }
   }
   if (outz==0){
    #then go back to ##***
    #when all your cs and ds are in lst2
    lst3=c(lst3,lst1);lst3=unique(lst3);
    if(length(lst2)>1){removes=unique(c(removes,lst2[2:length(lst2)]))}
   }
  }#if i is not in lst3
 }#for loop
 
 #now take the duplicates out of your fasta file
 badz=vector()
 old_fasta_file=read.csv(paste(fasta_file_csv),header=FALSE)# 
 #new_fasta_file=paste(fasta_file_csv,"NoDups.csv",sep="")
 for (i in 1:length(removes)){
  badz[i]=which(old_fasta_file$V1==removes[i])}
 new_fasta=old_fasta_file[-badz[!is.na(badz)],]
 write.csv(new_fasta,file=paste(fasta_file_csv,"NoDups.csv"))
 return(removes)
}#end of function

#HOW TO RUN THE FUNCTION
setwd("/path/to/data/")
RemoveDups("DATA_alignments.csv","DATA_output_haps.csv")
