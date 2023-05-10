##########################################################################
##########################################################################
# xHLA-R v.4.0
##########################################################################
##########################################################################


##########################################################################
## Libraries
##########################################################################
{
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

library(tidyverse)
library(jsonlite)
library(lpSolve)
library(IRanges)
}




##########################################################################
## Prep
##########################################################################
### Check installation of reference file
### If no reference file, create it.
if(!file.exists("chr6/chr6.fa.pac")){system("bwa index -a bwtsw chr6/chr6.fa")}

### Create Diamond file
system("diamond makedb --in data/hla.faa -d data/hla")

### Unzip any zipped fastq data
system("gunzip input/*.gz")

### Make all bin files executable
system("chmod +x bin/*")

##########################################################################
## Create a function to do HLA typing
##########################################################################

xhla<-function(sample, reads.to.use,seeder = sample(1:1e6, 1)){


  # calculate the file names

  readone<-paste(sample,"_L001_R1_001.fastq",sep="")
  readtwo<-paste(sample,"_L001_R2_001.fastq",sep="")

  message("Sampling ",reads.to.use," reads from read one")
  system(paste("seqtk sample -s", seeder ," input/", readone," ",reads.to.use,"  > input/input_tmp1.fastq",sep=""))

  message("Sampling ",reads.to.use," reads from read two")
  system(paste("seqtk sample -s", seeder ," input/", readtwo," ",reads.to.use,"  > input/input_tmp2.fastq",sep=""))

  message("Running bwa to compare chr6 to fastq")
  system(paste("bwa mem chr6/chr6.fa input/input_tmp1.fastq input/input_tmp2.fastq > input/",sample,".sam",sep=""))

  message("Running samtools to index and align reads")

  system(paste("samtools view -bT chr6/chr6.fa input/",sample,".sam > input/",sample,".bam",sep=""))
  system(paste("samtools sort input/",sample,".bam > input/",sample,".sorted.bam",sep=""))
  system(paste("samtools index input/",sample,".sorted.bam input/",sample,".sorted.bai",sep=""))

  system(paste("samtools view -u input/",sample,".sorted.bam chr6:29886751-33090696 | samtools view -L ./data/hla.bed - > input/",sample,"_extract.sam",sep=""))

  system(paste("perl bin/preprocess.pl input/",sample,"_extract.sam | gzip > input/",sample,".fq.gz",sep=""))

  #Aligning reads to IMGT database
  system(paste("bin/align.pl input/",sample,".fq.gz input/",sample,".tsv",sep=""))

  #performing typing
  system(paste("bin/typing.r input/",sample,".tsv input/",sample,".hla",sep=""))

  #echo "$BIN/report.py -in $OUT/${ID}.hla -out $OUT/${ID}.json -subject $ID -sample $ID"

  system(paste("python3 bin/report.py -in input/",sample,".hla -out input/",sample,".json -subject ",sample," -sample ",sample,sep=""))

  system(paste("bin/full.r input/",sample,".tsv input/",sample,".hla input/",sample,".hla.full",sep=""))

  message("running xHLA typer.sh")

  system(paste("perl bin/typer.sh input/",sample,".sorted.bam ",sample,sep=""))

  system("rm input/input_tmp*.fastq")
  system("rm input/*.sam")
  system("rm input/*.bam")
  system("rm input/*.bai")
  system("rm input/*.tsv")



  result = read_json(paste("input/",sample,".json",sep="")) %>%
    as_tibble

  result.report<-tibble(
    sample = result$subject_id,
    HLA_Ai = result$hla$alleles[[1]],
    HLA_Aii = result$hla$alleles[[2]],
    HLA_Bi = result$hla$alleles[[3]],
    HLA_Bii = result$hla$alleles[[4]],
    HLA_Ci = result$hla$alleles[[5]],
    HLA_Cii = result$hla$alleles[[6]],
    HLA_DRB1i = result$hla$alleles[[7]],
    HLA_DRB1ii = result$hla$alleles[[8]],
    HLA_DQB1i = result$hla$alleles[[9]],
    HLA_DQB1ii = result$hla$alleles[[10]],
    HLA_DPB1i = result$hla$alleles[[11]],
    HLA_DPB1ii = result$hla$alleles[[12]]
  )

  write_csv(result.report,file = paste("output/",sample,"_hla_typing_results.csv",sep=""))
  return(result.report)
}


##########################################################################
### Read list of samples
##########################################################################
a<-list.files(path = "input/",pattern = "fastq")
a<-as.data.frame(a)
b<-as.data.frame(regexpr(pattern = "*L001",text = a$a))
samples<-substr(a$a,start = 1,stop = (b[,1])-2)
samples<-as.data.frame(unique(samples))
names(samples)<-c("id")

head(samples)
rm(a,b)

##########################################################################
## Apply the function to a single fastq pair
##########################################################################
df<-xhla(sample = samples$id[1],reads.to.use = 7500)

##########################################################################
## Loop on remaining samples, binding each to df and saving individual reports
##########################################################################

for(i in 2:nrow(samples)){
                          message(paste("starting alignment",i,"of ",nrow(samples),sep = " "))
                          df<-bind_rows(df,xhla(sample = samples$id[i],reads.to.use = 7500))
                          }

write_csv(df,"output/full_typing_report.csv")
