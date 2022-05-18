if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("IRanges")

library(jsonlite)
library(lpSolve)
library(IRanges)


message("checking installation of reference file")
if(!file.exists("chr6/chr6.fa.pac")){system("bwa index -a bwtsw chr6/chr6.fa")}

message("creating diamond file")
system("diamond makedb --in data/hla.faa -d data/hla")

message("unzipping any gz files")
system("gunzip input/*.gz")
a<-list.files(path = "input/",pattern = "fastq")
a<-as.data.frame(a)
b<-as.data.frame(regexpr(pattern = "*L001",text = a$a))
samples<-substr(a$a,start = 1,stop = (b[,1])-2)
samples<-as.data.frame(unique(samples))
names(samples)<-c("id")


samples$readone<-paste(samples$id,"_L001_R1_001.fastq",sep="")
samples$readtwo<-paste(samples$id,"_L001_R2_001.fastq",sep="")
if(length(which(samples$id==""))>0){samples<-samples[-which(samples$id==""),]}



message("running HLA typing analysis on all specimens")

for (i in 1:dim(samples)[1])
{

  message(paste("starting analysis ",samples$readone[i],sep=""))
  message("getting sampled data set from fastq files")
  seeder<-sample(1:1e6, 1)

  system(paste("seqtk sample -s", seeder ," input/", samples$readone[i]," 15000 > input/input_tmp1.fastq",sep=""))
  system(paste("seqtk sample -s", seeder ," input/", samples$readtwo[i]," 15000 > input/input_tmp2.fastq",sep=""))

message("Running bwa to compare chr6 to fastq")
system(paste("bwa mem chr6/chr6.fa input/input_tmp1.fastq input/input_tmp2.fastq > input/",samples$id[i],".sam",sep=""))

message("Running samtools to index and align reads")

system(paste("samtools view -bT chr6/chr6.fa input/",samples$id[i],".sam > input/",samples$id[i],".bam",sep=""))
system(paste("samtools sort input/",samples$id[i],".bam > input/",samples$id[i],".sorted.bam",sep=""))
system(paste("samtools index input/",samples$id[i],".sorted.bam input/",samples$id[i],".sorted.bai",sep=""))

message("running xHLA typer.sh")

system(paste("perl bin/typer.sh input/",samples$id[i],".sorted.bam ",samples$id[i],sep=""))


message("reading results data")  
if (length(list.files(paste("hla-",samples$id[i],"/",sep=""),pattern=".json"))>0)
{  tmp1<-fromJSON(txt=paste("hla-",samples$id[i],"/",samples$id[i],".json",sep=""))
  write_json(x = tmp1,path = paste("output/",samples$id[i],".json",sep=""))
}


  #prepare format for first specimen
 if (samples$id[i] %in% gsub(".json","",(list.files("output/",pattern=".json"))))
{
  tmp<-fromJSON(txt=paste("output/",samples$id[i],".json",sep=""))
  sampleID<-tmp$subject_id
  outputreport<-as.data.frame(sampleID,stringsAsFactors = F)
  outputreport$analysis.date<-tmp$creation_time
  outputreport$report.version<-tmp$report_version
  outputreport$report_type<-tmp$report_type

  outputreport$hla.Ai<-tmp$hla$alleles[grep(pattern = "^A",x = tmp$hla$alleles)[1]]
  outputreport$hla.Aii<-tmp$hla$alleles[grep(pattern = "^A",x = tmp$hla$alleles)[2]]
  outputreport$hla.Bi<-tmp$hla$alleles[grep(pattern = "^B",x = tmp$hla$alleles)[1]]
  outputreport$hla.Bii<-tmp$hla$alleles[grep(pattern = "^B",x = tmp$hla$alleles)[2]]
  outputreport$hla.Ci<-tmp$hla$alleles[grep(pattern = "^C",x = tmp$hla$alleles)[1]]
  outputreport$hla.Cii<-tmp$hla$alleles[grep(pattern = "^C",x = tmp$hla$alleles)[2]]
  outputreport$hla.DPB1i<-tmp$hla$alleles[grep(pattern = "^DP",x = tmp$hla$alleles)[1]]
  outputreport$hla.DPB1ii<-tmp$hla$alleles[grep(pattern = "^DP",x = tmp$hla$alleles)[2]]
  outputreport$hla.DQB1i<-tmp$hla$alleles[grep(pattern = "^DQ",x = tmp$hla$alleles)[1]]
  outputreport$hla.DQB1ii<-tmp$hla$alleles[grep(pattern = "^DQ",x = tmp$hla$alleles)[2]]
  outputreport$hla.DRB1i<-tmp$hla$alleles[grep(pattern = "^DR",x = tmp$hla$alleles)[1]]
  outputreport$hla.DRB1ii<-tmp$hla$alleles[grep(pattern = "^DR",x = tmp$hla$alleles)[2]]
  #outputreport$crosschecked<-paste(tmp$hla$matches,collapse = ",")
  rm(tmp)
  data_out<-outputreport

  message("writing results to output/000_alldata_out.txt")
  if(file.exists("output/000_alldata_out.txt")){write.table(outputreport,append = T,file = "output/000_alldata_out.txt",col.names = F,row.names = F,sep = "\t",quote=F)}
  if(!file.exists("output/000_alldata_out.txt")){write.table(outputreport,append = F,file = "output/000_alldata_out.txt",col.names = T,row.names = F,sep = "\t",quote=F)}
}
	
  system("rm input/input_tmp*.fastq")
  system("rm input/*.sam")
  system("rm input/*.bam")
  system("rm input/*.bai")
  system(paste ("rm -rf hla-",samples$id[i],"/",sep="")) 
  message("done")
}

message("remove fastq files")
system("rm input/*fastq")
