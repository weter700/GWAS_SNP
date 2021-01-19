rm(list=ls())

library(snp.plotter)
library(SNPassoc)
DDD<-choose.files()
WORKING_FOLDER<-dirname(DDD)
BASE_NAME<-basename(WORKING_FOLDER)
NAME_FILE<-substr(BASE_NAME, 9, nchar(BASE_NAME))
DDD<-read.csv2(DDD)

# Changing TRUE into T
for(COUNT_1 in 13:ncol(DDD)){							
  if(DDD[1,COUNT_1]==TRUE){
    DDD[,COUNT_1]<-c("T")
  }
}

# Alleles to genotype
A<-seq(13,ncol(DDD), by = 2)
PAT<-c()
newMAT<-matrix(ncol=length(A),nrow=nrow(DDD))
COLNAMES_DDD<-c(colnames(DDD)[A])
colnames(newMAT)<-COLNAMES_DDD
for(x in 1:nrow(DDD)){
  for(y in A){
    PAT<-c(PAT,gsub("[^[:alnum:]]","", paste(DDD[x,y],DDD[x,y+1])))
  }
  newMAT[x,]<-PAT
  PAT<-c()
}
PHENO_GENO<-cbind(DDD[,c(1:12)],newMAT)

#Changing genotype 00 into NA
for(COUNT_62 in 13:ncol(PHENO_GENO)){
  for(COUNT_63 in 1:nrow(PHENO_GENO)){
    if(PHENO_GENO[COUNT_63,COUNT_62]=="00"){
      PHENO_GENO[COUNT_63,COUNT_62]<-NA
    }
  }
}

#Saving RAW data file
write.csv2(PHENO_GENO,paste0(WORKING_FOLDER,"/", NAME_FILE,"_PHENO_GENO_RAW.csv"), row.names = FALSE)

#Selection, only healthy patients with gender info
PHENO_GENO<-PHENO_GENO[PHENO_GENO$GENDER!=0,] 
PHENO_GENO<-PHENO_GENO[PHENO_GENO$DISEASE==1,] 

#Saving selected data file
write.csv2(PHENO_GENO,paste0(WORKING_FOLDER,"/", NAME_FILE,"_PHENO_GENO_SEL.csv"), row.names = FALSE)

