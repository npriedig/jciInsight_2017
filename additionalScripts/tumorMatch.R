# tumorMatch
# perform on two vcf's (either exome or RNA-seq derived) to get proportion of shared variants
# example script and workflow below function

tumorMatch <- function(vcf1, vcf2, id1=NULL, id2=NULL, c_ref1=NULL, c_ref2=NULL, c_geno1=NULL, c_geno2=NULL){
  if(is.null(id1)) id1=paste(vcf1[,1],vcf1[,2],sep="_")
  if(is.null(id2)) id2=paste(vcf2[,1],vcf2[,2],sep="_")
  if(is.null(c_ref1)) c_ref1=4
  if(is.null(c_ref2)) c_ref2=4
  if(is.null(c_geno1)) c_geno1=10
  if(is.null(c_geno2)) c_geno2=10
  
  ref1 <- vcf1[,c_ref1]
  geno1 <- vcf1[,c_geno1]
  ref2 <- vcf2[,c_ref2]
  geno2 <- vcf2[,c_geno2]  
  alt1 <- vcf1[,c_ref1+1]
  alt2 <- vcf2[,c_ref2+1]
  
  idx <- NULL
  idx <- na.omit(intersect(id1, id2))
  print(length(idx))
  idx1<-match(idx,id1)
  idx2<-match(idx,id2)
  ref_match <- rep(0,length(idx))
  for (i in 1:length(idx)){
    ref_match[i] <- as.character(ref1[idx1[i]]) == as.character(ref2[idx2[i]])
  }
  swap_idx <- which(ref_match=="0")
  print(length(swap_idx))
  geno1_1<-as.numeric(substr(geno1, 1, 1))
  geno1_2<-as.numeric(substr(geno1, 3, 3))
  geno2_1<-as.numeric(substr(geno2, 1, 1))
  geno2_2<-as.numeric(substr(geno2, 3, 3))
  swap_nomatch <- 0
  remove <-NULL
  for(i in swap_idx){
    if( (as.character(alt1[i]) != as.character(ref2[i])) |
          (as.character(alt2[i]) != as.character(ref1[i])) ){
      swap_nomatch<-swap_nomatch+1
      remove <- c(remove, i)
    }
  }
  swap_idx <- swap_idx[! swap_idx %in% remove]
  print(length(swap_idx))
  print(swap_nomatch)
  for(i in swap_idx){ #swap
    geno1_1[i]<-abs(geno1_1[i]-1)
    geno1_2[i]<-abs(geno1_2[i]-1)
  }
  if(substr(geno1, 2, 2)[1] == "/"){ #if unphased
    geno_match<-( ((geno1_1[match(idx, id1)]==geno2_1[match(idx, id2)]) & 
                     (geno1_2[match(idx, id1)]==geno2_2[match(idx, id2)])) |
                    ((geno1_1[match(idx, id1)]==geno2_2[match(idx, id2)]) & 
                       (geno1_2[match(idx, id1)]==geno2_1[match(idx, id2)])) ) 
    length(which(geno_match=="FALSE"))
  } else if (substr(geno1,2,2)[1]=="|") {
    geno_match<-( (geno1_1[match(idx, id1)]==geno2_1[match(idx, id2)]) & 
                    (geno1_2[match(idx, id1)]==geno2_2[match(idx, id2)]) )
  }
  score <- (length(which(geno_match)=="TRUE")-length(swap_nomatch))/length(idx)
  print(score)
  return(list(score=score, total=length(which(geno_match)=="TRUE")))
}

#################################################################
#################################################################
#################################################################

setwd("~/tumorMatchOutDirectory")
vcf.directory <- "~/hg19-vcf"

# read in all vcf files
vcf.file.list <- list.files(vcf.directory, recursive=TRUE, pattern = ".vcf")
sample.list <- substr(vcf.file.list,1,10) # get sample labels
setwd(vcf.directory)
i <- 1
for (file in vcf.file.list) { 
  assign(sample.list[i], read.table(file))
  i <- i+1
}

# create dummy dataframe where tumorMatch scores will populate
tumorMatch.scores.df <- data.frame(matrix(NA, nrow=length(sample.list), ncol=length(sample.list)))
rownames(tumorMatch.scores.df) <- sample.list
colnames(tumorMatch.scores.df) <- sample.list

# run tumorMatch
for (sample1 in sample.list) {
  for (sample2 in sample.list) {
    tumorMatch.out <- tumorMatch(get(sample1), get(sample2))
    tumorMatch.scores.df[sample1,sample2] <- tumorMatch.out$score
  }
}

save(tumorMatch.scores.df file = "tumorMatch.scores.Rda")

# create proportion of shared variants matrix plot
pdf(file = "bone-met-proportion-of-shared-variants-matrix.pdf", height = 8, width = 8)
corrplot(as.matrix(tumorMatch.scores.df), method = "square", is.corr = FALSE,tl.col = "black",tl.cex = 0.8)
dev.off()



