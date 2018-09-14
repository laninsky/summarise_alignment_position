summarise_alignment_position <- function(phylip_file) {
  temp <- readLines(phylip_file)
  refseq_name <- unlist(strsplit(temp[2],"[[:space:]]"))[1]
  refseq_seq <- unlist(strsplit(temp[2],"[[:space:]]"))[2]
  totz_length <- nchar(refseq_seq)
  outputmatrix <- matrix(c("Sample","Gapped start (bp)","Gapped end (bp)","Ungapped start (bp)","Ungapped end (bp)"),ncol=5,nrow=1)
  for (i in 3:length(temp)) {
    tempmatrix <- matrix(NA,ncol=5,nrow=1)
    tempmatrix[1,1] <- unlist(strsplit(temp[i],"[[:space:]]"))[1]
    temp_seq <- unlist(strsplit(temp[i],"[[:space:]]"))[2]
    tempmatrix[1,2] <- totz_length-nchar(gsub("^-*","",temp_seq))
    tempmatrix[1,3] <- totz_length-nchar(gsub("-*$","",temp_seq))
    tempmatrix[1,4] <- nchar(gsub("-","",substr(refseq_seq, 0, as.numeric(tempmatrix[1,2]))))
    tempmatrix[1,5] <- nchar(gsub("-","",substr(refseq_seq, 0, as.numeric(tempmatrix[1,3]))))
    outputmatrix <- rbind(outputmatrix,tempmatrix)
  }  
  write.table(outputmatrix,"summarize_alignment_position.txt",quote = FALSE,row.names=FALSE,col.names=FALSE,sep=",")
}
