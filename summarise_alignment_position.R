summarise_alignment_position <- function(phylip_file) {
  # Reading in the phylip file
  temp <- readLines(phylip_file)
  # Pulling out the reference name (should be in row 2 if it is the first sequence)
  refseq_name <- unlist(strsplit(temp[2],"[[:space:]]"))[1]
  # Pulling out the reference seq (should be in row 2 if it is the first sequence)
  refseq_seq <- unlist(strsplit(temp[2],"[[:space:]]"))[2]
  # Getting the total length of the alignment
  totz_length <- nchar(refseq_seq)
  # Setting up the output matrix
  outputmatrix <- matrix(c("Sample","Gapped start (bp)","Gapped end (bp)","Ungapped start (bp)","Ungapped end (bp)"),ncol=5,nrow=1)
  # Starting with the first non-reference sample and looping through
  for (i in 3:length(temp)) {
    # Creating a temporary matrix
    tempmatrix <- matrix(NA,ncol=5,nrow=1)
    # In which we store the sample's name
    tempmatrix[1,1] <- unlist(strsplit(temp[i],"[[:space:]]"))[1]
    # Extracting the samples sequence
    temp_seq <- unlist(strsplit(temp[i],"[[:space:]]"))[2]
    # We first calculate where it starts relative to the total alignment by stripping out all the starting gaps
    tempmatrix[1,2] <- totz_length-nchar(gsub("^-*","",temp_seq))
    # And then we do the same thing with the trailing gaps to find out where it ends
    tempmatrix[1,3] <- totz_length-nchar(gsub("-*$","",temp_seq))
    # We use these start and end alignment points to extract the reference sequence at these points and then strip out the gaps. This gives us our ungapped start and stop locations
    tempmatrix[1,4] <- nchar(gsub("-","",substr(refseq_seq, 0, as.numeric(tempmatrix[1,2]))))
    tempmatrix[1,5] <- nchar(gsub("-","",substr(refseq_seq, 0, as.numeric(tempmatrix[1,3]))))
    # We bind these to our outputmatrix
    outputmatrix <- rbind(outputmatrix,tempmatrix)
  }  
  # And then write the whole thing out
  write.table(outputmatrix,"summarize_alignment_position.txt",quote = FALSE,row.names=FALSE,col.names=FALSE,sep=",")
}
