require(metanetwork)

fileName <- as.character(commandArgs(TRUE)[[1]])
nodes <- as.numeric(commandArgs(TRUE)[[2]])
pathv <- as.character(commandArgs(TRUE)[[3]])
data <- data.matrix(read.csv(fileName,row.names=1))
sparrowMPI(data,nodes,pathv)
