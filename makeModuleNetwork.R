####function to take a network/module object, 
####and produce a weighted module network
makeModuleNetwork <- function(x,modKeep){
	network <- matrix(0,length(modKeep),length(modKeep))
	#restrict to necessary modules
	mods <- dplyr::filter(x$modules,modulelabels%in%modKeep)

	for (i in 1:(length(modKeep)-1)){
		for (j in (i+1):(length(modKeep))){
			ind1 <-  x$modules$GeneIDs[x$modules$modulelabels %in% modKeep[i]]
			ind2 <-  x$modules$GeneIDs[x$modules$modulelabels %in% modKeep[j]]			
			network[i,j] <- mean(x$network[ind1,ind2])
		}
	}
	network <- network + t(network)
	colnames(network) <- modKeep
	rownames(network) <- modKeep
	return(network)

}