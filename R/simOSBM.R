#' Simulates a matrix from the Stochastic Blockmodel.
#'
#' @param blockPMat The block probability matrix.
#' @param nMembers A vector specifying the number of nodes in each block.
#'
#' @export
#' @return The sparse representation of an adjacency matrix simulated
#' according to the Stochastic Blockmodel.
#'
#' @keywords stochastic blockmodel
simOSBM <- function(blockPMat, nMembers, mixed, thetaVec) {
    
    nBlocks = length(nMembers)
    nNodes = sum(nMembers)
    Z <- matrix(0,nNodes,kcluster)
    cum_nMembers = c(0, cumsum(nMembers))
    for( i in 1:kcluster){
        Z[(cum_nMembers[i]+1):cum_nMembers[i+1],i]=1
    }
    Z <- rbind(Z,mixed)
    adjMat = Z%*%blockPMat%*%t(Z)
    thetaVec <- matrix(thetaVec)
    adjMat = adjMat *(thetaVec%*%t(thetaVec))
    
    # success probabilities are not the same, thus simulate one element after another.
    nNodes <- nNodes+nrow(mixed)
    for ( i in 1: nNodes){
        block_i <- max(which(cum_nMembers <i))
        for ( j in i:nNodes){
            block_j <- max(which(cum_nMembers <j))
            adjMat[i,j] <- rbinom(n=1, size = 1, prob = adjMat[i,j])
        }
    }
    
    adjMat <- Matrix(adjMat)     ##
    return( forceSymmetric(adjMat, uplo = "U") ) 
}
