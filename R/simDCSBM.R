#' Simulates a matrix from the Degree-Corrected Stochastic Blockmodel.
#'
#' @param blockPMat The block probability matrix.
#' @param nMembers A vector specifying the number of nodes in each block.
#'
#' @export
#' @return The sparse representation of an adjacency matrix simulated
#' according to the Stochastic Blockmodel.
#'
#' @keywords degree corrected stochastic blockmodel
simDCSBM <- function(blockPMat, nMembers, thetaVec) {
    
    nBlocks = length(nMembers)
    nNodes = sum(nMembers)
    adjMat = matrix(runif(n= nNodes*nNodes), nNodes, nNodes)
    thetaM <- matrix(0, nNodes, nBlocks)
    csize <- c(0, cumsum(nMembers))
    for ( i in length(nMembers)){
        thetaM[(csize[i]+1):csize[i+1], i] <- thetaVec[(csize[i]+1):csize[i+1]]
    }
    
    adjMat <- (adjMat < thetaM %*% blockPMat %*% t(thetaM))
    adjMat <- Matrix(adjMat)     ##
    return( forceSymmetric(adjMat, uplo = "U") ) 
}


