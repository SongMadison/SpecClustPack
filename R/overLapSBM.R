#' Simulates a matrix from the Stochastic Blockmodel.
#'
#' @param blockPMat The block probability matrix.
#' @param memberZMat A matrix specifying the memberships of each node.
#'
#' @export
#' @return The sparse representation of an adjacency matrix simulated
#' according to the Stochastic Blockmodel.
#'
#' @keywords stochastic blockmodel
overlapSBM <- function(blockPMat, memberZMat, sparse = TRUE) {
    
    nBlocks = dim(memberZMat)[2]
    nNodes = dim(memberZMat)[1]
    populationAdjMat <- memberZMat%*% blockPMat %*% memberZMat
    AdjMat <- matrix(runif(nNodes*nNodes), nrow = nNodes, ncol = nNodes)
    AdjMat <- (AdjMat < populationAdjMat)
    return( forceSymmetric(triu(adjMat, k=1)) )
}

#' Simulates a Bernoulli random vector.
#'
#' @param nElem The number of elements in the vector.
#' @param p The probability of an element being equal to one.
#'
#' export
#' @return The sparse representation of a simulated Bernoulli vector.
#'
#' @keywords simulate sparse Bernoulli
simBernSparseVec <- function(nElem, p) {
    
    if(p == 0) {
        return( sparseVector(0, 1, length = nElem) )
    }
    
    # get expectation and standard deviation of ones in the vector
    expNumOnes = nElem*p
    sdNumOnes = sqrt(nElem*p*(1-p))
    
    # vector with intervals at which ones occur in the simulated vector
    # neg binom gives the # of failures before a success
    oneIntervals = rnbinom(expNumOnes + round(3*sdNumOnes), 1, p) + 1
    
    # take cumulative sum to get the index values for the ones
    oneIndices = cumsum(oneIntervals)
    
    # loop to ensure the indices cover the entire simulated vector
    while (max(oneIndices) < nElem) {
        oneIndices = c(oneIndices, max(oneIndices) +
                           rnbinom(1, 1, p) + 1)
    }
    
    # truncate the vector to remove access values
    oneIndices = oneIndices[oneIndices <= nElem]
    
    return( sparseVector(1, oneIndices, nElem) )
}
