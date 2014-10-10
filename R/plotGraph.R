#' Plots the adjacency matrix.
#'
#' @param adjMat The adjacency matrix to be plotted.
#' @param membership The node membership used to reorder the adjacency matrix
#' before plotting.
#'
#' @export
#'  
plotAdj <- function(adjMat, membership = NULL) {
    
    if(!is.null(membership)) {
        adjMat = adjMat[membership, membership]
    }

   # This function (Matrix::image) has a bug that cuts off the first and last
   # row of the matrix. To correct for this we specify ylim manually.
   image(adjMat, scales = list(draw=F), xlab = NULL, ylab = NULL, lwd = .1,
         sub = NULL, ylim = c(0, dim(adjMat)[1] + 1) )
}

#' Plots the Stochastic Blockmodel.
#'
#' @param blockPMat The block probability matrix.
#' @param nMembers A vector of the number of nodes in each block.
#'
#' @export 
#' 
plotSBM <- function(blockPMat, nMembers) {
    #adjust the block end points based on number of members
    blockGrid = c(0, cumsum(nMembers/sum(nMembers)))

    image(blockGrid, blockGrid, blockPMat, zlim = c(0,1), xlim = c(0,1),
          ylim = c(0,1), col = gray(0:100/100), axes = F, xlab = "",
          ylab = "", lwd = .1)
}