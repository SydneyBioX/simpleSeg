#' Calculates the tissue mask of a given image
#'
#' @param image A multiplexed image or list of multiplexed images, can be a CytoImageList
#' @param tissue_index a numerical index that corresponds to the tissue markers in the image object.
#' @param cores The number of parallel processing cores to be used
#' 
#' @return Tissue masks of given image/s (area of the image/s containing tissue)

#' @examples
#' 
#' image.mask <- calcTissueMaskParallel(image = image.list, tissue_index = c(2,3,4), cores = 40)
#' 
#' @export calcTissueMaskParallel
#' @rdname calcTissueMaskParallel
#' 
#' @importFrom BiocParallel bplapply
#' @importFrom EBImage normalize otsu bwlabel
#' @importFrom grDevices chull
#' @importFrom spatstat.geom owin
#' @importFrom spatstat.geom as.mask



############################## Tissue Mask Function ########################################
calcTissueMaskParallel <- function(
    image,
    tissue_index,
    size_selection = 10,
    cores = 50
){
    result <- BiocParallel::bplapply(image, calcTissueMask, tissue_index = tissue_index, size_selection = size_selection, BPPARAM  = BiocParallel::MulticoreParam(workers = cores))
    return(result)
}


calcTissueMask <- function(image, 
                             # an image object containing 1 or more markers
                             tissue_index,
                             # a numerical index that corresponds to the tissue markers in the image object.
                             size_selection = 10
                             # a size selection to remove extremely small cells/large noise that disrupts the calculation on tissue mask.
){
    
    tissue1 <- image[,,tissue_index[1]]
    
    tissue2 <- image[,,tissue_index[2]]
    
    tissue3 <- image[,,tissue_index[3]]
    
    
    # replace outside of circle with NA
    
    circle <- tissue1 + tissue2 + tissue3 # add intesities of these markers to try and highlight the circle structure.
    
    
    circle <- EBImage::normalize(log10(circle+0.01),
                                 ft = c(0,1))
    
    circle_otsu <- EBImage::otsu(circle,
                        range = c(0,1))
    
    circle <- circle > circle_otsu
    
    
    # Do a size selection 
    circle_label <- EBImage::bwlabel(circle) 
    tcircle <- table(circle_label)
    circle[circle_label%in%names(which(tcircle<=size_selection))] <- 0 # size selection of 10.
    
    # extract the convex hull
    nonZero <- as.data.frame(which(circle>0,2))
    colnames(nonZero) <- c("y", "x")
    
    ch <- grDevices::chull(nonZero[, c("x", "y")])
    poly <- nonZero[, c("x", "y")][rev(ch),]
    colnames(poly) <- c("x", "y")
    ow <- spatstat.geom::owin(xrange = range(nonZero$x), yrange = range(nonZero$y), poly = poly)
    tissueMask <- as.matrix(spatstat.geom::as.mask(ow, xy=list(y = 1:nrow(circle), x = 1:ncol(circle)))) # this is all the points inside the convex hull
    
    
    
    return(tissueMask)
    
    
    
}