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
    tissue_index = NULL,
    size_selection = 10,
    cores = 1
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
    
     # replace outside of tissue with NA
    
    if(is.null(tissue_index)) tissue_index <- seq_len(dim(image)[3])
  
    tissue <- apply(image[,,tissue_index], c(1,2), mean) # add intesities of these markers to try and highlight the tissue structure.
    
    
    tissue <- EBImage::normalize(log10(tissue+0.01),
                                 ft = c(0,1))
    
    tissue_otsu <- EBImage::otsu(tissue,
                        range = c(0,1))
    
    tissue <- tissue > tissue_otsu
    
    
    # Do a size selection 
    tissue_label <- EBImage::bwlabel(tissue) 
    ttissue <- table(tissue_label)
    tissue[tissue_label%in%names(which(ttissue<=size_selection))] <- 0 # size selection of 10.
    
    # extract the convex hull
    nonZero <- as.data.frame(which(tissue>0,2))
    colnames(nonZero) <- c("y", "x")
    
    ch <- grDevices::chull(nonZero[, c("x", "y")])
    poly <- nonZero[, c("x", "y")][rev(ch),]
    colnames(poly) <- c("x", "y")
    ow <- spatstat.geom::owin(xrange = range(nonZero$x), yrange = range(nonZero$y), poly = poly)
    tissueMask <- as.matrix(spatstat.geom::as.mask(ow, xy=list(y = 1:nrow(tissue), x = 1:ncol(tissue)))) # this is all the points inside the convex hull
    
    
    
    return(tissueMask)
    
    
    
}