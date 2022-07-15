#' Perform simple segmentation of multiplexed cellular images
#'
#' @param image An image or list of images, can be a CytoImageList
#' @param mask A mask or list of masks, can be a CytoImageList corresponding to the given image/s
#'
#' @return A matrix with individual cells as rows and features: mean marker intensities, x/y coordinates, area and image number as columns

#' @examples
#'
#' cell.matrix <- createCells(mask.list, image.list)
#'
#' @export createCells
#' @rdname createCells
#'
#' @importFrom cytomapper CytoImageList
#' @importFrom BiocParallel bpmapply
#' @importFrom EBImage computeFeatures.basic computeFeatures.moment computeFeatures.shape






################## Create cells ########################


createCells <- function(mask, image, cores = 1) {
    
    BPPARAM <- generateBPParam(cores)
    
    maskJ <- cytomapper::CytoImageList(mask)
    stackJ <- cytomapper::CytoImageList(image)
    
    imageID <- names(image)
    
    cells.list <- BiocParallel::bpmapply(calc_features,
                                         stackJ,
                                         maskJ,
                                         BPPARAM = BPPARAM)
    cells.list <- lapply(cells.list, as.data.frame)
    cells <- do.call("rbind", cells.list)
    
    cells
}

################################ calculate features
################################ #########################################

calc_features <- function(image, mask) {
    ncells <- length(names(table(mask))) -
        1  # number of cells in image. -1 since it includes 0.
    
    features <- matrix(nrow = ncells,
                       ncol = dim(image)[3],
                       byrow = TRUE)  # create a matrix where the mean intensities of each marker will be placed
    
    for (j in 1:dim(image)[3]) {
        features[, j] <-
            EBImage::computeFeatures.basic(mask, image[, , j])[1:ncells]
    }
    
    features <- cbind(
        features,
        as.numeric(EBImage::computeFeatures.moment(mask)[,1]),
        as.numeric(EBImage::computeFeatures.moment(mask)[,2])
    )
    features <- cbind(features, as.numeric(EBImage::computeFeatures.shape(mask)[,1]))
    
    return(features)
    
}
