#' Perform simple segmentation of multiplexed cellular images
#'
#' @param image An image
#' @param BPPARAM A BiocParallelParam object.
#' @param image An image or list of images or cytoimagelist to be read into the function.
#' @param cellBody method of cytoplasm identification. Can be 'dilate', 'diskModel' or the index / list of indexes of dedicated cytoplasm markers
#' @param nucleus the channel number or list of channel numbers corresponding to the nuclei marker/s
#' @param sizeSelection minimum pixels for an object to be recognised as signal and not noise
#' @param smooth the amount of smoothing to be applied to the nuclei marker channle
#' @param norm99perfrom 99th percentile transformation
#' @param maxThresh scale intensities between 0 and 1
#' @param watershed Method used to perform watershed.
#' @param tolerance The minimum height of the object in the units of image intensity between its highest point (seed) and the point where it contacts another object (checked for every contact pixel). If the height is smaller than the tolerance, the object will be combined with one of its neighbors, which is the highest. Tolerance should be chosen according to the range of x. Default value is 1, which is a reasonable value if x comes from distmap.
#' @param ext Radius of the neighborhood in pixels for the detection of neighboring objects. Higher value smoothes out small objects.
#' @param discSize size of dilation around nuclei to create cell disk #dilation size
#' @param minMax scale image channel intensities between 0 and 1
#' @param asinh perform asinh normalization on image channels
#' @param transform a list containing desired normalization/transformation methods to be performed prior to cytoplasm identification, accepted values are 'maxThresh', 'norm99', 'sqrt' and/or 'asinh'
#' @param cores The number or cores for parallel processing
#'
#' @return A list of image masks
#'
#' @examples
#'
#' 1+1
#'
#' @export simpleSeg
#' @rdname simpleSeg
#' @importFrom BiocParallel SerialParam bplapply MulticoreParam bpmapply
#' @importFrom EBImage gblur otsu bwlabel makeBrush filter2 watershed dilate distmap propagate Image
#' @importFrom terra predict
#' @import cytomapper
#' @importFrom stats prcomp quantile lm
#' @importFrom S4Vectors mcols
simpleSeg <- function(image,
                      nucleus = 1,
                      cellBody = "dilate",
                      sizeSelection = 10,
                      smooth = 1,
                      transform = NULL,
                      watershed = "combine",
                      tolerance = NULL,
                      ext = 1,
                      discSize = 3,
                      cores = 1) {
  imageClass <- class(image)
  if (!imageClass %in% c("list", "CytoImageList")) {
    image <- list(image)
    names(image) <- "image"
  }
  # do nmask (if cellBody is null return nuc mask)
  
  if (!(cellBody %in% c("none", "dilate", "discModel"))) {
    if (is.numeric(cellBody) ==
        FALSE) {
      stop(
        "cellBody must be one of the following: 'none', 'dilate', 'discModel' or the index of a cytoplasm marker channel"
      )
    }
  }
  wholeCell <- FALSE
  if (cellBody == "dilate") {
    wholeCell <- TRUE
  }
  
  BPPARAM <- generateBPParam(cores)
  
  nmask <- nucSegParallel(
    image,
    nucleus_index = nucleus,
    size_selection = sizeSelection,
    smooth = smooth,
    watershed = watershed,
    tolerance = tolerance,
    ext = ext,
    wholeCell = wholeCell,
    discSize = discSize,
    BPPARAM = BPPARAM
  )
  
  # if dilate or none
  if (cellBody %in% c("dilate", "none")) {
    nmask <- sapply(nmask, EBImage::Image, simplify = FALSE)
    cyto.nmask <- cytomapper::CytoImageList(nmask)
    S4Vectors::mcols(cyto.nmask) <-
      S4Vectors::DataFrame(imageID = names(cyto.nmask))
    objectNum <- as.list(sapply(cyto.nmask, max))
    mcols(cyto.nmask)$objectNum <- objectNum
    return(cyto.nmask)
  }
  
  if (cellBody == "discModel") {
    cells <- cytSegParallel(
      nmask,
      image,
      size_selection = sizeSelection,
      smooth = smooth,
      discSize = discSize,
      normalize = transform,
      BPPARAM = BPPARAM
    )
    #Converting from a tiff stack to individual images
    cellList <- NULL 
    for (i in 1:length(cells[1,1,])){ 
      cellList[[i]] <-as.Image(cells[,,i]) 
    }
    
    cyto.mask <- cytomapper::CytoImageList(cellList)
    S4Vectors::mcols(cyto.mask) <-
      S4Vectors::DataFrame(imageID = names(image))
    objectNum <- as.list(sapply(cyto.mask, max))
    mcols(cyto.mask)$objectNum <- objectNum
    
    return(cyto.mask)
  }
  
  if (any(cellBody %in% dimnames(image[[1]])[[3]])) {
    cells <- cytSeg2Parallel(
      nmask,
      image,
      channel = cellBody,
      size_selection = sizeSelection,
      smooth = smooth,
      normalize = transform,
      BPPARAM = BPPARAM
    )
    
    #Converting from a tiff stack to individual images
    cellList <- NULL 
    for (i in 1:length(cells[1,1,])){ 
      cellList[[i]] <-as.Image(cells[,,i]) 
      }
    
    cyto.mask <- cytomapper::CytoImageList(cellList)
    S4Vectors::mcols(cyto.mask) <-
      S4Vectors::DataFrame(imageID = names(image))
    objectNum <- as.list(sapply(cyto.mask, max))
    mcols(cyto.mask)$objectNum <- objectNum
    return(cyto.mask)
  }
}
