#' Perform simple segmentation of multiplexed cellular images
#'
#' @param image An image
#' @param BPPARAM A BiocParallelParam object.
#' @param image An image or list of images or cytoimagelist to be read into the function.
#' @param cellBody method of cytoplasm identification. Can be "dilate", "diskModel" or the index / list of indexes of dedicated cytoplasm markers
#' @param nucleus the channel number or list of channel numbers corresponding to the nuclei marker/s
#' @param sizeSelection minimum pixels for an object to be recognised as signal and not noise
#' @param smooth the amount of smoothing to be applied to the nuclei marker channle
#' @param norm99perfrom 99th percentile transformation
#' @param maxThresh scale intensities between 0 and 1

#' @param tolerance The minimum height of the object in the units of image intensity between its highest point (seed) and the point where it contacts another object (checked for every contact pixel). If the height is smaller than the tolerance, the object will be combined with one of its neighbors, which is the highest. Tolerance should be chosen according to the range of x. Default value is 1, which is a reasonable value if x comes from distmap.
#' @param ext Radius of the neighborhood in pixels for the detection of neighboring objects. Higher value smoothes out small objects.
#' @param discSize size of dilation around nuclei to create cell disk #dilation size

#' @param minMax scale image channel intensities between 0 and 1
#' @param asin perform asinh normalization on image channels
#' @param transform a list containing desired normalization/transformation methods to be performed prior to cytoplasm identification, accepted values are 'mixMax' and/or 'asin'


#' @param cores = 50 number or cores for paralell processing
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
#' @importFrom EBImage gblur otsu bwlabel makeBrush filter2 watershed dilate distmap propagate
#' @importFrom terra predict
#' @importFrom cytomapper CytoImageList
#' @importFrom stats prcomp quantile lm 
#' @importFrom BiocSingular RandomParam
#' @importFrom BiocNeighbors AnnoyParam



simpleSeg <- function(#nmask parameters
    image,
    nucleus = 1,
    
    cellBody = "dilate",
    
    sizeSelection = 10,
    smooth = 1,
    transform = c("norm99", "maxThresh", "asin"),
    tolerance = 0.01,
    ext = 1,
    discSize = 3,
    
    #minMax = FALSE,
    #asin = FALSE,
    #cytNormalize = c("minMax", "asin"),#use the nucNorm here Transform instead of normalize, always do minMax AFTER other transformations
    
    
    #cyt2 parameters
    #cyt_index=2,
    
    cores = 5
){
    # do nmask (if cellBody is null return nuc mask)
    
    if (!(cellBody %in% c("none", "dilate", "discModel"))){
        if (is.numeric(cellBody) == FALSE){
            stop("cellBody must be one of the following: 'none', 'dilate', 'discModel' or the index of a cytoplasm marker channel")
        }
    }
    whole_cell = FALSE
    if (cellBody == "dilate"){
        whole_cell = TRUE
    }
    
    nmask <- nucSegParalell( image,
                             nucleus_index = nucleus,
                             size_selection = sizeSelection,
                             smooth = smooth,
                             normalize = transform,
                             tolerance = tolerance,
                             ext = ext,
                             whole_cell = whole_cell,
                             discSize = discSize,
                             cores = cores)
    
    #if dilate
    if (cellBody == "dilate"){
        return(cytomapper::CytoImageList(nmask))
    }
    if (cellBody == "none"){
        return(cytomapper::CytoImageList(nmask))
    }
    
    if (cellBody == "discModel"){
        
        cells <- cytSegParalell (nmask,
                                 image,
                                 size_selection = sizeSelection,
                                 smooth = smooth,
                                 discSize = discSize,
                                 #minMax = minMax,
                                 #asin = asin,
                                 normalize = transform,
                                 cores = cores)
        cellList <- NULL
        for (i in 1:length(cells[1,1,])){
          cellList[[i]] <- as.Image(cells[,,i])
        }
        return(cytomapper::CytoImageList(cellList))
    }
    
    #if marker
    else if (is.numeric(cellBody)){
        
        cells <- cytSeg2Paralell(nmask,
                                 image,
                                 channel = cellBody,
                                 size_selection = sizeSelection,
                                 smooth = smooth,
                                 #minMax = minMax,
                                 #asin = asin
                                 normalize = transform,
                                 cores = cores)
        
        cellList <- NULL
        for (i in 1:length(cells[1,1,])){
          cellList[[i]] <- as.Image(cells[,,i])
        }
        
        return(cytomapper::CytoImageList(cellList))
    }
}







