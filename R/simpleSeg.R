#' Perform simple segmentation of multiplexed cellular images
#'
#' @param image An image
#' @param BPPARAM A BiocParallelParam object.
#' @param image An image or list of images or cytoimagelist to be read into the function.
#' @param cytIdentification method of cytoplasm identification. Can be "dilate", "diskModel" or "markerModel"
#' @param nucleus_index the channel number of the nuclei marker
#' @param size_selectionNuc minimum pixels for an object to be recognised as signal and not noise
#' @param smooth the amount of smoothing to be applied to the nuclei marker channle
#' @param norm99perfrom 99th percentile transformation
#' @param maxThresh scale intensities between 0 and 1
#' @param autoS dynamically scales smoothing based on signal to noise ratio of individual images
#' @param nucNormalize a list containing desired normalization/transformation methods to be performed prior to nucleus identification, accepted values are 'max Thresh' 'norm99perform' and/or 'autosmooth'
#' @param tolerance The minimum height of the object in the units of image intensity between its highest point (seed) and the point where it contacts another object (checked for every contact pixel). If the height is smaller than the tolerance, the object will be combined with one of its neighbors, which is the highest. Tolerance should be chosen according to the range of x. Default value is 1, which is a reasonable value if x comes from distmap.
#' @param ext Radius of the neighborhood in pixels for the detection of neighboring objects. Higher value smoothes out small objects.
#' @param discSize size of dilation around nuclei to create cell disk #dilation size

#' @param sizeSelectionCyt
#' @param minMax scale image channel intensities between 0 and 1
#' @param asin perform asinh normalization on image channels
#' @param cytNormalize a list containing desired normalization/transformation methods to be performed prior to cytoplasm identification, accepted values are 'mixMax' and/or 'asin'

#' @param cyt_index index of the cytoplasm marker channel. Use if cytidentification = markerModel

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
    cytIdentification = "dilate",
    nucleus_index = 1,
    size_selectionNuc = 10,
    smooth = 1,
    nucNormalize = c("norm99", "maxThresh", "autoS"),
    tolerance = 0.01,
    ext = 1,
    discSize = 3,
    
    #cyt1 parameters
    
    sizeSelectionCyt = 5,
    #minMax = FALSE,
    #asin = FALSE,
    cytNormalize = c("minMax", "asin"),
    
    #cyt2 parameters
    cyt_index=2,
    
    cores = 50
){
    # do nmask (if cytIdentification is null return nuc mask)
    
    if (!(cytIdentification %in% c("none", "dilate", "discModel", "markerModel"))){
        stop("cytIdentification must be one of the following: 'none', 'dilate', 'discModel', 'markerModel'")
    }
    whole_cell = FALSE
    if (cytIdentification == "dilate"){
        whole_cell = TRUE
    }
    
    nmask <- nucSegParalell( image,
                             nucleus_index = nucleus_index,
                             size_selection = size_selectionNuc,
                             smooth = smooth,
                             normalize = nucNormalize,
                             tolerance = tolerance,
                             ext = ext,
                             whole_cell = whole_cell,
                             discSize = discSize,
                             cores = cores)
    
    #if dilate
    if (cytIdentification == "dilate"){
        return(cytomapper::CytoImageList(nmask))
    }
    if (cytIdentification == "none"){
        return(cytomapper::CytoImageList(nmask))
    }
    
    if (cytIdentification == "discModel"){
        
        cells <- cytSegParalell (nmask,
                                 image,
                                 size_selection = sizeSelectionCyt,
                                 smooth = smooth,
                                 discSize = discSize,
                                 #minMax = minMax,
                                 #asin = asin,
                                 normalize = cytNormalize,
                                 cores = cores)
        
        return(cytomapper::CytoImageList(cells))
    }
    
    #if marker
    else if (cytIdentification == "markerModel"){
        
        cells <- cytSeg2Paralell(nmask,
                                 image,
                                 channel = cyt_index,
                                 size_selection = sizeSelectionCyt,
                                 smooth = smooth,
                                 #minMax = minMax,
                                 #asin = asin
                                 normalize = cytNormalize,
                                 cores = cores)
        
        return(cytomapper::CytoImageList(cells))
    }
}







