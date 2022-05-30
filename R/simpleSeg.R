#' Perform simple segmentation
#'
#' @param image An image
#' @param BPPARAM A BiocParallelParam object.
#'
#' @return A list of image masks
#'
#' @examples
#' 
#' 1+1
#'
#' @export
#' @rdname simpleSeg
#' @importFrom BiocParallel SerialParam bplapply

simpleSeg <- function(image,
                          nucleus_index = 1,
                          size_selection = 10,
                          smooth = 1,
                          tolerance = 0.01,
                          ext = 1){
    
    
    
    nuc <- image[,,nucleus_index]
    
    # Hotspot filtering. Intensities greater than the 99 percentile are changed to be exactly the 99th percentile intensity. This has the effect of removing outliers.
    nuc[nuc > quantile(nuc, 0.99,na.rm = TRUE)] <- quantile(nuc,
                                                            0.99,
                                                            na.rm = TRUE)
    nuc <- nuc/max(nuc,na.rm = TRUE)
    
    nuc1 <- nuc
    nuc1[is.na(nuc)] <- 0
    
    nuc1 <- gblur(sqrt(nuc1),
                  smooth)
    
    
    # Otsu thresholding. 
    nth <- otsu((nuc1),
                range = c(0,1)) # thresholding on the sqrt intensities works better. 
    nmask = nuc1 >nth # the threshold is squared to adjust of the sqrt previously.
    
    
    
    #### Watershed to segment nucleus.
    nuc1 <- nuc # make nuc a matrix
    
    nuc1[!nmask|is.na(nmask)] <- 0 #get the intensities values in the nmask, and set other values to 0 that are not part of nmask.
    
    
    # Size selection
    nMaskLabel <- bwlabel(nmask) 
    tnuc1 <- table(nMaskLabel)
    
    
    nuc1[nMaskLabel%in%names(which(tnuc1<=size_selection))] <- 0 # sizes less than 10 are set to zero.
    
    nmask[nMaskLabel%in%names(which(tnuc1<=size_selection))] <- 0 # sizes less than 10 are set to zero.
    
    
    
    kern = makeBrush(5, shape='disc')
    
    
    disk_blur <- filter2(sqrt(nuc1), kern)
    #disk_blur <- disk_blur/
    
    #creating a distance matrix
    #distNuc <- EBImage::distmap(nuc1)
    
    # water shed to segment
    #nuc1 <- (nuc1>0)*1
    #nuc1 <- EBImage::distmap(nuc1)
    nmask1 <- watershed(disk_blur * nmask,
                        tolerance = tolerance,
                        ext = ext)
    #nMaskLabel <- bwlabel(nmask1) 
    #nMaskLabel <- colorLabels(nmask1) 
    #display(nMaskLabel)
    kern = makeBrush(3, shape='disc')
    cell1 = dilate(nmask1, kern)
    disk1 = cell1-nmask1 >0
    disk1 <- watershed(disk1)
    output<-list(nmask1,nuc1, disk_blur, disk1, cell1)
    return(output)
    
    
}