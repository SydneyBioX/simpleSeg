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
    if(whole_cell){
        return(cell1)
    }
    else{
        return(nmask1)
    }
    
    
}
simpleSegParalell <- function(image,
                              nucleus_index = 1,
                              size_selection = 10,
                              smooth = 1,
                              tolerance = 0.01,
                              ext = 1,
                              whole_cell = TRUE,
                              cores = 50){
    output <- BiocParallel::bplapply(image, simpleSeg, nucleus_index = nucleus_index, tolerance = tolerance,  size_selection = size_selection, smooth = smooth, whole_cell = whole_cell,
                                       BPPARAM  = BiocParallel::MulticoreParam(workers = cores))
}

cytSeg <- function(nmask,
                           image,
                           size_selection = 5,
                           smooth = 1,
                           tolerance = 0.01,
                   kernSize = 3
){
    
    
    kern = makeBrush(kernSize, shape='disc')
    
    cell = dilate(nmask, kern)
    
    
    disk = cell - nmask > 0
    
    
    
    longImage_disk <- data.frame(apply(asinh(image),3, as.vector), disk = as.vector(disk))
    #longImage_background <- data.frame(apply(asinh(background),3, as.vector), disk = FALSE)
    #long_Image_all_disk <- longImage_disk %>% dplyr::filter(disk == TRUE)
    #long_Image_all_disk <- as.data.frame(long_Image_all_disk)
    
    
    #long_image_2 <- rbind(long_Image_all_disk, longImage_background)
    long_image_2 <- longImage_disk
    
    fit <- lm(disk ~ .-disk, data = long_image_2)
    
    cytpred <- nmask
    cytpred[] <- predict(fit, longImage_disk)
    cytpred <- cytpred - min(cytpred)
    cytpred <- cytpred/max(cytpred)
    
    cellTh <- otsu(cytpred,range = c(0,1))
    cell <- cytpred > cellTh
    
    cell <- cell + nmask > 0
    
    nuc_label <- bwlabel(nmask) 
    tnuc <- table(nuc_label)
    nmask[nuc_label%in%names(which(tnuc<=size_selection))] <- 0
    
    
    cmask4 <- propagate(cytpred, nmask, cell)
    justdisk <- propagate(disk, nmask, cell)
    
    #output<-list(cmask4,cmaskdisk, cell1, disk, justdisk)
    return(cmask4)
    
}
cytSegParalell <- function(nmask,
                           image,
                           size_selection = 5,
                           smooth = 1,
                           tolerance = 0.01,
                           kernSize = 3,
                           cores = 50){
    test.masks.cyt <- mcmapply(cytSeg, nmask, image, size_selection, mc.cores = 40)
}