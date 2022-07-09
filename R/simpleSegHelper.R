## Nuc Normalization ##

nucNormalize.helper <- function(image, nuc, normalize){
    #if ("autoS" %in% normalize){
    #    smooth <- autosmooth(image, smooth, 9, 4) # adjusting the smoothing parameter for low intensity images
    #}
    
    
    # Hotspot filtering. Intensities greater than the 99 percentile are changed to be exactly the 99th percentile intensity. This has the effect of removing outliers.
    
    if ("norm99" %in% normalize){
        nuc[nuc > quantile(nuc, 0.99,na.rm = TRUE)] <- quantile(nuc,
                                                                0.99,
                                                                na.rm = TRUE)
        
    if ("asin" %in% normalize) {
      nuc <- asinh(nuc)
    }   
    }
    if ("maxThresh" %in% normalize){
        nuc <- nuc/max(nuc,na.rm = TRUE)
    }
    #output <- list(smooth, nuc)
    return(nuc)
}


## Autosmooth ##

#autosmooth <- function(channel, smooth, threshold, adjustment){
#  signal <- mean(channel[,,9])
#  noise <- sd(channel[,,9])
#  SNR<-10*log10(signal/noise)
#  if (abs(SNR) > threshold){
#    smooth <- smooth + adjustment
#    return(smooth)
#  }
#  else{
#    return(smooth)
#  }
#}


## Segmentation functinos ##

## NucSeg ##
nucSeg <- function(image,
                   nucleus_index = 1,
                   size_selection = 10,
                   smooth = 1,
                   normalize = c("norm99", "maxThresh", "asin"),
                   tolerance = 0.01,
                   ext = 1,
                   discSize = 3,
                   whole_cell = TRUE){
  
  
  
  nuc <- image[,,nucleus_index]
  
  
  nuc <- nucNormalize.helper(image, nuc, normalize)
  
  #smooth <- nucNormRes[[1]]
  
  nuc1 <- nuc
  nuc1[is.na(nuc)] <- 0
  
  
  nuc1 <- EBImage::gblur(sqrt(nuc1),
                         smooth)
  
  # Otsu thresholding. 
  nth <- EBImage::otsu((nuc1),
                       range = c(0,1)) # thresholding on the sqrt intensities works better. 
  nmask = nuc1 >nth # the threshold is squared to adjust of the sqrt previously.
  
  
  #### Watershed to segment nucleus.
  nuc1 <- nuc # make nuc a matrix
  
  nuc1[!nmask|is.na(nmask)] <- 0 #get the intensities values in the nmask, and set other values to 0 that are not part of nmask.
  
  
  # Size selection
  nMaskLabel <- EBImage::bwlabel(nmask) 
  tnuc1 <- table(nMaskLabel)
  
  
  nuc1[nMaskLabel%in%names(which(tnuc1<=size_selection))] <- 0 # sizes less than 10 are set to zero.
  
  nmask[nMaskLabel%in%names(which(tnuc1<=size_selection))] <- 0 # sizes less than 10 are set to zero.
  
  
  
  kern = EBImage::makeBrush(5, shape='disc')
  
  
  disk_blur <- EBImage::filter2(sqrt(nuc1), kern) #a
  
  nmask1 <- EBImage::watershed(disk_blur * nmask,
                               tolerance = tolerance,
                               ext = ext)
  
  kern = EBImage::makeBrush(discSize, shape='disc')
  cell1 = EBImage::dilate(nmask1, kern)
  disk1 = cell1-nmask1 >0
  disk1 <- EBImage::watershed(disk1)
  if(whole_cell){
    return(cell1)
  }
  else{
    return(nmask1)
  }
  
  
}


## Nuc Seg Parallel ##

nucSegParalell <- function(image,
                           nucleus_index = 1,
                           size_selection = 10,
                           smooth = 1,
                           normalize = c("norm99", "maxThresh", "asin"),
                           tolerance = 0.01,
                           ext = 1,
                           discSize = 3,
                           whole_cell = TRUE,
                           cores = 5){
  output <- BiocParallel::bplapply(image, nucSeg, nucleus_index = nucleus_index, tolerance = tolerance, ext = ext, discSize = discSize, size_selection = size_selection, smooth = smooth, normalize = normalize,  whole_cell = whole_cell,
                                   BPPARAM  = BiocParallel::MulticoreParam(workers = cores))
}


## disc model ##

CytSeg <- function(nmask,
                   image,
                   size_selection = 5,
                   smooth = 1,#Does not appear to do anything here
                   discSize = 3,
                   #maxThresh = FALSE,
                   #asin = FALSE
                   normalize = c("maxThresh", "asin")){
  
  
  kern = EBImage::makeBrush(discSize, shape='disc')
  
  cell = EBImage::dilate(nmask, kern)
  
  
  disk = cell - nmask > 0
  
  #normalization
  image1 <- image
  test <- NULL
  
  if ("maxThresh" %in% normalize){
    for (i in 1:dim(image)[3]){
      image[,,i] <- image[,,i]/max(image[,,i])
    }
  }
  
  
  if ("asin" %in% normalize){
    image <- asinh(image)
  }
  if (is.null(smooth) == FALSE){
    image <- gblur(image, sigma = smooth)
  }
  
  
  
  longImage_disk <- data.frame(apply(asinh(image),3, as.vector), disk = as.vector(disk))
  #longImage_background <- data.frame(apply(asinh(background),3, as.vector), disk = FALSE)
  #long_Image_all_disk <- longImage_disk %>% dplyr::filter(disk == TRUE)
  #long_Image_all_disk <- as.data.frame(long_Image_all_disk)
  
  
  #long_image_2 <- rbind(long_Image_all_disk, longImage_background)
  long_image_2 <- longImage_disk
  
  fit <- lm(disk ~ .-disk, data = long_image_2)
  
  cytpred <- nmask
  cytpred[] <- terra::predict(fit, longImage_disk)
  cytpred <- cytpred - min(cytpred)
  cytpred <- cytpred/max(cytpred)
  
  cellTh <- EBImage::otsu(cytpred,range = c(0,1))
  cell <- cytpred > cellTh
  
  cell <- cell + nmask > 0
  
  nuc_label <- EBImage::bwlabel(nmask) 
  tnuc <- table(nuc_label)
  nmask[nuc_label%in%names(which(tnuc<=size_selection))] <- 0
  
  
  cmask4 <- EBImage::propagate(cytpred, nmask, cell)
  justdisk <- EBImage::propagate(disk, nmask, cell)
  
  #output<-list(cmask4,cmaskdisk, cell1, disk, justdisk)
  return(cmask4)
  
}


## Cyt seg parallel ##
cytSegParalell <- function(nmask,
                           image,
                           size_selection = 5,
                           smooth = 1,
                           discSize = 3,
                           #maxThresh = FALSE,
                           #asin = FALSE,
                           normalize = c("maxThresh", "asin"),
                           cores = 5){
  test.masks.cyt <- BiocParallel::bpmapply(CytSeg, nmask, image, size_selection = size_selection, smooth = smooth, discSize = discSize, normalize = normalize, BPPARAM  = BiocParallel::MulticoreParam(workers = cores))
}

list1 <- c(1,2,3,4,5)
list2 <- c(6,7,8,9,10)
test.function <- function(list1, list2){
  return(list1 + list2)
}

## Marker Model ##
## Cyt segmentation based on a specified cytoplasmic marker ##

CytSeg2 <- function(nmask,
                    image,
                    channel = 2,
                    size_selection = 5,
                    smooth = 1, #does not appear to do anything here
                    #maxThresh = FALSE,
                    #asin = FALSE
                    normalize = c("maxThresh", "asin")){
  
  CD44 <- asinh(image[,,channel])/asinh(max(image[,,channel])) #CD44 is the target protein for this channel
  
  if ("maxThresh" %in% normalize){
    CD44 <- CD44/max(CD44)
  }
  
  
  if ("asin" %in% normalize){
    CD44 <- asinh(CD44)
  }
  
  
  CD44smooth <- EBImage::gblur(CD44, sigma = smooth)
  
  
  longImage <- data.frame(apply(asinh(image),3, as.vector), CD44smooth = as.vector(CD44smooth))
  fit <- lm(CD44smooth ~ ., longImage) #using all the other variables (staining channels) to predict CD44
  
  CD44pred <- CD44
  CD44pred[] <- terra::predict(fit, longImage)
  CD44pred <- CD44pred - min(CD44pred)
  CD44pred <- CD44pred/max(CD44pred)
  
  cellTh <- EBImage::otsu(CD44pred,range = c(0,1))
  cell <- CD44pred > cellTh
  
  cell <- cell + nmask > 0
  
  nuc_label <- EBImage::bwlabel(nmask) 
  tnuc <- table(nuc_label)
  nmask[nuc_label%in%names(which(tnuc<=size_selection))] <- 0
  
  
  cmask4 <- EBImage::propagate(CD44pred, nmask, cell)
  
  return(cmask4)
}



## Marker model paralell ##
cytSeg2Paralell <- function(nmask,
                            image,
                            channel = 2,
                            size_selection = 5,
                            smooth = 1,
                            #maxThresh = FALSE,
                            #asin = FALSE,
                            normalize = c("maxThresh", "asin"),
                            cores = 5){
  test.masks.cyt <- BiocParallel::bpmapply(CytSeg2, nmask, image, channel = channel, size_selection = size_selection, smooth = smooth, normalize = normalize,  BPPARAM  = BiocParallel::MulticoreParam(workers = cores))
}
