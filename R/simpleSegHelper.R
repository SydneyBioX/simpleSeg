


## Autosmooth ##

# autosmooth <- function(channel, smooth, threshold, adjustment){ signal <-
# mean(channel[,,9]) noise <- sd(channel[,,9]) SNR<-10*log10(signal/noise) if
# (abs(SNR) > threshold){ smooth <- smooth + adjustment return(smooth) } else{
# return(smooth) } }


## Segmentation functinos ##

## NucSeg ##
#' @importFrom EBImage Image
nucSeg <- function(image,
                   nucleus_index = 1,
                   size_selection = 10,
                   smooth = 1,
                   tolerance = 0.01,
                   watershed = "combine",
                   ext = 1,
                   discSize = 3,
                   wholeCell = TRUE,
                   transform = NULL,
                   tissueIndex = NULL) {
  
  
  
  
  # Prepare matrix use to segment nuclei
  if ("tissueMask" %in% transform){ # calculate tissue mask
    tissueMask <- calcTissueMask(image,
                                   tissueIndex) # separate tissue from background
    
    image <- EBImage::Image(sweep(image, c(1,2), tissueMask, "*"))
    
  }

  nuc <- .prepNucSignal(image, nucleus_index, smooth)
  

  if (is.null(transform) == FALSE) nuc <- .Transform(nuc, transform)
  
  # Segment Nuclei
  nth <-
    EBImage::otsu(nuc, range = range(nuc))  # thresholding on the sqrt intensities works better.
  nMask <-
    nuc > nth  # the threshold is squared to adjust of the sqrt previously.
  
  
  # Size Selection
  nMaskLabel <- EBImage::bwlabel(nMask*1)
  tabNuc <- table(nMaskLabel)
  nMask[nMaskLabel %in% names(which(tabNuc <= size_selection))] <- 0  
  nMaskLabel[nMaskLabel %in% names(which(tabNuc <= size_selection))] <- 0 
  
  
  if(watershed == "distance"){
    if(is.null(tolerance)) tolerance <- 1
    dist <- EBImage::distmap(nMask)
    wMask <-
      EBImage::watershed(dist, tolerance = tolerance, ext = ext)
    
    if(wholeCell){
      kern <- EBImage::makeBrush(discSize, shape = "disc")
      wMask <- EBImage::dilate(wMask, kern)
    }
    
    return(wMask)
    
  }
  
  if(watershed == "combine"){
    
    # Scale cell intensities
    avg <- tapply(nuc, nMaskLabel, mean)
    AVG <- nMask
    AVG[] <- avg[as.character(nMaskLabel)]
    nuc <- (nuc/AVG)*nMask
    nuc <- nuc/median(avg[as.character(nMaskLabel)])*nMask
    #nuc <- nuc/sd(nuc[nuc!=0])/2
    
    dist <- EBImage::distmap(nMask)
    
    if(is.null(tolerance)){
      tolerance <- .estimateTolerance(dist*nuc, nMask)
    }
    
    wMask <-
      EBImage::watershed(dist*nuc, tolerance = tolerance, ext = ext)
    
    if(wholeCell){
      kern <- EBImage::makeBrush(discSize, shape = "disc")
      wMask <- EBImage::dilate(wMask, kern)
    }
    
    return(wMask)
  }
  
  
  # Add distance to nuc signal
  cellRadius <- 2*floor(sqrt(size_selection/pi)/2)+1
  nuc <- filter2(nuc, makeBrush(cellRadius, shape='disc'))
  nuc <- nuc * nMask
  
  if(is.null(tolerance)){
    tolerance <- .estimateTolerance(nuc, nMask)
  }
  
  wMask <-
    EBImage::watershed(nuc, tolerance = tolerance, ext = ext)
  
  # Size Selection
  tabNuc <- table(wMask)
  wMask[wMask %in% names(which(tabNuc <= size_selection))] <- 0  
  
  if(wholeCell){
    kern <- EBImage::makeBrush(discSize, shape = "disc")
    wMask <- EBImage::dilate(wMask, kern)
  }
  
  
  wMask
  
}


## Nuc Seg Parallel ##

nucSegParallel <- function(image,
                           nucleus_index = 1,
                           size_selection = 10,
                           smooth = 1,
                           tolerance = 0.01,
                           ext = 1,
                           discSize = 3,
                           wholeCell = TRUE,
                           watershed = "combine",
                           transform = NULL,
                           tissueIndex = NULL,
                           BPPARAM = BiocParallel::SerialParam()) {
  output <- BiocParallel::bplapply(
    image,
    nucSeg,
    nucleus_index = nucleus_index,
    tolerance = tolerance,
    watershed = watershed,
    ext = ext,
    discSize = discSize,
    size_selection = size_selection,
    smooth = smooth,
    wholeCell = wholeCell,
    transform = transform,
    tissueIndex = tissueIndex,
    BPPARAM = BPPARAM
  )
}



.prepNucSignal <- function(image, nucleus_index, smooth){
  
  #Default, assuming nucleus_index is an integer
  ind <- nucleus_index

  if("PCA" %in% nucleus_index){
    image <- apply(image, 3, function(x){
      x <- (x)
      EBImage::gblur(x, smooth)
    }, simplify = FALSE)
    
    image <- abind(image, along = 3)
    
    image.long <- apply(image,3, as.numeric)
    pca <- prcomp(image.long[, apply(image.long, 2, sd)>0])
    
    usePC <- 1
    if(any(nucleus_index%in%colnames(image.long))){
      ind <- intersect(nucleus_index, colnames(image.long))
      usePC <- which.max(abs(apply(pca$x, 2, cor, image.long[,nucleus_index[nucleus_index != "PCA"][1]])))
      
      PC <- pca$x[,usePC]
      PC <- PC*sign(cor(PC, image.long[,nucleus_index[nucleus_index != "PCA"][1]]))
    }
    else{
      PC <- pca$x[,usePC]
    }
  
    imagePC <- as.matrix(image[,,1])
    imagePC[] <- PC - min(PC)
    return(imagePC)
  }
  
  if(is(nucleus_index, "character"))
    ind <- intersect(nucleus_index, dimnames(image)[[3]])
  
  nuc <- image[, , ind]
  if(length(ind)>1) nuc <- apply(nuc, c(1,2), mean)
  nuc <- EBImage::gblur(nuc, smooth)
  
  nuc
  
}


.estimateTolerance <- function(input, nMask){
  y <- EBImage::distmap(nMask)
  fit <- lm(as.numeric(input[y>0]) ~ as.numeric(y[y>0])-1)
  tolerance <- coef(fit)[1]
  tolerance
  # tolerance <- sd(as.numeric(input[y>0]))/sd(as.numeric(y[y>0]))
  # tolerance
}

.Transform <- function(nuc, transform) {
  
  for (i in 1:length(transform)){
    nuc <- switch(transform[i],
                  "norm99" = .norm99(nuc),
                  "asinh" = asinh(nuc),
                  "maxThresh" = nuc / max(nuc, na.rm = TRUE),
                  "sqrt" = sqrt(nuc),
                  "tissueMask" = nuc
    )
  }
  
  return(nuc)
}

.norm99 <- function(nuc){
  nuc[nuc > quantile(nuc, 0.99, na.rm = TRUE)] <- quantile(nuc, 0.99, na.rm = TRUE)
  return(nuc)
}


## disc model ##

CytSeg <- function(nmask,
                   image,
                   size_selection = 5,
                   smooth = 1,
                   discSize = 3,
                   normalize = c("maxThresh", "asinh")) {
  kern <- EBImage::makeBrush(discSize, shape = "disc")
  
  cell <- EBImage::dilate(nmask, kern)
  
  
  disk <- cell - nmask > 0
  
  # normalization
  image1 <- image
  test <- NULL
  
  if ("maxThresh" %in% normalize) {
    for (i in 1:dim(image)[3]) {
      image[, , i] <- image[, , i] / max(image[, , i])
    }
  }
  if ("sqrt" %in% normalize) {
    image <- sqrt(image)
  }
  
  
  if ("asinh" %in% normalize) {
    image <- asinh(image)
  }
  if (is.null(smooth) ==
      FALSE) {
    image <- gblur(image, sigma = smooth)
  }
  
  
  
  longImage_disk <- data.frame(apply(asinh(image),
                                     3, as.vector),
                               disk = as.vector(disk))
  
  long_image_2 <- longImage_disk
  
  fit <- lm(disk ~ . - disk, data = long_image_2)
  
  cytpred <- nmask
  cytpred[] <- terra::predict(fit, longImage_disk)
  cytpred <- cytpred - min(cytpred)
  cytpred <- cytpred / max(cytpred)
  
  cellTh <- EBImage::otsu(cytpred, range = c(0, 1))
  cell <- cytpred > cellTh
  
  cell <- cell + nmask > 0
  
  nuc_label <- EBImage::bwlabel(nmask)
  tnuc <- table(nuc_label)
  nmask[nuc_label %in% names(which(tnuc <= size_selection))] <- 0
  
  
  cmask4 <- EBImage::propagate(cytpred, nmask, cell)
  justdisk <- EBImage::propagate(disk, nmask, cell)
  
  
  return(EBImage::Image(cmask4))
  
}


## Cyt seg parallel ##
cytSegParallel <- function(nmask,
                           image,
                           size_selection = 5,
                           smooth = 1,
                           discSize = 3,
                           normalize = c("maxThresh", "asinh"),
                           BPPARAM = BiocParallel::SerialParam()) {
  test.masks.cyt <- BiocParallel::bpmapply(
    CytSeg,
    nmask,
    image,
    MoreArgs = list(
      size_selection = size_selection,
      smooth = smooth,
      discSize = discSize,
      normalize = normalize
    ),
    BPPARAM = BPPARAM
  )
}


## Marker Model ## Cyt segmentation based on a specified cytoplasmic marker ##

CytSeg2 <- function(nmask,
                    image,
                    channel = 2,
                    size_selection = 5,
                    smooth = 1,
                    normalize = c("maxThresh", "asinh")) {
 
  
  cytpred <- EBImage::Image(apply(image[, , channel], c(1, 2),
                  mean))
  
  
  if ("maxThresh" %in% normalize) {
    cytpred <- cytpred / max(cytpred)
  }
  
  
  if ("asinh" %in% normalize) {
    cytpred <- asinh(cytpred)
  }
  if ("sqrt" %in% normalize) {
    cytpred <- sqrt(cytpred)
  }
  
  
  cytpredsmooth <- EBImage::gblur(cytpred, sigma = smooth)
  
  
  longImage <- data.frame(apply(asinh(image),
                                3, as.vector),
                          cytpredsmooth = as.vector(cytpredsmooth))
  fit <-
    lm(cytpredsmooth ~ ., longImage)  #using all the other variables (staining channels) to predict cytpred
  
  cytpredpred <- cytpred
  cytpredpred[] <- terra::predict(fit, longImage)
  cytpredpred <- cytpredpred - min(cytpredpred)
  cytpredpred <- cytpredpred / max(cytpredpred)
  
  cellTh <- EBImage::otsu(cytpredpred, range = c(0, 1))
  cell <- cytpredpred > cellTh
  
  cell <- cell + nmask > 0
  
  nuc_label <- EBImage::bwlabel(nmask)
  tnuc <- table(nuc_label)
  nmask[nuc_label %in% names(which(tnuc <= size_selection))] <- 0
  
  
  cmask4 <- EBImage::propagate(cytpredpred, nmask, cell)
  
  return(EBImage::Image(cmask4))
}



## Marker model Parallel ##
cytSeg2Parallel <- function(nmask,
                            image,
                            channel = 2,
                            size_selection = 5,
                            smooth = 1,
                            normalize = c("maxThresh", "asinh"),
                            BPPARAM = BiocParallel::SerialParam()) {
  test.masks.cyt <- BiocParallel::bpmapply(
    CytSeg2,
    nmask,
    image,
    MoreArgs = list(
      channel = channel,
      size_selection = size_selection,
      smooth = smooth,
      normalize = normalize
    ),
    BPPARAM = BPPARAM
  )
}
