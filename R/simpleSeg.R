#' Perform simple segmentation
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
#' @param autosmooth dynamically scales smoothing based on signal to noise ratio of individual images
#' @param tolerance
#' @param ext = 1,
#' @param kernSize size of dilation around nuclei to create cell disk

#' @param size_selectionCyt
#' @param minMax #scale image channel intensities between 0 and 1
# asin = FALSE #perform asinh normalization on image channels

#' @param cyt_index index of the cytoplasm marker channel. Use if cytidentification = markerModel

#' @param cores = 50 number or cores for paralell processing
#'
#' @return A list of image masks
#'
#' @examples
#' 
#' 1+1
#'
#' @export
#' @rdname simpleSeg
#' @importFrom BiocParallel SerialParam bplapply MulticoreParam
#' @importFrom EBImage gblur otsu bwlabel makeBrush filter2 watershed dilate distmap propagate
#' @importFrom parallel mcmapply
#' @importFrom terra predict
#' @importFrom cytomapper CytoImageList
#' @importFrom stats prcomp quantile lm 
#' @importFrom BiocSingular RandomParam
#' @importFrom BiocNeighbors AnnoyParam

autosmooth <- function(channel, smooth, threshold, adjustment){
    signal <- mean(channel[,,9])
    noise <- sd(channel[,,9])
    SNR<-10*log10(signal/noise)
    if (abs(SNR) > threshold){
        smooth <- smooth + adjustment
        return(smooth)
    }
    else{
        return(smooth)
    }
}

nucSeg <- function(image,
                          nucleus_index = 1,
                          size_selection = 10,
                          smooth = 1,
                          norm99 = TRUE,
                          maxThresh = TRUE,
                          autosmooth = TRUE,
                          tolerance = 0.01,
                          ext = 1,
                          kernSize = 3,
                          whole_cell = TRUE){
    
    
    
    nuc <- image[,,nucleus_index]
    
    if (autosmooth){
        smooth <- autosmooth(image, smooth, 9, 4) # adjusting the smoothing parameter for low intensity images
    }
    
    
    # Hotspot filtering. Intensities greater than the 99 percentile are changed to be exactly the 99th percentile intensity. This has the effect of removing outliers.
    
    if (norm99){
        nuc[nuc > quantile(nuc, 0.99,na.rm = TRUE)] <- quantile(nuc,
                                                                0.99,
                                                                na.rm = TRUE)
    }
    if (maxThresh){
        nuc <- nuc/max(nuc,na.rm = TRUE)
    }
    
    
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
    kern = makeBrush(kernSize, shape='disc')
    cell1 = dilate(nmask1, kern)
    disk1 = cell1-nmask1 >0
    disk1 <- watershed(disk1)
    #output<-list(nmask1,nuc1, disk_blur, disk1, cell1)
    #return(output)
    if(whole_cell){
        return(cell1)
    }
    else{
        return(nmask1)
    }
    
    
}
nucSegParalell <- function(image,
                           nucleus_index = 1,
                           size_selection = 10,
                           smooth = 1,
                           norm99 = TRUE,
                           maxThresh = TRUE,
                           autosmooth = TRUE,
                           tolerance = 0.01,
                           ext = 1,
                           kernSize = 3,
                           whole_cell = TRUE,
                           cores = 50){
    output <- BiocParallel::bplapply(image, nucSeg, nucleus_index = nucleus_index, tolerance = tolerance, ext = ext, kernSize = kernSize, size_selection = size_selection, smooth = smooth, norm99 = norm99,maxThresh = maxThresh, autosmooth = autosmooth,  whole_cell = whole_cell,
                                     BPPARAM  = BiocParallel::MulticoreParam(workers = cores))
}

## cytoplasm segmentation based on markers found in the nuclei disk

CytSeg <- function(nmask,
                   image,
                   size_selection = 5,
                   smooth = 1,#Does not appear to do anything here
                   kernSize = 3,
                   minMax = FALSE,
                   asin = FALSE){
    
    
    kern = makeBrush(kernSize, shape='disc')
    
    cell = dilate(nmask, kern)
    
    
    disk = cell - nmask > 0
    
    #normalization
    image1 <- image
    test <- NULL
    
    if (minMax){
        for (i in 1:dim(image)[3]){
            image[,,i] <- image[,,i]/max(image[,,i])
        }
    }
    
    
    if (asin){
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
                           kernSize = 3,
                           minMax = FALSE,
                           asin = FALSE,
                           cores = 50){
    test.masks.cyt <- mcmapply(CytSeg, nmask, image, size_selection = size_selection, smooth = smooth, kernSize = kernSize, minMax = minMax, asin = asin, mc.cores = cores)
}


## Cyt segmentation based on a specified cytoplasmic marker ##

CytSeg2 <- function(nmask,
                    image,
                    channel = 2,
                    size_selection = 5,
                    smooth = 1, #does not appear to do anything here
                    minMax = FALSE,
                    asin = FALSE){
    CD44 <- asinh(image[,,channel])/asinh(max(image[,,channel])) #CD44 is the target protein for this channel
    
    if (minMax){
        CD44 <- CD44/max(CD44)
    }
    
    
    if (asin){
        CD44 <- asinh(CD44)
    }
    
    
    CD44smooth <- gblur(CD44, sigma = smooth)
    
    
    longImage <- data.frame(apply(asinh(image),3, as.vector), CD44smooth = as.vector(CD44smooth))
    fit <- lm(CD44smooth ~ ., longImage) #using all the other variables (staining channels) to predict CD44
    
    CD44pred <- CD44
    CD44pred[] <- terra::predict(fit, longImage)
    CD44pred <- CD44pred - min(CD44pred)
    CD44pred <- CD44pred/max(CD44pred)
    
    cellTh <- otsu(CD44pred,range = c(0,1))
    cell <- CD44pred > cellTh
    
    cell <- cell + nmask > 0
    
    nuc_label <- bwlabel(nmask) 
    tnuc <- table(nuc_label)
    nmask[nuc_label%in%names(which(tnuc<=size_selection))] <- 0
    
    
    cmask4 <- propagate(CD44pred, nmask, cell)
    
    return(cmask4)
}



cytSeg2Paralell <- function(nmask,
                            image,
                            channel = 2,
                            size_selection = 5,
                            smooth = 1,
                            minMax = FALSE,
                            asin = FALSE,
                            cores = 50){
    test.masks.cyt <- mcmapply(CytSeg2, nmask, image, channel = channel, size_selection = size_selection, smooth = smooth,  mc.cores = 40)
}


### Simple seg
simpleSeg <- function(#nmask parameters
    image,
    cytIdentification = "dilate",
    nucleus_index = 1,
    size_selectionNuc = 10,
    smooth = 1,
    norm99 = TRUE,
    maxThresh = TRUE,
    autosmooth = TRUE,
    tolerance = 0.01,
    ext = 1,
    kernSize = 3,
    
    #cyt1 parameters
    #nmask,
    
    size_selectionCyt = 5,
    minMax = FALSE,
    asin = FALSE,
    
    #cyt2 parameters
    #nmask,
    cyt_index=2,
    
    cores = 50
){
    # do nmask (if cytIdentification is null return nuc mask)
    whole_cell = FALSE
    if (cytIdentification == "dilate"){
        whole_cell = TRUE
    }
    
    nmask <- nucSegParalell( image,
                             nucleus_index = nucleus_index,
                             size_selection = size_selectionNuc,
                             smooth = smooth,
                             norm99 = norm99,
                             maxThresh = maxThresh,
                             autosmooth = autosmooth,
                             tolerance = tolerance,
                             ext = ext,
                             whole_cell = whole_cell,
                             kernSize = kernSize,
                             cores = cores)
    
    #if dilate
    if (cytIdentification == "dilate"){
        return(CytoImageList(nmask))
    }
    
    if (cytIdentification == "diskModel"){
        
        cells <- cytSegParalell (nmask,
                                 image,
                                 size_selection = size_selectionCyt,
                                 smooth = smooth,
                                 kernSize = kernSize,
                                 minMax = minMax,
                                 asin = asin,
                                 cores = cores)
        
        return(CytoImageList(cells))
    }
    
    #if marker
    else if (cytIdentification == "markerModel"){
        
        cells <- cytSeg2Paralell(nmask,
                                 image,
                                 channel = cyt_index,
                                 size_selection = size_selectionCyt,
                                 smooth = smooth,
                                 minMax = minMax,
                                 asin = asin,
                                 cores = cores)
        
        return(CytoImageList(cells))
    }
}






normalize.cells <- function(cells,
                            markers,
                            isSCE = FALSE,
                            assayName = NULL,
                            imageNb = NULL, #requires the image numbering varlable from sce list data to stratify the cells on
                            transformation = NULL,
                            method = NULL){
    
    #handeling sce and se
    if (isSCE == TRUE){
        if (is.null(assayName)){
            cellsdf <- data.frame(t(assay(cells)))
        }
        else{
            cellsdf <- data.frame(t(cells@assays@data@listData[[assayName]]))
        }
        cellsdf$imageID <- cells@colData@listData[[imageNb]]
        cells <- cellsdf
    }
    
    #Transformations
    if (transformation == "asinh"){
        cells[,markers] <- asinh(cells[,markers])
    }
    if (transformation == "sqrt"){
        cells[,markers] <- sqrt(cells[,markers])
    }
    if (transformation == "log"){
        cells[,markers] <- log10(cells[,markers])
    }
    #Methods
    if (method == "meandiv"){
        for(i in 1:length(unique(cells$imageID))){
            cells[cells$imageID==i,markers] <- sweep(cells[cells$imageID==i,markers], 2, apply(cells[cells$imageID==i,markers],2,mean, 0.2), "/")
        }
    }
    if (method == "99perc"){
        cells[,markers] <- data.frame(apply(cells[,markers],
                                            2, 
                                            function(x){
                                                q <- quantile(x,0.99)
                                                pmin(x,q)/q
                                            }))
    }
    if (method == "1stPC"){
        pca <- prcomp(cells[,markers])
        
        PC1 <- pca$x[,"PC1"]
        
        
        for(i in markers){
            y <- cells[,i]
            fit <- lm(y~PC1)
            int <- fit$coefficients["(Intercept)"]
            cells[,i] <- pmax(resid(fit) + int,0)
            q <- quantile(cells[,i],0.99)
            cells[,i] <- pmin(cells[,i],q)/q
        }
        
        
    }
    if (method == "scMerge"){
        #sc merge
        ncores <- 64
        use_bpparam <- BiocParallel::MulticoreParam(workers = ncores)
        use_bsparam <- BiocSingular::RandomParam()
        use_bnparam <- BiocNeighbors::AnnoyParam()
        # dat_sub <- dat[sample(nrow(dat), 500000), ]
        # ctl_genes <- rownames(sce)
        exprsMat <- t(cells[, markers])
        colnames(exprsMat) <- seq_len(ncol(exprsMat))
        scMerge_res <- scMerge2(exprsMat = exprsMat, #the exprs matrix to be normalised
                                batch = cells$imageID, # batch labels
                                cellTypes = NULL, # set NULL clustering will be performed within scMerge2... can also try the published cell type labels, which will match between the cell types
                                use_bpparam = use_bpparam,
                                use_bsparam = use_bsparam,
                                use_bnparam = use_bnparam,
                                ruvK = 2, # Number of unwanted variation to be removed
                                ctl = markers, # negative control genes
                                k_psuedoBulk = 5, # Number of pseudo bulk to be created for each cell type each batch
                                k_celltype = 20, # Number of neighbours when perform graph clustering
                                pseudoBulk_fn = create_pseudoBulk, # ways of constructing pseudo bulk
                                ncores = ncores,
                                chosen.hvg = markers, #Highly variable genes to be used to identify pseudo-replicates... since IMC has very few features, using all features.
                                cosineNorm = F,
                                return_subset = FALSE,
                                normalised = T)
        dat_norm <- as.data.frame(scMerge_res$newY)
        cells[, markers] <- dat_norm
    }
    return(cells)
}



