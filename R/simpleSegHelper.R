#' @importFrom EBImage Image abind
#' @importFrom methods is
#' @importFrom stats coef cor median resid runif sd
.nucSeg <- function(image,
                   nucleusIndex = 1,
                   sizeSelection = 10,
                   smooth = 1,
                   tolerance = NULL,
                   watershed = "combine",
                   ext = 1,
                   discSize = 3,
                   wholeCell = TRUE,
                   transform = NULL,
                   tissueIndex = NULL,
                   pca = FALSE) {

  ## Prepare matrix use to segment nuclei
  if ("tissueMask" %in% transform) { ## calculate tissue mask
    tissueMask <- .calcTissueMask(
      image,
      tissueIndex
    ) ## separate tissue from background

    image <- EBImage::Image(sweep(image, c(1, 2), tissueMask, "*"))
  }

  

  if (is.null(transform) == FALSE) {
    image <- .Transform(image, transform, isNuc = FALSE)
  }
  #nmask
  nuc <- .prepNucSignal(image, nucleusIndex, smooth, pca)

  ## Segment Nuclei
  nth <- EBImage::otsu(nuc, range = range(nuc))
  ## threshold on the sqrt intensities works better.

  nMask <- nuc > nth
  ## the threshold is squared to adjust of the sqrt previously.

  ## Size Selection
  nMaskLabel <- EBImage::bwlabel(nMask * 1)
  tabNuc <- table(nMaskLabel)
  nMask[nMaskLabel %in% names(which(tabNuc <= sizeSelection))] <- 0
  nMaskLabel[nMaskLabel %in% names(which(tabNuc <= sizeSelection))] <- 0

  if (watershed == "distance") {
    if (is.null(tolerance)) tolerance <- 1
    dist <- EBImage::distmap(nMask)
    wMask <-
      EBImage::watershed(dist, tolerance = tolerance, ext = ext)

    if (wholeCell) {
      kern <- EBImage::makeBrush(discSize, shape = "disc")
      wMask <- EBImage::dilate(wMask, kern)
    }
    return(wMask)
  }

  if (watershed == "combine") {
    ## Scale cell intensities
    avg <- tapply(nuc, nMaskLabel, mean)
    AVG <- nMask
    AVG[] <- avg[as.character(nMaskLabel)]
    nuc <- (nuc / AVG) * nMask
    nuc <- nuc / median(avg[as.character(nMaskLabel)]) * nMask

    dist <- EBImage::distmap(nMask)

    if (is.null(tolerance)) {
      tolerance <- .estimateTolerance(dist * nuc, nMask)
    }

    wMask <- EBImage::watershed(dist * nuc, tolerance = tolerance, ext = ext)

    if (wholeCell) {
      kern <- EBImage::makeBrush(discSize, shape = "disc")
      wMask <- EBImage::dilate(wMask, kern)
    }

    return(wMask)
  }

  ## Add distance to nuc signal
  cellRadius <- 2 * floor(sqrt(sizeSelection / pi) / 2) + 1
  nuc <- filter2(nuc, makeBrush(cellRadius, shape = "disc"))
  nuc <- nuc * nMask

  if (is.null(tolerance)) {
    tolerance <- .estimateTolerance(nuc, nMask)
  }

  wMask <- EBImage::watershed(nuc, tolerance = tolerance, ext = ext)

  ## Size Selection
  tabNuc <- table(wMask)
  wMask[wMask %in% names(which(tabNuc <= sizeSelection))] <- 0

  if (wholeCell) {
    kern <- EBImage::makeBrush(discSize, shape = "disc")
    wMask <- EBImage::dilate(wMask, kern)
  }

  wMask
}

.nucSegParallel <- function(image,
                           nucleusIndex = 1,
                           sizeSelection = 10,
                           smooth = 1,
                           tolerance = 0.01,
                           ext = 1,
                           discSize = 3,
                           wholeCell = TRUE,
                           watershed = "combine",
                           transform = NULL,
                           tissueIndex = NULL,
                           pca = FALSE,
                           BPPARAM = BiocParallel::SerialParam()) {
  output <- BiocParallel::bplapply(image,
    .nucSeg,
    nucleusIndex = nucleusIndex,
    tolerance = tolerance,
    watershed = watershed,
    ext = ext,
    discSize = discSize,
    sizeSelection = sizeSelection,
    smooth = smooth,
    wholeCell = wholeCell,
    transform = transform,
    tissueIndex = tissueIndex,
    pca = pca,
    BPPARAM = BPPARAM
  )
}

.prepNucSignal <- function(image, nucleusIndex, smooth, pca) {

  #TODO: PCA fails without smoothing, why?
  image <- apply(image, 3, function(x) {
    x <- (x)
    EBImage::gblur(x, smooth)
  }, simplify = FALSE)

  if (pca) {
    image <- EBImage::abind(image, along = 3)
    image.long <- apply(image, 3, as.numeric)

    #TODO: This should be somewhere in input validation, not here.
    ind <- intersect(nucleusIndex, colnames(image.long))

    image_nucleus <- image[, , ind]
    # if there is more than one nucluear marker, average them
    if (length(ind) > 1) image_nucleus <- apply(image_nucleus, c(1, 2), mean)
    image_nucleus <- EBImage::gblur(image_nucleus, smooth)

    otsu_thresh <- EBImage::otsu(image_nucleus, range = range(image_nucleus))

    nucleus_mask <- image_nucleus > otsu_thresh

    #TODO: figure out if 3 is a sesible default in this case
    kern <- EBImage::makeBrush(3, shape = "disc")
    cell <- EBImage::dilate(nucleus_mask, kern)

    use <- as.vector(cell)
    pca <- prcomp(image.long[use, apply(image.long, 2, sd) > 0], scale = TRUE)

    #TODO: how to select PC?
    # currently, PC is selected by highest correlation with the fist marker.
    # this is arbitrary.
    usePC <- which.max(abs(
      apply(
        pca$x, 2, cor,
        image.long[use, nucleusIndex[1]]
    )))
    PC <- pca$x[, usePC]
    PC <- PC * sign(cor(
      PC, image.long[use, nucleusIndex[1]]
    ))

    imagePC <- as.matrix(image[, , 1])
    imagePC[] <- 0
    imagePC[cell] <- PC - min(PC)

    return(imagePC)
  }

  if (is(nucleusIndex, "character")) {
    ind <- intersect(nucleusIndex, dimnames(image)[[3]])
  }

  nuc <- image[, , ind]
  if (length(ind) > 1) nuc <- apply(nuc, c(1, 2), mean)
  nuc <- EBImage::gblur(nuc, smooth)

  return(nuc - min(nuc))
}

.estimateTolerance <- function(input, nMask) {
  y <- EBImage::distmap(nMask)
  fit <- lm(as.numeric(input[y > 0]) ~ as.numeric(y[y > 0]) - 1)
  tolerance <- coef(fit)[1]
  tolerance
}

.Transform <- function(image, transform, isNuc) {
  if (isNuc) {
    for (i in seq_along(transform)) {
      image <- switch(transform[i],
        "norm99" = .norm99(image),
        "asinh" = asinh(image),
        "maxThresh" = image / max(image, na.rm = TRUE),
        "sqrt" = sqrt(image),
        "tissueMask" = image
      )
    }
    return(image)
  } else {
    for (i in seq_along(transform)) {
      image <- switch(transform[i],
        "norm99" = .norm99(image),
        "asinh" = asinh(image),
        "maxThresh" = .maxThresh(image),
        "sqrt" = sqrt(image),
        "tissueMask" = image
      )
    }

    return(image)
  }
}

.norm99 <- function(nuc) {
  nuc[nuc > quantile(nuc, 0.99, na.rm = TRUE)] <- quantile(nuc, 0.99, na.rm = TRUE)
  return(nuc)
}

.maxThresh <- function(image) {
  for (i in 1:dim(image)[3]) {
    if (max(image[, , i]) > 0) {
      image[, , i] <- image[, , i] / max(image[, , i])
    }
  }
  return(image)
}

## disc model ##
.CytSeg <- function(nmask,
                   image,
                   sizeSelection = 5,
                   smooth = 1,
                   discSize = 3,
                   transform = NULL) {
  kern <- EBImage::makeBrush(discSize, shape = "disc")

  cell <- EBImage::dilate(nmask, kern)

  disk <- cell - nmask > 0

  ## normalization
  if (is.null(transform) == FALSE) {
    image <- .Transform(image, transform, isNuc = FALSE)
  }

  longImage_disk <- data.frame(apply(
    asinh(image),
    3, as.vector
  ),
  disk = as.vector(disk)
  )

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
  nmask[nuc_label %in% names(which(tnuc <= sizeSelection))] <- 0


  cmask4 <- EBImage::propagate(cytpred, nmask, cell)
  justdisk <- EBImage::propagate(disk, nmask, cell)


  return(EBImage::Image(cmask4))
}

## Cyt seg parallel ##
.cytSegParallel <- function(nmask,
                           image,
                           sizeSelection = 5,
                           smooth = 1,
                           discSize = 3,
                           transform = NULL,
                           BPPARAM = BiocParallel::SerialParam()) {
  test.masks.cyt <- BiocParallel::bpmapply(.CytSeg,
    nmask,
    image,
    MoreArgs = list(
      sizeSelection = sizeSelection,
      smooth = smooth,
      discSize = discSize,
      transform = transform
    ),
    BPPARAM = BPPARAM
  )
}

## Marker Model ## Cyt segmentation based on a specified cytoplasmic marker ##
.CytSeg2 <- function(nmask,
                    image,
                    channel = 2,
                    sizeSelection = 5,
                    smooth = 1,
                    transform = c("maxThresh", "asinh")) {
  cytpred <- EBImage::Image(apply(image[, , channel], c(1, 2), mean))

  if (is.null(transform) == FALSE) {
    image <- .Transform(image,
      transform,
      isNuc = FALSE
    )
  }


  cytpredsmooth <- EBImage::gblur(cytpred, sigma = smooth)


  longImage <- data.frame(apply(asinh(image), 3, as.vector),
    cytpredsmooth = as.vector(cytpredsmooth)
  )
  fit <- lm(cytpredsmooth ~ ., longImage)
  ## using all the other variables (staining channels) to predict cytpred

  cytpredpred <- cytpred
  cytpredpred[] <- terra::predict(fit, longImage)
  cytpredpred <- cytpredpred - min(cytpredpred)
  cytpredpred <- cytpredpred / max(cytpredpred)

  cellTh <- EBImage::otsu(cytpredpred, range = c(0, 1))
  cell <- cytpredpred > cellTh

  cell <- cell + nmask > 0

  nuc_label <- EBImage::bwlabel(nmask)
  tnuc <- table(nuc_label)
  nmask[nuc_label %in% names(which(tnuc <= sizeSelection))] <- 0


  cmask4 <- EBImage::propagate(cytpredpred, nmask, cell)

  return(EBImage::Image(cmask4))
}

## Marker model Parallel ##
.cytSeg2Parallel <- function(nmask,
                            image,
                            channel = 2,
                            sizeSelection = 5,
                            smooth = 1,
                            transform = c("maxThresh", "asinh"),
                            BPPARAM = BiocParallel::SerialParam()) {
  test.masks.cyt <-
    BiocParallel::bpmapply(.CytSeg2,
      nmask,
      image,
      MoreArgs = list(
        channel = channel,
        sizeSelection = sizeSelection,
        smooth = smooth,
        transform = transform
      ),
      BPPARAM = BPPARAM
    )
}
