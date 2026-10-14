#' @importFrom EBImage Image
#' @importFrom methods is
#' @importFrom stats coef cor median resid runif sd
NULL

## Segment one image: nuclei, then (optionally) the cell bodies.
.segmentImage <- function(image,
                          name = "image",
                          nucleusIndex = 1,
                          cellBody = "dilate",
                          useCellMarkers = FALSE,
                          sizeSelection = 10,
                          smooth = 1,
                          tolerance = NULL,
                          watershed = "intensity",
                          ext = 1,
                          discSize = 3,
                          transform = NULL,
                          tissueIndex = NULL,
                          pca = FALSE) {
  ## Plain arrays: subsetting and arithmetic on Image objects copy the whole
  ## image through S4 methods.
  image <- as.array(image)
  cellModel <- useCellMarkers || cellBody == "discModel"
  ## The tissue mask only applies to the nucleus signal; the cell body models
  ## use the transformed image without it.
  if ("tissueMask" %in% transform) {
    tissueMask <- .calcTissueMask(image, tissueIndex, sizeSelection)
    nucImage <- .Transform(image * as.vector(tissueMask), transform)
    timage <- if (cellModel) .Transform(image, transform)
  } else {
    timage <- if (is.null(transform)) image else .Transform(image, transform)
    nucImage <- timage
  }
  if (!cellModel) rm(timage)
  if (!useCellMarkers) rm(image)

  nmask <- .nucSeg(nucImage,
    name = name,
    nucleusIndex = nucleusIndex,
    sizeSelection = sizeSelection,
    smooth = smooth,
    tolerance = tolerance,
    watershed = watershed,
    ext = ext,
    discSize = discSize,
    wholeCell = !useCellMarkers && cellBody == "dilate",
    pca = pca
  )
  rm(nucImage)
  if (!any(nmask > 0)) {
    return(EBImage::Image(nmask))
  }

  if (useCellMarkers) {
    return(.CytSeg2(nmask, image, timage,
      channel = cellBody,
      sizeSelection = sizeSelection,
      smooth = smooth
    ))
  }
  if (cellBody == "discModel") {
    return(.CytSeg(nmask, timage,
      sizeSelection = sizeSelection,
      discSize = discSize
    ))
  }
  EBImage::Image(nmask)
}

.nucSeg <- function(image,
                    name = "image",
                    nucleusIndex = 1,
                    sizeSelection = 10,
                    smooth = 1,
                    tolerance = NULL,
                    watershed = "intensity",
                    ext = 1,
                    discSize = 3,
                    wholeCell = TRUE,
                    pca = FALSE) {
  nuc <- .prepNucSignal(image, nucleusIndex, smooth, pca, discSize)
  rm(image)

  ## Segment Nuclei
  nuc[nuc < 0.0001] <- NA
  if (all(is.na(nuc))) {
    return(.noNuclei(nuc, name))
  }
  nth <- EBImage::otsu(nuc, range = range(nuc, na.rm = TRUE))
  nuc[is.na(nuc)] <- 0

  nMask <- nuc > nth

  ## Size Selection
  nMaskLabel <- EBImage::bwlabel(nMask * 1)
  small <- .smallObjects(nMaskLabel, sizeSelection)
  nMask[small] <- 0
  nMaskLabel[small] <- 0

  if (!any(nMask > 0)) {
    return(.noNuclei(nuc, name))
  }

  if (watershed == "distance") {
    if (is.null(tolerance)) tolerance <- 1
    dist <- EBImage::distmap(nMask)
    wMask <-
      .watershed(dist, tolerance = tolerance, ext = ext)

    if (wholeCell) {
      kern <- EBImage::makeBrush(discSize, shape = "disc")
      wMask <- EBImage::dilate(wMask, kern)
    }
    return(wMask)
  }

  if (watershed == "combine") {
    ## Scale cell intensities
    AVG <- .labelMeans(nuc, nMaskLabel)
    nuc <- (nuc / AVG) * nMask
    nuc <- nuc / median(AVG) * nMask

    dist <- EBImage::distmap(nMask)

    if (is.null(tolerance)) {
      tolerance <- .estimateTolerance(dist * nuc, dist, discSize = discSize)
    }

    wMask <- .watershed(dist * nuc, tolerance = tolerance, ext = ext)

    if (wholeCell) {
      kern <- EBImage::makeBrush(discSize, shape = "disc")
      wMask <- EBImage::dilate(wMask, kern)
    }

    return(wMask)
  }

  # intensity watershedding

  if (is.null(tolerance)) {
    tolerance <- .estimateTolerance(nuc, EBImage::distmap(nMask), discSize)
  }

  wMask <- .watershed(nuc * nMask, tolerance = tolerance, ext = ext)

  ## Size Selection
  wMask[.smallObjects(wMask, sizeSelection)] <- 0

  if (wholeCell) {
    kern <- EBImage::makeBrush(discSize, shape = "disc")
    wMask <- EBImage::dilate(wMask, kern)
  }

  wMask
}

## EBImage::watershed() after separating pixels of equal height by a tiny
## amount that depends on their position in the whole image (x0, y0: position of
## the first pixel), so that the result does not depend on how the watershed
## orders ties, nor on whether the image is processed whole or in tiles.
.watershed <- function(h, tolerance, ext, x0 = 1L, y0 = 1L) {
  h <- as.array(h)
  storage.mode(h) <- "double"
  h <- .breakTiesCpp(h, as.integer(x0), as.integer(y0), 1e-6 * tolerance)
  EBImage::watershed(h, tolerance = tolerance, ext = ext)
}

## An empty mask, with a warning, for an image where no nuclei pass the
## threshold and size selection.
.noNuclei <- function(nuc, name) {
  warning(
    "No nuclei found in image '", name, "'; its mask is empty. ",
    "If a few very bright pixels dominate the nucleus marker, try ",
    "transform = \"sqrt\", \"asinh\" or \"norm99\".",
    call. = FALSE
  )
  matrix(0, nrow(nuc), ncol(nuc))
}

## Pixels of a label image whose object has at most `size` pixels.
.smallObjects <- function(label, size) {
  lab <- as.integer(label)
  counts <- tabulate(lab + 1L, nbins = max(lab, 0L) + 1L)
  (counts <= size)[lab + 1L]
}

## Each pixel's mean intensity over its label (background included as label 0).
.labelMeans <- function(x, label) {
  lab <- as.integer(label)
  avg <- vapply(split(as.numeric(x), lab), mean, numeric(1))
  full <- numeric(max(lab, 0L) + 1L)
  full[as.integer(names(avg)) + 1L] <- avg
  out <- label
  out[] <- full[lab + 1L]
  out
}

.prepNucSignal <- function(image, nucleusIndex, smooth, pca, discSize) {
  if (pca) {
    ## Smoothed pixels x channels matrix, built in place.
    d <- dim(image)
    image.long <- matrix(0, d[1] * d[2], d[3],
      dimnames = list(NULL, dimnames(image)[[3]])
    )
    for (j in seq_len(d[3])) {
      image.long[, j] <- .gblur(image[, , j], smooth)
    }
    rm(image)

    ind <- if (is.character(nucleusIndex)) {
      intersect(nucleusIndex, colnames(image.long))
    } else {
      nucleusIndex
    }

    # if there is more than one nuclear marker, average them
    image_nucleus <- if (length(ind) > 1) {
      rowMeans(image.long[, ind, drop = FALSE])
    } else {
      image.long[, ind]
    }
    dim(image_nucleus) <- d[1:2]

    otsu_thresh <- EBImage::otsu(image_nucleus, range = range(image_nucleus))

    nucleus_mask <- image_nucleus > otsu_thresh

    kern <- EBImage::makeBrush(discSize, shape = "disc")
    cell <- EBImage::dilate(nucleus_mask, kern)

    use <- as.vector(cell)

    # subsample: use at most 100k pixels for the PCA
    if (length(which(use)) > 1e5L) use <- sample(which(use), min(1e5L, length(use)))

    ## markers that vary over the pixels used for the PCA
    useMarker <- vapply(seq_len(d[3]), function(j) sd(image.long[use, j]), numeric(1)) > 0
    pca <- prcomp(image.long[use, useMarker], scale = TRUE)

    usePC <- which.max(abs(
      apply(
        pca$x, 2, cor,
        image_nucleus[use]
      )
    ))
    PC <- pca$x[, usePC]
    PC_sign <- sign(cor(
      PC, image_nucleus[use]
    ))

    if (!all(useMarker)) image.long <- image.long[, useMarker]
    imagePC <- image.long %*% (pca$rotation[, usePC] * PC_sign)
    imagePC <- imagePC - min(imagePC)
    dim(imagePC) <- d[1:2]

    return(imagePC)
  }

  nuc <- if (length(nucleusIndex) > 1) {
    rowMeans(image[, , nucleusIndex, drop = FALSE], dims = 2)
  } else {
    image[, , nucleusIndex]
  }
  nuc <- .gblur(nuc, smooth)

  return(nuc - min(nuc))
}

## Slope of `input` on the distance to the nucleus edge, near the edge.
## Same QR fit as lm(), without building a model frame.
.estimateTolerance <- function(input, dist, discSize) {
  max_tresh <- max(3, 3 * discSize)
  near <- dist > 0 & dist < max_tresh
  y <- as.numeric(input[near])
  x <- as.numeric(dist[near])
  fit <- .lm.fit(cbind(1, x), y)
  fit$coefficients[2]
}

.Transform <- function(image, transform) {
  for (i in seq_along(transform)) {
    image <- switch(transform[i],
      "norm99" = .norm99(image),
      "asinh" = asinh(image),
      "maxThresh" = .maxThresh(image),
      "sqrt" = sqrt(image),
      "tissueMask" = image
    )
  }
  image
}

.norm99 <- function(nuc) {
  q <- quantile(nuc, 0.99, na.rm = TRUE)
  nuc[nuc > q] <- q
  return(nuc)
}

.maxThresh <- function(image) {
  m <- apply(image, 3, max)
  m[!(m > 0)] <- 1
  image / rep(m, each = prod(dim(image)[1:2]))
}

## Fitted values of a linear model of `y` on the asinh of every channel, the
## same numbers as predict(lm(y ~ ., data)) without the model frame copies.
.channelModel <- function(image, y) {
  d <- dim(image)
  X <- asinh(image)
  dim(X) <- c(d[1] * d[2], d[3])
  y <- as.numeric(y)
  if (anyNA(X) || anyNA(y)) {
    ## lm() drops incomplete pixels; keep its behaviour exactly.
    data <- data.frame(X, y = y)
    return(stats::predict(stats::lm(y ~ ., data = data), data))
  }
  X <- cbind(1, X)
  fit <- .lm.fit(X, y)
  r <- seq_len(fit$rank)
  piv <- fit$pivot[r]
  drop(X[, piv, drop = FALSE] %*% fit$coefficients[r])
}

## disc model ##
.CytSeg <- function(nmask,
                    image,
                    sizeSelection = 5,
                    discSize = 3) {
  kern <- EBImage::makeBrush(discSize, shape = "disc")

  cell <- EBImage::dilate(nmask, kern)

  disk <- cell - nmask > 0

  cytpred <- nmask
  cytpred[] <- .channelModel(image, disk)
  rm(image)
  cytpred <- cytpred - min(cytpred)
  cytpred <- cytpred / max(cytpred)

  cellTh <- EBImage::otsu(cytpred, range = c(0, 1))
  cell <- cytpred > cellTh

  cell <- cell + nmask > 0

  nmask[.smallObjects(EBImage::bwlabel(nmask), sizeSelection)] <- 0

  cmask4 <- EBImage::propagate(cytpred, nmask, cell)

  return(EBImage::Image(cmask4))
}

## Marker Model ## Cyt segmentation based on a specified cytoplasmic marker ##
## `image` is the raw image (the marker signal); `timage` the transformed one.
.CytSeg2 <- function(nmask,
                     image,
                     timage,
                     channel = 2,
                     sizeSelection = 5,
                     smooth = 1) {
  cytpred <- EBImage::Image(
    rowMeans(image[, , channel, drop = FALSE], dims = 2)
  )
  rm(image)

  cytpredsmooth <- .gblur(cytpred, smooth)

  ## using all the staining channels to predict the smoothed marker
  cytpredpred <- cytpred
  cytpredpred[] <- .channelModel(timage, cytpredsmooth)
  rm(timage)
  cytpredpred <- cytpredpred - min(cytpredpred)
  cytpredpred <- cytpredpred / max(cytpredpred)

  cellTh <- EBImage::otsu(cytpredpred, range = c(0, 1))
  cell <- cytpredpred > cellTh
  cell <- cell + nmask > 0

  nmask[.smallObjects(EBImage::bwlabel(nmask), sizeSelection)] <- 0

  cmask4 <- EBImage::propagate(cytpredpred, nmask, cell)

  return(EBImage::Image(cmask4))
}
