## Cells of a label mask: area, centroid (pixels, from 1) and the sum of each
## channel, as the tiles report them.
.maskCells <- function(mask, image) {
  lab <- matrix(as.integer(mask), nrow(mask))
  n <- max(lab, 0L)
  m <- .labelMomentsCpp(lab, n, 1, 1)
  ids <- which(m[, 1] > 0)
  sums <- matrix(0, length(ids), dim(image)[3])
  for (k in seq_len(dim(image)[3])) {
    sums[, k] <- .labelSumsCpp(lab, as.numeric(image[, , k]), n)[ids]
  }
  list(
    sums = sums, area = m[ids, 1], x = m[ids, 2] / m[ids, 1], y = m[ids, 3] / m[ids, 1],
    object_id = ids
  )
}

## One table of cells from the tiles of an image.
.combineCells <- function(cells, img) {
  out <- list(
    sums = do.call(rbind, lapply(cells, `[[`, "sums")),
    area = unlist(lapply(cells, `[[`, "area")),
    x = unlist(lapply(cells, `[[`, "x")),
    y = unlist(lapply(cells, `[[`, "y")),
    object_id = unlist(lapply(cells, `[[`, "object_id"))
  )
  if (is.null(out$sums)) out <- .noCells(img)
  out
}

.noCells <- function(img) {
  list(
    sums = matrix(0, 0, length(img$channels)), area = numeric(0), x = numeric(0),
    y = numeric(0), object_id = integer(0)
  )
}

## A SpatialExperiment of the cells of all images: mean intensity of each
## channel (assay "counts"), image, object number (the cell's label in its
## mask) and area in the colData, and centroids as spatialCoords, named as in
## cytomapper::measureObjects().
#' @importFrom S4Vectors DataFrame
.cellsSPE <- function(cells, channels) {
  if (!requireNamespace("SpatialExperiment", quietly = TRUE)) {
    stop("output = \"cells\" needs the SpatialExperiment package.")
  }
  imageID <- rep(names(cells), vapply(cells, function(x) length(x$area), integer(1)))
  means <- do.call(rbind, lapply(cells, function(x) x$sums / x$area))
  if (is.null(means)) means <- matrix(0, 0, length(channels))
  objectID <- as.integer(unlist(lapply(cells, `[[`, "object_id")))
  coords <- cbind(
    m.cx = as.numeric(unlist(lapply(cells, `[[`, "x"))),
    m.cy = as.numeric(unlist(lapply(cells, `[[`, "y")))
  )
  counts <- t(means)
  dimnames(counts) <- list(channels, paste(imageID, objectID, sep = "_"))
  spe <- SpatialExperiment::SpatialExperiment(
    assays = list(counts = counts),
    colData = DataFrame(
      imageID = imageID,
      object_id = objectID,
      s.area = as.numeric(unlist(lapply(cells, `[[`, "area")))
    ),
    spatialCoords = coords
  )
  colnames(spe) <- colnames(counts)
  spe
}
