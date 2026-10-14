#' Perform simple segmentation of multiplexed cellular images
#'
#' @param image An image, a list of images or a CytoImageList; or paths to
#'              image files or folders (see Details).
#' @param nucleus The marker or list of markers corresponding to the nuclei.
#' @param cellBody Method of cytoplasm identification. Can be 'none', dilate',
#'                 'discModel' or the name of a dedicated cytoplasm marker
#' @param sizeSelection Minimum pixels for an object to be recognized as a cell
#'                      and not noise.
#' @param smooth The amount of Gaussian smoothing to be applied to the image/s
#' @param transform A transformation or list of transformations and
#'                  normalizations to be performed prior to nuclei or cytoplasm
#'                  identification. Accepted vales: "sqrt", "asinh", "norm99",
#'                  "maxThresh" and "tissueMask". Tissue mask may be used when
#'                  the sample does not take up the entirety of the image
#'                  (typically a circular sample inside the image. When tissue
#'                  mask is specified the background noise present outside the
#'                  sample area is removed).
#' @param watershed Method used to perform watersheding. Accepted values:
#'                  "intensity", "distance" or "combine".
#' @param tolerance The minimum height of the object in the units of image
#'                  intensity between its highest point (seed) and the point
#'                  where it contacts another object (checked for every contact
#'                  pixel). If the height is smaller than the tolerance, the
#'                  object will be combined with one of its neighbors, which is
#'                  the highest. Tolerance should be chosen according to the
#'                  range of x. Default value is 1, which is a reasonable value
#'                  if x comes from distmap.
#' @param ext Radius of the neighborhood in pixels for the detection of
#'            neighboring objects. Higher value smooths out small objects.
#' @param discSize The size of dilation around nuclei to create cell disc or
#'                 capture cytoplasm
#' @param tissue Channels to be used to create the tissue mask if specified
#'               in transforms.
#' @param pca Whether to run PCA on aggregated nucleus markers in order to
#'            detect the cellular nucclei.
#' @param cores The number or cores for parallel processing or a BPPARAM object
#' @param output What to return: "masks" (a CytoImageList of cell masks),
#'               "cells" (a SpatialExperiment of the cells, see Value) or
#'               "none" (only write the masks to `maskDir`). Defaults to
#'               "masks" for images and "cells" for paths.
#' @param maskDir Optional folder in which to save each image's mask as a TIFF
#'                file named after the image.
#' @param tileSize Images read from files that are wider or taller than
#'                 `tileSize` pixels are segmented in tiles of this size,
#'                 without reading the whole image into memory.
#' @param overlap Pixels added on each side of a tile. It should be larger
#'                than the largest cell.
#'
#' @details
#' `image` may be given as paths:
#' * a multi-channel image file (OME-TIFF, QPTIFF, or a TIFF with one page or
#'   band per channel) is one image; channel names are read from the OME-XML or
#'   QPTIFF metadata;
#' * a folder of single-channel TIFF files is one image split by marker, with
#'   the file names (without extension) as channel names; if the files have
#'   more than three channels each, each file is instead treated as an image;
#' * a folder of such folders holds one image per folder.
#'
#' Images that fit within `tileSize` are read into memory and segmented as
#' images in R are. Larger images are segmented tile by tile: the threshold,
#' transformation constants and watershed tolerance are first calculated over
#' the whole image, each tile is segmented with `overlap` pixels around it, and
#' each cell belongs to the tile containing its centroid. The cells are those
#' found when segmenting the whole image. Tiled images support `cellBody`
#' "dilate" and "none", and all transformations except "tissueMask";
#' "norm99" needs integer pixel values.
#'
#' Before the watershed, pixels of equal height are separated by a tiny amount
#' (a millionth of the tolerance) that depends on their position, so that ties
#' are resolved the same way however the image is processed.
#'
#' @return With output = "masks", a CytoImageList of masks. With output =
#'   "cells", a SpatialExperiment with one column per cell: the mean intensity
#'   of each channel (assay "counts"), the image (`imageID`), the cell's label
#'   in its mask (`object_id`) and its area in pixels (`s.area`) in the
#'   colData, and its centroid (`m.cx`, `m.cy`, in pixels) as spatialCoords.
#'   With output = "none", the paths of the mask files, invisibly.
#'
#' @examples
#'
#' library(cytomapper)
#' data("pancreasImages")
#' masks <- simpleSeg(pancreasImages,
#'   nucleus = "H3",
#'   cellBody = "discModel",
#'   sizeSelection = 8,
#'   smooth = 1.2,
#'   transform = "sqrt",
#'   watershed = "combine",
#'   tolerance = 1, ext = 1,
#'   discSize = 3,
#'   cores = 1
#' )
#'
#' # images stored as one folder of single-channel files per image
#' path <- system.file("extdata", package = "simpleSeg")
#' cells <- simpleSeg(path, nucleus = "HH3", transform = "sqrt")
#'
#' @export simpleSeg
#' @rdname simpleSeg
#' @useDynLib simpleSeg, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom EBImage otsu bwlabel makeBrush watershed dilate distmap propagate Image
#' @importFrom cytomapper CytoImageList
#' @importFrom stats prcomp quantile median sd .lm.fit
#' @importFrom S4Vectors mcols mcols<-
#' @importFrom methods is
simpleSeg <- function(image,
                      nucleus,
                      cellBody = "dilate",
                      sizeSelection = 10,
                      smooth = 1,
                      transform = NULL,
                      watershed = "intensity",
                      tolerance = NULL,
                      ext = 1,
                      discSize = 3,
                      tissue = NULL,
                      pca = FALSE,
                      cores = 1,
                      output = NULL,
                      maskDir = NULL,
                      tileSize = 4096,
                      overlap = 64) {
  fromFiles <- is.character(image)
  if (fromFiles) {
    sources <- .imageSources(image)
    channels <- sources[[1]]$channels
    for (src in sources) {
      if (!identical(src$channels, channels)) {
        stop(
          "All images must have the same channels; '", src$id,
          "' differs from '", sources[[1]]$id, "'."
        )
      }
    }
    nChannels <- length(channels)
    ids <- names(sources)
  } else {
    if (!is(image, "list") && !is(image, "CytoImageList")) {
      image <- list(image = image)
    }
    if (length(image) == 0) stop("`image` contains no images.")
    if (is.null(names(image))) {
      names(image) <- as.character(seq_along(image))
    }
    channels <- dimnames(image[[1]])[[3]]
    nChannels <- dim(image[[1]])[3]
    if (is.na(nChannels)) {
      stop("Each image must have a third dimension holding its channels.")
    }
    ids <- names(image)
  }

  if (is.null(output)) output <- if (fromFiles) "cells" else "masks"
  if (!(length(output) == 1 && output %in% c("masks", "cells", "none"))) {
    stop("`output` must be \"masks\", \"cells\" or \"none\".")
  }
  if (output == "none" && is.null(maskDir)) {
    stop("output = \"none\" needs `maskDir`, the folder to save the masks in.")
  }

  if ("PCA" %in% nucleus) {
    if (length(nucleus) == 1) {
      stop(
        "simpleSeg no longer supports running PCA on all available markers.\n",
        "Please input a list of nuclear markers."
      )
    } else {
      warning(
        "PCA is depreciated as a nuclear marker.\n",
        "Please set simpleSeg's pca parameter to TRUE instead"
      )
    }
    pca <- TRUE
    # remove PCA from nucleus markers
    nucleus <- nucleus[nucleus != "PCA"]
  }
  nucleus <- .checkChannels(nucleus, channels, nChannels, "nucleus")
  if (!is.null(tissue)) {
    tissue <- .checkChannels(tissue, channels, nChannels, "tissue")
  }

  .checkPositive(sizeSelection, "sizeSelection")
  .checkPositive(smooth, "smooth")
  .checkPositive(ext, "ext")
  .checkPositive(discSize, "discSize")
  .checkPositive(tileSize, "tileSize")
  if (!is.null(tolerance)) .checkPositive(tolerance, "tolerance")
  if (!(isTRUE(pca) || isFALSE(pca))) stop("`pca` must be TRUE or FALSE.")
  if (!(length(overlap) == 1 && is.numeric(overlap) && overlap >= 0)) {
    stop("`overlap` must be a number of pixels.")
  }

  valid_watershed <- c("intensity", "distance", "combine")
  if (!(length(watershed) == 1 && watershed %in% valid_watershed)) {
    stop(
      sprintf("Invalid watershed: '%s'. Choose one of: ", paste(watershed, collapse = ", ")),
      paste(valid_watershed, collapse = ", "), "."
    )
  }

  # cellBody input validation
  useCellMarkers <- FALSE
  valid_cellbody <- c("none", "dilate", "discModel")
  if (length(cellBody) == 1 && cellBody %in% valid_cellbody) {
    # an approximation method
  } else if (length(cellBody) > 0 && all(cellBody %in% channels)) {
    useCellMarkers <- TRUE
  } else {
    missing <- setdiff(cellBody, channels)
    stop(
      "Invalid cellBody: '", paste(missing, collapse = "', '"), "'. ",
      "Choose one of the methods ", paste(valid_cellbody, collapse = ", "),
      ", or one or more cell body markers from: ",
      paste(channels, collapse = ", "), "."
    )
  }

  # transform input validation
  valid_transforms <- c("sqrt", "asinh", "norm99", "maxThresh", "tissueMask")
  bad <- setdiff(transform, valid_transforms)
  if (length(bad) > 0) {
    stop(
      sprintf("Invalid transform: '%s'. ", paste(bad, collapse = "', '")),
      "Choose from ", paste(valid_transforms, collapse = ", "), "."
    )
  }

  .nWorkers(cores)
  if (!is.null(maskDir)) dir.create(maskDir, showWarnings = FALSE, recursive = TRUE)
  maskFile <- function(id) if (!is.null(maskDir)) file.path(maskDir, paste0(id, ".tif"))

  ## Each image has its own random number stream (used by the PCA subsample).
  x <- runif(1) # nolint
  streams <- .rngStreams(length(ids), .Random.seed[1])

  ## Segment one image held in memory, then measure and/or save its mask.
  segment <- function(i, img) {
    mask <- .withStream(streams[[i]], .segmentImage(img,
      name = ids[i],
      nucleusIndex = nucleus,
      cellBody = cellBody,
      useCellMarkers = useCellMarkers,
      sizeSelection = sizeSelection,
      smooth = smooth,
      tolerance = tolerance,
      watershed = watershed,
      ext = ext,
      discSize = discSize,
      transform = transform,
      tissueIndex = tissue,
      pca = pca
    ))
    if (!is.null(maskDir)) .writeMask(mask, maskFile(ids[i]))
    list(
      mask = if (output == "masks") mask,
      cells = if (output == "cells") .maskCells(mask, img)
    )
  }

  if (!fromFiles) {
    res <- .parLapply(seq_along(image), function(i) segment(i, image[[i]]), cores = cores)
  } else {
    big <- vapply(sources, function(s) s$nx > tileSize || s$ny > tileSize, logical(1))
    if (any(big) && (!cellBody %in% c("dilate", "none") || pca || "tissueMask" %in% transform)) {
      stop(
        "Images larger than `tileSize` (", paste(ids[big], collapse = ", "),
        ") are segmented in tiles, which supports cellBody \"dilate\" or \"none\", ",
        "without pca or \"tissueMask\". Increase `tileSize` to read them whole."
      )
    }
    res <- vector("list", length(sources))
    small <- which(!big)
    res[small] <- .parLapply(small, function(i) segment(i, .readImage(sources[[i]])), cores = cores)
    par <- list(
      nucleus = nucleus, smooth = smooth, sizeSelection = sizeSelection,
      watershed = watershed, tolerance = tolerance, ext = ext, discSize = discSize,
      wholeCell = cellBody == "dilate", transform = transform
    )
    for (i in which(big)) {
      r <- .segmentTiled(sources[[i]], par, tileSize, overlap, cores,
        maskFile = maskFile(ids[i]), keepMask = output == "masks"
      )
      res[[i]] <- list(mask = r$mask, cells = r$cells)
    }
  }
  names(res) <- ids

  if (output == "none") return(invisible(vapply(ids, maskFile, "")))
  if (output == "cells") return(.cellsSPE(lapply(res, `[[`, "cells"), channels))

  masks <- lapply(res, function(r) EBImage::Image(r$mask))
  cyto.mask <- cytomapper::CytoImageList(masks)
  if (is(image, "CytoImageList")) {
    mcols(cyto.mask) <- mcols(image)
  }
  S4Vectors::mcols(cyto.mask)$imageID <- ids
  objectNum <- as.list(vapply(masks, max, numeric(1)))
  mcols(cyto.mask)$objectNum <- objectNum
  cyto.mask
}

## Channel names (or positions) that must exist in the images.
.checkChannels <- function(x, channels, nChannels, arg) {
  if (is.numeric(x)) {
    if (length(x) == 0 || any(is.na(x)) || any(x < 1 | x > nChannels | x != round(x))) {
      stop("`", arg, "` channel positions must be between 1 and ", nChannels, ".")
    }
    return(if (is.null(channels)) x else channels[x])
  }
  missing <- setdiff(x, channels)
  if (length(x) == 0 || length(missing) > 0) {
    stop(
      "`", arg, "` marker", if (length(missing) > 1) "s", " '",
      paste(missing, collapse = "', '"), "' not found in the images. ",
      "Channels are: ", paste(channels, collapse = ", "), "."
    )
  }
  x
}

.checkPositive <- function(x, arg) {
  if (!(is.numeric(x) && length(x) == 1 && !is.na(x) && x > 0)) {
    stop(sprintf("Invalid %s: '%s'. %s must be a positive number.", arg, paste(x, collapse = ", "), arg))
  }
}
