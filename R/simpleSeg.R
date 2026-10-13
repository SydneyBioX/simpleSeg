#' Perform simple segmentation of multiplexed cellular images
#'
#' @param image An image or list of images or CytoImageList to be read into the
#'              function.
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
#'
#' @return A list of image masks
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
#'   cores = 5
#' )
#'
#' @export simpleSeg
#' @rdname simpleSeg
#' @importFrom EBImage gblur otsu bwlabel makeBrush watershed dilate distmap propagate Image
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
                      cores = 1) {
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
  if (!is.null(tolerance)) .checkPositive(tolerance, "tolerance")
  if (!(isTRUE(pca) || isFALSE(pca))) stop("`pca` must be TRUE or FALSE.")

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

  ## Each image has its own random number stream (used by the PCA subsample).
  x <- runif(1) # nolint
  streams <- .rngStreams(length(image), .Random.seed[1])

  masks <- .parLapply(seq_along(image), function(i) {
    .withStream(streams[[i]], .segmentImage(image[[i]],
      name = names(image)[i],
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
  }, cores = cores)
  names(masks) <- names(image)

  cyto.mask <- cytomapper::CytoImageList(masks)

  if (is(image, "CytoImageList")) {
    mcols(cyto.mask) <- mcols(image)
  }

  S4Vectors::mcols(cyto.mask)$imageID <- names(image)

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
