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
#' @importFrom BiocParallel SerialParam bplapply MulticoreParam bpmapply
#' @importFrom EBImage gblur otsu bwlabel makeBrush filter2 watershed dilate distmap propagate Image as.Image
#' @importFrom terra predict
#' @importFrom cytomapper CytoImageList
#' @importFrom stats prcomp quantile lm coef cor median resid runif sd
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

  #sizeSelection validation
  if (!(sizeSelection > 0)) {
    stop(
      paste0(
        sprintf(
          "Invalid sizeSelection: '%s'. sizeSelection must be positive",
          sizeSelection
        )
      )
    )
  }
  if (!(smooth > 0)) {
    stop(
      paste0(
        sprintf(
          "Invalid smooth: '%s'. smooth must be positive",
          smooth
        )
      )
    )
  }

  if (!(ext > 0)) {
    stop(
      paste0(
        sprintf(
          "Invalid ext: '%s'. ext must be positive",
          ext
        )
      )
    )
  }
  if (!(discSize > 0)) {
    stop(
      paste0(
        sprintf(
          "Invalid discSize: '%s'. discSize must be positive",
          discSize
        )
      )
    )
  }
  # cellBody input validation
  valid_cellbody <- c("none", "dilate", "discModel")
  if (!(cellBody %in% valid_cellbody)) {
    # Throw informative error message.
    stop(
      paste0(
        sprintf(
          "Invalid method of cytoplasm identification: '%s'. Must be ",
          cellBody
        ),
        paste(valid_cellbody, collapse = ", "),
        "."
      )
    )
  }

  # transform input validation
  valid_tansforms <- c("sqrt", "asinh", "norm99", "maxThresh", "tissueMask")
  if (is.null(transform)) {

  } else if (is.vector(transform)) {
    dummy <- TRUE
    for (element in transform) {
      dummy <- (element %in% valid_tansforms) && dummy
      if (!dummy) {
        # Throw informative error message.
        stop(
          paste0(
            sprintf(
              "Transform list contains invalid transform: '%s'. ",
              element
            ),
            paste(valid_tansforms, collapse = ", "),
            "."
          )
        )
      }
    }
  } else {
    if (!(transform %in% c())) {
      stop(
        paste0(
          sprintf("Invalid transform: '%s'. Choose from ", transform),
          paste(valid_tansforms, collapse = ", "),
          "."
        )
      )
    }
  }

  imageClass <- class(image)

  if (!imageClass %in% c("list", "CytoImageList")) {
    image <- list(image)
    names(image) <- "image"
  }

  wholeCell <- FALSE
  if (cellBody == "dilate") {
    wholeCell <- TRUE
  }

  x <- runif(1) # nolint

  if (!is(cores, 'BPParam')) {
    BPPARAM <- generateBPParam(cores)
  } else {
    BPPARAM <- cores
  }

  nmask <- .nucSegParallel(image,
    nucleusIndex = nucleus,
    sizeSelection = sizeSelection,
    smooth = smooth,
    watershed = watershed,
    tolerance = tolerance,
    ext = ext,
    wholeCell = wholeCell,
    discSize = discSize,
    transform = transform,
    tissueIndex = tissue,
    pca = pca,
    BPPARAM = BPPARAM
  )

  # if dilate or none
  if (cellBody %in% c("dilate", "none")) {
    nmask <- lapply(nmask, EBImage::Image)
    cyto.nmask <- cytomapper::CytoImageList(nmask)

    if (is.null(names(cyto.nmask))) {
      names(cyto.nmask) <- c(seq_along(cyto.nmask))
    }

    if (is(image, "CytoImageList")) {
      mcols(cyto.nmask) <- mcols(image)
    }

    S4Vectors::mcols(cyto.nmask)$imageID <- names(image)

    objectNum <- as.list(vapply(cyto.nmask, max, numeric(1)))

    mcols(cyto.nmask)$objectNum <- objectNum
    return(cyto.nmask)
  }

  if (cellBody == "discModel") {
    cells <- .cytSegParallel(nmask,
      image,
      sizeSelection = sizeSelection,
      smooth = smooth,
      discSize = discSize,
      transform = transform,
      BPPARAM = BPPARAM
    )

    # Converting from a tiff stack to individual images
    cellList <- NULL

    for (i in seq_along(cells[1, 1, ])) {
      cellList[[i]] <- as.Image(cells[, , i])
    }

    cyto.mask <- cytomapper::CytoImageList(cellList)

    if (is.null(names(image))) {
      names(image) <- c(seq_along(image))
    }

    if (is(image, "CytoImageList")) {
      mcols(cyto.mask) <- mcols(image)
    }

    S4Vectors::mcols(cyto.mask)$imageID <- names(image)

    objectNum <- as.list(vapply(cyto.mask, max, numeric(1)))

    mcols(cyto.mask)$objectNum <- objectNum
    return(cyto.mask)
  }

  if (any(cellBody %in% dimnames(image[[1]])[[3]])) {
    cells <- .cytSeg2Parallel(nmask,
      image,
      channel = cellBody,
      sizeSelection = sizeSelection,
      smooth = smooth,
      transform = transform,
      BPPARAM = BPPARAM
    )

    # Converting from a tiff stack to individual images
    cellList <- NULL
    for (i in seq_along(cells[1, 1, ])) {
      cellList[[i]] <- as.Image(cells[, , i])
    }

    cyto.mask <- cytomapper::CytoImageList(cellList)

    if (is.null(names(image))) {
      names(image) <- c(seq_along(image))
    }

    if (is(image, "CytoImageList")) {
      mcols(cyto.mask) <- mcols(image)
    }

    S4Vectors::mcols(cyto.mask)$imageID <- names(image)

    objectNum <- as.list(vapply(cyto.mask, max, numeric(1)))

    mcols(cyto.mask)$objectNum <- objectNum
    return(cyto.mask)
  }
}
