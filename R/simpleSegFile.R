#' Segment a large image file tile by tile
#'
#' Segments a multiplexed image stored as a (pyramidal, tiled) TIFF file, such
#' as an OME-TIFF from CODEX or a PhenoCycler QPTIFF, without reading the whole
#' image into memory. The image is processed in tiles with an overlap; the
#' quantities that `simpleSeg()` computes over the whole image (Otsu threshold,
#' transformation constants, watershed tolerance) are accumulated over all tiles
#' first, so the segmentation follows `simpleSeg()`. Each cell is kept by the
#' tile containing its centroid, and the mean intensity of every channel in each
#' cell is measured from the raw image as the tiles are segmented.
#'
#' @param file Path to a TIFF file with one page per channel (OME-TIFF, QPTIFF
#'   or a multi-page TIFF).
#' @param nucleus The marker or markers corresponding to the nuclei.
#' @param cellBody Method of cytoplasm identification: 'dilate' or 'none'.
#' @param sizeSelection,smooth,watershed,tolerance,ext,discSize As in
#'   [simpleSeg()].
#' @param transform A transformation or list of transformations: "sqrt",
#'   "asinh", "norm99" and "maxThresh", as in [simpleSeg()]. "norm99" needs an
#'   image with integer pixel values.
#' @param channels Optional channel names, in page order. By default the names
#'   are read from the OME-XML or QPTIFF metadata.
#' @param tileSize Width and height of the tiles, in pixels.
#' @param overlap Pixels added on each side of a tile. It should be larger than
#'   the largest cell, so that every cell lies wholly within the tile that keeps
#'   it.
#' @param imageID Image name stored with the cells; defaults to the file name.
#' @param cores The number of cores or a BiocParallel param object; tiles are
#'   processed in parallel.
#'
#' @return A SpatialExperiment with one column per cell: the mean intensity of
#'   each channel (assay "intensities"), the cell area in pixels and the image
#'   name in the colData, and the centroid (pixels, x then y) as
#'   spatialCoords. The pixel size, if known, is stored in the metadata.
#'
#' @examples
#' \dontrun{
#' cells <- simpleSegFile("slide.ome.tif",
#'   nucleus = "DAPI",
#'   transform = "sqrt", cores = 8
#' )
#' }
#'
#' @export simpleSegFile
#' @importFrom graphics hist.default
simpleSegFile <- function(file,
                          nucleus,
                          cellBody = "dilate",
                          sizeSelection = 10,
                          smooth = 1,
                          transform = NULL,
                          watershed = "intensity",
                          tolerance = NULL,
                          ext = 1,
                          discSize = 3,
                          channels = NULL,
                          tileSize = 2048,
                          overlap = 64,
                          imageID = NULL,
                          cores = 1) {
  for (pkg in c("terra", "SpatialExperiment")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("simpleSegFile() needs the ", pkg, " package.")
    }
  }
  img <- .openImageFile(file, channels)
  if (is.null(imageID)) imageID <- sub("\\.(ome\\.)?(tiff?|qptiff)$", "", basename(file))

  nucleus <- .checkChannels(nucleus, img$channels, length(img$channels), "nucleus")
  if (is.numeric(nucleus)) nucleus <- img$channels[nucleus]
  .checkPositive(sizeSelection, "sizeSelection")
  .checkPositive(smooth, "smooth")
  .checkPositive(ext, "ext")
  .checkPositive(discSize, "discSize")
  .checkPositive(tileSize, "tileSize")
  if (!is.null(tolerance)) .checkPositive(tolerance, "tolerance")
  if (!(length(overlap) == 1 && is.numeric(overlap) && overlap >= 0)) {
    stop("`overlap` must be a number of pixels.")
  }
  valid_watershed <- c("intensity", "distance", "combine")
  if (!(length(watershed) == 1 && watershed %in% valid_watershed)) {
    stop("Invalid watershed. Choose one of: ", paste(valid_watershed, collapse = ", "), ".")
  }
  if (!(length(cellBody) == 1 && cellBody %in% c("dilate", "none"))) {
    stop("simpleSegFile() supports cellBody = \"dilate\" or \"none\".")
  }
  valid_transforms <- c("sqrt", "asinh", "norm99", "maxThresh")
  bad <- setdiff(transform, valid_transforms)
  if (length(bad) > 0) {
    stop(
      "Invalid transform for simpleSegFile(): '", paste(bad, collapse = "', '"),
      "'. Choose from ", paste(valid_transforms, collapse = ", "), "."
    )
  }
  .nWorkers(cores)

  tiles <- .tileGrid(img$nx, img$ny, as.integer(tileSize), as.integer(overlap))
  par <- list(
    img = img, nucleus = nucleus, smooth = smooth,
    sizeSelection = sizeSelection, watershed = watershed, ext = ext,
    discSize = discSize, wholeCell = cellBody == "dilate"
  )

  ## Constants of the transformations, from the whole image.
  par$plan <- .transformPlan(img, transform, nucleus, tiles, cores)

  ## The nucleus signal is shifted by its minimum over the whole image, and
  ## thresholded by Otsu's method over the whole image.
  r <- .parLapply(tiles, function(tile) {
    range(.tileSignal(par, tile)[tile$cx, tile$cy])
  }, cores = cores)
  par$shift <- min(vapply(r, `[`, numeric(1), 1))

  r <- .parLapply(tiles, function(tile) {
    v <- .tileSignal(par, tile)[tile$cx, tile$cy] - par$shift
    v <- v[v >= 0.0001]
    if (length(v)) range(v) else c(Inf, -Inf)
  }, cores = cores)
  otsuRange <- c(min(vapply(r, `[`, numeric(1), 1)), max(vapply(r, `[`, numeric(1), 2)))
  if (!all(is.finite(otsuRange))) {
    warning("No nucleus signal found in '", file, "'.", call. = FALSE)
    return(.cellsToSPE(NULL, img, imageID))
  }
  breaks <- seq(otsuRange[1], otsuRange[2], length.out = 257)
  counts <- .parLapply(tiles, function(tile) {
    v <- .tileSignal(par, tile)[tile$cx, tile$cy] - par$shift
    graphics::hist.default(v[v >= 0.0001], breaks = breaks, plot = FALSE)$counts
  }, cores = cores)
  par$nth <- .otsuCounts(Reduce(`+`, counts), breaks)

  ## Whole-image statistics of the thresholded nuclei: the median of the
  ## nucleus means (combine) and the watershed tolerance.
  if (watershed == "combine" || (is.null(tolerance) && watershed == "intensity")) {
    stats <- .parLapply(tiles, function(tile) .tileNucStats(par, tile), cores = cores)
    if (watershed == "combine") {
      ## background pixels all take the background mean over the whole image
      bgN <- sum(vapply(stats, `[[`, numeric(1), "bgN"))
      bgMean <- sum(vapply(stats, `[[`, numeric(1), "bgSum")) / bgN
      par$median <- .weightedMedian(
        c(unlist(lapply(stats, `[[`, "value")), if (bgN > 0) bgMean),
        c(unlist(lapply(stats, `[[`, "weight")), if (bgN > 0) bgN)
      )
    }
    if (is.null(tolerance)) {
      fit <- Reduce(.mergeMoments, lapply(stats, `[[`, "moments"))
      if (is.null(fit) || fit$n < 2) {
        warning("No nuclei found in '", file, "'.", call. = FALSE)
        return(.cellsToSPE(NULL, img, imageID))
      }
      tolerance <- fit$cxy / fit$cxx
      if (watershed == "combine") tolerance <- tolerance / par$median
    }
  }
  if (is.null(tolerance)) tolerance <- 1
  par$tolerance <- tolerance

  cells <- .parLapply(tiles, function(tile) .tileCells(par, tile), cores = cores)
  edge <- sum(vapply(cells, function(x) x$edge, numeric(1)))
  if (edge > 0) {
    warning(edge, " cells touch the edge of their tile's overlap; ",
      "consider a larger `overlap`.",
      call. = FALSE
    )
  }
  .cellsToSPE(cells, img, imageID)
}

## ---- reading -------------------------------------------------------------------------------

## Pages, channel names and size of a TIFF with one page per channel.
.openImageFile <- function(file, channels = NULL) {
  file <- normalizePath(file, mustWork = TRUE)
  info <- terra::describe(file)
  desc <- grep("^  TIFFTAG_IMAGEDESCRIPTION=", info, value = TRUE)
  desc <- if (length(desc)) sub("^  TIFFTAG_IMAGEDESCRIPTION=", "", desc[1]) else ""
  pageDesc <- grep("^  SUBDATASET_[0-9]+_DESC=", info, value = TRUE)
  size <- regmatches(pageDesc, regexec("\\(([0-9]+)P x ([0-9]+)L", pageDesc))
  if (length(size)) {
    nx <- as.integer(vapply(size, `[`, "", 2))
    ny <- as.integer(vapply(size, `[`, "", 3))
  } else {
    r <- .rast(file)
    nx <- terra::ncol(r)
    ny <- terra::nrow(r)
  }
  pixelSize <- NA_real_

  if (grepl("<OME", desc)) {
    names <- .xmlAttr(desc, "Channel", "Name")
    firstC <- as.integer(.xmlAttr(desc, "TiffData", "FirstC"))
    ifd <- as.integer(.xmlAttr(desc, "TiffData", "IFD"))
    pages <- if (length(ifd) == length(names)) ifd[order(firstC)] + 1L else seq_along(names)
    pixelSize <- suppressWarnings(as.numeric(.xmlAttr(desc, "Pixels", "PhysicalSizeX")[1]))
  } else {
    ## Full-resolution pages; QPTIFF names each channel in its page description.
    pages <- if (length(nx)) which(nx == nx[1] & ny == ny[1]) else 1L
    names <- vapply(pages, function(p) {
      d <- grep("TIFFTAG_IMAGEDESCRIPTION=", terra::describe(.pageName(file, p)), value = TRUE)
      m <- regmatches(d, regexpr("<(Biomarker|Name)>[^<]*</", d))
      if (length(m)) sub("^<(Biomarker|Name)>(.*)</$", "\\2", m[1]) else NA_character_
    }, character(1))
    if (anyNA(names)) names <- paste0("channel", seq_along(pages))
  }
  if (!is.null(channels)) {
    if (length(channels) != length(pages)) {
      stop("`channels` has ", length(channels), " names but the file has ", length(pages), " channels.")
    }
    names <- channels
  }
  list(
    file = file, nx = nx[1], ny = ny[1], channels = names, pages = pages,
    paged = length(size) > 0, pixelSize = pixelSize
  )
}

.xmlAttr <- function(xml, tag, attr) {
  m <- regmatches(xml, gregexpr(paste0("<", tag, " [^>]*>"), xml))[[1]]
  v <- regmatches(m, regexpr(paste0(" ", attr, "=\"[^\"]*\""), m))
  sub(paste0(" ", attr, "=\"([^\"]*)\""), "\\1", v)
}

.pageName <- function(file, page) sprintf("GTIFF_DIR:%d:%s", page, file)

.rast <- function(x) suppressWarnings(terra::rast(x))

## One channel of a window of the image, as an x by y matrix.
.readWindow <- function(img, channel, win) {
  k <- match(channel, img$channels)
  r <- if (img$paged) .rast(.pageName(img$file, img$pages[k])) else .rast(img$file)[[k]]
  terra::readStart(r)
  on.exit(terra::readStop(r))
  v <- terra::readValues(r,
    row = win$y0, nrows = win$h, col = win$x0, ncols = win$w,
    mat = FALSE
  )
  matrix(as.numeric(v), win$w, win$h)
}

## Tiles of the image: `core` pixels are owned by the tile, `win` is the window
## read (the core plus `overlap` on each side), `cx`/`cy` the core within it.
.tileGrid <- function(nx, ny, size, overlap) {
  g <- expand.grid(x0 = seq(1L, nx, by = size), y0 = seq(1L, ny, by = size))
  lapply(seq_len(nrow(g)), function(i) {
    x0 <- g$x0[i]
    y0 <- g$y0[i]
    x1 <- min(x0 + size - 1L, nx)
    y1 <- min(y0 + size - 1L, ny)
    wx0 <- max(1L, x0 - overlap)
    wy0 <- max(1L, y0 - overlap)
    wx1 <- min(nx, x1 + overlap)
    wy1 <- min(ny, y1 + overlap)
    list(
      win = list(x0 = wx0, y0 = wy0, w = wx1 - wx0 + 1L, h = wy1 - wy0 + 1L),
      cx = (x0 - wx0 + 1L):(x1 - wx0 + 1L),
      cy = (y0 - wy0 + 1L):(y1 - wy0 + 1L),
      ## window edges that are inside the image
      inner = c(wx0 > 1L, wx1 < nx, wy0 > 1L, wy1 < ny)
    )
  })
}

## ---- transformations -----------------------------------------------------------------------

## The constants each transformation needs from the whole image, found by
## running the transformations on each channel's histogram of pixel values.
.transformPlan <- function(img, transform, nucleus, tiles, cores) {
  if (!any(c("norm99", "maxThresh") %in% transform)) {
    return(lapply(transform, function(t) list(type = t)))
  }
  used <- if ("norm99" %in% transform) img$channels else nucleus
  hists <- .parLapply(tiles, function(tile) {
    core <- list(
      x0 = tile$win$x0 + tile$cx[1] - 1L, y0 = tile$win$y0 + tile$cy[1] - 1L,
      w = length(tile$cx), h = length(tile$cy)
    )
    lapply(used, function(ch) {
      v <- .readWindow(img, ch, core)
      if (any(v != round(v)) || min(v) < 0 || max(v) > 65535) {
        stop("transform \"norm99\" and \"maxThresh\" need non-negative integer pixel values.")
      }
      n <- tabulate(as.integer(v) + 1L, nbins = 65536L)
      list(at = which(n > 0), n = n[n > 0])
    })
  }, cores = cores)
  counts <- lapply(seq_along(used), function(j) {
    n <- numeric(65536L)
    for (h in hists) n[h[[j]]$at] <- n[h[[j]]$at] + h[[j]]$n
    n
  })
  names(counts) <- used
  levels <- lapply(counts, function(n) as.numeric(which(n > 0) - 1L))
  counts <- lapply(counts, function(n) n[n > 0])

  plan <- list()
  for (t in transform) {
    step <- list(type = t)
    if (t == "maxThresh") {
      step$max <- vapply(levels, max, numeric(1))
      step$max[!(step$max > 0)] <- 1
      levels <- Map(`/`, levels, step$max)
    } else if (t == "norm99") {
      step$q <- .weightedQuantile(unlist(levels), unlist(counts), 0.99)
      levels <- lapply(levels, pmin, step$q)
    } else {
      levels <- lapply(levels, match.fun(t))
    }
    plan[[length(plan) + 1]] <- step
  }
  plan
}

.applyPlan <- function(x, plan, channel) {
  for (step in plan) {
    x <- switch(step$type,
      "sqrt" = sqrt(x),
      "asinh" = asinh(x),
      "maxThresh" = x / step$max[[channel]],
      "norm99" = {
        x[x > step$q] <- step$q
        x
      }
    )
  }
  x
}

## quantile(type = 7) of `x` repeated `w` times.
.weightedQuantile <- function(x, w, p) {
  o <- order(x)
  x <- x[o]
  cw <- cumsum(as.numeric(w[o]))
  index <- 1 + (cw[length(cw)] - 1) * p
  lo <- floor(index)
  hi <- ceiling(index)
  xlo <- x[findInterval(lo - 1, cw) + 1L]
  xhi <- x[findInterval(hi - 1, cw) + 1L]
  h <- index - lo
  if (index > lo && xhi != xlo) (1 - h) * xlo + h * xhi else xlo
}

## median() of `x` repeated `w` times.
.weightedMedian <- function(x, w) {
  o <- order(x)
  x <- x[o]
  cw <- cumsum(as.numeric(w[o]))
  n <- cw[length(cw)]
  at <- function(k) x[findInterval(k - 1, cw) + 1L]
  half <- (n + 1) %/% 2
  if (n %% 2 == 1) at(half) else mean(c(at(half), at(half + 1)))
}

## EBImage::otsu() from histogram counts.
.otsuCounts <- function(counts, breaks) {
  counts <- as.double(counts)
  mids <- 0.5 * (breaks[-1L] + breaks[-length(breaks)])
  len <- length(counts)
  w1 <- cumsum(counts)
  w2 <- w1[len] + counts - w1
  cm <- counts * mids
  m1 <- cumsum(cm)
  m2 <- m1[len] + cm - m1
  var <- w1 * w2 * (m2 / w2 - m1 / w1)^2
  maxi <- which(var == max(var, na.rm = TRUE))
  (mids[maxi[1]] + mids[maxi[length(maxi)]]) / 2
}

## ---- segmentation of a tile ----------------------------------------------------------------

## Smoothed nucleus signal of a tile's window (before the whole-image shift).
## The window is read with a margin of the blur radius that wraps around the
## image edges, so the blur gives exactly the values of the whole image.
.tileSignal <- function(par, tile) {
  r <- ceiling(3 * par$smooth)
  xs <- (tile$win$x0 - r - 1L + seq_len(tile$win$w + 2L * r) - 1L) %% par$img$nx + 1L
  ys <- (tile$win$y0 - r - 1L + seq_len(tile$win$h + 2L * r) - 1L) %% par$img$ny + 1L
  sig <- array(0, c(length(xs), length(ys), length(par$nucleus)))
  for (k in seq_along(par$nucleus)) {
    ch <- par$nucleus[k]
    sig[, , k] <- .applyPlan(.readIndexed(par$img, ch, xs, ys), par$plan, ch)
  }
  sig <- if (length(par$nucleus) > 1) rowMeans(sig, dims = 2) else sig[, , 1]
  .gaussianBlurCpp(sig, par$smooth, FALSE)
}

## One channel at the given x and y pixel positions (runs of consecutive
## positions are read as blocks).
.readIndexed <- function(img, channel, xs, ys) {
  runs <- function(v) split(seq_along(v), cumsum(c(1L, diff(v) != 1L)))
  out <- matrix(0, length(xs), length(ys))
  for (rx in runs(xs)) {
    for (ry in runs(ys)) {
      out[rx, ry] <- .readWindow(img, channel, list(
        x0 = xs[rx[1]], y0 = ys[ry[1]], w = length(rx), h = length(ry)
      ))
    }
  }
  out
}

## Thresholded, size-selected nuclei of a tile's window.
.tileMask <- function(par, tile) {
  nuc <- .tileSignal(par, tile) - par$shift
  nuc[nuc < 0.0001] <- 0
  nMask <- nuc > par$nth
  nMaskLabel <- EBImage::bwlabel(nMask * 1)
  small <- .smallObjects(nMaskLabel, par$sizeSelection)
  nMask[small] <- 0
  nMaskLabel[small] <- 0
  list(nuc = nuc, nMask = nMask, label = nMaskLabel)
}

## Per-label centroids (window coordinates) and whether they fall in the core.
.ownedLabels <- function(label, tile) {
  lab <- matrix(as.integer(label), nrow(label))
  n <- max(lab, 0L)
  m <- .labelMomentsCpp(lab, n)
  present <- m[, 1] > 0
  cx <- m[, 2] / m[, 1]
  cy <- m[, 3] / m[, 1]
  own <- present & cx >= tile$cx[1] - 0.5 & cx < tail(tile$cx, 1) + 0.5 &
    cy >= tile$cy[1] - 0.5 & cy < tail(tile$cy, 1) + 0.5
  list(labels = which(own), cx = cx[own], cy = cy[own], n = m[own, 1], max = n)
}

## Nucleus means and sizes (for the combine median) and the moments of the
## tolerance regression, over the pixels a tile owns.
.tileNucStats <- function(par, tile) {
  m <- .tileMask(par, tile)
  core <- matrix(FALSE, nrow(m$nuc), ncol(m$nuc))
  core[tile$cx, tile$cy] <- TRUE
  out <- list(value = numeric(0), weight = numeric(0))

  dist <- EBImage::distmap(m$nMask)
  if (par$watershed == "combine") {
    own <- .ownedLabels(m$label, tile)
    AVG <- .labelMeans(m$nuc, m$label)
    lab <- as.integer(m$label)
    first <- match(own$labels, lab)
    bg <- core & m$label == 0
    out$value <- AVG[first]
    out$weight <- own$n
    out$bgSum <- sum(m$nuc[bg])
    out$bgN <- sum(bg)
    input <- dist * ((m$nuc / AVG) * m$nMask)
  } else {
    input <- m$nuc
  }
  maxTresh <- max(3, 3 * par$discSize)
  near <- core & dist > 0 & dist < maxTresh
  out$moments <- .moments(as.numeric(dist[near]), as.numeric(input[near]))
  out
}

## Count, means and centred cross-products of x and y, merged across tiles.
.moments <- function(x, y) {
  n <- length(x)
  if (n == 0) return(NULL)
  mx <- mean(x)
  my <- mean(y)
  list(n = n, mx = mx, my = my, cxx = sum((x - mx)^2), cxy = sum((x - mx) * (y - my)))
}

.mergeMoments <- function(a, b) {
  if (is.null(a)) return(b)
  if (is.null(b)) return(a)
  n <- a$n + b$n
  dx <- b$mx - a$mx
  dy <- b$my - a$my
  list(
    n = n, mx = a$mx + dx * b$n / n, my = a$my + dy * b$n / n,
    cxx = a$cxx + b$cxx + dx^2 * a$n * b$n / n,
    cxy = a$cxy + b$cxy + dx * dy * a$n * b$n / n
  )
}

## Segment a tile and measure the cells it owns.
.tileCells <- function(par, tile) {
  m <- .tileMask(par, tile)
  if (par$watershed == "distance") {
    wMask <- EBImage::watershed(EBImage::distmap(m$nMask), tolerance = par$tolerance, ext = par$ext)
  } else if (par$watershed == "combine") {
    AVG <- .labelMeans(m$nuc, m$label)
    nuc <- (m$nuc / AVG) * m$nMask
    nuc <- nuc / par$median * m$nMask
    dist <- EBImage::distmap(m$nMask)
    wMask <- EBImage::watershed(dist * nuc, tolerance = par$tolerance, ext = par$ext)
  } else {
    wMask <- EBImage::watershed(m$nuc * m$nMask, tolerance = par$tolerance, ext = par$ext)
    wMask[.smallObjects(wMask, par$sizeSelection)] <- 0
  }
  if (par$wholeCell) {
    wMask <- EBImage::dilate(wMask, EBImage::makeBrush(par$discSize, shape = "disc"))
  }

  own <- .ownedLabels(wMask, tile)
  lab <- as.integer(wMask)

  ## cells reaching a window edge inside the image may be cut off
  edgeLab <- c(
    if (tile$inner[1]) wMask[1, ], if (tile$inner[2]) wMask[nrow(wMask), ],
    if (tile$inner[3]) wMask[, 1], if (tile$inner[4]) wMask[, ncol(wMask)]
  )
  edge <- sum(own$labels %in% edgeLab)

  means <- matrix(0, length(own$labels), length(par$img$channels),
    dimnames = list(NULL, par$img$channels)
  )
  if (length(own$labels)) {
    for (ch in par$img$channels) {
      v <- .readWindow(par$img, ch, tile$win)
      means[, ch] <- .labelSumsCpp(lab, v, own$max)[own$labels]
    }
    means <- means / own$n
  }
  list(
    means = means,
    area = own$n,
    x = own$cx + tile$win$x0 - 1,
    y = own$cy + tile$win$y0 - 1,
    edge = edge
  )
}

## ---- output --------------------------------------------------------------------------------

#' @importFrom S4Vectors DataFrame
.cellsToSPE <- function(cells, img, imageID) {
  means <- do.call(rbind, lapply(cells, `[[`, "means"))
  if (is.null(means)) {
    means <- matrix(0, 0, length(img$channels), dimnames = list(NULL, img$channels))
  }
  n <- nrow(means)
  coords <- cbind(
    m.cx = unlist(lapply(cells, `[[`, "x")),
    m.cy = unlist(lapply(cells, `[[`, "y"))
  )
  if (is.null(coords) || ncol(coords) < 2) coords <- matrix(numeric(0), 0, 2, dimnames = list(NULL, c("m.cx", "m.cy")))
  colnames(means) <- img$channels
  rownames(means) <- NULL
  spe <- SpatialExperiment::SpatialExperiment(
    assays = list(intensities = t(means)),
    colData = DataFrame(
      imageID = rep(imageID, n),
      objectNum = seq_len(n),
      area = as.numeric(unlist(lapply(cells, `[[`, "area")))
    ),
    spatialCoords = coords
  )
  colnames(spe) <- paste(imageID, seq_len(n), sep = "_")
  S4Vectors::metadata(spe)$pixelSize <- img$pixelSize
  spe
}
