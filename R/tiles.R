## Segmentation of a large image file tile by tile.
##
## The quantities simpleSeg() computes over the whole image (the transformation
## constants, the shift and Otsu threshold of the nucleus signal, the combine
## median and the watershed tolerance) are first accumulated over all tiles.
## Each tile is then segmented with an overlap around it, and keeps the cells
## whose centroid lies in its core. The blur reads a margin that wraps around
## the image edges and the watershed breaks ties by position in the whole
## image, so the cells are those of simpleSeg() on the whole image.

#' @importFrom graphics hist.default
.segmentTiled <- function(img, par, tileSize, overlap, cores, maskFile = NULL, keepMask = FALSE) {
  tiles <- .tileGrid(img$nx, img$ny, as.integer(tileSize), as.integer(overlap))
  par$img <- img
  empty <- list(cells = .noCells(img), mask = if (keepMask) matrix(0L, img$nx, img$ny))

  par$plan <- .transformPlan(img, par$transform, par$nucleus, tiles, cores)

  ## shift and Otsu threshold of the nucleus signal
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
    warning("No nuclei found in image '", img$id, "'.", call. = FALSE)
    return(empty)
  }
  breaks <- seq(otsuRange[1], otsuRange[2], length.out = 257)
  counts <- .parLapply(tiles, function(tile) {
    v <- .tileSignal(par, tile)[tile$cx, tile$cy] - par$shift
    graphics::hist.default(v[v >= 0.0001], breaks = breaks, plot = FALSE)$counts
  }, cores = cores)
  par$nth <- .otsuCounts(Reduce(`+`, counts), breaks)

  ## combine median and watershed tolerance
  tolerance <- par$tolerance
  if (par$watershed == "combine" || (is.null(tolerance) && par$watershed == "intensity")) {
    stats <- .parLapply(tiles, function(tile) .tileNucStats(par, tile), cores = cores)
    if (par$watershed == "combine") {
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
        warning("No nuclei found in image '", img$id, "'.", call. = FALSE)
        return(empty)
      }
      tolerance <- fit$cxy / fit$cxx
      if (par$watershed == "combine") tolerance <- tolerance / par$median
    }
  }
  if (is.null(tolerance)) tolerance <- 1
  par$tolerance <- tolerance

  ## segment and measure; keep each tile's labels for the mask
  saveLabels <- keepMask || !is.null(maskFile)
  labelDir <- if (saveLabels) tempfile("simpleSegTiles")
  if (saveLabels) {
    dir.create(labelDir)
    on.exit(unlink(labelDir, recursive = TRUE), add = TRUE)
  }
  cells <- .parLapply(seq_along(tiles), function(i) {
    .tileCells(par, tiles[[i]], if (saveLabels) file.path(labelDir, paste0(i, ".rds")))
  }, cores = cores)
  edge <- sum(vapply(cells, `[[`, numeric(1), "edge"))
  if (edge > 0) {
    warning(edge, " cells in image '", img$id, "' touch the edge of their tile's overlap; ",
      "consider a larger `overlap`.",
      call. = FALSE
    )
  }
  n <- vapply(cells, function(x) length(x$area), integer(1))
  ids <- split(seq_len(sum(n)), rep(seq_along(cells), n))
  for (i in seq_along(cells)) cells[[i]]$object_id <- if (n[i]) ids[[as.character(i)]] else integer(0)
  out <- list(cells = .combineCells(cells, img))

  if (saveLabels) {
    out$mask <- .stitchMask(img, tiles, labelDir, c(0L, cumsum(n)), cores, maskFile, keepMask)
  }
  out
}

## Tiles of the image: `win` is the window read (the core plus `overlap` on each
## side), `cx`/`cy` the core within it.
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
      core = c(x0, x1, y0, y1),
      win = list(x0 = wx0, y0 = wy0, w = wx1 - wx0 + 1L, h = wy1 - wy0 + 1L),
      cx = (x0 - wx0 + 1L):(x1 - wx0 + 1L),
      cy = (y0 - wy0 + 1L):(y1 - wy0 + 1L),
      ## window edges that are inside the image
      inner = c(wx0 > 1L, wx1 < nx, wy0 > 1L, wy1 < ny)
    )
  })
}

## ---- transformations ---------------------------------------------------------------------------

## The constants each transformation needs from the whole image, found by
## running the transformations on each channel's histogram of pixel values.
## norm99 needs the histograms (integer pixel values); maxThresh only each
## channel's maximum, as the transformations keep the order of values.
.transformPlan <- function(img, transform, nucleus, tiles, cores) {
  if (!any(c("norm99", "maxThresh") %in% transform)) {
    return(lapply(transform, function(t) list(type = t)))
  }
  full <- "norm99" %in% transform
  used <- if (full) img$channels else nucleus
  hists <- .parLapply(tiles, function(tile) {
    core <- list(
      x0 = tile$core[1], y0 = tile$core[3],
      w = length(tile$cx), h = length(tile$cy)
    )
    lapply(used, function(ch) {
      v <- .readWindow(img, ch, core)
      if (!full) return(list(at = max(v), n = 1))
      if (any(v != round(v)) || min(v) < 0 || max(v) > 65535) {
        stop("transform \"norm99\" on images larger than `tileSize` needs non-negative integer pixel values.")
      }
      n <- tabulate(as.integer(v) + 1L, nbins = 65536L)
      list(at = which(n > 0), n = n[n > 0])
    })
  }, cores = cores)
  if (full) {
    counts <- lapply(seq_along(used), function(j) {
      n <- numeric(65536L)
      for (h in hists) n[h[[j]]$at] <- n[h[[j]]$at] + h[[j]]$n
      n
    })
    levels <- lapply(counts, function(n) as.numeric(which(n > 0) - 1L))
    counts <- lapply(counts, function(n) n[n > 0])
  } else {
    levels <- lapply(seq_along(used), function(j) max(vapply(hists, function(h) h[[j]]$at, numeric(1))))
    counts <- lapply(levels, function(l) 1)
  }
  names(levels) <- names(counts) <- used

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

## ---- segmentation of a tile --------------------------------------------------------------------

## Smoothed nucleus signal of a tile's window (before the whole-image shift),
## read with a margin of the blur radius that wraps around the image edges.
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

## One channel at the given x and y pixel positions, reading runs of
## consecutive positions as blocks.
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

## Pixel counts and coordinate sums (whole-image positions) of every label in
## a tile's window, and which labels the tile owns (centroid in its core).
.labelTable <- function(label, tile) {
  lab <- matrix(as.integer(label), nrow(label))
  n <- max(lab, 0L)
  m <- .labelMomentsCpp(lab, n, tile$win$x0, tile$win$y0)
  cx <- m[, 2] / m[, 1]
  cy <- m[, 3] / m[, 1]
  own <- m[, 1] > 0 & cx >= tile$core[1] - 0.5 & cx < tile$core[2] + 0.5 &
    cy >= tile$core[3] - 0.5 & cy < tile$core[4] + 0.5
  list(lab = lab, m = m, own = which(own), cx = cx, cy = cy)
}

## Nucleus means and sizes (for the combine median) and the moments of the
## tolerance regression, over the pixels a tile owns.
.tileNucStats <- function(par, tile) {
  m <- .tileMask(par, tile)
  core <- matrix(FALSE, nrow(m$nuc), ncol(m$nuc))
  core[tile$cx, tile$cy] <- TRUE
  out <- list()

  dist <- EBImage::distmap(m$nMask)
  if (par$watershed == "combine") {
    tab <- .labelTable(m$label, tile)
    AVG <- .labelMeans(m$nuc, m$label)
    first <- match(tab$own, tab$lab)
    bg <- core & m$label == 0
    out$value <- AVG[first]
    out$weight <- tab$m[tab$own, 1]
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

## Segment a tile and measure the cells it owns; optionally save the labels of
## its core for the mask.
.tileCells <- function(par, tile, labelFile = NULL) {
  m <- .tileMask(par, tile)
  x0 <- tile$win$x0
  y0 <- tile$win$y0
  if (par$watershed == "distance") {
    wMask <- .watershed(EBImage::distmap(m$nMask), par$tolerance, par$ext, x0, y0)
  } else if (par$watershed == "combine") {
    AVG <- .labelMeans(m$nuc, m$label)
    nuc <- (m$nuc / AVG) * m$nMask
    nuc <- nuc / par$median * m$nMask
    wMask <- .watershed(EBImage::distmap(m$nMask) * nuc, par$tolerance, par$ext, x0, y0)
  } else {
    wMask <- .watershed(m$nuc * m$nMask, par$tolerance, par$ext, x0, y0)
    wMask[.smallObjects(wMask, par$sizeSelection)] <- 0
  }
  if (par$wholeCell) {
    wMask <- EBImage::dilate(wMask, EBImage::makeBrush(par$discSize, shape = "disc"))
  }
  tab <- .labelTable(wMask, tile)
  own <- tab$own

  ## cells reaching a window edge inside the image may be cut off
  edgeLab <- c(
    if (tile$inner[1]) tab$lab[1, ], if (tile$inner[2]) tab$lab[nrow(tab$lab), ],
    if (tile$inner[3]) tab$lab[, 1], if (tile$inner[4]) tab$lab[, ncol(tab$lab)]
  )

  if (!is.null(labelFile)) {
    ## the window's pixels of the cells this tile owns, numbered 1, 2, ...
    lookup <- integer(nrow(tab$m))
    lookup[own] <- seq_along(own)
    lab <- tab$lab
    lab[lab > 0] <- lookup[lab[lab > 0]]
    saveRDS(lab, labelFile, compress = FALSE)
  }

  sums <- matrix(0, length(own), length(par$img$channels))
  if (length(own)) {
    for (k in seq_along(par$img$channels)) {
      v <- .readWindow(par$img, par$img$channels[k], tile$win)
      sums[, k] <- .labelSumsCpp(tab$lab, v, nrow(tab$m))[own]
    }
  }
  list(
    sums = sums, area = tab$m[own, 1], x = tab$cx[own], y = tab$cy[own],
    edge = sum(own %in% edgeLab)
  )
}

## The whole-image label mask from the tiles' saved labels: every tile writes
## the pixels of the cells it owns (cell k of tile i is cell offset[i] + k of
## the image), including pixels beyond its core. Rows of tiles are assembled in
## turn, and written to `maskFile` and/or kept.
.stitchMask <- function(img, tiles, labelDir, offset, cores, maskFile, keepMask) {
  if (!is.null(maskFile)) {
    r <- .rast(terra::rast(nrows = img$ny, ncols = img$nx, xmin = 0, xmax = img$nx, ymin = 0, ymax = img$ny))
    type <- if (offset[length(offset)] < 65536) "INT2U" else "INT4U"
    suppressWarnings(terra::writeStart(r, maskFile,
      overwrite = TRUE,
      wopt = list(datatype = type, gdal = c("COMPRESS=DEFLATE", "TILED=YES"))
    ))
    on.exit(terra::writeStop(r), add = TRUE)
  }
  mask <- if (keepMask) matrix(0L, img$nx, img$ny)
  rows <- vapply(tiles, function(t) t$core[3], numeric(1))
  for (y0 in sort(unique(rows))) {
    t1 <- tiles[[which(rows == y0)[1]]]
    rowsOut <- t1$core[3]:t1$core[4]
    ## tiles whose window reaches these rows
    from <- which(vapply(tiles, function(t) {
      t$win$y0 <= max(rowsOut) && t$win$y0 + t$win$h - 1L >= min(rowsOut)
    }, logical(1)))
    parts <- .parLapply(from, function(i) {
      t <- tiles[[i]]
      lab <- readRDS(file.path(labelDir, paste0(i, ".rds")))
      ys <- t$win$y0 - 1L + seq_len(t$win$h)
      keep <- ys %in% rowsOut
      lab <- lab[, keep, drop = FALSE]
      pix <- which(lab > 0)
      xy <- arrayInd(pix, dim(lab))
      list(
        index = (t$win$x0 - 1L + xy[, 1]) + (ys[keep][xy[, 2]] - min(rowsOut)) * img$nx,
        id = lab[pix] + offset[i]
      )
    }, cores = cores)
    full <- integer(img$nx * length(rowsOut))
    for (p in parts) full[p$index] <- pmax(full[p$index], p$id)
    if (!is.null(maskFile)) terra::writeValues(r, full, y0, length(rowsOut))
    if (keepMask) mask[, rowsOut] <- full
  }
  mask
}
