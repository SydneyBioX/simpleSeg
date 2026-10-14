## Segmentation of a large image file tile by tile.
##
## Quantities that simpleSeg() computes over the whole image are found first:
## exactly where they come from histograms or sums over pixels (transformation
## constants, tissue mask and nucleus thresholds, combine median, watershed
## tolerance), and from a random sample of pixels across the tiles for the PCA
## and the cell body models. Each tile is then segmented with an overlap
## around it, and keeps the cells whose centroid lies in its core. The blur
## reads a margin that wraps around the image edges and the watershed breaks
## ties by position in the whole image, so without PCA or a cell body model
## the cells are those of simpleSeg() on the whole image.

#' @importFrom graphics hist.default
.segmentTiled <- function(img, par, tileSize, overlap, cores, maskFile = NULL,
                          keepMask = FALSE, salt = 0) {
  tiles <- .tileGrid(img$nx, img$ny, as.integer(tileSize), as.integer(overlap))
  par$img <- img
  par$salt <- salt
  par$cache <- tempfile("simpleSegTiles")
  dir.create(par$cache)
  on.exit(unlink(par$cache, recursive = TRUE), add = TRUE)
  empty <- list(cells = .noCells(img), mask = if (keepMask) matrix(0L, img$nx, img$ny))
  noNuclei <- function() {
    warning("No nuclei found in image '", img$id, "'.", call. = FALSE)
    empty
  }

  ## tissue mask and transformation constants: the nucleus signal uses the
  ## image within the tissue, the cell body models the whole image
  if ("tissueMask" %in% par$transform) par$tissue <- .tissueSetup(par, tiles, cores)
  allChannels <- par$pca || par$cellModel != "none" || "norm99" %in% par$transform
  used <- if (allChannels) img$channels else par$nucleus
  par$plan <- .transformPlan(par, used, tiles, cores, masked = !is.null(par$tissue))
  par$cellPlan <- if (par$cellModel != "none" && !is.null(par$tissue)) {
    .transformPlan(par, img$channels, tiles, cores, masked = FALSE)
  } else {
    par$plan
  }

  if (par$pca) {
    par <- .pcaSetup(par, tiles, cores)
    if (is.null(par)) return(noNuclei())
  }

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
  if (!all(is.finite(otsuRange))) return(noNuclei())
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
      if (is.null(fit) || fit$n < 2) return(noNuclei())
      tolerance <- fit$cxy / fit$cxx
      if (par$watershed == "combine") tolerance <- tolerance / par$median
    }
  }
  if (is.null(tolerance)) tolerance <- 1
  par$tolerance <- tolerance

  if (par$cellModel != "none") par <- .cellModelSetup(par, tiles, cores)

  ## segment and measure; keep each tile's labels for the mask
  saveLabels <- keepMask || !is.null(maskFile)
  cells <- .parLapply(seq_along(tiles), function(i) {
    .tileCells(par, tiles[[i]], i, saveLabels)
  }, cores = cores)
  edge <- sum(vapply(cells, `[[`, numeric(1), "edge"))
  if (edge > 0) {
    warning(edge, " cells in image '", img$id, "' touch the edge of their tile's overlap; ",
      "consider a larger `overlap`.",
      call. = FALSE
    )
  }
  n <- vapply(cells, function(x) length(x$area), integer(1))
  offset <- c(0L, cumsum(n))
  for (i in seq_along(cells)) cells[[i]]$object_id <- offset[i] + seq_len(n[i])
  out <- list(cells = .combineCells(cells, img))
  if (saveLabels) {
    out$mask <- .stitchMask(img, tiles, par$cache, offset, cores, maskFile, keepMask)
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
## `masked`: the pixels outside the tissue are 0.
.transformPlan <- function(par, used, tiles, cores, masked = FALSE) {
  transform <- setdiff(par$transform, "tissueMask")
  if (!any(c("norm99", "maxThresh") %in% transform)) {
    return(lapply(transform, function(t) list(type = t)))
  }
  img <- par$img
  full <- "norm99" %in% transform
  hists <- .parLapply(tiles, function(tile) {
    xs <- tile$core[1]:tile$core[2]
    ys <- tile$core[3]:tile$core[4]
    inside <- if (masked) .tileTissue(par, xs, ys)
    lapply(used, function(ch) {
      v <- .readIndexed(img, ch, xs, ys)
      if (masked) v <- v * inside
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

## A tile window's x and y positions, with `margin` pixels around it that
## wrap around the image edges.
.windowIndex <- function(par, tile, margin = 0L) {
  list(
    xs = (tile$win$x0 - margin - 1L + seq_len(tile$win$w + 2L * margin) - 1L) %% par$img$nx + 1L,
    ys = (tile$win$y0 - margin - 1L + seq_len(tile$win$h + 2L * margin) - 1L) %% par$img$ny + 1L
  )
}

## One transformed channel at positions xs, ys; within the tissue if `masked`.
.tileChannel <- function(par, ch, xs, ys, plan, masked = FALSE, inside = NULL) {
  v <- .readIndexed(par$img, ch, xs, ys)
  if (masked && !is.null(par$tissue)) {
    if (is.null(inside)) inside <- .tileTissue(par, xs, ys)
    v <- v * inside
  }
  .applyPlan(v, plan, ch)
}

## Smoothed nucleus signal of a tile's window (before the whole-image shift):
## the mean of the nucleus channels, blurred with a margin of the blur radius
## that wraps around the image edges; or the cached PCA projection.
.tileSignal <- function(par, tile) {
  if (par$pca) return(readRDS(.cacheFile(par, "pca", tile)))
  r <- ceiling(3 * par$smooth)
  w <- .windowIndex(par, tile, r)
  inside <- if (!is.null(par$tissue)) .tileTissue(par, w$xs, w$ys)
  sig <- array(0, c(length(w$xs), length(w$ys), length(par$nucleus)))
  for (k in seq_along(par$nucleus)) {
    sig[, , k] <- .tileChannel(par, par$nucleus[k], w$xs, w$ys, par$plan, TRUE, inside)
  }
  sig <- if (length(par$nucleus) > 1) rowMeans(sig, dims = 2) else sig[, , 1]
  .gaussianBlurCpp(sig, par$smooth, FALSE)
}

.cacheFile <- function(par, what, tile) {
  file.path(par$cache, sprintf("%s_%d_%d.rds", what, tile$core[1], tile$core[3]))
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

## Nuclei of a tile's window, labelled by the watershed (cached when a cell
## body model needs them twice).
.tileNuclei <- function(par, tile) {
  f <- .cacheFile(par, "nuclei", tile)
  if (file.exists(f)) return(readRDS(f))
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
  wMask
}

## Cell bodies grown from the nuclei of a tile with the cell body model.
.tileCellBody <- function(par, tile, nmask) {
  w <- .windowIndex(par, tile)
  model <- par$model
  pred <- matrix(model$beta[1], length(w$xs), length(w$ys))
  for (k in seq_along(par$img$channels)) {
    b <- model$beta[k + 1]
    if (b != 0) {
      pred <- pred + b * asinh(.tileChannel(par, par$img$channels[k], w$xs, w$ys, par$cellPlan))
    }
  }
  pred <- (pred - model$lo) / (model$hi - model$lo)
  cell <- pred > model$th
  cell <- cell + nmask > 0
  nmask[.smallObjects(EBImage::bwlabel(nmask), par$sizeSelection)] <- 0
  EBImage::propagate(pred, nmask, cell)
}

## Segment a tile and measure the cells it owns; optionally save the labels of
## its cells for the mask.
.tileCells <- function(par, tile, i, saveLabels = FALSE) {
  nmask <- .tileNuclei(par, tile)
  wMask <- if (par$cellModel != "none") {
    .tileCellBody(par, tile, nmask)
  } else if (par$wholeCell) {
    EBImage::dilate(nmask, EBImage::makeBrush(par$discSize, shape = "disc"))
  } else {
    nmask
  }
  tab <- .labelTable(wMask, tile)
  own <- tab$own

  ## cells reaching a window edge inside the image may be cut off
  edgeLab <- c(
    if (tile$inner[1]) tab$lab[1, ], if (tile$inner[2]) tab$lab[nrow(tab$lab), ],
    if (tile$inner[3]) tab$lab[, 1], if (tile$inner[4]) tab$lab[, ncol(tab$lab)]
  )

  if (saveLabels) {
    ## the window's pixels of the cells this tile owns, numbered 1, 2, ...
    lookup <- integer(nrow(tab$m))
    lookup[own] <- seq_along(own)
    lab <- tab$lab
    lab[lab > 0] <- lookup[lab[lab > 0]]
    saveRDS(lab, file.path(par$cache, paste0(i, ".rds")), compress = FALSE)
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

## ---- whole-image models from samples of pixels ---------------------------------------------

## Tissue mask: the log mean of the tissue channels, scaled to [0, 1] over the
## whole image and thresholded by Otsu's method over the whole image.
.tissueSetup <- function(par, tiles, cores) {
  channels <- if (is.null(par$tissueChannels)) par$img$channels else par$tissueChannels
  logMean <- function(tile) {
    xs <- tile$core[1]:tile$core[2]
    ys <- tile$core[3]:tile$core[4]
    m <- 0
    for (ch in channels) m <- m + .readIndexed(par$img, ch, xs, ys)
    log10(m / length(channels) + 0.01)
  }
  r <- .parLapply(tiles, function(tile) range(logMean(tile)), cores = cores)
  lo <- min(vapply(r, `[`, numeric(1), 1))
  hi <- max(vapply(r, `[`, numeric(1), 2))
  if (!(hi > lo)) hi <- lo + 1
  breaks <- seq(0, 1, length.out = 257)
  counts <- .parLapply(tiles, function(tile) {
    graphics::hist.default((logMean(tile) - lo) / (hi - lo), breaks = breaks, plot = FALSE)$counts
  }, cores = cores)
  list(channels = channels, lo = lo, hi = hi, th = .otsuCounts(Reduce(`+`, counts), breaks))
}

## The tissue mask at positions xs, ys, without tissue regions of
## `sizeSelection` pixels or fewer.
.tileTissue <- function(par, xs, ys) {
  t <- par$tissue
  m <- 0
  for (ch in t$channels) m <- m + .readIndexed(par$img, ch, xs, ys)
  inside <- (log10(m / length(t$channels) + 0.01) - t$lo) / (t$hi - t$lo) > t$th
  inside[.smallObjects(EBImage::bwlabel(inside * 1), par$sizeSelection)] <- FALSE
  inside
}

## PCA of the blurred channels: the nucleus signal is the principal component
## most correlated with the nucleus markers, fitted on up to 100,000 pixels
## sampled from around the nuclei across the whole image. The projection of
## each tile's window is cached.
.pcaSetup <- function(par, tiles, cores) {
  r <- ceiling(3 * par$smooth)
  blurred <- function(ch, w, inside) {
    .gaussianBlurCpp(.tileChannel(par, ch, w$xs, w$ys, par$plan, TRUE, inside), par$smooth, FALSE)
  }
  window <- function(tile) {
    w <- .windowIndex(par, tile, r)
    w$inside <- if (!is.null(par$tissue)) .tileTissue(par, w$xs, w$ys)
    w
  }
  nucMean <- function(w) {
    Reduce(`+`, lapply(par$nucleus, blurred, w = w, inside = w$inside)) / length(par$nucleus)
  }
  kern <- EBImage::makeBrush(par$discSize, shape = "disc")
  useMask <- function(nm) EBImage::dilate(nm > par$pcaTh, kern) > 0

  rr <- .parLapply(tiles, function(tile) range(nucMean(window(tile))[tile$cx, tile$cy]), cores = cores)
  breaks <- seq(min(vapply(rr, `[`, numeric(1), 1)), max(vapply(rr, `[`, numeric(1), 2)), length.out = 257)
  counts <- .parLapply(tiles, function(tile) {
    graphics::hist.default(nucMean(window(tile))[tile$cx, tile$cy], breaks = breaks, plot = FALSE)$counts
  }, cores = cores)
  par$pcaTh <- .otsuCounts(Reduce(`+`, counts), breaks)

  n <- .parLapply(tiles, function(tile) sum(useMask(nucMean(window(tile)))[tile$cx, tile$cy]), cores = cores)
  n <- sum(unlist(n))
  if (n == 0) return(NULL)
  p <- min(1, 1e5 / n)
  samples <- .parLapply(tiles, function(tile) {
    w <- window(tile)
    nm <- nucMean(w)
    idx <- .samplePixelsCpp(
      useMask(nm), tile$win$x0, tile$win$y0, tile$cx[1], tile$cx[length(tile$cx)],
      tile$cy[1], tile$cy[length(tile$cy)], p, par$salt
    )
    X <- vapply(par$img$channels, function(ch) blurred(ch, w, w$inside)[idx], numeric(length(idx)))
    list(X = matrix(X, length(idx), length(par$img$channels)), nuc = nm[idx])
  }, cores = cores)
  X <- do.call(rbind, lapply(samples, `[[`, "X"))
  nuc <- unlist(lapply(samples, `[[`, "nuc"))

  useMarker <- apply(X, 2, sd) > 0
  pca <- prcomp(X[, useMarker], scale = TRUE)
  usePC <- which.max(abs(apply(pca$x, 2, cor, nuc)))
  rot <- numeric(ncol(X))
  rot[useMarker] <- pca$rotation[, usePC] * sign(cor(pca$x[, usePC], nuc))

  .parLapply(tiles, function(tile) {
    w <- window(tile)
    proj <- 0
    for (k in which(rot != 0)) proj <- proj + rot[k] * blurred(par$img$channels[k], w, w$inside)
    saveRDS(proj, .cacheFile(par, "pca", tile), compress = FALSE)
    NULL
  }, cores = cores)
  par
}

## Cell body model: a linear model of the disc around each nucleus (discModel)
## or of the smoothed cell body markers on the asinh of every channel, fitted on
## about 1,000,000 pixels sampled across the whole image. Its predictions are
## scaled to [0, 1] and thresholded by Otsu's method over the sample.
.cellModelSetup <- function(par, tiles, cores) {
  p <- min(1, 1e6 / (as.numeric(par$img$nx) * par$img$ny))
  r <- ceiling(3 * par$smooth)
  samples <- .parLapply(tiles, function(tile) {
    nmask <- .tileNuclei(par, tile)
    saveRDS(nmask, .cacheFile(par, "nuclei", tile), compress = FALSE)
    if (par$cellModel == "discModel") {
      cell <- EBImage::dilate(nmask, EBImage::makeBrush(par$discSize, shape = "disc"))
      y <- (cell - nmask > 0) * 1
    } else {
      wm <- .windowIndex(par, tile, r)
      raw <- 0
      for (ch in par$markers) raw <- raw + .readIndexed(par$img, ch, wm$xs, wm$ys)
      y <- .gaussianBlurCpp(raw / length(par$markers), par$smooth, FALSE)
    }
    all <- matrix(TRUE, tile$win$w, tile$win$h)
    idx <- .samplePixelsCpp(
      all, tile$win$x0, tile$win$y0, tile$cx[1], tile$cx[length(tile$cx)],
      tile$cy[1], tile$cy[length(tile$cy)], p, par$salt + 1
    )
    w <- .windowIndex(par, tile)
    X <- vapply(par$img$channels, function(ch) {
      asinh(.tileChannel(par, ch, w$xs, w$ys, par$cellPlan))[idx]
    }, numeric(length(idx)))
    list(X = matrix(X, length(idx), length(par$img$channels)), y = y[idx])
  }, cores = cores)
  X <- cbind(1, do.call(rbind, lapply(samples, `[[`, "X")))
  y <- unlist(lapply(samples, `[[`, "y"))

  fit <- .lm.fit(X, y)
  r <- seq_len(fit$rank)
  beta <- numeric(ncol(X))
  beta[fit$pivot[r]] <- fit$coefficients[r]
  pred <- drop(X %*% beta)
  lo <- min(pred)
  hi <- max(pred)
  if (!(hi > lo)) hi <- lo + 1
  breaks <- seq(0, 1, length.out = 257)
  h <- graphics::hist.default((pred - lo) / (hi - lo), breaks = breaks, plot = FALSE)$counts
  par$model <- list(beta = beta, lo = lo, hi = hi, th = .otsuCounts(h, breaks))
  par
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
