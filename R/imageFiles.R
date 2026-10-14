## Reading images from files, a window at a time, with terra (GDAL).
##
## An image is described by its size, channel names and, for each channel, the
## file and page (or band) holding it. Multi-channel files (OME-TIFF, QPTIFF,
## multi-page or multi-band TIFF) hold every channel; an image split by marker
## has one file per channel.

.imageExt <- "\\.(ome\\.)?(tiff?|qptiff)$"

## Images found at `paths`: files, folders of single-channel files (one image
## split by marker) or folders of such folders.
.imageSources <- function(paths) {
  if (!requireNamespace("terra", quietly = TRUE)) {
    stop("Reading images from files needs the terra package.")
  }
  out <- list()
  for (p in paths) {
    if (dir.exists(p)) {
      files <- .imageFiles(p)
      if (length(files)) {
        out <- c(out, .folderImages(p, files))
      } else {
        subdirs <- list.dirs(p, recursive = FALSE)
        subdirs <- subdirs[vapply(subdirs, function(d) length(.imageFiles(d)) > 0, logical(1))]
        if (!length(subdirs)) stop("No TIFF images found in '", p, "'.")
        for (d in subdirs) out <- c(out, .folderImages(d, .imageFiles(d)))
      }
    } else if (file.exists(p)) {
      out <- c(out, list(.openImageFile(p)))
    } else {
      stop("'", p, "' does not exist.")
    }
  }
  ids <- make.unique(vapply(out, `[[`, "", "id"))
  for (i in seq_along(out)) out[[i]]$id <- ids[i]
  names(out) <- ids
  out
}

.imageFiles <- function(dir) {
  sort(list.files(dir, .imageExt, full.names = TRUE, ignore.case = TRUE))
}

.stripExt <- function(file) sub(.imageExt, "", basename(file), ignore.case = TRUE)

## The files of a folder as one image split by marker, or, if they hold more
## than three channels each, as one image per file.
.folderImages <- function(dir, files) {
  nChannels <- vapply(files, function(f) length(.openImageFile(f)$channels), integer(1))
  if (any(nChannels > 3)) {
    message(
      "The images in '", dir, "' have more than three channels, so each file ",
      "is treated as a separate image."
    )
    return(lapply(files, .openImageFile))
  }
  if (any(nChannels > 1)) {
    message("Some images in '", dir, "' have ", max(nChannels), " channels; the first is used.")
  }
  info <- lapply(files, .openImageFile)
  nx <- vapply(info, `[[`, numeric(1), "nx")
  ny <- vapply(info, `[[`, numeric(1), "ny")
  if (any(nx != nx[1] | ny != ny[1])) {
    stop("The channel images in '", dir, "' have different sizes.")
  }
  list(list(
    id = basename(normalizePath(dir)), nx = nx[1], ny = ny[1],
    channels = .stripExt(files),
    src = do.call(rbind, lapply(info, function(x) x$src[1, ])),
    pixelSize = info[[1]]$pixelSize
  ))
}

## Size, channels and pages of one image file.
.openImageFile <- function(file) {
  file <- normalizePath(file, mustWork = TRUE)
  info <- terra::describe(file)
  desc <- grep("^  TIFFTAG_IMAGEDESCRIPTION=", info, value = TRUE)
  desc <- if (length(desc)) sub("^  TIFFTAG_IMAGEDESCRIPTION=", "", desc[1]) else ""
  pageDesc <- grep("^  SUBDATASET_[0-9]+_DESC=", info, value = TRUE)
  size <- regmatches(pageDesc, regexec("\\(([0-9]+)P x ([0-9]+)L", pageDesc))
  pixelSize <- NA_real_

  if (length(size)) {
    ## one channel per page; lower-resolution pages of a pyramid are skipped
    nx <- as.integer(vapply(size, `[`, "", 2))
    ny <- as.integer(vapply(size, `[`, "", 3))
    if (grepl("<OME", desc)) {
      names <- .xmlAttr(desc, "Channel", "Name")
      firstC <- as.integer(.xmlAttr(desc, "TiffData", "FirstC"))
      ifd <- as.integer(.xmlAttr(desc, "TiffData", "IFD"))
      pages <- if (length(ifd) == length(names)) ifd[order(firstC)] + 1L else seq_along(names)
      pixelSize <- suppressWarnings(as.numeric(.xmlAttr(desc, "Pixels", "PhysicalSizeX")[1]))
    } else {
      pages <- which(nx == nx[1] & ny == ny[1])
      ## QPTIFF names each channel in its page description
      names <- vapply(pages, function(p) {
        d <- grep("TIFFTAG_IMAGEDESCRIPTION=", terra::describe(.pageName(file, p)), value = TRUE)
        m <- regmatches(d, regexpr("<(Biomarker|Name)>[^<]*</", d))
        if (length(m)) sub("^<(Biomarker|Name)>(.*)</$", "\\2", m[1]) else NA_character_
      }, character(1))
    }
    src <- data.frame(file = file, page = pages, band = 1L)
    nx <- nx[pages[1]]
    ny <- ny[pages[1]]
  } else {
    ## one page with one band per channel
    r <- .rast(file)
    nx <- terra::ncol(r)
    ny <- terra::nrow(r)
    ## band names, unless they are terra's defaults
    names <- names(r)
    stem <- .stripExt(file)
    default <- grepl("^lyr\\.[0-9]+$", names) | names == stem |
      names == paste0(stem, "_", seq_along(names))
    if (any(default) || anyDuplicated(names)) names <- rep(NA_character_, terra::nlyr(r))
    src <- data.frame(file = file, page = 0L, band = seq_len(terra::nlyr(r)))
  }
  if (anyNA(names) || length(names) != nrow(src)) {
    names <- if (nrow(src) == 1) .stripExt(file) else paste0("channel", seq_len(nrow(src)))
  }
  list(
    id = .stripExt(file), nx = as.integer(nx), ny = as.integer(ny), channels = names,
    src = src, pixelSize = pixelSize
  )
}

.xmlAttr <- function(xml, tag, attr) {
  m <- regmatches(xml, gregexpr(paste0("<", tag, " [^>]*>"), xml))[[1]]
  v <- regmatches(m, regexpr(paste0(" ", attr, "=\"[^\"]*\""), m))
  sub(paste0(" ", attr, "=\"([^\"]*)\""), "\\1", v)
}

.pageName <- function(file, page) sprintf("GTIFF_DIR:%d:%s", page, file)

.rast <- function(x) suppressWarnings(terra::rast(x))

## Rasters opened for reading, kept open for the rest of a simpleSeg() call
## (opening a file costs more than reading a tile). Forked workers open their
## own: the process id is part of the key.
.rasters <- new.env(parent = emptyenv())

.channelRaster <- function(img, channel) {
  s <- img$src[match(channel, img$channels), ]
  key <- paste(Sys.getpid(), s$file, s$page, s$band)
  r <- .rasters[[key]]
  if (is.null(r)) {
    r <- if (s$page > 0) .rast(.pageName(s$file, s$page)) else .rast(s$file)[[s$band]]
    terra::readStart(r)
    assign(key, r, envir = .rasters)
  }
  r
}

## Close the rasters this process opened.
.closeRasters <- function() {
  mine <- grep(paste0("^", Sys.getpid(), " "), ls(.rasters), value = TRUE)
  for (key in mine) try(terra::readStop(.rasters[[key]]), silent = TRUE)
  rm(list = mine, envir = .rasters)
}

## One channel of a window of the image, as an x by y matrix.
.readWindow <- function(img, channel, win) {
  v <- terra::readValues(.channelRaster(img, channel),
    row = win$y0, nrows = win$h, col = win$x0, ncols = win$w,
    mat = FALSE
  )
  matrix(as.numeric(v), win$w, win$h)
}

## The whole image as an x by y by channel array.
.readImage <- function(img) {
  win <- list(x0 = 1L, y0 = 1L, w = img$nx, h = img$ny)
  out <- array(0, c(img$nx, img$ny, length(img$channels)),
    dimnames = list(NULL, NULL, img$channels)
  )
  for (k in seq_along(img$channels)) out[, , k] <- .readWindow(img, img$channels[k], win)
  out
}

## A label mask written as a compressed TIFF (16-bit when the labels fit).
.writeMask <- function(mask, file) {
  r <- .rast(t(matrix(as.integer(mask), nrow(mask))))
  type <- if (max(mask, 0) < 65536) "INT2U" else "INT4U"
  suppressWarnings(terra::writeRaster(r, file,
    overwrite = TRUE, datatype = type,
    gdal = c("COMPRESS=DEFLATE", "TILED=YES")
  ))
  file
}
