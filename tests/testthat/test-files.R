## An example image cropped to 300 x 250 pixels, written as a multi-band TIFF.
write_example <- function(dir = tempfile()) {
  dir.create(dir, showWarnings = FALSE)
  path <- system.file("extdata", package = "simpleSeg")
  folder <- dir(path, "Point", full.names = TRUE)[1]
  img <- EBImage::readImage(list.files(folder, "tif", full.names = TRUE), as.is = TRUE)
  img <- as.array(img)[201:500, 151:400, ]
  file <- file.path(dir, "example.tif")
  suppressWarnings(terra::writeRaster(terra::rast(aperm(img, c(2, 1, 3))), file, datatype = "FLT8S"))
  list(file = file, image = img, channels = dimnames(img)[[3]])
}

order_cells <- function(spe) {
  xy <- SpatialExperiment::spatialCoords(spe)
  spe[, order(spe$imageID, xy[, 1], xy[, 2])]
}

expect_same_cells <- function(a, b) {
  a <- order_cells(a)
  b <- order_cells(b)
  expect_equal(ncol(a), ncol(b))
  expect_equal(unname(SpatialExperiment::spatialCoords(a)), unname(SpatialExperiment::spatialCoords(b)))
  expect_equal(a$s.area, b$s.area)
  expect_equal(unname(SummarizedExperiment::assay(a)), unname(SummarizedExperiment::assay(b)))
}

## Two label masks with the same cells, whatever their numbers.
expect_same_masks <- function(a, b) {
  a <- as.integer(a)
  b <- as.integer(b)
  expect_identical(a > 0, b > 0)
  pairs <- unique(cbind(a, b)[a > 0, , drop = FALSE])
  expect_false(anyDuplicated(pairs[, 1]) > 0)
  expect_false(anyDuplicated(pairs[, 2]) > 0)
}

test_that("tiled segmentation gives the cells and masks of the whole image", {
  skip_if_not_installed("terra")
  skip_if_not_installed("SpatialExperiment")
  ex <- write_example()
  setups <- list(
    list(ws = "intensity", tr = "sqrt"),
    list(ws = "combine", tr = c("asinh", "maxThresh")),
    list(ws = "distance", tr = "sqrt")
  )
  for (setup in setups) {
    args <- list(nucleus = "channel25", transform = setup$tr, watershed = setup$ws)
    whole <- do.call(simpleSeg, c(list(ex$file), args))
    tiled <- do.call(simpleSeg, c(list(ex$file, tileSize = 110, overlap = 40), args))
    expect_same_cells(whole, tiled)
    wholeMask <- do.call(simpleSeg, c(list(ex$file, output = "masks"), args))
    tiledMask <- do.call(simpleSeg, c(list(ex$file, output = "masks", tileSize = 110, overlap = 40), args))
    expect_same_masks(wholeMask[[1]], tiledMask[[1]])
  }
})

test_that("cells from files match cells from images in R", {
  skip_if_not_installed("terra")
  skip_if_not_installed("SpatialExperiment")
  ex <- write_example()
  fromFile <- simpleSeg(ex$file, nucleus = 25, transform = "sqrt")
  image <- ex$image
  dimnames(image)[[3]] <- paste0("channel", seq_len(dim(image)[3]))
  inR <- simpleSeg(list(example = image), nucleus = 25, transform = "sqrt", output = "cells")
  expect_same_cells(fromFile, inR)
  ## the means are those of cytomapper::measureObjects()
  masks <- simpleSeg(list(example = image), nucleus = 25, transform = "sqrt")
  images <- cytomapper::CytoImageList(list(example = EBImage::Image(image)))
  S4Vectors::mcols(images)$imageID <- "example"
  sce <- cytomapper::measureObjects(masks, images, img_id = "imageID")
  o <- match(sce$object_id, inR$object_id)
  expect_equal(unname(SummarizedExperiment::assay(inR)[, o]), unname(SummarizedExperiment::assay(sce, "counts")))
  expect_equal(unname(SpatialExperiment::spatialCoords(inR)[o, ]), unname(as.matrix(SummarizedExperiment::colData(sce)[, c("m.cx", "m.cy")])))
})

test_that("masks are saved to maskDir", {
  skip_if_not_installed("terra")
  ex <- write_example()
  dir <- tempfile()
  files <- simpleSeg(ex$file, nucleus = 25, transform = "sqrt", output = "none", maskDir = dir, tileSize = 110)
  masks <- simpleSeg(ex$file, nucleus = 25, transform = "sqrt", output = "masks", tileSize = 110)
  saved <- terra::values(terra::rast(files[[1]]), mat = FALSE)
  expect_equal(as.numeric(as.array(masks[[1]])), as.numeric(saved))
})

test_that("folders of single-channel images are read as images split by marker", {
  skip_if_not_installed("terra")
  skip_if_not_installed("SpatialExperiment")
  path <- system.file("extdata", package = "simpleSeg")
  cells <- simpleSeg(path, nucleus = "HH3", transform = "sqrt")
  images <- load_data()
  inR <- simpleSeg(images, nucleus = "HH3", transform = "sqrt", output = "cells")
  expect_equal(sort(unique(cells$imageID)), sort(names(images)))
  expect_setequal(rownames(cells), rownames(inR))
  expect_same_cells(cells, inR[rownames(cells), ])
})

test_that("a folder of multi-channel images gives one image per file", {
  skip_if_not_installed("terra")
  skip_if_not_installed("SpatialExperiment")
  dir <- tempfile()
  ex <- write_example(dir)
  file.copy(ex$file, file.path(dir, "copy.tif"))
  expect_message(
    cells <- simpleSeg(dir, nucleus = 25, transform = "sqrt"),
    "separate image"
  )
  expect_setequal(unique(cells$imageID), c("example", "copy"))
})
