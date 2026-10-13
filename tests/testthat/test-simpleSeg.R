# helper function to load data for testing.
load_data <- function() {
  # Get path to image directory
  path_to_images <- system.file("extdata", package = "simpleSeg")

  # Get directories of images

  image_dirs <- dir(path_to_images, "Point", full.names = TRUE)
  names(image_dirs) <- dir(path_to_images, "Point", full.names = FALSE)

  # Get files in each directory
  files <- lapply(
    image_dirs,
    list.files,
    pattern = "tif",
    full.names = TRUE
  )

  # Read files with readImage from EBImage
  images <- lapply(files, EBImage::readImage, as.is = TRUE)

  # Convert to cytoImageList
  images <- cytomapper::CytoImageList(images)
  mcols(images)$imageID <- names(images)

  return(images)
}

# test that simpleSeg is running as expected
# test_that("Test masks are the same as the saved ones.", {
#   # load saved masks
#   saved_masks <- readRDS("saved_masks.rds")

#   # compute masks on example data
#   images <- load_data()

#   masks <- simpleSeg::simpleSeg(images,
#     nucleus = "HH3",
#     transform = "sqrt"
#   )

#   # compare the computed data in one of the masks to the saved correct results
#   expect_equal(
#     saved_masks$Point2203_pt1072_31606@.Data,
#     masks$Point2203_pt1072_31606@.Data
#   )
# })

# test transform cellBody parameter
test_that("Test if cellBody parameter is valid.", {
  # load images
  images <- load_data()

  expect_silent(simpleSeg(images,
    nucleus = "HH3",
    cellBody = "none"
  ))

  # TODO: Test expected error for incorrect input.
})

# test transform parameter
test_that("Test if transform parameter is valid.", {
  # load images
  images <- load_data()

  # code runs as expected for valid input
  expect_silent(simpleSeg(images,
    nucleus = "HH3",
    transform = "sqrt"
  ))
  expect_silent(simpleSeg(images,
    nucleus = "HH3",
    transform = c("sqrt", "norm99")
  ))
  expect_silent(simpleSeg(images,
    nucleus = "HH3",
    transform = c()
  ))

  # error on invalid input
  expect_error(simpleSeg(images,
    nucleus = "HH3",
    transform = c("sqrt", "bad")
  ))
  expect_error(simpleSeg(images,
    nucleus = "HH3",
    transform = "fake"
  ))
})

# small crops of the example images keep the tests fast
crop_data <- function() {
  images <- load_data()
  cropped <- cytomapper::CytoImageList(lapply(
    as.list(images),
    function(x) x[1:300, 1:300, ]
  ))
  mcols(cropped)$imageID <- names(cropped)
  cropped
}

test_that("cell body models return one mask per image", {
  images <- crop_data()
  for (cellBody in list("discModel", "PanKRT", c("PanKRT", "VIM"))) {
    masks <- simpleSeg(images,
      nucleus = "HH3", transform = "sqrt",
      cellBody = cellBody
    )
    expect_s4_class(masks, "CytoImageList")
    expect_equal(names(masks), names(images))
    expect_equal(dim(masks[[1]]), c(300, 300))
    expect_equal(mcols(masks)$imageID, names(images))
  }
})

test_that("a single image works with a cell body marker", {
  images <- crop_data()
  masks <- simpleSeg(images[[1]],
    nucleus = "HH3", transform = "sqrt",
    cellBody = "PanKRT"
  )
  expect_length(masks, 1)
})

test_that("nucleus can be given by position", {
  images <- crop_data()
  pos <- match("HH3", channelNames(images))
  expect_identical(
    as.list(simpleSeg(images, nucleus = pos)),
    as.list(simpleSeg(images, nucleus = "HH3"))
  )
})

test_that("invalid arguments give informative errors", {
  images <- crop_data()
  expect_error(simpleSeg(images, nucleus = "HH33"), "not found")
  expect_error(simpleSeg(images, nucleus = "HH3", tissue = "nope", transform = "tissueMask"), "not found")
  expect_error(simpleSeg(images, nucleus = "HH3", watershed = "intensty"), "Invalid watershed")
  expect_error(simpleSeg(images, nucleus = "HH3", cellBody = "dilat"), "Invalid cellBody")
  expect_error(simpleSeg(images, nucleus = "HH3", cores = "two"), "cores")
})

test_that("results do not depend on the number of cores", {
  skip_on_os("windows")
  images <- crop_data()
  serial <- simpleSeg(images,
    nucleus = c("HH3", "Nuc"), transform = "sqrt",
    pca = TRUE, cellBody = "discModel"
  )
  parallel <- simpleSeg(images,
    nucleus = c("HH3", "Nuc"), transform = "sqrt",
    pca = TRUE, cellBody = "discModel", cores = 2
  )
  expect_identical(as.list(serial), as.list(parallel))
})

test_that("an image without nuclei gives a warning and an empty mask", {
  images <- crop_data()[1]
  img <- as.array(images[[1]])
  img[1, 1, "HH3"] <- 1e6
  images[[1]] <- EBImage::Image(img)
  expect_warning(
    masks <- simpleSeg(images, nucleus = "HH3", cellBody = "discModel"),
    "No nuclei"
  )
  expect_equal(max(masks[[1]]), 0)
})

test_that("normalizeCells checks imageID and handles one marker", {
  data("pancreasSCE", package = "cytomapper")
  expect_error(
    normalizeCells(pancreasSCE,
      assayIn = "counts", imageID = "nope",
      method = "minMax"
    ),
    "imageID"
  )
  one <- normalizeCells(pancreasSCE,
    markers = "CD99", assayIn = "counts",
    imageID = "ImageNb", method = c("trim99", "minMax")
  )
  norm <- assay(one, "norm")["CD99", ]
  expect_true(all(tapply(norm, pancreasSCE$ImageNb, max) == 1))
  other <- rownames(one) != "CD99"
  expect_equal(assay(one, "norm")[other, ], assay(one, "counts")[other, ])
})

test_that("normalizeCells gives the same values for data frames", {
  data("pancreasSCE", package = "cytomapper")
  sce <- normalizeCells(pancreasSCE,
    assayIn = "counts", imageID = "ImageNb",
    transformation = "asinh", method = c("trim99", "mean")
  )
  df <- data.frame(t(assay(pancreasSCE, "counts")),
    ImageNb = pancreasSCE$ImageNb, check.names = FALSE
  )
  df <- normalizeCells(df,
    imageID = "ImageNb",
    transformation = "asinh", method = c("trim99", "mean")
  )
  expect_equal(
    unname(t(assay(sce, "norm"))),
    unname(as.matrix(df[, rownames(pancreasSCE)]))
  )
})
