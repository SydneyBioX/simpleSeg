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
