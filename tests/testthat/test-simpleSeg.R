
load_data <- function() {
    # Get path to image directory
    pathToImages <- system.file("extdata", package = "simpleSeg")

    # Get directories of images

    imageDirs <- dir(pathToImages, "Point", full.names = TRUE)
    names(imageDirs) <- dir(pathToImages, "Point", full.names = FALSE)

    # Get files in each directory
    files <- sapply(imageDirs, list.files, pattern = "tif", full.names = TRUE, simplify = FALSE)

    # Read files with readImage from EBImage
    images <- lapply(files, EBImage::readImage, as.is = TRUE)

    # Convert to cytoImageList
    images <- cytomapper::CytoImageList(images)
    mcols(images)$imageID <- names(images)

    return(images)
}

test_that("Test masks are the sames as the saved ones.", {
    # load saved masks
    saved_masks <- readRDS("saved_masks.rds")

    # compute masks on example data
    images <- load_data()

    masks <- simpleSeg::simpleSeg(images,
                                 nucleus = "HH3",
                                 transform = "sqrt")

    expect_equal(saved_masks, masks)
})

test_that("Test if cellBody parameter is valid.", {
    # load images
    images <- load_data()

    expect_silent(simpleSeg(images, 
                            nucleus = "HH3",
                            cellBody='none'))

    expect_error(simpleSeg(images, 
                           nucleus = "HH3",
                           cellBody='large'))
})

test_that("Test if transform parameter is valid.", {
    # load images
    images <- load_data()

    # code runs as expected for valid input
    expect_silent(simpleSeg(images, 
                            nucleus = "HH3",
                            transform = 'sqrt'))
    expect_silent(simpleSeg(images, 
                            nucleus = "HH3",
                            transform = c('sqrt','norm99')))
    expect_silent(simpleSeg(images, 
                            nucleus = "HH3",
                            transform = c()))

    # error on invalid input
    expect_error(simpleSeg(images, 
                           nucleus = "HH3",
                           transform = c('sqrt','bad')))
    expect_error(simpleSeg(images, 
                           nucleus = "HH3",
                           transform = 'fake'))

})