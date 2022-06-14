nucNormalize.helper <- function(image, nuc, smooth, normalize){
    if ("autoS" %in% normalize){
        smooth <- autosmooth(image, smooth, 9, 4) # adjusting the smoothing parameter for low intensity images
    }
    
    
    # Hotspot filtering. Intensities greater than the 99 percentile are changed to be exactly the 99th percentile intensity. This has the effect of removing outliers.
    
    if ("norm99" %in% normalize){
        nuc[nuc > quantile(nuc, 0.99,na.rm = TRUE)] <- quantile(nuc,
                                                                0.99,
                                                                na.rm = TRUE)
    }
    if ("maxThresh" %in% normalize){
        nuc <- nuc/max(nuc,na.rm = TRUE)
    }
    output <- list(smooth, nuc)
    return(output)
}