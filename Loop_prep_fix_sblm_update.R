library(raster)     
library(terra)
library(zoo)         
library(signal)
library(entropy)

setwd("C:/Users/rakhm/OneDrive/Documents/Mapping_Jatim/New_Model_UGM_Comm/Data")

rasterOptions(progress = "Text", tmpdir = paste0(getwd(), "/tmp"))

# data chirps dan terrain
ch <- stack('chirps_idn_monthly_202311_202410.tif')
ter <- raster('ElevationExport.tif')

# list data raster sen2-sen1 yang akan diolah
setwd("C:/Users/rakhm/OneDrive/Documents/Mapping_Jatim/New_Model_UGM_Comm/Data/Loop_1/")
raster_files <- list.files(getwd(), ".tif", full.names = T)
raster_files <- raster_files[1:15]

# Define function to smooth and interpolate time series
smooth_and_interpolate <- function(ts) {
  ts_filled <- na.approx(ts, na.rm = FALSE)  
  ts_smoothed <- sgolayfilt(ts_filled, p = 3, n = 5) 
  return(ts_smoothed)
}

compute_stats <- function(x) {
  c(mean = mean(x, na.rm = TRUE),
    median = median(x, na.rm = TRUE),
    sd = sd(x, na.rm = TRUE),
    cov = ifelse(mean(x, na.rm = TRUE) != 0, sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE), NA),
    min = min(x, na.rm = TRUE),
    max = max(x, na.rm = TRUE),
    range = diff(range(x, na.rm = TRUE)),
    iqr = IQR(x, na.rm = TRUE),
    skewness = e1071::skewness(x, na.rm = TRUE),
    kurtosis = e1071::kurtosis(x, na.rm = TRUE))
}

calculate_entropy <- function(x) {
  entropy(x, method = "ML")
}

# reorder
for (file_index in seq_along(raster_files)) {
  file <- raster_files[file_index]
  cat(sprintf("\nProcessing file %d of %d: %s\n", file_index, length(raster_files), basename(file)))
  
  # Record start time
  start_time <- Sys.time()
  cat(sprintf("Start time: %s\n", start_time))
  
  # Load the raster stack
  cat("- Loading raster stack...\n")
  raster_stack <- stack(file)
  
  # Get the number of bands per month and months
  num_bands <- 10     
  num_months <- 12
  
  # Rearrange the raster stack: group by band first (B2, B2, ..., B3, B3, ..., etc.)
  cat("- Rearranging bands by time-series order...\n")
  order1 <- stack(raster_stack[[1:10]])
  order2 <- stack(raster_stack[[11:20]])
  order11 <- stack(raster_stack[[21:30]])
  order12 <- stack(raster_stack[[31:40]])
  order3 <- stack(raster_stack[[41:50]])
  order4 <- stack(raster_stack[[51:60]])
  order5 <- stack(raster_stack[[61:70]])
  order6 <- stack(raster_stack[[71:80]])
  order7 <- stack(raster_stack[[81:90]])
  order8 <- stack(raster_stack[[91:100]])
  order9 <- stack(raster_stack[[101:110]])
  order10 <- stack(raster_stack[[111:120]])
  
  B2 <- stack(order1[[1]], order2[[1]], order3[[1]], order4[[1]], order5[[1]], order6[[1]], order7[[1]], order8[[1]], order9[[1]], order10[[1]], order11[[1]], order12[[1]])
  B3 <- stack(order1[[2]], order2[[2]], order3[[2]], order4[[2]], order5[[2]], order6[[2]], order7[[2]], order8[[2]], order9[[2]], order10[[2]], order11[[2]], order12[[2]])
  B4 <- stack(order1[[3]], order2[[3]], order3[[3]], order4[[3]], order5[[3]], order6[[3]], order7[[3]], order8[[3]], order9[[3]], order10[[3]], order11[[3]], order12[[3]])
  B5 <- stack(order1[[4]], order2[[4]], order3[[4]], order4[[4]], order5[[4]], order6[[4]], order7[[4]], order8[[4]], order9[[4]], order10[[4]], order11[[4]], order12[[4]])
  B6 <- stack(order1[[5]], order2[[5]], order3[[5]], order4[[5]], order5[[5]], order6[[5]], order7[[5]], order8[[5]], order9[[5]], order10[[5]], order11[[5]], order12[[5]])
  B7 <- stack(order1[[6]], order2[[6]], order3[[6]], order4[[6]], order5[[6]], order6[[6]], order7[[6]], order8[[6]], order9[[6]], order10[[6]], order11[[6]], order12[[6]])
  B8 <- stack(order1[[7]], order2[[7]], order3[[7]], order4[[7]], order5[[7]], order6[[7]], order7[[7]], order8[[7]], order9[[7]], order10[[7]], order11[[7]], order12[[7]])
  B8A <- stack(order1[[8]], order2[[8]], order3[[8]], order4[[8]], order5[[8]], order6[[8]], order7[[8]], order8[[8]], order9[[8]], order10[[8]], order11[[8]], order12[[8]])
  VV <- stack(order1[[9]], order2[[9]], order3[[9]], order4[[9]], order5[[9]], order6[[9]], order7[[9]], order8[[9]], order9[[9]], order10[[9]], order11[[9]], order12[[9]])
  VH <- stack(order1[[10]], order2[[10]], order3[[10]], order4[[10]], order5[[10]], order6[[10]], order7[[10]], order8[[10]], order9[[10]], order10[[10]], order11[[10]], order12[[10]])
  
  cat("- Index Calculation...\n")
  
  ndvi_stack <- (B8 - B4)/(B8 + B4) #NIR | red
  ndwi_stack <- (B8 - B3)/(B8 + B3) #green | NIR
  
  cat("- smoothing...\n")
  B2_smoothed <- calc(B2, fun = smooth_and_interpolate)
  B3_smoothed <- calc(B3, fun = smooth_and_interpolate)
  B4_smoothed <- calc(B4, fun = smooth_and_interpolate)
  B5_smoothed <- calc(B5, fun = smooth_and_interpolate)
  B6_smoothed <- calc(B6, fun = smooth_and_interpolate)
  B7_smoothed <- calc(B7, fun = smooth_and_interpolate)
  B8_smoothed <- calc(B8, fun = smooth_and_interpolate)
  B8A_smoothed <- calc(B8A, fun = smooth_and_interpolate)
  ndvi_smoothed <- calc(ndvi_stack, fun = smooth_and_interpolate)
  ndwi_smoothed <- calc(ndwi_stack, fun = smooth_and_interpolate)
  
  # Add smoothed NDVI and NDWI to the smoothed stack
  cat("- Creating the final stack...\n")
  
  # Calculate Statistics VV, VG and Entropy
  vh_st <- stack(app(rast(VH), fun = compute_stats))
  vv_en <- stack(app(rast(VV), fun = calculate_entropy))
  vv_st <- stack(app(rast(VV), fun = compute_stats))
  vh_en <- stack(app(rast(VH), fun = calculate_entropy))
  
  # Resample CH + Terain
  ch_crop <- raster::resample(ch, vh_st, 'bilinear')
  ter_crop <- raster::resample(ter, vh_st, 'bilinear')
  slope_cr <- terrain(ter_crop, opt = 'slope', unit = 'degrees')
  aspect_cr <- terrain(ter_crop, opt = 'aspect', unit = 'degrees')
  ter_crop <- stack(ter_crop, slope_cr, aspect_cr)
  
  # Combine NDVI and NDWI with the rearranged stack
  cat("- Combining Sentinel-2, Sentinel-1, NDVI, and NDWI layers...\n")
  combined_stack <- stack(B2_smoothed, B3_smoothed, B4_smoothed, B5_smoothed, B6_smoothed, B7_smoothed, B8_smoothed, B8A_smoothed, ndvi_smoothed, ndwi_smoothed, VV, VH, ch_crop, ter_crop, vv_st, vh_st, vv_en, vh_en)
  
  # Save the smoothed stack
  output_path <- getwd()
  output_file <- paste0(output_path, tools::file_path_sans_ext(basename(file)), "_processed.tif")
  cat(sprintf("- Saving processed raster to %s...\n", output_file))
  writeRaster(combined_stack, output_file, format = "GTiff", overwrite = TRUE)
  
  # Record end time
  end_time <- Sys.time()
  cat(sprintf("End time: %s\n", end_time))
  
  # Calculate processing duration
  processing_duration <- difftime(end_time, start_time, units = "secs")
  cat(sprintf("Processing time for %s: %f seconds\n", basename(file), processing_duration))
  
  
  cat(sprintf("- Finished processing %s\n", basename(file)))
}
