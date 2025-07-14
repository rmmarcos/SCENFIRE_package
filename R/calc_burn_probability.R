#' Convert Fire Length Points (FLP) data to a raster layer
#'
#' This function reads a CSV file containing Fire Length Points (FLP) data,
#' processes the `FIL1` to `FIL6` columns to determine an estimated flame length (FL)
#' based on predefined intervals (e.g., FIL1 = 0-0.6m, FIL2 = 0.6-1.2m, etc.).
#' It then converts these points into a `terra` SpatVector and rasterizes them
#' onto a reference raster grid, aggregating values by mean for overlapping points.
#'
#' @param n Integer index indicating which FLP file from `flp.files` to process.
#' @param r A `terra::SpatRaster` object that serves as the reference grid
#'   for rasterization. The extent and resolution of `r` will define the output raster.
#' @return A `terra::SpatRaster` object where each cell contains the mean
#'   flame length (FL) value derived from the FLP data.
#' @importFrom utils read.csv
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr filter mutate select case_when
#' @importFrom terra vect rasterize
#' @export
#' @examples
#' \dontrun{
#' # Assuming 'flp.files' is a character vector of paths to FLP CSVs
#' # and 'reference_raster' is a pre-existing terra SpatRaster object.
#' #
#' # # Example reference raster (replace with your actual raster)
#' # library(terra)
#' # r_ref_example <- rast(nrows=10, ncols=10, xmin=0, xmax=100, ymin=0, ymax=100,
#' #                       crs="EPSG:25830")
#' #
#' # # Example flp.files (replace with your actual file paths)
#' # # Create a dummy CSV for demonstration
#' # dummy_flp_path <- tempfile(fileext = ".csv")
#' # write.csv(data.frame(XPos=c(10,20,30), YPos=c(10,20,30),
#' #                      FIL1=c(1,0,0), FIL2=c(0,1,0), FIL3=c(0,0,1),
#' #                      FIL4=0, FIL5=0, FIL6=0),
#' #           dummy_flp_path, row.names = FALSE)
#' # flp.files <- c(dummy_flp_path)
#' #
#' # fl_raster <- flp_to_raster(n = 1, r = r_ref_example)
#' # plot(fl_raster, main = "Flame Length Raster")
#' # unlink(dummy_flp_path) # Clean up dummy file
#' }
flp_to_raster <- function(n,r){

  foo <- read.csv(flp.files[n]) |>
    pivot_longer(cols = FIL1:FIL6) |>
    filter(value>0) |>
    mutate(FL = case_when(

      value == 1 & name == "FIL1" ~ 0.3,
      value == 1 & name == "FIL2" ~ 0.6,
      value == 1 & name == "FIL3" ~ 1.2,
      value == 1 & name == "FIL4" ~ 1.8,
      value == 1 & name == "FIL5" ~ 2.4,
      value == 1 & name == "FIL6" ~ 3.7,
    )) |>
    select(XPos,YPos,FL)

  # fl.r <- rast(foo,"xyz") # Original comment for reference

  xyz_points <- vect(foo, geom = c("XPos", "YPos"))
  fl.r <- rasterize(xyz_points, r, field = "FL", fun = "mean") 	# or fun = "first", etc.

  return(fl.r)

}


#' Convert 20-Category Fire Length Points (FLP20) to Raster
#'
#' This function reads a specific type of Fire Length Points (FLP20) CSV file,
#' which contains 20 categories of flame length. It processes these categories
#' to calculate a central flame length (`llama_m`) for each point and then
#' rasterizes these points onto a provided reference raster grid.
#'
#' @param v A data frame or tibble that contains the `path` to the FLP20 CSV files.
#'   This parameter defaults to `candidate_surfaces_cleansed`, implying it's often
#'   used within a larger workflow where this variable is available.
#' @param n Integer index indicating which file path from `v` to process.
#' @param r_ref A `terra::SpatRaster` object that serves as the reference grid
#'   for rasterization. Its CRS will be used for the new `SpatVector`.
#' @return A `terra::SpatRaster` object where each cell contains the mean
#'   flame length (`llama_m`) value, rasterized to the extent and resolution of `r_ref`.
#' @importFrom readr read_csv
#' @importFrom terra vect rasterize crs
#' @export
#' @examples
#' \dontrun{
#' # Assuming 'candidate_surfaces_cleansed' is a data frame with a 'path' column
#' # and 'reference_raster' is a pre-existing terra SpatRaster object.
#' #
#' # # Example reference raster
#' # library(terra)
#' # r_ref_example <- rast(nrows=10, ncols=10, xmin=0, xmax=100, ymin=0, ymax=100,
#' #                       crs="EPSG:25830")
#' #
#' # # Create a dummy directory and FLP20 file for demonstration
#' # dummy_dir <- tempdir()
#' # dummy_flp20_path <- file.path(dummy_dir, "Run_FLP_filtered.csv")
#' # write.csv(data.frame(XPos=c(10,20), YPos=c(10,20),
#' #                      FIL1=c(1,0), FIL2=c(0,1),
#' #                      FIL3=0, FIL4=0, FIL5=0, FIL6=0, FIL7=0, FIL8=0, FIL9=0, FIL10=0,
#' #                      FIL11=0, FIL12=0, FIL13=0, FIL14=0, FIL15=0, FIL16=0, FIL17=0,
#' #                      FIL18=0, FIL19=0, FIL20=0),
#' #           dummy_flp20_path, row.names = FALSE)
#' #
#' # # Dummy 'v' data frame
#' # dummy_v <- data.frame(path = dummy_dir)
#' #
#' # flp20_raster_out <- flp20_to_raster(v = dummy_v, n = 1, r_ref = r_ref_example)
#' # plot(flp20_raster_out, main = "FLP20 Flame Length Raster")
#' # unlink(dummy_dir, recursive = TRUE) # Clean up dummy directory
#' }
flp20_to_raster <- function(v=candidate_surfaces_cleansed,n,r_ref){

  # === 2. Read the CSV ===
  df <- read_csv(paste0(v[n,]$path_rel,'/Run_FLP_filtered.csv'))

  # === 3. Identify flame length (center of the interval)
  fil_cols <- paste0("FIL", 1:20)
  intervalos <- seq(0.25, 9.75, by = 0.5) 	# center of each interval

  # For each row, find which FIL has value 1
  df$llama_m <- apply(df[, fil_cols], 1, function(row) {
    idx <- which(row == 1)
    if (length(idx) == 1) {
      return(intervalos[idx])
    } else {
      return(NA) 	# handle cases without 1s or multiple 1s
    }
  })

  # === 4. Create a raster with the values ===
  # Convert to SpatVector
  v <- vect(df[, c("XPos", "YPos", "llama_m")],
            geom = c("XPos", "YPos"),
            crs = crs(r_ref))

  # Rasterize to the reference raster
  r_llama <- rasterize(v, r_ref, field = "llama_m", fun = "mean")

  # === 5. Save output raster ===
  # writeRaster(r_llama, paste0(v[n,]$path,'/FL.tif'), overwrite = TRUE) # Original comment
  return(r_llama)
}


#' Calculate Conditional Burned Probability (CBP)
#'
#' This function calculates the conditional burned probability from a stack of
#' flame length raster layers. For each pixel, it counts how many layers have
#' a flame length value greater than a specified threshold (`FL`) and then
#' calculates this as a proportion of the total number of layers.
#' Non-numeric (NA) values in the raster are treated as 0 (not burned).
#' This function is identical in its implementation to `calc_bp` but is named
#' `calc_cbp` to reflect its specific use case for conditional probability based on flame length.
#'
#' @param r_stack A `terra::SpatRaster` object containing multiple layers.
#'   Each layer typically represents fire presence/intensity (e.g., flame length).
#' @param FL Numeric value. The flame length threshold. Pixels with values
#'   greater than `FL` are considered burned (default is 0).
#' @return A numeric value representing the overall conditional burned proportion
#'   across the raster stack.
#' @importFrom terra app nlyr
#' @export
#' @examples
#' # Create a dummy raster stack for demonstration
#' library(terra)
#' r1 <- rast(nrows=10, ncols=10, vals=sample(0:5, 100, replace=TRUE))
#' r2 <- rast(nrows=10, ncols=10, vals=sample(0:5, 100, replace=TRUE))
#' r_stack_example <- c(r1, r2)
#'
#' # Calculate conditional burned probability with default FL=0
#' cbp_result_0 <- calc_cbp(r_stack = r_stack_example, FL = 0)
#' print(paste("Conditional Burned Probability (FL=0):", cbp_result_0))
#'
#' # Calculate conditional burned probability with FL=2
#' cbp_result_2 <- calc_cbp(r_stack = r_stack_example, FL = 2)
#' print(paste("Conditional Burned Probability (FL=2):", cbp_result_2))
calc_cbp <- function(r_stack,FL=0){

  r <- app(r_stack, fun = function(x) ifelse(x >= FL, 1, 0))
  r[is.na(r)] <- 0

  cbp <- sum(r)/nlyr(r)

  return(cbp)
}

#' Calculate Burned Probability (BP) from Perimeters
#'
#' This function calculates the burned probability raster from a set of
#' simulated fire perimeters. It first converts the perimeters into a `terra`
#' vector object, ensures their validity, and then rasterizes them onto a
#' template raster. The resulting raster indicates the number of times each
#' pixel has been burned, which is then normalized by a scaling factor derived
#' from the total simulated area and a reference surface.
#'
#' @param template_raster A `terra::SpatRaster` object that defines the
#'   extent, resolution, and CRS for the output burned probability raster.
#' @param candidate_surfaces An `sf` object (simple features) representing
#'   the simulated fire perimeters. It is assumed to have a column named `size`
#'   containing the area of each perimeter.
#' @param reference_surface Numeric value. A reference total surface area used
#'   for normalization. This typically represents the target total burned area
#'   (e.g., historical average annual burned area).
#' @return A `terra::SpatRaster` object representing the burned probability.
#'   Each cell value indicates the proportion of times that pixel was burned,
#'   normalized by the ratio of total simulated area to the reference surface.
#' @importFrom terra vect makeValid rasterize
#' @importFrom sf st_geometry
#' @export
#' @examples
#' \dontrun{
#' # Create a dummy template raster
#' library(terra)
#' template_raster_ex <- rast(nrows=10, ncols=10, xmin=0, xmax=100, ymin=0, ymax=100,
#'                            crs="EPSG:25830")
#'
#' # Create dummy candidate_surfaces (sf object with polygons)
#' # Ensure 'size' column exists
#' library(sf)
#' poly1 <- st_polygon(list(cbind(c(10,10,30,30,10), c(10,30,30,10,10))))
#' poly2 <- st_polygon(list(cbind(c(50,50,70,70,50), c(50,70,70,50,50))))
#' poly3 <- st_polygon(list(cbind(c(20,20,40,40,20), c(20,40,40,20,20))))
#'
#' candidate_surfaces_ex <- st_sf(
#'   id = 1:3,
#'   size = c(400, 400, 400), # Dummy sizes for example
#'   geometry = st_sfc(poly1, poly2, poly3),
#'   crs = 25830
#' )
#'
#' # Example reference surface
#' reference_surface_ex <- 1000
#'
#' # Calculate burned probability
#' bp_raster <- calc_bp(template_raster = template_raster_ex,
#'                      candidate_surfaces = candidate_surfaces_ex,
#'                      reference_surface = reference_surface_ex)
#'
#' plot(bp_raster, main = "Burned Probability Raster")
#' }
calc_bp <- function(template_raster,candidate_surfaces,reference_surface) {

  # Convert sf object to terra SpatVector
  polygons <- vect(candidate_surfaces)
  # Ensure geometries are valid
  polygons <- makeValid(polygons)

  # Rasterize the number of times a pixel has been burned
  # 'field = NA' means it will count the occurrences of polygons
  # 'fun = "count"' counts how many polygons overlap each cell
  # 'background = 0' sets cells with no overlap to 0
  bp <- rasterize(polygons, template_raster, field = NA, fun = "count", background = 0)
  # Normalize the burned count by the ratio of total simulated area to reference surface
  # This scales the BP to represent probability relative to the target area
  bp <- bp/(sum(candidate_surfaces$size)/reference_surface)

  return(bp)

}
