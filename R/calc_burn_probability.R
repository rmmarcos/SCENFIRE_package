#' @title Convert FLP20 file to a data frame
#'
#' @description
#' Reads a .csv file, typically an FLP20 output from FConstMTT, for single fire event simulation,
#' and reshapes the data to long format, calculating flame length (FL) from `FIL` columns.
#' `PBurn` is expected to be 1 for individual fires, and `FIL` columns are binary (1 or 0).
#'
#' @param file A character string specifying the path to the FLP20 .csv file.
#' @return A data frame with reshaped data, including `name` (original FIL column index)
#'   and `FL` (calculated flame length).
#' @importFrom readr read_csv
#' @importFrom dplyr filter mutate
#' @importFrom tidyr pivot_longer starts_with
#' @export
#'
#' @examples
#' # Create a dummy FLP20 CSV file representing individual fires
#' # Each row is a unique fire event at a specific XPos, YPos.
#' # PBurn is 1, and only one FIL column is 1 for the flame length category.
#' temp_flp20_file <- tempfile(fileext = ".csv")
#' dummy_data <- data.frame(
#'   XPos = c(10, 20, 10, 30),
#'   YPos = c(10, 20, 20, 10),
#'   PBurn = c(1, 1, 1, 1), # PBurn set to 1 for individual fires
#'   FIL1 = c(1, 0, 0, 0),  # Fire 1 is in FL bin 1 (0.25)
#'   FIL2 = c(0, 1, 0, 0),  # Fire 2 is in FL bin 2 (0.75)
#'   FIL3 = c(0, 0, 1, 0),  # Fire 3 is in FL bin 3 (1.25)
#'   FIL4 = c(0, 0, 0, 1)   # Fire 4 is in FL bin 4 (1.75)
#'   # FIL5 to FIL20 would be 0 for these examples
#' )
#' write.csv(dummy_data, temp_flp20_file, row.names = FALSE)
#'
#' # Process the dummy file
#' df_result <- flp20_to_df(temp_flp20_file)
#' print(head(df_result))
#'
#' # Expected output for first few rows after processing:
#' # XPos YPos PBurn name    value   FL
#' # 1   10   10     1    1       1 0.25
#' # 2   20   20     1    2       1 0.75
#' # 3   10   20     1    3       1 1.25
#' # 4   30   10     1    4       1 1.75
#'
#' # Clean up the dummy file
#' unlink(temp_flp20_file)
flp20_to_df <- function(file){

  df <- read.csv(file)

  # Ensure FIL columns are read, even if not all 20 are present in smaller files
  # This makes the code more robust if file only has FIL1, FIL2, FIL3 etc.
  fil_cols_present <- names(df)[grepl("^FIL[0-9]+$", names(df))]

  df <- df |>
    filter(PBurn > 0) |> # Still filters for positive PBurn (e.g., PBurn = 1)
    mutate(path=tools::file_path_sans_ext(basename(file))) |>
    pivot_longer(cols = all_of(fil_cols_present), names_to = "name", values_to = "value") |>
    # Filter out entries where FILX was 0 (i.e., this fire was not in this FL bin)
    filter(value == 1) |>
    mutate(name = sub("FIL", "", name), # Extracts the number from FIL1, FIL2, etc.
           FL = ((as.numeric(name) * 0.5) - 0.25))
}


 #' @title Process multiple FLP20 files into a combined data frame
 #'
 #' @description
 #' Applies `flp20_to_df` to a list of file paths and combines the resulting data frames
 #' into a single, large data frame. Each file represents individual fire events.
 #'
 #' @param files A character vector of paths to FLP20 .csv files.
 #' @return A single, combined data frame containing processed data from all input files.
 #' @importFrom dplyr bind_rows
 #' @export
 #'
 #' @examples
 #' \dontrun{
 #' # In a real scenario, 'files' would be paths to your actual FLP20 CSV files.
 #' # For this example, we'll create a few dummy files representing individual fires.
 #'
 #' temp_dir <- tempdir()
 #' file1 <- file.path(temp_dir, "flp20_fire1.csv")
 #' file2 <- file.path(temp_dir, "flp20_fire2.csv")
 #'
 #' # Create dummy data for file1 (one fire, FL bin 10)
 #' dummy_data1 <- data.frame(
 #'   XPos = 50, YPos = 50, PBurn = 1,
 #'   FIL1 = 0, FIL2 = 0, FIL3 = 0, FIL4 = 0,
 #'   FIL5 = 0, FIL6 = 0, FIL7 = 0, FIL8 = 0,
 #'   FIL9 = 0,
 #'   FIL10 = 1, # This fire is in FL bin 10
 #'   FIL11 = 0, FIL12 = 0, FIL13 = 0, FIL14 = 0,
 #'   FIL15 = 0, FIL16 = 0, FIL17 = 0, FIL18 = 0,
 #'   FIL19 = 0, FIL20 = 0
 #' )
 #' write.csv(dummy_data1, file1, row.names = FALSE)
 #'
 #' # Create dummy data for file2 (one fire, FL bin 15)
 #' dummy_data2 <- data.frame(
 #'   XPos = 70, YPos = 80, PBurn = 1,
 #'   FIL1 = 0, FIL2 = 0, FIL3 = 0, FIL4 = 0,
 #'   FIL5 = 0, FIL6 = 0, FIL7 = 0, FIL8 = 0,
 #'   FIL9 = 0,
 #'   FIL10 = 0, FIL11 = 0, FIL12 = 0, FIL13 = 0,
 #'   FIL14 = 0, FIL15 = 1, # This fire is in FL bin 15
 #'   FIL16 = 0, FIL17 = 0, FIL18 = 0, FIL19 = 0, FIL20 = 0
 #' )
 #' write.csv(dummy_data2, file2, row.names = FALSE)
 #'
 #' # List of files to process
 #' my_files <- c(file1, file2)
 #'
 #' # Process the files
 #' combined_df <- flp20_to_bp_df(my_files)
 #' print(head(combined_df))
 #'
 #' # Clean up dummy files
 #' unlink(my_files)
 #' }
 flp20_to_bp_df <- function(files){

   message('Processing files') # Changed print to message, common for package functions
   # Using `files` parameter directly
   fl.dfs <- lapply(files, function(x) flp20_to_df(x))
   foo <- do.call(bind_rows, fl.dfs)

   return(foo)
  }

#' @title Convert processed FLP data to raster format and calculate conditional annual burn probability
#'
#' @description
#' Filters a data frame of fire data based on a flame length threshold,
#' groups by spatial position (`XPos`, `YPos`), and summarizes it into
#' conditional annual burning probability (CBP) estimates and
#' mean flame length (FL_mean). These summarized points are then rasterized
#' to a reference raster. `PBurn` is treated as a count of individual fires at a location.
#'
#' @param df A data frame processed by `flp20_to_df` or `flp20_to_bp_df`,
#'   containing `XPos`, `YPos`, `FL`, and `PBurn` (expected to be 1 for individual fires) columns.
#' @param fl_threshold A numeric value representing the minimum flame length (FL)
#'   to include in the summary.
#' @param selected_surf A data frame or list containing a `size` column, used to
#'   calculate the total selected surface.
#' @param reference_surface A numeric value representing a reference surface area,
#'   used in the BP calculation.
#' @param r_ref A terra SpatRaster object that serves as the reference raster for
#'   rasterization.
#' @return A list containing two raster objects: `CBP` (Burning Probability),
#'    `CFL` (Mean Flame Length) and `ID_fires` (list of perimeter IDs).
#' @importFrom dplyr filter group_by summarise
#' @importFrom sf st_as_sf
#' @importFrom terra rast rasterize
#' @export
#'
#' @examples
#' \dontrun{
#' # This example requires 'terra' and 'sf' packages.
#'
#' # Create dummy data frame 'df_example' to simulate output from flp20_to_df
#' # Each row here represents an individual fire at a specific location, with its
#' # derived FL bin value. PBurn is 1.
#' df_example <- data.frame(
#'   XPos = c(10, 10, 20, 20, 30, 30),
#'   YPos = c(10, 10, 10, 20, 20, 10),
#'   PBurn = c(1, 1, 1, 1, 1, 1), # All individual fires, PBurn = 1
#'   # FL values derived from FIL bins (e.g., 0.25 for FIL1, 0.75 for FIL2, etc.)
#'   FL = c(0.25, 0.75, 0.25, 1.25, 0.75, 0.75)
#' )
#'
#' # Create dummy 'selected_surf' data
#' selected_surf_example <- data.frame(size = c(1000, 2000, 500, 1500))
#'
#' # Define a reference surface
#' reference_surface_example <- 10000
#'
#' # Create a dummy SpatRaster for r_ref
#' # In a real scenario, this would be your actual reference raster.
#' # Make sure the extent covers the XPos, YPos range of dummy_df
#' r_ref_example <- terra::rast(
#'   nrows = 3, ncols = 3,
#'   xmin = 0, xmax = 40, ymin = 0, ymax = 30,
#'   crs = "EPSG:4326" # Example CRS
#' )
#'
#' # Run the function with example data
#' raster_results <- fpl_to_raster(
#'   df = df_example,
#'   fl_threshold = 0.5, # Only include fires with FL >= 0.5 (e.g., FL 0.75 and 1.25)
#'   selected_surf = selected_surf_example,
#'   reference_surface = reference_surface_example,
#'   r_ref = r_ref_example
#' )
#'
#' # Print and plot the results (optional)
#' print(raster_results$CBP)
#' plot(raster_results$CBP, main = "CBP Raster Example")
#' plot(raster_results$CFL, main = "CFL Raster Example")
#' }
flp20_to_raster <- function(df, fl_threshold, selected_surf, reference_surface, r_ref){

  foo <- df |>
    filter(FL >= fl_threshold) |>
    mutate(ID_fire = str_split_i(path, '_', 6))

  list.fires <- unique(foo$ID_fire)

  foo <- foo |> # Filter value > 0 is handled in flp20_to_df now
    group_by(XPos, YPos) |>
    summarise(
      # sum(PBurn) will now effectively count fires at this location
      BP = sum(PBurn) / ((sum(selected_surf$size) / reference_surface)),
      # weighted.mean(FL, w=PBurn) becomes simple mean(FL) if PBurn is always 1
      FL_mean = weighted.mean(FL, w = PBurn)
    )
  # Rasterize to the reference raster
  points_sf <- st_as_sf(foo, coords = c("XPos", "YPos"), crs = NA)
  bp <- rasterize(points_sf, r_ref, field = "BP", fun = "mean")
  fl <- rasterize(points_sf, r_ref, field = "FL_mean", fun = "mean")

  return(list(CBP = bp, CFL = fl, ID_fires = list.fires))
  }

#' Calculate annual burned probability (BP) from perimeters
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
calc_abp <- function(template_raster,candidate_surfaces,reference_surface) {

  # Convert sf object to terra SpatVector
  polygons <- vect(candidate_surfaces)
  # Ensure geometries are valid
  polygons <- makeValid(polygons)

  # Rasterize the number of times a pixel has been burned
  # 'field = NA' means it will count the occurrences of polygons
  # 'fun = "count"' counts how many polygons overlap each cell
  # 'background = 0' sets cells with no overlap to 0
  bp <- rasterize(polygons, template_raster, fun = "sum", background = 0)
  # Normalize the burned count by the ratio of total simulated area to reference surface
  # This scales the BP to represent probability relative to the target area
  bp <- bp/(sum(candidate_surfaces$size)/reference_surface)

  return(bp)

}
