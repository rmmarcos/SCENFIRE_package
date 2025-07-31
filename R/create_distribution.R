#' Function to build power law distribution to fit a theortical target histogram
#'
#' Generates random numbers that follow a truncated power-law distribution.
#' This function uses the inverse cumulative distribution function (CDF)
#' to sample from the specified power-law within given minimum and maximum bounds.
#'
#' @param n The number of random variates (points) to generate (default 100).
#' @param alpha The exponent (alpha parameter) of the power-law distribution (default 1.5).
#' @param xmin The lower bound (minimum value) of the truncated distribution.
#' @param xmax The upper bound (maximum value) of the truncated distribution.
#' @return A numeric vector of `n` random variates sampled from the truncated power-law distribution.
#' @importFrom stats runif
#' @export
#' @examples
#' # Generate 100 random numbers from a power-law distribution
#' # with alpha=2.5, xmin=10, xmax=1000
#' set.seed(123)
#' powerlaw_samples <- fit_powerlaw(n = 100, alpha = 2.5, xmin = 10, xmax = 1000)
#' # hist(powerlaw_samples, main = "Power-Law Distribution Samples")
fit_powerlaw <- function(n=100, alpha=1.5, xmin, xmax) {
  # Inverse CDF of the truncated distribution
  r <- runif(n)
  C <- (xmax^(1 - alpha) - xmin^(1 - alpha))
  return((r * C + xmin^(1 - alpha))^(1 / (1 - alpha)))
}

#' Function to obtain the target histogram from a vector of fire size data
#'
#' This function builds a histogram (specifically, its density) from a set of
#' fire sizes, which is intended to represent a target distribution (e.g., from
#' historical fire events). It allows for an optional logarithmic transformation
#' of the fire sizes before binning and determines bin breakpoints dynamically.
#'
#' @param num_bins Integer. The desired number of bins for the histogram (default: 20).
#' @param logaritmic Logical. If `TRUE`, a logarithmic transformation (`log(sizes + 1e-6)`)
#'   is applied to `sizes`. If `FALSE`, the original scale is used. This parameter
#'   influences both the distribution transformation and the binning.
#' @param sizes Numeric vector. The fire sizes (e.g., area in hectares) from the
#'   data used to build the target distribution.
#' @param event_surfaces Numeric vector. The fire sizes (e.g., area in hectares) from the
#'   simulated fire perimeters.
#' @return A list containing two elements:
#'   \describe{
#'     \item{target_hist}{A numeric vector representing the density values of the target histogram.}
#'     \item{bins}{A numeric vector specifying the breakpoints (edges) of the bins used in the histogram.}
#'   }
#' @export
#' @examples
#' # Dummy data for demonstration (replace with your actual data)
#' set.seed(123)
#' historical_data_for_target <- floor(fit_powerlaw(n = 500, alpha = 2, xmin = 10, xmax = 10000))
#' event_surfaces <- fit_powerlaw(n = 10000, alpha = 2, xmin = 10, xmax = 10000)
#'
#' # Discard simulated fires that are too large (below 110% max historical size)
#' event_surfaces <- event_surfaces[event_surfaces < max(historical_data_for_target) * 1.1]
#'
#' # Default configuration: logaritmic transformation on fire size
#' target_info_example <- build_target_hist(num_bins = 10, logaritmic = TRUE,
#'                                          sizes = historical_data_for_target,
#'                                          event_surfaces = event_surfaces)
#'
#' # Print results
#' target_hist <- target_info_example$target_hist
#' print(target_hist)
#' bins <- target_info_example$bins
#' print(bins)
#'
#' # Alternate configuration: original frequency distribution on fire size
#' target_info_example <- build_target_hist(num_bins = 10, logaritmic = FALSE,
#'                                          sizes = historical_data_for_target,
#'                                          event_surfaces = event_surfaces)
#'
#' # Print results
#' target_hist <- target_info_example$target_hist
#' print(target_hist)
#' bins <- target_info_example$bins
#' print(bins)
build_target_hist <- function(num_bins=20, logaritmic=T, sizes, event_surfaces){

  if(logaritmic == T){
    # If logarithm is used, the distribution of historical fire sizes is transformed
    target_distribution <- log(sizes + 1e-6) 	# A small value is added to avoid log(0)

    # Bin limits are adjusted based on the logarithmic distribution
    bins <- seq(min(c(log(event_surfaces + 1e-6), target_distribution, na.rm = TRUE)),
                max(c(log(event_surfaces + 1e-6), target_distribution, na.rm = TRUE)),
                length.out = num_bins + 1)
  } else {
    # If logarithm is not used, the original scale is maintained
    target_distribution <- (sizes)

    # Bin limits are adjusted based on the original distribution
    bins <- seq(min(c((event_surfaces), target_distribution, na.rm = TRUE)),
                max(c((event_surfaces), target_distribution, na.rm = TRUE)),
                length.out = num_bins + 1)
  }

  # The density of the historical fire distribution is calculated
  target_hist <- hist(target_distribution, breaks = bins, plot = FALSE)$density

  return(list(target_hist=target_hist, bins=bins))


}


#' Function to calculate the relative distance between the target distribution and the selected perimeters
#'
#' Calculates the relative absolute discrepancy between the density histogram of
#' a set of `selected_surfaces` and a `target_hist`. This metric quantifies
#' how well the distribution of the selected surfaces matches the target distribution.
#' It can optionally apply a logarithmic transformation to the selected surfaces.
#'
#' @param selected_surfaces Numeric vector of surface values for the selected events.
#' @param target_hist Numeric vector representing the density values of the target histogram.
#' @param bins Numeric vector of bin breakpoints used for both histograms.
#' @param logaritmic Logical. If `TRUE`, a logarithmic transformation (`log(selected_surfaces + 1e-6)`)
#'   is applied to `selected_surfaces` before calculating its histogram density.
#' @return A numeric value representing the relative discrepancy. A value of 0 indicates a perfect match.
#' @export
#' @examples
#' # Example target histogram and bins (from build_target_hist)
#' historical_sizes_ex <- c(10, 50, 100, 200, 500, 1000)
#' # Dummy 'event_surfaces' for example context (replace with actual data)
#' event_surfaces <- c(5, 15, 25, 35, 45, 55)
#' target_info_ex <- build_target_hist(num_bins = 5, logaritmic = TRUE, sizes = historical_sizes_ex)
#' target_hist_ex <- target_info_ex$target_hist
#' bins_ex <- target_info_ex$bins
#'
#' # Example selected surfaces
#' simulated_surfaces_ex <- c(20, 80, 450, 900)
#'
#' # Calculate discrepancy with logarithmic transformation
#' discrepancy_log <- calculate_discrepancy(selected_surfaces = simulated_surfaces_ex,
#'                                          target_hist = target_hist_ex,
#'                                          bins = bins_ex,
#'                                          logaritmic = TRUE)
#' print(paste("Discrepancy (log):", discrepancy_log))
#'
#' # Example selected surfaces for linear transformation
#' simulated_surfaces_linear_ex <- c(15, 60, 110, 210, 480, 950)
#' target_info_linear_ex <- build_target_hist(num_bins = 5, logaritmic = FALSE, sizes = historical_sizes_ex)
#' target_hist_linear_ex <- target_info_linear_ex$target_hist
#' bins_linear_ex <- target_info_linear_ex$bins
#'
#' # Calculate discrepancy with linear transformation
#' discrepancy_linear <- calculate_discrepancy(selected_surfaces = simulated_surfaces_linear_ex,
#'                                             target_hist = target_hist_linear_ex,
#'                                             bins = bins_linear_ex,
#'                                             logaritmic = FALSE)
#' print(paste("Discrepancy (linear):", discrepancy_linear))
calculate_discrepancy <- function(selected_surfaces, target_hist, bins, logaritmic=T) {
  # Calculate the histogram density of the selected set
  if(logaritmic==T){
    hist_selected <- hist(log(selected_surfaces + 1e-6), breaks = bins, plot = FALSE)$density
  }else{
    hist_selected <- hist((selected_surfaces), breaks = bins, plot = FALSE)$density
  }


  # Ensure that densities have the same length as the target histogram
  if (length(hist_selected) != length(target_hist)) {
    stop("Selected histogram densities do not match target histogram length.")
  }

  # Calculate the absolute difference between densities
  abs_diff <- abs(hist_selected - target_hist)

  # Calculate discrepancy as relative area
  total_area_target <- sum(target_hist) 	# Total area of the target histogram
  if (total_area_target == 0) {
    stop("The total area of the target histogram is zero.")
  }

  relative_discrepancy <- sum(abs_diff) / total_area_target 	# Relative discrepancy

  return(relative_discrepancy)
}


#' Function to select simulated perimeters
#'
#' Selects a subset of simulated fire perimeters by matching their surface
#' distribution to a predefined target histogram. The selection is iterative
#' and probabilistic, aiming to minimize the discrepancy while accumulating
#' a total surface area above a certain threshold (e.g., mean annul burned area).
#'
#' @param event_sizes Numeric vector of surface values for all available
#'   simulated events.
#' @param event_probabilities Numeric vector of probabilities corresponding to each event in `event_surfaces`.
#'   These probabilities are used to influence the selection of events within each bin.
#' @param target_hist Numeric vector representing the density of the target histogram distribution.
#' @param bins Numeric vector of bin breakpoints used for classifying event surfaces and calculating histograms.
#' @param reference_surface Numeric value representing the total target surface area that the selected events should approximate.
#' @param surface_threshold Numeric value between 0 and 1. The selection process
#'   continues until the cumulative surface area of selected events is at least
#'   `reference_surface * surface_threshold`.
#' @param tolerance Numeric value. A tolerance level for the final discrepancy. (Note:
#'   the current implementation finds the best discrepancy, not necessarily stopping
#'   once tolerance is met, but aims for the minimum).
#' @param max_it Integer for the maximum number of iterations for the inner loop (default: 100).
#' @param iter_limit Integer for the maximum number of iterations in selection (default: 100000).
#' @param logaritmic Logical. If `TRUE`, a logarithmic transformation is applied
#'   to `event_surfaces` before binning and to `selected_surfaces` for histogram
#'   calculations (default: `TRUE`).
#' @return A list containing the best selection found across all iterations:
#'   \describe{
#'     \item{selected_surfaces}{Numeric vector of surface values of the events in the best selection.}
#'     \item{surface_index}{Integer vector of the original indices (from `event_surfaces`)
#'       of the events in the best selection.}
#'     \item{total_surface}{Numeric value. The sum of surface areas of the events in the best selection.}
#'     \item{final_discrepancy}{Numeric value. The relative discrepancy between the
#'       selected events' distribution and the target histogram for the best selection.}
#'   }
#'   Returns `NULL` if no valid selection could be made (e.g., no valid results after iterations).
#' @importFrom foreach foreach %dopar%
#' @export
#' @examples
#' \dontrun{
#' # This example requires the 'foreach' package and a parallel backend (e.g., 'doParallel')
#' # to be set up.
#' # library(doParallel)
#'
#' # Dummy data for demonstration (replace with your actual data)
#' set.seed(123)
#' historical_data_for_target <- floor(fit_powerlaw(n = 500, alpha = 2, xmin = 10, xmax = 10000))
#' event_surfaces <- fit_powerlaw(n = 10000, alpha = 2, xmin = 10, xmax = 10000)
#' # Discard simulated fires that are too large (below 110% max historical size)
#' event_surfaces <- event_surfaces[event_surfaces<max(historical_data_for_target)*1.1]
#' event_probabilities <- rnorm(length(event_surfaces))
#' event_probabilities <- (event_probabilities-min(event_probabilities))/(max(event_probabilities)-min(event_probabilities))
#'
#' y <- 100 #number of years spanning historical fire data
#' check_fire_data(fires_hist_size = historical_data_for_target,
#'                 sim_perimeters_size = event_surfaces,
#'                 n_years = y)
#'
#' reference_surface_example <- sum(historical_data_for_target)/y
#' surface_threshold_example <- check_fire_data(fires_hist_size = historical_data_for_target,
#'                                              sim_perimeters_size = event_surfaces,
#'                                              n_years = 10)
#' tolerance_example <- 0.1
#'
#' # Create a dummy target histogram (assuming 'event_surfaces' from historical data)
#' # For a real scenario, 'event_surfaces' here would be your historical fire sizes.
#'
#' target_info_example <- build_target_hist(num_bins = 10, logaritmic = TRUE,
#'                                          sizes = historical_data_for_target,
#'                                          event_surfaces = event_surfaces)
#' target_hist <- target_info_example$target_hist
#' bins <- target_info_example$bins
#'
#' # Run the selection process
#' selected_events_result <- select_events(
#'   event_sizes = event_surfaces,
#'   event_probabilities = event_probabilities,
#'   target_hist = target_hist,
#'   bins = bins,
#'   reference_surface = reference_surface_example,
#'   surface_threshold = surface_threshold_example,
#'   tolerance = tolerance_example,
#'   max_it = 2 # Reduced iterations for example
#' )
#'
#' # Stop the parallel cluster when done
#' # stopImplicitCluster()
#' }
select_events <- function(event_sizes, event_probabilities, target_hist, bins,
                          reference_surface, surface_threshold, tolerance, max_it = 5, iter_limit=100000,logaritmic = T) {

  num_cores <- max_it  # Número de núcleos a usar (idealmente todos menos uno)
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)

  event_surfaces <- event_sizes

  # Exporta variables y funciones necesarias al clúster paralelo
  clusterExport(cl, list("calculate_discrepancy",
                         "event_surfaces",
                         "event_probabilities",
                         "bins",
                         "target_hist"))

  max_iterations <- max_it

  # Classify events into bins
  if(logaritmic == T) {
    bin_indices <- cut(log(event_surfaces + 1e-6), breaks = bins, include.lowest = TRUE, labels = FALSE)
  } else {
    bin_indices <- cut(event_surfaces, breaks = bins, include.lowest = TRUE, labels = FALSE)
  }

  results <- ({
    foreach(block = 1:max_iterations,
            .combine = list,
            .multicombine = TRUE,
            .packages = c("stats"),
            .export = "calculate_discrepancy") %dopar% {

              best_selection <- NULL
              best_discrepancy <- Inf

              for (i in 1:max_it) {
                n_iter <- 0
                # Use tryCatch to handle potential errors
                result <- tryCatch({
                  selected_surfaces <- c()
                  surface_index <- c()
                  total_surface <- 0
                  current_hist <- rep(0, length(target_hist)) 	# Initialize histogram of the current selection

                  # Calculate the frequency of the bins with the contingency table
                  bin_counts <- table(bin_indices)
                  bin_priorities <- rep(0, length(bins) - 1) 	# Initialize priorities for each bin

                  for (bin in 1:(length(bins) - 1)) {
                    if (bin %in% names(bin_counts)) {
                      bin_discrepancy <- abs(target_hist[bin] - current_hist[bin])
                      bin_priorities[bin] <- bin_discrepancy
                    }
                  }

                  # Keep a record of selected indices
                  selected_indices_global <- c()

                  while ((total_surface < reference_surface * surface_threshold)) {

                    if(n_iter>iter_limit){
                      cat("Maximum iterations where reached. Execution interrumpted. Increase iter_limit.")
                      break
                    }

                    # Probabilistic selection based on bin priority
                    if (any(bin_priorities > 0)) { 	# Only select if there is any bin with positive priority
                      selected_bin <- sample(1:length(bin_priorities), 1, prob = bin_priorities + 1e-6)

                      # Filter eligible indices excluding previously selected ones
                      eligible_indices <- unique(which(bin_indices == selected_bin)) 	# Ensure uniqueness
                      eligible_indices <- setdiff(eligible_indices, selected_indices_global) 	# Exclude already selected

                      while (length(eligible_indices) > 0) {

                        if(n_iter>iter_limit){
                          cat("Maximum iterations where reached. Execution interrumpted. Increase iter_limit.")
                          break
                        }

                        if (length(eligible_indices) == 1) {
                          selected_index <- eligible_indices[1] 	# Direct selection if only one eligible index
                        } else if (length(eligible_indices) > 1) {
                          selected_index <- sample(eligible_indices, 1, prob = event_probabilities[eligible_indices])
                        } else {
                          break 	# Exit if no eligible indices
                        }

                        # Check if the index has already been selected
                        if (!(selected_index %in% selected_indices_global)) {
                          # If not repeated, process the selection
                          selected_surface <- event_surfaces[selected_index]

                          # Update selection
                          selected_surfaces <- c(selected_surfaces, selected_surface)
                          surface_index <- c(surface_index, selected_index)
                          total_surface <- sum(selected_surfaces)

                          # Register the selected index globally
                          selected_indices_global <- c(selected_indices_global, selected_index)

                          # Update histogram of current selection without logarithmic transformation
                          if(logaritmic == T) {
                            current_hist <- hist(log(selected_surfaces + 1e-6), breaks = bins, plot = FALSE)$density
                          } else {
                            current_hist <- hist(selected_surfaces, breaks = bins, plot = FALSE)$density
                          }

                          break
                        } else {
                          # If repeated, remove it from eligible indices and continue searching
                          eligible_indices <- setdiff(eligible_indices, selected_indices_global)

                          if (length(eligible_indices) == 0) {
                            break
                          }
                        }
                      }
                    } else {
                      break 	# Exit if no bins with priority
                    }

                    temp_discrepancy <- calculate_discrepancy(selected_surfaces, target_hist, bins, logaritmic)
                  }

                  # Calculate the discrepancy
                  final_discrepancy <- calculate_discrepancy(selected_surfaces, target_hist, bins, logaritmic)

                  # Check if it is a valid result
                  if (final_discrepancy < best_discrepancy) {
                    best_selection <- list(
                      selected_surfaces = selected_surfaces,
                      surface_index = surface_index,
                      total_surface = total_surface,
                      final_discrepancy = final_discrepancy
                    )
                    best_discrepancy <- final_discrepancy
                  }

                  best_selection

                }, error = function(e) {
                  cat("Error in iteration", i, ":", conditionMessage(e), "\n")
                  next
                })

                # Store the result of the iteration
                if (!is.null(result)) {
                  return(result)
                }
              }

              # Return the best selection from this block
              best_selection
            }
  })

  # Filter null results and select the best
  valid_results <- Filter(Negate(is.null), results)

  if (length(valid_results) == 0) {
    return(NULL) 	# If no valid results
  }

  # Select the best result
  best_result <- valid_results[[which.min(sapply(valid_results, function(x) x$final_discrepancy))]]

  # Check and remove duplicates
  if (!is.null(best_result)) {
    # Find unique indices keeping the original order
    unique_indices <- !duplicated(best_result$surface_index)

    # Update all components of the result
    best_result$surface_index <- best_result$surface_index[unique_indices]
    best_result$selected_surfaces <- best_result$selected_surfaces[unique_indices]
    best_result$total_surface <- sum(best_result$selected_surfaces)

    # Recalculate the final discrepancy with the unique events
    best_result$final_discrepancy <- calculate_discrepancy(
      best_result$selected_surfaces,
      target_hist,
      bins,
      logaritmic
    )
  }

  return(best_result)
  stopCluster(cl)

}

#' Function to remove duplicated perimeters
#'
#' This function iteratively checks for and removes duplicated geometries from an
#' `sf` object containing fire perimeters. Duplicates are identified using
#' `sf::st_equals()`. The process repeats until no more geometrically identical
#' polygons are found.
#'
#' @param candidates An `sf` object (simple features) containing fire perimeters.
#' @return An `sf` object with all geometrically duplicated perimeters removed.
#' @importFrom sf st_equals
#' @importFrom methods is
#' @export
#' @examples
#' \dontrun{
#' # Create a dummy sf object with some duplicated geometries for demonstration
#' library(sf)
#' library(dplyr)
#' polygon1 <- st_polygon(list(cbind(c(0,0,1,1,0), c(0,1,1,0,0))))
#' polygon2 <- st_polygon(list(cbind(c(1,1,2,2,1), c(1,2,2,1,1))))
#' # Duplicate of polygon1
#' polygon1_dup <- st_polygon(list(cbind(c(0,0,1,1,0), c(0,1,1,0,0))))
#'
#' dummy_perimeters <- st_sf(
#'   id = 1:4,
#'   geometry = st_sfc(polygon1, polygon2, polygon1_dup, st_polygon(list(cbind(c(3,3,4,4,3), c(3,4,4,3,3))))),
#'   crs = 25830
#' )
#'
#' print("Original perimeters with potential duplicates:")
#' print(dummy_perimeters)
#'
#' cleaned_perimeters <- cleanse_duplicates(candidates = dummy_perimeters)
#'
#' print("Perimeters after cleansing duplicates:")
#' print(cleaned_perimeters)
#' }
cleanse_duplicates <- function(candidates){

  candidate_surfaces <- candidates

  duplicate_indices <- list(-1)

  while(length(duplicate_indices)>0){
    # Comprobar si las geometrías están repetidas
    duplicates <- st_equals(candidate_surfaces)

    # Mostrar índices de geometrías duplicadas
    duplicate_indices <- which(sapply(duplicates, length) > 1)

    # Imprimir los resultados
    if (length(duplicate_indices) > 0) {
      print("Se encontraron polígonos duplicados en los siguientes índices:\n")
      print(duplicate_indices)
      candidate_surfaces <- candidate_surfaces[-duplicate_indices,]
    } else {
      print("No se encontraron polígonos duplicados.\n")
    }

  }

  return(candidate_surfaces)
}

#' Visualize selected fire perimeter distribution
#'
#' This function generates a histogram comparing the distribution of selected
#' fire surfaces (e.g., areas) with a predefined target distribution. It also
#' prints key statistics about the selection, such as the number of events,
#' total surface area, and the discrepancy with the target. The visualization
#' automatically adapts to whether a logarithmic transformation was used.
#'
#' @param result A list object, typically the output from the `select_events` function.
#'   It is expected to contain at least the following elements:
#'   \describe{
#'     \item{selected_surfaces}{Numeric vector of surface values for the selected events.}
#'     \item{total_surface}{Numeric value representing the sum of surface areas of the selected events.}
#'     \item{final_discrepancy}{Numeric value representing the final discrepancy metric for the selection.}
#'   }
#' @param logaritmic Logical (assumed to be available in the calling environment).
#'   If `TRUE`, it indicates that a logarithmic transformation was applied to the
#'   fire sizes during the selection process and should be used for visualization.
#'   If `FALSE`, the original scale is used.
#' @param bins Numeric vector (assumed to be available in the calling environment).
#'   The breakpoints for the histogram bins. This should be consistent with the
#'   `bins` output from `build_target_hist`.
#' @param target_hist Numeric vector (assumed to be available in the calling environment).
#'   The density values of the target histogram, typically obtained from `build_target_hist`.
#' @return This function does not return a value. It produces a plot as a side effect
#'   and prints information to the console.
#' @importFrom graphics hist lines legend par
#' @export
#' @examples
#' \dontrun{
#' # This example requires the 'foreach' package and a parallel backend (e.g., 'doParallel')
#' # to be set up.
#' # library(doParallel)
#'
#' # Dummy data for demonstration (replace with your actual data)
#' set.seed(123)
#' historical_data_for_target <- floor(fit_powerlaw(n = 500, alpha = 2, xmin = 10, xmax = 10000))
#' event_surfaces <- fit_powerlaw(n = 10000, alpha = 2, xmin = 10, xmax = 10000)
#' # Discard simulated fires that are too large (below 110% max historical size)
#' event_surfaces <- event_surfaces[event_surfaces<max(historical_data_for_target)*1.1]
#' event_probabilities <- rnorm(length(event_surfaces))
#' event_probabilities <- (event_probabilities-min(event_probabilities))/(max(event_probabilities)-min(event_probabilities))
#'
#' y <- 100 #number of years spanning historical fire data
#' check_fire_data(fires_hist_size = historical_data_for_target,
#'                 sim_perimeters_size = event_surfaces,
#'                 n_years = y)
#'
#' reference_surface_example <- sum(historical_data_for_target)/y
#' surface_threshold_example <- check_fire_data(fires_hist_size = historical_data_for_target,
#'                                              sim_perimeters_size = event_surfaces,
#'                                              n_years = 10)
#' tolerance_example <- 0.1
#'
#' # Create a dummy target histogram (assuming 'event_surfaces' from historical data)
#' # For a real scenario, 'event_surfaces' here would be your historical fire sizes.
#'
#' target_info_example <- build_target_hist(num_bins = 10, logaritmic = TRUE,
#'                                          sizes = historical_data_for_target,
#'                                          event_surfaces = event_surfaces)
#' target_hist <- target_info_example$target_hist
#' bins <- target_info_example$bins
#'
#' # Run the selection process
#' selected_events_result <- select_events(
#'   event_sizes = event_surfaces,
#'   event_probabilities = event_probabilities,
#'   target_hist = target_hist,
#'   bins = bins,
#'   reference_surface = reference_surface_example,
#'   surface_threshold = surface_threshold_example,
#'   tolerance = tolerance_example,
#'   max_it = 2 # Reduced iterations for example
#' )
#'
#' visualize_selected_dist(result = selected_events_result,
#'   logaritmic = T,
#'   target_hist = target_hist,
#'   bins = bins)
#'
#' # Stop the parallel cluster when done
#' # stopImplicitCluster()
#'
#' visualize_selected_dist(result = selected_events_result)
#' }
visualize_selected_dist <- function(result = result, logaritmic = TRUE, target_hist = target_hist, bins = bins) {

  # Calculate bin midpoints for the target histogram
  bin_mids <- bins[-length(bins)] + (bins[2] - bins[1]) / 2

  # Check if there are results to visualize
  if (!is.null(result$selected_surfaces)) {
    selected_surfaces <- result$selected_surfaces
    total_surface_selected <- result$total_surface
    final_discrepancy <- result$final_discrepancy

    # Display statistics about the selected events
    cat("Number of selected events:", length(selected_surfaces), "\n")
    cat("Total surface:", total_surface_selected, "\n")
    cat("Discrepancy:", final_discrepancy, "\n")

    # Prepare data for ggplot2
    # Create a dataframe for selected surfaces
    df_selected <- data.frame(
      surface = if(logaritmic) log(selected_surfaces + 1e-6) else selected_surfaces
    )

    # Create a dataframe for the target density
    df_target <- data.frame(
      x = bin_mids,
      density = target_hist
    )

    # Build the plot with ggplot2
    p <- ggplot(df_selected, aes(x = surface)) +
      # Layer for the histogram of selected surfaces
      geom_histogram(aes(y = after_stat(density), fill = "Selected"),
                     breaks = bins,
                     color = "steelblue4", # White border between bars
                     alpha = 0.8) +
      # Layer for the target distribution line
      geom_line(data = df_target, aes(x = x, y = density, color = "Target"),
                size = 1.2) + # 'size' in ggplot2 is 'lwd' in base R
      # Define colors for legend labels
      scale_fill_manual(name = "Distribution",
                        values = c("Selected" = "#799fbf"),
                        labels = c("Selected" = "Selected")) + # Label for fill legend
      scale_color_manual(name = "Distribution",
                         values = c("Target" = "#a65455"),
                         labels = c("Target" = "Target")) + # Label for color legend
      # Titles and labels
      labs(title = if(logaritmic) "Selected vs. Target Distribution (log-transformed)" else "Selected vs. Target Distribution",
           x = if(logaritmic) "Fire Size (log)" else "Fire Size",
           y = "Density") +
      # Visual theme (optional, can be customized)
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"), # Center and bold title
        legend.position = "bottom", # Legend position
        legend.title = element_blank() # Remove legend title if not needed
      )

    # Display the plot
    print(p)
    return(p)

  } else {
    cat("Could not find a solution for the desired distribution.\n")
  }
}

#' Check fire data for sufficiency
#'
#' This function assesses whether the available simulated fire perimeters are
#' sufficient in terms of maximum fire size and total burned area, relative
#' to historical fire data. It provides messages to the user if deficiencies
#' are found and, if data is sufficient, calculates and recommends a surface
#' threshold for selecting events.
#'
#' @param fires_hist_size Numeric vector. A vector of historical fire sizes (e.g., areas in hectares).
#' @param sim_perimeters_size Numeric vector. A vector of simulated fire perimeter sizes (e.g., areas in hectares).
#' @param n_years Integer. The number of years represented by the historical fire data.
#' @return If simulated data is sufficient, returns an integer representing a
#'   recommended surface threshold (typically 10% of the maximum possible threshold).
#'   Otherwise, returns `NULL` (invisibly, as it only prints messages) if
#'   simulated data is deemed insufficient.
#' @seealso \code{\link{get_select_params}}
#' @export
#' @examples
#' # Example 1: Simulated fires are too small and not enough area
#' check_fire_data(fires_hist_size = c(100, 500, 1000, 2000),
#'                 sim_perimeters_size = c(10, 50, 80, 120),
#'                 n_years = 1)
#'
#' # Example 2: Simulated fires are sufficient
#' check_fire_data(fires_hist_size = c(100, 500, 1000, 2000),
#'                 sim_perimeters_size = c(10, 50, 80, 120, 500, 1500, 2500, 5000),
#'                 n_years = 1)
#'
#' # Example 3: Insufficient area, but max size is okay
#' check_fire_data(fires_hist_size = c(100, 500, 1000, 2000),
#'                 sim_perimeters_size = c(10, 50, 80, 120, 500, 1500),
#'                 n_years = 10)
#' # Example 4: Correct data
#'check_fire_data(fires_hist_size = c(1, 5, 10, 600),
#'                sim_perimeters_size = c(10, 50, 80, 120, 500, 2500),
#'                n_years = 10)
check_fire_data <- function(fires_hist_size, sim_perimeters_size, n_years) {

  size_check <- max(fires_hist_size) > max(sim_perimeters_size)
  area_check <- sum(fires_hist_size) > sum(sim_perimeters_size)

  if(size_check) {print("Simulated fires are too small. Consider running extra simulations or trim historical fires to match.")}
  if(area_check) {print("Not enough simulated fires. Extra simulations are needed.")}

  if(!size_check & !area_check){

    st <- get_select_params(fires_hist_size, sim_perimeters_size, n_years)
    num_seasons <- floor(st*0.1)

    if(num_seasons>0){

      cat("Sufficient simulated perimeters and burned area. Maximum surface threshold: ",
          st,".\n Recommended surface threshold: ",floor(st*0.1), "\n")

      return(num_seasons)

    }else{"Not enough simulated fires. Extra simulations are needed."}

  }

}

#' Get Selection Parameters
#'
#' Calculates a parameter for event selection based on the total area of
#' historical fires (normalized by number of years) and the total area of
#' simulated fire perimeters. This parameter can be used to determine a
#' maximum surface threshold for selecting a subset of simulated events.
#'
#' @param fires_hist_size Numeric vector. A vector of historical fire sizes.
#' @param sim_perimeters_size Numeric vector. A vector of simulated fire perimeter sizes.
#' @param n_years Integer. The number of years covered by the historical fire data.
#' @return An integer representing the ratio of total simulated area to the
#'   annual average historical area. This value can be interpreted as the
#'   maximum number of historical fire "equivalents" that can be formed from
#'   the simulated data.
#' @seealso \code{\link{check_fire_data}}
#' @export
#' @examples
#' # Calculate selection parameter for a scenario where simulated data
#' # covers 10 times the annual historical area.
#' get_select_params(fires_hist_size = c(100, 500, 1000),
#'                   sim_perimeters_size = c(1000, 2000, 3000, 4000, 5000),
#'                   n_years = 1)
#'
#' # If historical data represents 10 years
#' get_select_params(fires_hist_size = c(100, 500, 1000),
#'                   sim_perimeters_size = c(1000, 2000, 3000, 4000, 5000),
#'                   n_years = 10)
get_select_params <- function(fires_hist_size, sim_perimeters_size, n_years) {

  area_hist <- sum(fires_hist_size) / n_years
  area_sim <- sum(sim_perimeters_size)
  return(floor(area_sim/area_hist))

}
