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
#' target_info_ex <- build_target_hist(event_surfaces = event_surfaces,
#'                                     num_bins = 5,
#'                                     logaritmic = TRUE,
#'                                     sizes = historical_sizes_ex)
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
#' target_info_linear_ex <- build_target_hist(event_surfaces = event_surfaces,
#'                                            num_bins = 5,
#'                                            logaritmic = FALSE,
#'                                            sizes = historical_sizes_ex)
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
#' event_probabilities <- (event_probabilities-min(event_probabilities))/
#'                         (max(event_probabilities)-min(event_probabilities))
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
                          reference_surface, surface_threshold, tolerance, max_it = 5,
                          iter_limit=100000,logaritmic = T) {

  num_cores <- max_it  # N�mero de n�cleos a usar (idealmente todos menos uno)
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)

  event_surfaces <- event_sizes

  # Exporta variables y funciones necesarias al cl�ster paralelo
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
              n_iter <- 0
              status <- "Completed correctly"

              for (i in 1:max_it) {

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
                      status <- "Interrupted early"
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
                          status <- "Interrupted early"
                          cat("Maximum iterations where reached. Execution interrumpted. Increase iter_limit.")
                          break
                        }

                        if (length(eligible_indices) == 1) {
                          selected_index <- eligible_indices[1] 	# Direct selection if only one eligible index
                        } else if (length(eligible_indices) > 1) {
                          selected_index <- sample(eligible_indices, 1, prob = event_probabilities[eligible_indices])
                          n_iter <- n_iter + 1
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

                        # n_iter <- n_iter + 1
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
                      final_discrepancy = final_discrepancy,
                      run_status = status
                    )
                    best_discrepancy <- final_discrepancy
                  }

                  best_selection

                }, error = function(e) {
                  cat("Error in iteration", i, ":", conditionMessage(e), "\n")

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

#' Visualize Selected vs. Target Fire Size Distribution
#'
#' Esta función genera un gráfico comparando la distribución de los
#' tamaños de incendio seleccionados con una distribución objetivo.
#' Permite aplicar transformación logarítmica a los tamaños seleccionados
#' y visualizar ambas distribuciones de manera conjunta.
#'
#' @param result Lista devuelta por un modelo o proceso de selección
#'   que debe contener al menos:
#'   \describe{
#'     \item{selected_surfaces}{Vector numérico con los tamaños de superficie seleccionados.}
#'     \item{total_surface}{Valor numérico con la superficie total seleccionada.}
#'     \item{final_discrepancy}{Valor numérico con la discrepancia final respecto a la distribución objetivo.}
#'   }
#' @param logaritmic Lógico. Si \code{TRUE}, los tamaños seleccionados se transforman logarítmicamente
#'   antes de graficar. Por defecto \code{TRUE}.
#' @param target_hist Vector numérico con las densidades de la distribución objetivo,
#'   correspondiente a los puntos medios de los intervalos definidos en \code{bins}.
#' @param bins Vector numérico con los límites de los intervalos (breaks) usados
#'   para construir el histograma.
#'
#' @return Un objeto \code{ggplot} que representa la comparación entre la distribución
#' seleccionada y la distribución objetivo. Además, la función imprime en consola:
#' \itemize{
#'   \item Número de eventos seleccionados.
#'   \item Superficie total seleccionada.
#'   \item Discrepancia respecto a la distribución objetivo.
#' }
#' Si no se encuentra una solución, devuelve un mensaje y no genera gráfico.
#'
#' @examples
#' \dontrun{
#' # Use example for select_events
#' visualize_selected_dist(result, logaritmic = TRUE,
#'                         target_hist = target_hist, bins = bins)
#' }
#'
#' @import ggplot2
#' @export
visualize_selected_dist <- function(result = result, logaritmic = TRUE,
                                    target_hist = target_hist, bins = bins) {

  bin_mids <- bins[-length(bins)] + (bins[2] - bins[1])/2

  if (!is.null(result$selected_surfaces)) {
    selected_surfaces <- result$selected_surfaces
    total_surface_selected <- result$total_surface
    final_discrepancy <- result$final_discrepancy
    cat("Number of selected events:", length(selected_surfaces),
        "\n")
    cat("Total surface:", total_surface_selected, "\n")
    cat("Discrepancy:", final_discrepancy, "\n")
    df_selected <- data.frame(surface = if (logaritmic)
      log(selected_surfaces + 1e-06)
      else selected_surfaces)
    df_target <- data.frame(x = bin_mids, density = target_hist)
    p <- ggplot(df_selected, aes(x = surface)) +
      geom_histogram(aes(y = after_stat(density),
                         fill = "Selected"), breaks = bins, color = "steelblue4", alpha = 0.8) +
      geom_line(data = df_target, aes(x = x, y = density, color = "Target"), size = 1.2) +
      scale_fill_manual(name = "Distribution", values = c(Selected = "#799fbf"),
                        labels = c(Selected = "Selected")) +
      scale_color_manual(name = "Distribution",
                         values = c(Target = "#a65455"),
                         labels = c(Target = "Target")) +
      labs(title = if (logaritmic)
        "Selected vs. Target Distribution (log-transformed)"
        else "Selected vs. Target Distribution", x = if (logaritmic)
          "Fire Size (log)"
        else "Fire Size", y = "Density") + theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"),
            legend.position = "bottom", legend.title = element_blank())
    print(p)
    return(p)
  }
  else {
    cat("Could not find a solution for the desired distribution.\n")
  }
}

#' Get Selection Parameters for Fire Simulations
#'
#' This function calculates a selection parameter that represents how many times
#' the simulated burned area covers, on average, the historical annual burned area.
#'
#' @param fires_hist_size Numeric vector of historical fire sizes.
#' @param sim_perimeters_size Numeric vector of simulated fire sizes.
#' @param n_years Integer. Number of years considered in the historical dataset.
#'
#' @return An integer value representing the ratio between the total simulated area
#' and the average annual historical burned area.
#' Larger values indicate that the simulation covers several historical fire seasons.
#'
#' @details
#' The parameter is calculated as:
#' \deqn{ \text{param} = \left\lfloor \frac{\sum(\text{sim\_perimeters\_size})}{\sum(\text{fires\_hist\_size}) / n\_years} \right\rfloor }
#'
#' Where:
#' \itemize{
#'   \item \eqn{\sum(\text{fires\_hist\_size}) / n\_years} = average historical burned area per year.
#'   \item \eqn{\sum(\text{sim\_perimeters\_size})} = total simulated burned area.
#' }
#'
#' @examples
#' \dontrun{
#' # Example with toy data
#' hist_sizes <- c(100, 200, 150, 300)
#' sim_sizes  <- c(80, 120, 200, 250, 300)
#'
#' get_select_params(fires_hist_size = hist_sizes,
#'                   sim_perimeters_size = sim_sizes,
#'                   n_years = 4)
#' }
#'
#' @seealso [check_fire_data]
#'
#' @export
get_select_params <- function(fires_hist_size, sim_perimeters_size, n_years) {
  area_hist <- sum(fires_hist_size) / n_years
  area_sim <- sum(sim_perimeters_size)
  return(floor(area_sim / area_hist))
}


#' Check Consistency of Historical vs. Simulated Fire Data
#'
#' This function compares historical fire sizes with simulated fire sizes
#' to assess whether the simulations are sufficient in terms of
#' maximum burned area and total burned area.
#' If the simulated fires are not adequate, the function prints guidance
#' messages suggesting additional simulations or adjustments.
#'
#' @param fires_hist_size Numeric vector of historical fire sizes.
#' @param sim_perimeters_size Numeric vector of simulated fire sizes.
#' @param n_years Integer. Number of years represented in the simulation.
#'
#' @return
#' If the simulated fires are sufficient, returns an integer representing
#' the recommended surface threshold (in number of seasons).
#' If the simulations are not sufficient, prints diagnostic messages
#' and returns nothing.
#'
#' @details
#' The function checks:
#' \itemize{
#'   \item Maximum fire size in historical vs. simulated data.
#'   \item Total burned area in historical vs. simulated data.
#' }
#'
#' If both criteria are satisfied, it calculates a recommended threshold
#' using \code{\link{get_select_params}}.
#'
#' @examples
#' \dontrun{
#' # Example with toy data
#' set.seed(123)
#' hist_sizes <- rgamma(100, shape = 2, scale = 5)
#' sim_sizes  <- rgamma(200, shape = 2, scale = 4)
#'
#' check_fire_data(fires_hist_size = hist_sizes,
#'                 sim_perimeters_size = sim_sizes,
#'                 n_years = 10)
#' }
#'
#' @seealso \code{\link{get_select_params}}
#'
#' @export
check_fire_data <- function(fires_hist_size, sim_perimeters_size, n_years) {
  size_check <- max(fires_hist_size) > max(sim_perimeters_size)
  area_check <- sum(fires_hist_size) > sum(sim_perimeters_size)

  if (size_check) {
    print("Simulated fires are too small. Consider running extra simulations or trim historical fires to match.")
  }
  if (area_check) {
    print("Not enough simulated fires. Extra simulations are needed.")
  }
  if (!size_check & !area_check) {
    st <- get_select_params(fires_hist_size, sim_perimeters_size, n_years)
    num_seasons <- floor(st * 0.25)

    if (num_seasons > 0) {
      cat("Sufficient simulated perimeters and burned area. Maximum surface threshold: ",
          st, ".\n Recommended surface threshold: ",
          floor(st * 0.1), "\n")
      return(num_seasons)
    } else {
      "Not enough simulated fires. Extra simulations are needed."
    }
  }
}

#' Remove Duplicate Polygons from Candidate Set
#'
#' This function checks a set of candidate polygons and removes duplicates
#' based on spatial equality. It iteratively removes duplicated polygons
#' until none remain.
#'
#' @param candidates An \code{sf} object containing candidate polygons.
#'
#' @return An \code{sf} object with duplicate polygons removed.
#' The function also prints messages to the console indicating whether
#' duplicates were found and their indices.
#'
#' @details
#' Duplicates are detected using \code{sf::st_equals()}, which checks
#' for geometric equality between polygons.
#' The function continues removing duplicates until no further matches are found.
#'
#' @examples
#' \dontrun{
#' library(sf)
#'
#' # Create a simple sf object with duplicate polygons
#' poly1 <- st_polygon(list(rbind(c(0,0), c(1,0), c(1,1), c(0,1), c(0,0))))
#' poly2 <- st_polygon(list(rbind(c(0,0), c(1,0), c(1,1), c(0,1), c(0,0)))) # duplicate
#' poly3 <- st_polygon(list(rbind(c(2,2), c(3,2), c(3,3), c(2,3), c(2,2))))
#'
#' candidates <- st_sf(geometry = st_sfc(poly1, poly2, poly3))
#'
#' cleansed <- cleanse_duplicates(candidates)
#' }
#'
#' @seealso \code{\link[sf]{st_equals}}
#'
#' @import sf
#' @export
cleanse_duplicates <- function(candidates) {
  candidate_surfaces <- candidates
  duplicate_indices <- list(-1)

  while (length(duplicate_indices) > 0) {
    duplicates <- st_equals(candidate_surfaces)
    duplicate_indices <- which(sapply(duplicates, length) > 1)

    if (length(duplicate_indices) > 0) {
      print("Duplicate polygons found at the following indices:\n")
      print(duplicate_indices)
      candidate_surfaces <- candidate_surfaces[-duplicate_indices, ]
    } else {
      print("No duplicate polygons found.\n")
    }
  }
  return(candidate_surfaces)
}

