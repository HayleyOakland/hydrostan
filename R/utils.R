# Data Processing Functions ####

#' Query data{} block requirements for a given stan file
#' 
#' @param stan_file_name Single item character vector that defines the file path + file name + .stan extension of the stan file to query data structure for
#' @param file_path Default is NULL, which applies when the stan file is saved in the working directory. Otherwise, specify the single-item character vector defining the file path that contains the stan file
#' @export
get_stan_data_reqs <- function(stan_file_name, file_path=NULL) {
  # Read the stan file
  stan_code <- readLines(paste0(file_path, stan_file_name, ".stan"))
  
  # Find the data block
  data_block_start <- which(grepl("^\\s*data\\s*\\{", stan_code))
  data_block_end <- which(grepl("^\\s*\\}", stan_code))
  data_block_end <- data_block_end[data_block_end > data_block_start][1]
  
  # Extract data block lines
  data_lines <- stan_code[(data_block_start+1):(data_block_end-1)]
  
  # Remove comments
  data_lines <- gsub("//.*$", "", data_lines)
  data_lines <- trimws(data_lines)
  data_lines <- data_lines[nzchar(data_lines)]
  
  # Parse variable declarations
  result <- data.frame(variable = character(), type = character(), 
                       lower = character(), upper = character(),
                       dimensions = character(), stringsAsFactors = FALSE)
  
  for (line in data_lines) {
    if (!grepl(";", line)) next # Skip if not a variable declaration
    
    # Split off the comment if any
    var_decl <- trimws(gsub(";.*$", "", line))
    
    # Extract variable type and name with dimensions
    pattern <- "^([^\\s]+)\\s+(.+)$"
    full_type <- gsub(pattern, "\\1", var_decl)
    var_with_dim <- trimws(gsub(pattern, "\\2", var_decl))
    
    # Extract constraints from the type
    lower_bound <- NA_character_
    upper_bound <- NA_character_
    
    if (grepl("<.*>", full_type)) {
      constraints <- regmatches(full_type, regexpr("<.*>", full_type))
      constraints <- gsub("<|>", "", constraints)
      
      # Parse lower bound
      if (grepl("lower\\s*=", constraints)) {
        lower_match <- regmatches(constraints, regexpr("lower\\s*=\\s*[^,>]+", constraints))
        if (length(lower_match) > 0) {
          lower_bound <- trimws(gsub("lower\\s*=\\s*", "", lower_match))
        }
      }
      
      # Parse upper bound
      if (grepl("upper\\s*=", constraints)) {
        upper_match <- regmatches(constraints, regexpr("upper\\s*=\\s*[^,>]+", constraints))
        if (length(upper_match) > 0) {
          upper_bound <- trimws(gsub("upper\\s*=\\s*", "", upper_match))
        }
      }
    }
    
    # Clean up the type - remove constraints and extract base type
    base_type <- gsub("<.*?>", "", full_type)
    
    # Check if dimensions are specified after the type (like vector[N])
    type_has_dims <- grepl("\\[", full_type)
    var_has_dims <- grepl("\\[", var_with_dim)
    
    if (type_has_dims) {
      # Extract dimension from type (e.g., vector[N])
      bracket_start <- regexpr("\\[", full_type)
      bracket_end <- regexpr("\\]", full_type)
      if (bracket_start > 0 && bracket_end > bracket_start) {
        dim_desc <- substr(full_type, bracket_start + 1, bracket_end - 1)
      } else {
        dim_desc <- "1"
      }
      var_name <- var_with_dim
      
      # For types with dimensions, they're still scalar variables
      type_desc <- paste(base_type, "scalar", sep = ", ")
      
    } else if (var_has_dims) {
      # Extract dimension from variable name (e.g., real var[N])
      var_name <- gsub("([^\\[]*)\\[.*", "\\1", var_with_dim)
      
      # Extract ALL dimensions, including nested brackets
      # Find the first '[' and then find the matching ']' (accounting for nested brackets)
      bracket_start <- regexpr("\\[", var_with_dim)
      if (bracket_start > 0) {
        # Extract everything after the first '['
        after_bracket <- substr(var_with_dim, bracket_start + 1, nchar(var_with_dim))
        
        # Find the matching closing bracket by counting brackets
        bracket_count <- 0
        end_pos <- 0
        for (i in 1:nchar(after_bracket)) {
          char <- substr(after_bracket, i, i)
          if (char == "[") {
            bracket_count <- bracket_count + 1
          } else if (char == "]") {
            if (bracket_count == 0) {
              end_pos <- i
              break
            } else {
              bracket_count <- bracket_count - 1
            }
          }
        }
        
        if (end_pos > 0) {
          dim_desc <- substr(after_bracket, 1, end_pos - 1)
        } else {
          dim_desc <- "1"
        }
      } else {
        dim_desc <- "1"
      }
      
      # Determine array type
      if (base_type %in% c("vector", "simplex", "ordered", "positive_ordered")) {
        type_desc <- paste(base_type, "vector", sep = ", ")
      } else {
        type_desc <- paste(base_type, "array", sep = ", ")
      }
      
    } else {
      # No dimensions specified anywhere
      var_name <- var_with_dim
      
      if (base_type %in% c("vector", "simplex", "ordered", "positive_ordered")) {
        type_desc <- paste(base_type, "vector", sep = ", ")
        dim_desc <- "1"
      } else {
        type_desc <- paste(base_type, "scalar", sep = ", ")
        dim_desc <- "1"
      }
    }
    
    # Ensure all variables have values
    if (is.null(var_name) || length(var_name) == 0 || var_name == "") {
      next  # Skip this iteration if variable name is empty
    }
    
    result <- rbind(result, data.frame(
      variable = trimws(var_name),
      type = type_desc,
      lower = lower_bound,
      upper = upper_bound,
      dimensions = dim_desc,
      stringsAsFactors = FALSE
    ))
  }
  
  return(result)
}

#' Make list of data required for stan run
#' 
#' @param concData Data frame containing concentration by time data per trial (trialIdx x C x time). 
#' @param Cvals Data frame containing starting and ending concentration values per trial.
#' @param flumeInfo Data frame containing information about the flumes used in the trials, including `flumeIdx`, `phase` (1=beginnning of experiment, 2=end of experiment), and `finalCaddisN` (the final abundance of caddisflies in each flume)
#' @param avgKprimeVals Data frame containing avgKprime values per trial.
#' @param fracExTau Vector of length 1 containing numeric value of time to find the fraction of the hyporheic zone that has exchanged by this given time (e.g., "30", if concentration DF times are in minutes, will find the fraction of the hyporheic zone exchanged by 30 minutes); equivalent to I(tau)
#' @param model The model to run in stan. Options include "wellMixedModel" for a basic exponential decay model (assuming a well-mixed hyporheic zone), "powerLawRTD" for the power law functional form of an RTD model of the hyporheic zone, or "exponentialRTD" for the exponential function form of the RTD model of the hyporheic zone
#' @param resultsMatrix Needed if model != wellMixedModel. Matrix of concentrations from the hydrogeom model for interpolation, relative to values of the curvature parameter (alpha for PLRTD or sigma for ExpRTD), volume ratio, and time
#' @param include_drift Logical statement, whether to include drift as a parameter in the stan model. 
#' @param include_qChange Logical statement, whether to include a change in q across trials (e.g., if there is more than one release per flume, are we testing the statistical model that q_down changed as a function of caddisfly density)
#' @export
makeStanDataList <- function(concData, 
                             Cvals, 
                             flumeInfo, 
                             avgKprimeVals, 
                             fracExTau,
                             model, 
                             resultsMatrix,
                             include_drift,
                             include_qChange) {
  # data required regardless of model
  N = nrow(concData)
  Trl = length(unique(concData$trialIdx))
  if(any(grepl("newTrialIdx", names(concData)))) {
    trialIdx = concData$newTrialIdx
  }
  else {trialIdx = concData$trialIdx}
  time = concData$time
  conc = concData$C
  C_pre = Cvals$C_pre          #C_pre = labFlumeData |> distinct(flumeDeploymentIdx,C_star) |> select(C_star) |> as_vector() |> unname()/
  C_max_t = Cvals$C_max_t      #C_max_t = labFlumeData |> distinct(flumeDeploymentIdx,C_max_t) |> select(C_max_t) |> as_vector() |> unname(),
  max_t = Cvals$max_t          #max_t = labFlumeData |> distinct(flumeDeploymentIdx,timeAtC_max_t) |> select(timeAtC_max_t) |> as_vector() |> unname(),
  avgKprime = avgKprimeVals
  fracExTau = fracExTau
  
  if(include_qChange == T) {
    calc_qChange = 1
    Flm = length(unique(flumeInfo$flumeIdx))
    flumeIdx = flumeInfo$flumeIdx
    phase = flumeInfo$phase
    density = flumeInfo |>
      dplyr::distinct(flumeIdx, finalCaddisN) |>
      dplyr::arrange(flumeIdx) |>
      dplyr::mutate(finalCaddisDens_sqm = finalCaddisN / 0.13) |>
      dplyr::pull(finalCaddisDens_sqm) #convert final caddis N to final caddis density by dividing by surface area of flume stream bed (0.13m^2; to get n per m^2)
  }
  else {
    calc_qChange = 0
    Flm = 1
    flumeIdx = rep(1, times=Trl)
    phase = rep(1, times=Trl) 
    density = array(1, dim=1) #make sure it's a vector, not scalar
  }

  if(include_drift == T) {
    calc_drift = 1
  }
  else {
    calc_drift = 0
  }
  
  #model-specific data
  if(model == "powerLawRTD") {
    use_hydrogeom_model = 1
    is_powerLaw = 1
    N_pars = c(length(unique(as.numeric(dimnames(resultsMatrix)$alpha))), length(unique(as.numeric(dimnames(resultsMatrix)$V_h))), length(unique(as.numeric(dimnames(resultsMatrix)$time))))
    curvParVals = unique(as.numeric(dimnames(resultsMatrix)$alpha))
    V_ratios = unique(as.numeric(dimnames(resultsMatrix)$V_h))
    times = unique(as.numeric(dimnames(resultsMatrix)$time))
    m3d = resultsMatrix
  } 
  else if (model == "exponentialRTD") {
    use_hydrogeom_model = 1
    is_powerLaw = 0
    N_pars = c(length(unique(as.numeric(dimnames(resultsMatrix)$sigma))), length(unique(as.numeric(dimnames(resultsMatrix)$V_h))), length(unique(as.numeric(dimnames(resultsMatrix)$time))))
    curvParVals = unique(as.numeric(dimnames(resultsMatrix)$sigma))
    V_ratios = unique(as.numeric(dimnames(resultsMatrix)$V_h))
    times = unique(as.numeric(dimnames(resultsMatrix)$time))
    m3d = resultsMatrix
  }
  else if (model == "wellMixedModel") {
    use_hydrogeom_model = 0
    is_powerLaw = 0
    N_pars = c(1,1,1)
    curvParVals = 0
    V_ratios = 0
    times = 0
    m3d = 0
  }
  stanDataList <- list(
    N = N,
    Trl = Trl,
    Flm = Flm,
    trialIdx = trialIdx,
    flumeIdx = flumeIdx,
    phase = phase,
    time = time, 
    conc = conc,
    C_pre = C_pre,
    C_max_t = C_max_t,
    max_t = max_t,
    fracExTau = fracExTau,
    density = density,
    avgKprime = avgKprime,
    use_hydrogeom_model = use_hydrogeom_model,
    is_powerLaw = is_powerLaw,
    calc_drift = calc_drift,
    calc_qChange = calc_qChange,
    N_pars = N_pars,
    curvParVals = curvParVals,
    V_ratios = V_ratios,
    times = times,
    m3d = m3d
  )
  return(stanDataList)
}


#' Make 3D matrix of channel concentration values.
#' 
#' @param result_tbl A tibble of concentration values relative to the three dimensions of the eventual matrix (created using hydrogeom and flumeTracer)
#' @param target A character vector of length 1, default is "C_c" (channel concentration) 
#' 
#' @description
#' Make 3D concentration matrix relative to the maximum residence time (tau_n), 
#' the ratio of volume between the hyporheic zone : total (V_ratio) 
#' and the curvature parameter in the residence time distribution 
#' (alpha for power law or sigma for exponential functional forms)
#' 
#' @return 3D concentration matrix with dimensionless concentration values 
#' (decaying from 1 to 0, e.g., scaled relative to initial concentration and 
#' final concentration); and relative to dimensionless time values 
#' (e.g., increasing from 0 to 1 as a fraction of maximum residence time, 
#' time 0 is when "C_c" = 1.0, time 1.0 is when "C_c" = 0)
#' @rdname make_matrix
#' @export
make_matrix <- function(results_tbl, target = "C_c") {
  #extract the names of the axes of the tibble (aside from the target)
  axes_names <- names(results_tbl)[names(results_tbl) != target]
  #extract the range of values for each of the axes
  axis_values <- sapply(axes_names, \(x) sort(unique(results_tbl[[x]])), simplify = F)
  #create unique index for the axis values
  axis_idx <- sapply(axis_values, \(x) 1:length(x), simplify = F)
  #make tibbles of the axis values / indices
  axis_tbls <- Map(\(i,val) dplyr::tibble(idx = i, value = val), axis_idx, axis_values)
  #update names of axis tibbles
  axis_tbls <- Map(\(nm, tbl) {
    names(tbl)[1] <- paste0(nm, "_idx")
    names(tbl)[2] <- nm
    tbl
  }, 
  axes_names,
  axis_tbls)
  # join axis tibbles
  joined_result_tbl <- results_tbl
  for(join_tbl in axis_tbls) {
    joined_result_tbl <- inner_join(joined_result_tbl, join_tbl)
  }
  #query variables of axes
  vars <- paste0(rev(axes_names), "_idx")
  #arrange values 
  joined_result_tbl <- dplyr::arrange(joined_result_tbl, across(all_of(vars)))
  #create 3D matrix of values
  array(
    data = as.vector(unlist(joined_result_tbl[,target])),
    dim = sapply(axis_tbls, nrow),
    dimnames = lapply(axis_values, as.character))
}

#' Dedensify data using a log10 spacing
#'
#' Sample data using evenly spaced orders of magnitude.
#'
#' @param x vector of data values in x dimension or a data.frame
#' @param y vector of data values in y dimension (equal in length to x)
#' @param n number of samples to extract
#' @param x_offset value added to x before log10 is applied; large values
#'   provide more evenly spaced points.  Small values crowd points toward the
#'   left (smaller values of x).
#' @param x_col,y_col name or number of data.frame column containing x dimension
#'   and y dimension
#' @param approx logical (default=T) argument. If TRUE, dedensify() uses approx() 
#'    to find y values associated with de-densified x values;
#'    If FALSE, dedensify finds indices of original data x values closest to the de-densified
#'    x values, then uses those indices to filter y values (from the original data, 
#'    e.g., instead of approximating y values)
#' @param ... other values passed on to dedensify.
#' @return a tibble of n {x,y} pairs
#' @export
#' @rdname dedensify
#' @export
dedensify <- function(x, y, n, x_offset = 1, decimals = 10, approx=T) {

  #make sure x and y are ordered
  y = y[order(x)]
  x = sort(x)
  
  #add offset (see documentation - controls spacing of de-densified pts)
  x_off = x + x_offset
  
  #find log of offset x values
  log_x = log10(x_off)
  #de-densified x is n-values from min to max of logged values
  log_deden_x = seq(min(log_x), max(log_x), length.out = n)
  #exponentiate back to original x scale
  deden_x = round(10^log_deden_x - x_offset, decimals)
  #find y-values associated with de-densified x values using approx, and return tibble
  if(approx==T){
    approx(x, y, deden_x) |>
    dplyr::as_tibble() %>% dplyr::rename(time=x, C=y)
  }
  else {
    closest_indices <- sapply(deden_x, function(d) {which.min(abs(x-d))})
    x = x[closest_indices]
    y = y[closest_indices]
    dplyr::bind_cols(time=x,C=y)
  }
}

#' @rdname dedensify
#' @export
dedensify_df <- function(x, x_col, y_col, n, ...) {
  #if x is a dataframe (instead of individual vectors for x and y), use the column names provided to run dedensify()
  dedensify(x[[x_col]], x[[y_col]], n, ...)
}

# Post-processing stanfit functions ####

#' Extract posterior distributions for parameter estimates
#' @param posterior is the posterior distributions from rstan::extract() of the stanFit object
#' @param param is the desired parameter to extract (e.g., "tau_n", "predicted_concentration"). If param = "predicted_concentration", dataDF cannot be NA
#' @param dataDF (only required if `param` = "predicted_concentration"; default=NA) is the original dataset dataframe fit by stan (e.g., nrow(dataDF) must be equal to the N in the list of data sent to the stan model), which must include columns `trialIdx` identifying the trials (with values from 1:trialN), `time` and `C` (concentration). The default is NA, as this is not required when param != "predicted_concentration"
#' @export
paramPostDists <- function(posterior, param, dataDF=NA){
  postDistDF <- as.data.frame(posterior[[param]])
  if(param == "predicted_concentration") {
    postDistDF <- postDistDF |>
      dplyr::mutate(iter=1:n()) |>
      tidyr::pivot_longer(cols = dplyr::contains("V"), 
                   names_to = "obs_idx", 
                   names_prefix = "V",
                   values_to = "pred_conc") |>
      dplyr::mutate(obs_idx = as.numeric(obs_idx)) |>
      dplyr::left_join(data.frame(obs_idx = 1:nrow(dataDF), 
                           trial = dataDF$trialIdx, 
                           time = dataDF$time, 
                           conc = dataDF$C), 
                by = "obs_idx")
  }
  else{
    trialN = ncol(postDistDF)
    postDistDF <- postDistDF |>
      `names<-`(paste0("trial", 1:trialN)) |>
      dplyr::mutate(iter = 1:dplyr::n()) |>
      tidyr::pivot_longer( 
        cols = dplyr::contains("trial"),
        names_to = "trial", 
        names_prefix = "trial",
        values_to = paste0(param, "_estimate")
      ) 
  }
  return(postDistDF)
}

#' Summary of predicted concentration lines across iterations. 
#' @description
#' Returns list with two dataframes, `predSummary` and `iterRMSE`. `predSummary` contains a summary of the predicted concentrations across iterations, either based on median/quantile of iteration predictions or based on a random sample of predicted concentration iterations. 
#' @param posterior Posterior distribution of stan fit (result of rstan::extract(stanFitObject)).
#' @param dataDF Original dataset dataframe fit by stan (e.g., nrow(dataDF) must be equal to the N in the list of data sent to the stan model), which must include columns `trialIdx` identifying the trials (with values from 1:trialN), `time` and `C` (concentration)
#' @param summary Desired summary type. Either "quantile" (based on quantiles of predicted concentration), or "sample" (random sample of iterations)
#' @param summaryArgs If `summary`="quantile", a two-item vector of upper and lower quantile limits to summarize predicted concentration across iterations (e.g., c(0.0275, 0.975)). If `summary`="sample", a single-item numeric vector defining the number of iterations to sample in the predicted concentration summary table.  
#' @param predType Single item character vector defining the prediction type, typically either "exponential RTD" or "power law RTD"
#' @export
predConcSummary <- function(posterior, dataDF,
                            summary,
                            summaryArgs, 
                            predType) {
  
  predConc <- paramPostDists(posterior = posterior, 
                             param = "predicted_concentration",
                             dataDF = dataDF)
  
  iterN <- nrow(posterior[["predicted_concentration"]])
  obvN <- predConc |>
    dplyr::group_by(trial) |>
    dplyr::summarise(nObv = n_distinct(obs_idx)) |>
    dplyr::ungroup() %>% 
    dplyr::distinct(nObv) |> 
    dplyr::pull()
  
  iterRMSE <- predConc |>
    dplyr::group_by(trial) |>
    dplyr::mutate(iter = rep(1:iterN, each = n_distinct(obs_idx))) |>
    dplyr::ungroup() |>
    dplyr::mutate(sqErr = (pred_conc - conc)^2) |>
    dplyr::group_by(trial, iter) |>
    dplyr::summarise(rmse = sqrt(mean(sqErr))) |>
    dplyr::ungroup() |>
    dplyr::mutate(predType = predType)
  
  if(summary == "quantile"){
    predSummary <- predConc |>
      dplyr::group_by(trial, time) %>%
      dplyr::summarize(
        pred_median = median(pred_conc),
        pred_low = quantile(pred_conc, summaryArgs[1]),
        pred_high = quantile(pred_conc, summaryArgs[2]),
        observed_conc = first(conc)
      ) %>%
      dplyr::mutate(predType = predType)
  }
  if(summary == "sample") {
    predSummary <- predConc %>%
      filter(iter %in% sample(iterN, summaryArgs, replace=F))
  }
  return(list(predSummary = predSummary, iterRMSE = iterRMSE))
}

#' Summarize parameter estimates from stan fit.
#' @description Returns a dataframe containing the mean, median sd and se of the parameter estimates per trial, summarized across iterations of the stan model fit for all parameters included in `params`.
#' @param posterior Posterior distribution of stan fit, created using rstan::extract(stanFitObject)
#' @param params Character vector of parameter names for those parameters estimates to be summarized. Options can be queried using `names(posterior)`. `params` can only contain the name(s) of one of the main parameters estimated per trial (e.g., not "sigma") or transformed parameters estimated per trial (e.g., not "predicted_concentration" or "lp__", see predConcSumm() to summarize predicted concentration).
#' @param predType Single-item character vector with the type of prediction used in the stan fit (typically either "power law RTD" or "exponential RTD") 
#' @export
estimateSummary <- function(posterior,
                            params, 
                            predType) {
  estSumm <- lapply(params, function(x) posterior[[x]] |> 
                            dplyr::as_tibble() |>
                            dplyr::mutate(iter=1:n()) |>
                            tidyr::pivot_longer(cols=dplyr::contains("V"),
                                                names_to = "trial",
                                                names_prefix = "V",
                                                values_to = "value") |>
                            dplyr::mutate(estimate=paste(x))) |>
    dplyr::bind_rows() |>
    dplyr::group_by(trial, estimate) |>
    dplyr::reframe(
      mean = mean(value, na.rm=T),
      median = median(value, na.rm=T),
      sd = sd(value, na.rm=T),
      se = sd/sqrt(n())
    ) %>%
    mutate(predType = predType)
  
  return(estSumm)
}

#' Compute fraction of hyporheic zone that has exchanged by time t
#' @param shape is a character vector of length 1, the shape of the residence time distribution (options: "exponent" or "powerLaw")
#' @param t is a numeric vector of length 1, the time of interest (same units as tau_0 and tau_n) to compute the fraction of the hyporheic zone that has exchanged with the channel
#' @param tau_0 is a numeric vector of length 1, the minimum residence time (can be 0 if shape="exponent", must be >0 if shape="powerLaw") (same time units as t and tau_n)
#' @param tau_n is a numeric vector of length 1, the maximum residence time (same time units as t and tau_0)
#' @param curvParVal is a numeric vector of length 1, the curvature parameter in the residence time distribution function (sigma if shape = "exponent", alpha if shape = "powerLaw")
#' 
#' @description
#' The fraction of the hyporheic zone that has exchanged by a given time `t`
#' can be computed as the integral from tau_0 to t of the washout function
#' divided by the integral of the washout function.
#' 
#' @return 
#' Numeric vector of length 1, the value of the fraction of the 
#' hyporheic zone that has exchanged with the channel by given 
#' time, t, relative to residence time distribution parameters 
#' (tau_0, tau_n and the curvature parameter, sigma or alpha).
#' 
#' @rdname fracExchanged
#' @export
fracExchanged <- function(shape, t, tau_0, tau_n, curvParVal) {
  if(t > tau_n) return(1) #by definition, after tau_n, 100% of the hyporheic zone has exchanged with the channel
  if(shape=="exponent") {
  integrate(
    function(tau) hydrogeom::exponentCCDF(tau, tau_0, 
                                          tau_n, curvParVal) / 
      hydrogeom::exponentIntCCDF(tau_0, tau_n,
                                 tau_0, tau_n, 
                                 sigma = curvParVal),
    tau_0,
    t)$value
  }
  if(shape=="powerLaw") {
    integrate(
      function(tau) hydrogeom::powerLawCCDF(tau, tau_0, 
                                            tau_n, curvParVal) / 
        hydrogeom::powerLawIntCCDF(tau_0, tau_n,
                                   tau_0, tau_n, 
                                   alpha = curvParVal),
      tau_0,
      t)$value
  }
  else {return("Options for shape are 'exponent' or 'powerLaw'")}
}


