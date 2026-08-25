# ==============================================================================
# Free Terminal Time Optimal Control for Climate Temperature Overshoot
#
# Part of: MACROM: An Optimal Control Model for Balancing Climate Change
# Abatement and Damage Trade-offs
# Authors: Nina Rynne, Michael Bode, Melanie Roberts, Ryan Heneghan
# Institution: Griffith University
#
# Copyright (c) 2026 Nina Rynne
# Licensed under CC-BY-4.0 - see LICENSE file for details
#
# Citation: If you use this code, please cite:
#   Rynne, N., Bode, M., Roberts, M., & Heneghan, R. (2025).
#   MACROM: An Optimal Control Model for Balancing Climate Change
#   Abatement and Damage Trade-offs.
#   https://doi.org/10.5281/zenodo.18463951
#
# Description:
# Converts MACROM's fixed-endpoint, fixed-time optimal control problem
# (x(T) = target_emissions at T = 2100) into a free-terminal-time problem:
# find the cost-optimal year T* at which cumulative emissions first return to
# the 1.5°C-consistent target, for K/r combinations that cannot reach it by
# 2100 at all (see cdr_scale_free_time_sensitivity.R for the stage-1 gate that
# decides which cells reach this solver).
#
# WHAT T* MEANS: T* is the year the historical backlog of overshoot emissions
# (accumulated while CDR capacity was still ramping up) is first cleared. It
# does NOT mean emissions or CDR deployment stop at T* -- cumulative emissions
# keep evolving after T*, same as always. What stops is the SCOPE of this
# optimization: the model only prices the cost of clearing the backlog, and is
# silent on the (comparatively simple) maintenance problem of keeping sources
# and sinks balanced once the target is first reached. This mirrors the
# classical free-terminal-time convention (Lenhart & Workman, "Optimal Control
# Applied to Biological Models", Ch. 20 -- e.g. disease treatment, where the
# state is "cured" at T*): the optimization's horizon ends at the terminal
# event even though the underlying state continues to evolve after it in
# reality.
#
# WHY THIS ISN'T AN EQUALITY-CONSTRAINED PROBLEM: an earlier version of this
# solver searched for a year T where x(T) = target held exactly, using the
# transversality condition H(T) = dV/dT = 0 (or its boundary form) from
# Lenhart & Workman Ch. 20. That formulation is degenerate here because, for
# strong-CDR cells, x(T) = target can hold at MULTIPLE T's: a genuine descent
# through the target, followed by overshoot below it, followed by drifting
# back UP to touch the target again later "from below" as CDR eases off. A
# per-candidate-year require_return_from_above filter was added to reject the
# "from below" touches, but it only had years actually sampled by a coarse
# pre-scan grid to work with -- and for strong enough CDR, the genuine
# descending crossing fell entirely inside a single unsampled gap in that
# grid, so the search found no valid year anywhere and wrongly reported the
# cell infeasible.
#
# THE FIX: reframe the question as minimum-cost-to-first-reach-a-target-SET
# (x(T) <= target, an inequality) rather than minimum-cost-to-hit-an-exact-
# VALUE (x(T) = target, an equality). Two facts make this well-posed and
# simple to solve:
#   1. EMPIRICAL FINDING (Phase 2 validation, July 2026): for the MACROM cost
#      structure (positive discounted damages every year, CDR capacity-
#      limited), the minimized cost of returning by year T is monotone
#      INCREASING in T -- delaying the return is never cost-optimal. So the
#      cost-optimal T is always the EARLIEST year the target is achievable,
#      never a later one.
#   2. At the TRUE earliest achievable year, arrival is GUARANTEED to be a
#      genuine descent from above: an instant earlier the target wasn't
#      achievable at all (cumulative emissions were still above budget, with
#      no control path able to bring them down in time), so the earliest
#      achievable year cannot be a later, already-past, drifted-back-up
#      touch -- that would contradict it being the earliest. No per-year
#      "arrived from above" filter is needed; it is a structural consequence
#      of asking for the EARLIEST reachable year specifically.
#
# Because "is the target achievable by year T" (does the inner fixed-endpoint
# shooting solve converge within tolerance) is a monotonic property of T --
# false for every T before the true crossing, true for every T at or after it
# -- the earliest such T is found by plain integer BISECTION, no coarse
# pre-scan grid required and no risk of stepping over a narrow window:
# typically 7-8 inner solves regardless of how wide [t_min, t_max] is.
#
# Numerical method: nested search. The OUTER loop is integer bisection on
# reachability, driving the bracket [not-yet-achievable, achievable] down to
# adjacent years; the INNER loop is the existing, unmodified
# optimal_control_shooting() from optimal_control_capacity.R, driving the
# emissions gap -> 0 for each candidate T. Because the model uses an annual
# grid (dt = 1), T is effectively discrete: candidate years are rounded to
# integers and reachability is memoized per year so no year is ever re-solved.
#
# compute_terminal_hamiltonian() (Section 1) is retained and still reported
# per solution (hamiltonian_at_T_star) as a diagnostic -- it is the marginal
# cost, in trillion $/yr, of delaying the return one year beyond T* -- but it
# no longer drives the search; H(T) > 0 everywhere is now a documented
# empirical property used to justify the bisection shortcut, not something
# recomputed and sign-tested at every candidate year.
#
# DEPENDENCIES: optimal_control_capacity.R must be sourced in the session
# before calling optimal_control_free_terminal_time(). Do NOT source
# optimal_control_core.R alongside it — both define same-named functions, and
# the capacity-aware version is the one used throughout this analysis chain.
# The extended-horizon emissions/economic data frames should be built with
# interpolate_ssp_emissions()/interpolate_ssp_economic() using an end_year
# that covers the full [t_min, t_max] search range.
#
# Version: 2.0.0
# Last updated: July 2026
# ==============================================================================


# ==============================================================================
# Section 1: Terminal Hamiltonian
# ==============================================================================

#' @title Compute Terminal Hamiltonian of a Solved Optimal Control Solution
#' @description
#' Evaluates the Hamiltonian H = f + lambda * g at the final year of a solution
#' returned by optimal_control_solve() / optimal_control_shooting(). Pure
#' post-processing: uses only fields already present in the solution list.
#'
#' f(T) = total_costs_annual[T] is the discounted running cost and
#' g(T) = baseline_annual_emissions[T] - qty_mitig[T] - qty_remov[T] is the
#' state dynamics (annual net emissions). adjoint_var is already the correct
#' present-value adjoint for this Hamiltonian (verified against the control
#' characterisation and adjoint derivative in optimal_control_solve()).
#'
#' At the cost-optimal free terminal time, H(T*) = 0. For T below the optimum
#' more time still reduces total cost (expected H < 0); above it, more time
#' adds cost through accumulated discounted damage (expected H > 0).
#'
#' @param solution Solution list from optimal_control_solve() or
#'   optimal_control_shooting()
#'
#' @return Numeric scalar: the Hamiltonian at the final year of the solution
compute_terminal_hamiltonian <- function(solution) {
  n <- solution$n_years
  solution$total_costs_annual[n] + solution$adjoint_var[n] *
    (solution$baseline_annual_emissions[n] - solution$qty_mitig[n] -
       solution$qty_remov[n])
}


# ==============================================================================
# Section 2: Data Truncation Helper
# ==============================================================================

#' @title Truncate Scenario Data to a Candidate Terminal Year
#' @description
#' Filters the emissions and economic data frames to Year <= t_max so the
#' fixed-time solver can be run on a shortened horizon. Both data frames are
#' expected to already extend to the full search horizon (e.g. 2200).
#'
#' @param emissions_df Data frame with emissions data (long format, Year column)
#' @param economic_df Data frame with economic data (long format, Year column)
#' @param t_max Numeric last year to retain
#'
#' @return List with elements emissions_df and economic_df, each truncated
truncate_scenario_data <- function(emissions_df, economic_df, t_max) {
  list(
    emissions_df = emissions_df %>% dplyr::filter(Year <= t_max),
    economic_df  = economic_df  %>% dplyr::filter(Year <= t_max)
  )
}


# ==============================================================================
# Section 3: Free Terminal Time Solver
# ==============================================================================

#' @title Free Terminal Time Optimal Control via Bisection on Reachability
#' @description
#' Finds T*, the earliest year at which cumulative emissions can first be
#' brought down to target_emissions via a genuine descent from overshoot
#' (see file header for why this is both the cost-optimal year and
#' guaranteed to arrive "from above", with no separate filter required).
#'
#' Each candidate year T is evaluated by truncating the scenario data to that
#' horizon and solving the fixed-endpoint problem with the EXISTING
#' optimal_control_shooting() (inner loop, unmodified): T is "reachable" if
#' that inner solve converges and hits target_emissions within
#' inner_gap_tolerance. Reachability is monotone in T (false, false, ...,
#' false, true, true, ..., true), so the earliest reachable year is found by
#' plain integer bisection between a known-unreachable year (t_min) and a
#' known-reachable one (t_max) -- no coarse pre-scan grid, no risk of
#' stepping over a narrow feasible window.
#'
#' @param parameter_df Single-row data frame containing all model parameters
#' @param emissions_df Data frame with emissions data extending to t_max
#'   (build with interpolate_ssp_emissions(..., end_year = t_max))
#' @param economic_df Data frame with economic data extending to t_max
#' @param scenario Scenario name to filter data (e.g. "SSP2-Baseline")
#' @param target_emissions Target cumulative emissions (GtCO2). Default NULL
#'   uses co2_target_2100 from parameter_df.
#' @param t_min Earliest candidate terminal year -- must be genuinely
#'   unreachable (a buffer past deployment start, not the true crossing).
#'   Default NULL uses start_year + max(mitigation_delay_years,
#'   cdr_delay_years) + 5.
#' @param t_max Latest candidate terminal year (default: 2200). Must be
#'   reachable, or the cell is reported infeasible.
#' @param mitigation_delay_years Years to delay mitigation deployment (default: 0)
#' @param cdr_delay_years Years to delay CDR deployment (default: 0)
#' @param use_mitigation_capacity_limit Passed through to the inner solver
#'   (default: FALSE)
#' @param mitigation_capacity_function Passed through to the inner solver
#'   (default: NULL)
#' @param use_cdr_capacity_limit Passed through to the inner solver
#'   (default: FALSE)
#' @param cdr_capacity_function Passed through to the inner solver (default: NULL)
#' @param inner_gap_tolerance Max |emission_gap| (GtCO2) for an inner solve to
#'   count as "reachable" (default: 1.0, matching the inner shooting_tolerance)
#' @param max_bisection_iterations Cap on bisection iterations, purely a
#'   safety net -- an integer bracket of width W converges in ceiling(log2(W))
#'   steps, e.g. ~8 for a 200-year search range (default: 40)
#' @param verbose Print progress information (default: TRUE)
#'
#' @return List containing:
#'   - T_star: earliest cost-optimal terminal year (NA if infeasible)
#'   - optimum_type: "earliest_reachable" if converged, else NA
#'   - converged: TRUE when bisection found the earliest reachable year
#'   - feasible: same convention as cdr_scale_sensitivity.R -- TRUE when a
#'     genuine return is achievable within [t_min, t_max]
#'   - hamiltonian_at_T_star: H(T*), trillion $/yr -- diagnostic only (the
#'     marginal cost of delaying the return one further year); does not
#'     drive T* selection (see file header)
#'   - solution: full inner solution list at T_star (states, controls, costs
#'     time series) for plotting; best-effort closest attempt if not converged
#'   - evaluations: data frame of (T, reachable, net_emissions_at_T,
#'     emission_gap, final_temperature) for every year actually evaluated
#'     during bisection, for diagnostics
#'   - outer_iterations: number of bisection steps used
#'   - t_min, t_max, target_emissions, scenario
#'   - infeasible_reason plus best-effort diagnostics (final temperature at
#'     t_max) when no solution is found
#'
#' @examples
#' \dontrun{
#' result <- optimal_control_free_terminal_time(
#'   parameter_df           = parameter_df[1, ],
#'   emissions_df           = emissions_df_ext,   # end_year = 2300
#'   economic_df            = economic_df_ext,    # end_year = 2300
#'   scenario               = "SSP2-Baseline",
#'   use_cdr_capacity_limit = TRUE,
#'   cdr_capacity_function  = make_logistic_from_zero(2, 100, 0.1, 2025)
#' )
#' result$T_star
#' }
optimal_control_free_terminal_time <- function(parameter_df,
                                               emissions_df,
                                               economic_df,
                                               scenario,
                                               target_emissions              = NULL,
                                               t_min                         = NULL,
                                               t_max                         = 2200,
                                               mitigation_delay_years        = 0,
                                               cdr_delay_years               = 0,
                                               use_mitigation_capacity_limit = FALSE,
                                               mitigation_capacity_function  = NULL,
                                               use_cdr_capacity_limit        = FALSE,
                                               cdr_capacity_function         = NULL,
                                               inner_gap_tolerance           = 1.0,
                                               max_bisection_iterations      = 40,
                                               verbose                       = TRUE) {

  # ============================================================================
  # Input validation and setup
  # ============================================================================

  if (!is.data.frame(parameter_df) || nrow(parameter_df) != 1) {
    stop("parameter_df must be a single-row data frame")
  }

  if (is.null(target_emissions)) {
    if ("co2_target_2100" %in% names(parameter_df)) {
      target_emissions <- parameter_df$co2_target_2100
    } else {
      stop("target_emissions must be specified or parameter_df must contain co2_target_2100")
    }
  }

  # Filter to the requested scenario once; per-candidate truncation below.
  # optimal_control_shooting() filters by scenario itself, so the Scenario
  # column is retained.
  emissions_scenario <- emissions_df %>%
    dplyr::filter(Scenario == scenario) %>%
    dplyr::arrange(Year)
  economic_scenario <- economic_df %>%
    dplyr::filter(Scenario == scenario) %>%
    dplyr::arrange(Year)

  if (nrow(emissions_scenario) == 0 || nrow(economic_scenario) == 0) {
    stop("No data found for scenario: ", scenario)
  }

  start_year    <- min(emissions_scenario$Year)
  data_end_year <- max(emissions_scenario$Year)

  if (data_end_year < t_max) {
    stop("Scenario data end at ", data_end_year, " but t_max = ", t_max,
         ". Rebuild emissions_df/economic_df with end_year = ", t_max,
         " (interpolate_ssp_emissions / interpolate_ssp_economic).")
  }

  # T_min: never evaluate a horizon at or before deployment start (must be a
  # genuinely unreachable anchor for the bisection, not the true crossing)
  if (is.null(t_min)) {
    t_min <- start_year + max(mitigation_delay_years, cdr_delay_years) + 5
  }
  t_min <- as.integer(round(t_min))
  t_max <- as.integer(round(t_max))

  if (t_min >= t_max) {
    stop("t_min (", t_min, ") must be less than t_max (", t_max, ")")
  }

  if (verbose) {
    cat("=== FREE TERMINAL TIME OPTIMAL CONTROL (earliest-reachable search) ===\n")
    cat("Scenario:            ", scenario, "\n")
    cat("Target emissions:    ", target_emissions, "GtCO2\n")
    cat("Terminal year range: [", t_min, ",", t_max, "]\n\n")
  }

  # ============================================================================
  # Memoized reachability evaluation via the existing inner shooting solver
  # ============================================================================
  # A candidate year T is "reachable" when the fixed-endpoint problem
  # truncated to that horizon converges and hits target_emissions within
  # inner_gap_tolerance. Each year is solved at most once (memoized).

  memo <- new.env(parent = emptyenv())

  evaluate_T <- function(T_cand) {
    T_cand <- as.integer(round(T_cand))
    key <- as.character(T_cand)
    if (!is.null(memo[[key]])) return(memo[[key]])

    truncated <- truncate_scenario_data(emissions_scenario, economic_scenario,
                                        t_max = T_cand)

    result <- tryCatch(
      optimal_control_shooting(
        parameter_df                  = parameter_df,
        emissions_df                  = truncated$emissions_df,
        economic_df                   = truncated$economic_df,
        scenario                      = scenario,
        target_emissions              = target_emissions,
        mitigation_delay_years        = mitigation_delay_years,
        cdr_delay_years               = cdr_delay_years,
        use_mitigation_capacity_limit = use_mitigation_capacity_limit,
        mitigation_capacity_function  = mitigation_capacity_function,
        use_cdr_capacity_limit        = use_cdr_capacity_limit,
        cdr_capacity_function         = cdr_capacity_function,
        verbose                       = FALSE
      ),
      error = function(e) NULL
    )

    reachable <- !is.null(result) &&
      isTRUE(result$converged) &&
      is.finite(result$emission_gap) &&
      abs(result$emission_gap) <= inner_gap_tolerance

    net_emissions_at_T <- if (!is.null(result)) {
      n <- result$n_years
      result$baseline_annual_emissions[n] - result$qty_mitig[n] -
        result$qty_remov[n]
    } else NA_real_

    out <- list(T = T_cand, reachable = reachable, net_emis = net_emissions_at_T,
               sol = result)
    memo[[key]] <- out

    if (verbose) {
      cat(sprintf("  T = %d: %s\n", T_cand,
                  if (reachable) sprintf("reachable (net emissions %+.2f)", net_emissions_at_T)
                  else "not reachable"))
    }
    out
  }

  # Collect every memoized evaluation into a diagnostic data frame
  collect_evaluations <- function() {
    evals <- lapply(ls(memo), function(k) {
      ev <- memo[[k]]
      data.frame(
        T                  = ev$T,
        reachable          = ev$reachable,
        net_emissions_at_T = ev$net_emis,
        emission_gap       = if (!is.null(ev$sol)) ev$sol$emission_gap else NA_real_,
        final_temperature  = if (!is.null(ev$sol)) ev$sol$final_temperature else NA_real_
      )
    })
    evals <- dplyr::bind_rows(evals)
    evals[order(evals$T), , drop = FALSE]
  }

  # ============================================================================
  # Shared result assembly
  # ============================================================================

  build_result <- function(ev, converged, infeasible_reason = NA_character_,
                           bisection_iterations = 0L) {

    if (converged) {
      H_val <- compute_terminal_hamiltonian(ev$sol)
      if (verbose) {
        cat(sprintf("\nEarliest reachable year: T* = %d (H(T*) = %+.4f, %d bisection steps)\n",
                    ev$T, H_val, bisection_iterations))
      }
      # Sanity check: the monotonicity argument in the file header guarantees
      # arrival from above at the TRUE earliest reachable year. If this ever
      # fires, something about this cell doesn't match that assumption and
      # the result deserves a closer look.
      if (!isTRUE(ev$net_emis < 0)) {
        warning("T* = ", ev$T, " for scenario '", scenario, "' does not arrive ",
                "from above (net_emissions_at_T = ", round(ev$net_emis, 3),
                "). This contradicts the monotonicity argument the earliest-",
                "reachable search relies on (see file header) -- treat this ",
                "cell's result with caution.")
      }
    } else {
      H_val <- NA_real_
    }

    list(
      T_star                = if (converged) ev$T else NA_real_,
      optimum_type          = if (converged) "earliest_reachable" else NA_character_,
      converged             = converged,
      feasible              = converged && isTRUE(ev$reachable),
      hamiltonian_at_T_star = H_val,
      solution              = if (converged) ev$sol else NULL,
      evaluations           = collect_evaluations(),
      outer_iterations      = bisection_iterations,
      infeasible_reason     = infeasible_reason,
      t_min                 = t_min,
      t_max                 = t_max,
      target_emissions      = target_emissions,
      scenario              = scenario
    )
  }

  # ============================================================================
  # Step 1: is the target reachable anywhere in [t_min, t_max] at all?
  # ============================================================================

  ev_max <- evaluate_T(t_max)

  if (!ev_max$reachable) {
    if (verbose) {
      cat("\nINFEASIBLE: target not reachable even at t_max =", t_max, "\n")
    }
    res <- build_result(
      ev = list(T = NA_real_, reachable = FALSE, net_emis = NA_real_, sol = NULL),
      converged = FALSE,
      infeasible_reason = paste0("target not reachable by t_max = ", t_max)
    )
    res$final_temperature_at_t_max <- if (!is.null(ev_max$sol)) {
      ev_max$sol$final_temperature
    } else NA_real_
    return(res)
  }

  # ============================================================================
  # Step 2: edge case -- is t_min itself already reachable?
  # ============================================================================

  ev_min <- evaluate_T(t_min)

  if (ev_min$reachable) {
    if (verbose) {
      cat("\nt_min itself is reachable; T* = t_min",
          "(earlier returns may exist below the search range)\n")
    }
    return(build_result(ev_min, converged = TRUE, bisection_iterations = 0L))
  }

  # ============================================================================
  # Step 3: bisect on reachability to find the true earliest reachable year
  # ============================================================================
  # Reachability is monotone in T (false, ..., false, true, ..., true; see
  # file header for why), so standard integer bisection between the known-
  # unreachable t_min and known-reachable t_max converges on the true
  # earliest reachable year in ceiling(log2(t_max - t_min)) steps -- no
  # coarse grid, no risk of stepping over a narrow feasible window.

  lo <- t_min                     # known not reachable
  hi <- t_max; ev_hi <- ev_max     # known reachable
  iterations <- 0L

  while (hi - lo > 1L && iterations < max_bisection_iterations) {
    iterations <- iterations + 1L
    mid <- as.integer(floor((lo + hi) / 2))
    ev_mid <- evaluate_T(mid)
    if (ev_mid$reachable) {
      hi <- mid; ev_hi <- ev_mid
    } else {
      lo <- mid
    }
  }

  build_result(ev_hi, converged = TRUE, bisection_iterations = iterations)
}
