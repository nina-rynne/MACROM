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
# find the cost-optimal year T* at which the return to the 1.5°C target is
# achieved. Reference: Lenhart & Workman, "Optimal Control Applied to
# Biological Models", Ch. 20 (Free Terminal Time) and Ch. 21.4/21.5 (Adapted
# Forward-Backward Sweep).
#
# EMPIRICAL FINDING (Phase 2 validation, July 2026): for the MACROM cost
# structure (positive discounted damages every year, CDR capacity-limited),
# the minimized cost V(T) of returning by year T is monotone INCREASING in T
# across the entire feasible range: H(T) = dV/dT > 0 everywhere (verified
# numerically: computed H matches finite-difference dV/dT). Delaying the
# return is never cost-optimal, so the optimum usually sits at the BOUNDARY
# of the feasible set — the earliest year a genuine return is achievable —
# rather than at an interior H(T*) = 0 crossing. This solver therefore
# reports three optimum types:
#   - "interior":          H(T) sign change found and refined per the book
#   - "boundary_earliest": all valid H(T) > 0; T* = earliest feasible return
#                          year (feasibility edge refined to 1-year resolution)
#   - "boundary_latest":   all valid H(T) < 0; T* = latest feasible year
#                          (would indicate delay always pays; not observed)
# H(T*) is reported in all cases; for boundary_earliest it is positive and
# measures the cost (trillion $/yr) of each year of delay beyond T*.
#
# RETURN-FROM-ABOVE FILTER: the raw endpoint constraint x(T) = target admits
# degenerate early solutions that end exactly as cumulative emissions FIRST
# reach the budget (temperature still rising through 1.5°C, near-zero CDR,
# lambda(T) ~ 0). These say nothing about overshoot recovery, so by default
# a candidate year is only valid when annual net emissions at T are negative
# — cumulative emissions falling back through the target from above,
# i.e. temperature returning to 1.5°C after overshoot. Disable with
# require_return_from_above = FALSE to reproduce the raw constraint.
#
# The necessary conditions are unchanged from the fixed-time case (control
# characterisation, adjoint ODE, shooting on lambda(T) to hit the emissions
# target) PLUS one new transversality condition at the optimal terminal time:
#
#   H(T*) = f(T*) + lambda(T*) * g(T*) = 0
#
# where f is the discounted running cost (total_costs_annual) and
# g = baseline_emissions - u_m - u_r is the state dynamics. This condition is
# informative here only because the problem is non-autonomous (discounting and
# time-varying baseline data).
#
# Numerical method: nested shooting. The OUTER loop is a secant search on the
# terminal year T driving H(T) -> 0; the INNER loop is the existing,
# unmodified optimal_control_shooting() from optimal_control_capacity.R
# driving the emissions gap -> 0 for each candidate T. Because the model uses
# an annual grid (dt = 1), T is effectively discrete: candidate years are
# rounded to integers, H(T) is memoized per year, and outer convergence is
# declared when the sign-change bracket collapses to adjacent years.
#
# This machinery was validated against the analytic free-terminal-time
# solution of Lenhart & Workman Example 20.1 in
# scratch/free_terminal_time_validation.R before being applied here.
#
# DEPENDENCIES: optimal_control_capacity.R must be sourced in the session
# before calling optimal_control_free_terminal_time(). Do NOT source
# optimal_control_core.R alongside it — both define same-named functions, and
# the capacity-aware version is the one used throughout this analysis chain.
# The extended-horizon emissions/economic data frames should be built with
# interpolate_ssp_emissions()/interpolate_ssp_economic() using end_year = 2200
# (baseline data beyond 2100 are flat-lined at their 2100 values via
# approx(..., rule = 2)).
#
# Version: 1.0.0
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

#' @title Free Terminal Time Optimal Control via Nested Shooting
#' @description
#' Finds the cost-optimal terminal year T* at which cumulative emissions reach
#' target_emissions via a genuine return from overshoot. The outer search
#' evaluates the terminal Hamiltonian H(T) = dV/dT over integer candidate
#' years; each candidate is evaluated by truncating the scenario data to that
#' horizon and solving the fixed-endpoint problem with the EXISTING
#' optimal_control_shooting() (inner loop, unmodified).
#'
#' If H(T) changes sign, the interior optimum is refined by secant (per
#' Lenhart & Workman). If H(T) > 0 at every valid candidate — the typical
#' MACROM case, see file header — the optimum is the boundary of the feasible
#' set: T* = earliest year a return is achievable, located by integer
#' bisection of the feasibility edge. optimum_type in the returned list says
#' which case occurred.
#'
#' Robustness measures (the outer search is the hardest-converging piece of
#' the method, per Lenhart & Workman):
#' - T candidates never fall at or before deployment start (t_min default is
#'   start_year + max(delays) + 5-year buffer)
#' - A coarse pre-scan across [t_min, t_max] first maps out H(T), confirms
#'   whether a sign change exists, and supplies the secant bracket or the
#'   feasibility-edge bisection interval; with no valid candidate anywhere
#'   the combination is reported infeasible rather than searched blindly
#' - Candidate years where the inner shooting fails to converge, misses the
#'   target, or (by default) arrives at the target from below are treated as
#'   unusable and skipped, not treated as valid H(T) evaluations
#' - Secant steps that leave the bracket (or revisit an evaluated year) fall
#'   back to bisection; H(T) is memoized per integer year so no year is ever
#'   re-solved
#' - A running best (minimum |H|) result is kept so a non-converged search
#'   still returns its closest attempt with diagnostics
#'
#' @param parameter_df Single-row data frame containing all model parameters
#' @param emissions_df Data frame with emissions data extending to t_max
#'   (build with interpolate_ssp_emissions(..., end_year = t_max))
#' @param economic_df Data frame with economic data extending to t_max
#' @param scenario Scenario name to filter data (e.g. "SSP2-Baseline")
#' @param target_emissions Target cumulative emissions (GtCO2). Default NULL
#'   uses co2_target_2100 from parameter_df.
#' @param t_min Earliest candidate terminal year. Default NULL uses
#'   start_year + max(mitigation_delay_years, cdr_delay_years) + 5.
#' @param t_max Latest candidate terminal year (default: 2200)
#' @param mitigation_delay_years Years to delay mitigation deployment (default: 0)
#' @param cdr_delay_years Years to delay CDR deployment (default: 0)
#' @param use_mitigation_capacity_limit Passed through to the inner solver
#'   (default: FALSE)
#' @param mitigation_capacity_function Passed through to the inner solver
#'   (default: NULL)
#' @param use_cdr_capacity_limit Passed through to the inner solver
#'   (default: FALSE)
#' @param cdr_capacity_function Passed through to the inner solver (default: NULL)
#' @param prescan_step Coarse pre-scan spacing in years (default: 15)
#' @param inner_gap_tolerance Max |emission_gap| (GtCO2) for an inner solve to
#'   count as a valid H(T) evaluation (default: 1.0, matching the inner
#'   shooting_tolerance)
#' @param require_return_from_above Logical (default TRUE): a candidate year
#'   is only valid when annual net emissions at T are negative, i.e. the
#'   trajectory arrives at the emissions target from above (temperature
#'   falling back through 1.5°C after overshoot). Set FALSE to accept any
#'   solution with x(T) = target, including degenerate endings where the
#'   budget is first exhausted with temperature still rising.
#' @param max_outer_iterations Cap on outer secant/bisection iterations
#'   (default: 40)
#' @param outer_tolerance |H| below which the outer search stops early
#'   (default: 0.01, in the same units as total_costs_annual — trillion
#'   dollars/year). With integer years the usual stopping condition is the
#'   bracket collapsing to adjacent years, not this tolerance.
#' @param verbose Print progress information (default: TRUE)
#'
#' @return List containing:
#'   - T_star: cost-optimal terminal year (NA if infeasible)
#'   - optimum_type: "interior" (H = 0 crossing refined by secant),
#'     "boundary_earliest" (all valid H > 0; T* = earliest feasible return
#'     year), "boundary_latest" (all valid H < 0; T* = latest feasible year),
#'     or NA when infeasible
#'   - converged: TRUE when the search completed (interior bracket narrowed
#'     to adjacent years / |H| <= outer_tolerance, or boundary edge refined
#'     to 1-year resolution)
#'   - feasible: same convention as cdr_scale_sensitivity.R — TRUE only when
#'     the search converged and the solution at T_star is a valid return
#'   - hamiltonian_at_T_star: H(T*) (NA if infeasible). For boundary_earliest
#'     this is positive: the marginal cost (trillion $/yr) of delaying the
#'     return one year beyond T*.
#'   - solution: full inner solution list at T_star (states, controls, costs
#'     time series) for plotting; best-effort closest attempt if not converged
#'   - prescan_trace: data frame of (T, H, valid, inner_converged,
#'     emission_gap, net_emissions_at_T, final_temperature) from the coarse
#'     pre-scan, for diagnostic plots of H(T)
#'   - evaluations: same columns for every year evaluated (pre-scan + secant
#'     + edge refinement), sorted by T
#'   - outer_iterations, best_abs_H, t_min, t_max, target_emissions, scenario
#'   - infeasible_reason plus best-effort diagnostics (min |H| observed,
#'     final temperature at t_max) when no solution is found
#'
#' @examples
#' \dontrun{
#' result <- optimal_control_free_terminal_time(
#'   parameter_df           = parameter_df[1, ],
#'   emissions_df           = emissions_df_ext,   # end_year = 2200
#'   economic_df            = economic_df_ext,    # end_year = 2200
#'   scenario               = "SSP2-Baseline",
#'   use_cdr_capacity_limit = TRUE,
#'   cdr_capacity_function  = make_logistic_from_zero(2, 100, 0.1, 2025)
#' )
#' result$T_star
#' plot(result$prescan_trace$T, result$prescan_trace$H)
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
                                               prescan_step                  = 15,
                                               inner_gap_tolerance           = 1.0,
                                               require_return_from_above     = TRUE,
                                               max_outer_iterations          = 40,
                                               outer_tolerance               = 0.01,
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

  # T_min: never evaluate a horizon at or before deployment start (initial
  # search nodes must sit strictly inside the controllable window)
  if (is.null(t_min)) {
    t_min <- start_year + max(mitigation_delay_years, cdr_delay_years) + 5
  }
  t_min <- as.integer(round(t_min))
  t_max <- as.integer(round(t_max))

  if (t_min >= t_max) {
    stop("t_min (", t_min, ") must be less than t_max (", t_max, ")")
  }

  if (verbose) {
    cat("=== FREE TERMINAL TIME OPTIMAL CONTROL ===\n")
    cat("Scenario:            ", scenario, "\n")
    cat("Target emissions:    ", target_emissions, "GtCO2\n")
    cat("Terminal year range: [", t_min, ",", t_max, "]\n")
    cat("Pre-scan step:       ", prescan_step, "years\n\n")
  }

  # ============================================================================
  # Memoized H(T) evaluation via the existing inner shooting solver
  # ============================================================================
  # Each candidate year is solved at most once. An evaluation is valid only
  # when (a) the inner shooting converged and hit the emissions target within
  # inner_gap_tolerance, and (b) if require_return_from_above, the trajectory
  # arrives at the target from above (negative annual net emissions at T —
  # a genuine post-overshoot return, not the budget being first exhausted
  # with temperature still rising). Invalid years are skipped by the outer
  # search.

  memo <- new.env(parent = emptyenv())

  evaluate_H <- function(T_cand) {
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

    inner_ok <- !is.null(result) &&
      isTRUE(result$converged) &&
      is.finite(result$emission_gap) &&
      abs(result$emission_gap) <= inner_gap_tolerance

    net_emissions_at_T <- if (!is.null(result)) {
      n <- result$n_years
      result$baseline_annual_emissions[n] - result$qty_mitig[n] -
        result$qty_remov[n]
    } else NA_real_

    ok <- inner_ok &&
      (!require_return_from_above || net_emissions_at_T < 0)

    out <- list(
      T        = T_cand,
      ok       = ok,
      inner_ok = inner_ok,
      net_emis = net_emissions_at_T,
      H        = if (ok) compute_terminal_hamiltonian(result) else NA_real_,
      sol      = result
    )
    memo[[key]] <- out
    out
  }

  # Collect every memoized evaluation into a diagnostic data frame
  collect_evaluations <- function() {
    evals <- lapply(ls(memo), function(k) {
      ev <- memo[[k]]
      data.frame(
        T                  = ev$T,
        H                  = ev$H,
        valid              = ev$ok,
        inner_converged    = ev$inner_ok,
        emission_gap       = if (!is.null(ev$sol)) ev$sol$emission_gap else NA_real_,
        net_emissions_at_T = ev$net_emis,
        final_temperature  = if (!is.null(ev$sol)) ev$sol$final_temperature else NA_real_
      )
    })
    evals <- dplyr::bind_rows(evals)
    evals[order(evals$T), , drop = FALSE]
  }

  # ============================================================================
  # Coarse pre-scan: map H(T) across [t_min, t_max]
  # ============================================================================

  scan_years <- as.integer(round(unique(c(seq(t_min, t_max, by = prescan_step),
                                          t_max))))

  if (verbose) cat("Pre-scanning", length(scan_years), "candidate years...\n")

  for (Ts in scan_years) {
    ev <- evaluate_H(Ts)
    if (verbose) {
      if (ev$ok) {
        cat(sprintf("  T = %d: H = %+.4f\n", Ts, ev$H))
      } else if (ev$inner_ok) {
        cat(sprintf("  T = %d: target hit from below (net emissions %+.2f) — not a return, skipped\n",
                    Ts, ev$net_emis))
      } else {
        cat(sprintf("  T = %d: inner shooting failed (skipped)\n", Ts))
      }
    }
  }

  prescan_trace <- collect_evaluations()

  # ============================================================================
  # Shared result assembly
  # ============================================================================

  outer_iterations <- 0L

  build_result <- function(ev, converged, optimum_type,
                           infeasible_reason = NA_character_) {
    if (verbose && converged) {
      cat(sprintf("\nConverged (%s): T* = %d, H(T*) = %+.4f (%d outer iterations)\n",
                  optimum_type, ev$T, ev$H, outer_iterations))
    }
    evaluations <- collect_evaluations()
    list(
      T_star                = if (converged) ev$T else NA_real_,
      optimum_type          = if (converged) optimum_type else NA_character_,
      converged             = converged,
      feasible              = converged && isTRUE(ev$ok),
      hamiltonian_at_T_star = if (converged) ev$H else NA_real_,
      solution              = ev$sol,
      prescan_trace         = prescan_trace,
      evaluations           = evaluations,
      outer_iterations      = outer_iterations,
      best_abs_H            = if (any(evaluations$valid, na.rm = TRUE)) {
        min(abs(evaluations$H[evaluations$valid]), na.rm = TRUE)
      } else NA_real_,
      infeasible_reason     = infeasible_reason,
      t_min                 = t_min,
      t_max                 = t_max,
      target_emissions      = target_emissions,
      scenario              = scenario
    )
  }

  # ============================================================================
  # Infeasible: no valid return year anywhere in [t_min, t_max]
  # ============================================================================

  scan_valid <- prescan_trace[prescan_trace$valid & is.finite(prescan_trace$H),
                              , drop = FALSE]

  if (nrow(scan_valid) == 0) {
    reason <- if (any(prescan_trace$inner_converged)) {
      "target only reachable from below (no post-overshoot return) at every scanned year"
    } else {
      "inner shooting failed to converge at every scanned terminal year"
    }
    if (verbose) cat("\nINFEASIBLE:", reason, "\n")

    ev_tmax <- memo[[as.character(t_max)]]

    res <- build_result(
      ev = list(T = NA_real_, ok = FALSE, H = NA_real_, sol = NULL),
      converged = FALSE, optimum_type = NA_character_,
      infeasible_reason = reason
    )
    res$final_temperature_at_t_max <-
      if (!is.null(ev_tmax) && !is.null(ev_tmax$sol)) {
        ev_tmax$sol$final_temperature
      } else NA_real_
    return(res)
  }

  # ============================================================================
  # Interior optimum: secant on T within a sign-change bracket
  # ============================================================================
  # Runs when the pre-scan (or boundary-edge refinement below) finds
  # consecutive valid years with opposite H signs — the textbook case.

  interior_search <- function(T_lo, H_lo, T_hi, H_hi) {

    if (verbose) {
      cat(sprintf("\nBracket found: T in [%d, %d], H in [%+.4f, %+.4f]\n",
                  T_lo, T_hi, H_lo, H_hi))
    }

    best_ev    <- NULL
    best_abs_H <- Inf
    for (Tb in c(T_lo, T_hi)) {
      ev <- evaluate_H(Tb)  # memoized, no re-solve
      if (abs(ev$H) < best_abs_H) {
        best_ev <- ev
        best_abs_H <- abs(ev$H)
      }
    }

    for (iteration in seq_len(max_outer_iterations)) {
      outer_iterations <<- outer_iterations + 1L

      # Bracket collapsed to adjacent years: the discrete optimum is the
      # endpoint with the smaller |H|
      if ((T_hi - T_lo) <= 1L) {
        ev_lo <- evaluate_H(T_lo)
        ev_hi <- evaluate_H(T_hi)
        ev <- if (abs(ev_lo$H) <= abs(ev_hi$H)) ev_lo else ev_hi
        return(build_result(ev, converged = TRUE, optimum_type = "interior"))
      }

      # Secant proposal, rounded to an integer year; bisection fallback when
      # the proposal is non-finite, leaves the open bracket, or lands on an
      # already-evaluated endpoint
      T_new <- as.integer(round(T_lo - H_lo * (T_hi - T_lo) / (H_hi - H_lo)))
      if (!is.finite(T_new) || T_new <= T_lo || T_new >= T_hi) {
        T_new <- as.integer(floor((T_lo + T_hi) / 2))
      }

      ev <- evaluate_H(T_new)

      if (!ev$ok) {
        # Invalid year inside the bracket (should be rare). Try the bisection
        # midpoint instead; if that also fails, return the best attempt.
        T_new <- as.integer(floor((T_lo + T_hi) / 2))
        ev <- evaluate_H(T_new)
        if (!ev$ok) {
          if (verbose) {
            cat("Invalid candidate years inside the bracket; returning best attempt\n")
          }
          return(build_result(
            if (!is.null(best_ev)) best_ev else ev, converged = FALSE,
            optimum_type = NA_character_,
            infeasible_reason = "invalid candidate years inside the secant bracket"
          ))
        }
      }

      if (abs(ev$H) < best_abs_H) {
        best_ev <- ev
        best_abs_H <- abs(ev$H)
      }

      if (verbose) {
        cat(sprintf("  Outer iter %2d: T = %d, H(T) = %+.4f, bracket [%d, %d]\n",
                    iteration, ev$T, ev$H, T_lo, T_hi))
      }

      if (abs(ev$H) <= outer_tolerance) {
        return(build_result(ev, converged = TRUE, optimum_type = "interior"))
      }

      # Narrow the bracket, preserving the sign change
      if (sign(ev$H) == sign(H_lo)) {
        T_lo <- ev$T; H_lo <- ev$H
      } else {
        T_hi <- ev$T; H_hi <- ev$H
      }
    }

    if (verbose) {
      cat("Outer search reached max iterations; returning best attempt\n")
    }
    build_result(
      best_ev, converged = FALSE, optimum_type = NA_character_,
      infeasible_reason = "outer search reached max iterations without collapsing the bracket"
    )
  }

  # Look for a sign change between consecutive VALID scan points
  if (nrow(scan_valid) >= 2) {
    for (k in 2:nrow(scan_valid)) {
      if (sign(scan_valid$H[k - 1]) != sign(scan_valid$H[k])) {
        return(interior_search(scan_valid$T[k - 1], scan_valid$H[k - 1],
                               scan_valid$T[k],     scan_valid$H[k]))
      }
    }
  }

  # ============================================================================
  # Boundary optimum: all valid H(T) share one sign
  # ============================================================================
  # H = dV/dT (verified numerically for this model). All H > 0 means the cost
  # of returning by T rises with T everywhere it is achievable, so the
  # cheapest achievable return is the EARLIEST feasible year — a boundary
  # minimum. Refine the feasibility edge to 1-year resolution by integer
  # bisection between the last invalid and first valid known years. (All
  # H < 0 — never observed — would symmetrically give the latest year.)

  if (all(scan_valid$H > 0)) {

    T_first_valid <- min(scan_valid$T)
    below <- scan_years[scan_years < T_first_valid]

    if (length(below) == 0) {
      # Even the first candidate year is a valid return: t_min is binding
      if (verbose) {
        cat("\nAll valid H(T) > 0 and t_min itself is feasible;",
            "T* = t_min (earlier returns may exist below the search range)\n")
      }
      T_edge <- T_first_valid
    } else {
      lo <- max(below)          # invalid (or unusable) scanned year
      hi <- T_first_valid       # valid year
      if (verbose) {
        cat(sprintf("\nAll valid H(T) > 0: refining earliest feasible return year in (%d, %d]...\n",
                    lo, hi))
      }
      while (hi - lo > 1L) {
        outer_iterations <- outer_iterations + 1L
        mid <- as.integer(floor((lo + hi) / 2))
        ev_mid <- evaluate_H(mid)
        if (verbose) {
          cat(sprintf("  T = %d: %s\n", mid,
                      if (ev_mid$ok) sprintf("valid, H = %+.4f", ev_mid$H)
                      else "not a valid return"))
        }
        if (ev_mid$ok) hi <- mid else lo <- mid
      }
      T_edge <- hi
    }

    ev_edge <- evaluate_H(T_edge)

    # If the refined edge year turns out to have H < 0, an interior crossing
    # exists between the edge and the first coarse valid point after all —
    # hand over to the secant
    if (ev_edge$H < 0) {
      T_above <- min(scan_valid$T[scan_valid$T > T_edge])
      ev_above <- evaluate_H(T_above)
      return(interior_search(T_edge, ev_edge$H, T_above, ev_above$H))
    }

    return(build_result(ev_edge, converged = TRUE,
                        optimum_type = "boundary_earliest"))
  }

  # All valid H < 0: cost falls with T everywhere achievable; the optimum is
  # the latest feasible year (t_max-limited if the last scan point is valid)
  T_last_valid <- max(scan_valid$T)
  above <- scan_years[scan_years > T_last_valid]

  if (length(above) == 0) {
    T_edge <- T_last_valid
    if (verbose) {
      cat("\nAll valid H(T) < 0; T* = t_max — optimum truncated by the search horizon\n")
    }
  } else {
    lo <- T_last_valid
    hi <- min(above)
    if (verbose) {
      cat(sprintf("\nAll valid H(T) < 0: refining latest feasible year in [%d, %d)...\n",
                  lo, hi))
    }
    while (hi - lo > 1L) {
      outer_iterations <- outer_iterations + 1L
      mid <- as.integer(floor((lo + hi) / 2))
      ev_mid <- evaluate_H(mid)
      if (ev_mid$ok) lo <- mid else hi <- mid
    }
    T_edge <- lo
  }

  ev_edge <- evaluate_H(T_edge)

  if (ev_edge$H > 0) {
    # Symmetric hand-over: crossing between the last coarse valid point and
    # the refined edge
    T_below <- max(scan_valid$T[scan_valid$T < T_edge])
    ev_below <- evaluate_H(T_below)
    return(interior_search(T_below, ev_below$H, T_edge, ev_edge$H))
  }

  build_result(ev_edge, converged = TRUE, optimum_type = "boundary_latest")
}
