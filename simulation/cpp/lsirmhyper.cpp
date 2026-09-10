#include <RcppArmadillo.h>
#include <cmath>
#include <limits>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

#include "hyperboloid_utils.h"

namespace {

inline double log1pexp(double x) {
  return x > 0.0 ? x + std::log1p(std::exp(-x)) : std::log1p(std::exp(x));
}

inline double bernoulli_loglik(double y, double eta) {
  return y == 1.0 ? -log1pexp(-eta) : -log1pexp(eta);
}

inline bool observed_binary(double y) { return y == 0.0 || y == 1.0; }

inline double log_halfnormal_sd(double sigma, double scale) {
  return std::log(2.0) + R::dnorm4(sigma, 0.0, scale, 1);
}

inline double log_invgamma_variance_as_sd(double sigma, double shape,
                                          double rate) {
  const double variance = sigma * sigma;
  return shape * std::log(rate) - R::lgammafn(shape) -
         (shape + 1.0) * std::log(variance) - rate / variance +
         std::log(2.0 * sigma);
}

template <typename VectorX, typename VectorY>
inline double euclidean_distance(const VectorX &x, const VectorY &y) {
  double sum_of_squares = 0.0;
  for (arma::uword q = 0; q < x.n_elem; ++q) {
    const double difference = x(q) - y(q);
    sum_of_squares += difference * difference;
  }
  return std::sqrt(sum_of_squares);
}

inline double latent_radius(const arma::rowvec &x, bool hyperbolic) {
  return hyperbolic ? hyperboloid_utils::radius_hyperboloid(x)
                    : arma::norm(x, 2);
}

} // namespace

// Paper-production sampler.  The statistical target is
//   logit P(Y_ij=1) = theta_i + beta_j - d(Z_i,W_j),
// with gamma fixed to one and both main-effect prior means fixed to zero.
// H^d uses ndim=d+1 ambient coordinates; E^d uses ndim=d coordinates.
// [[Rcpp::export]]
Rcpp::List lsirm1pl_cpp(
    arma::mat data, const int ndim, const int niter, const int nburn,
    const int nthin, const double jump_beta, const double jump_theta,
    const double jump_z, const double jump_w, double pr_sd_beta,
    const double pr_a_beta, const double pr_b_beta, const double pr_a_theta,
    const double pr_b_theta, double pr_sd_theta, double pr_sd_z, double pr_sd_w,
    const bool hyperbolic, const bool estimate_beta_sd = false,
    const bool adaptation = true, const int adapt_interval = 25,
    const double target_accept_1d = 0.44, const double adapt_exponent = 0.6,
    Rcpp::Nullable<Rcpp::List> radial_block_control = R_NilValue,
    const bool stable_intrinsic_proposal = false,
    const double max_hyperbolic_radius = 13.0,
    const bool polar_split = false) {
  const int nrow = data.n_rows;
  const int ncol = data.n_cols;
  const int intrinsic_dim = hyperbolic ? ndim - 1 : ndim;

  if (nrow < 1 || ncol < 1) {
    Rcpp::stop("data must be nonempty.");
  }
  if (niter <= 0 || nburn < 0 || nburn >= niter) {
    Rcpp::stop("Require niter > 0 and 0 <= nburn < niter.");
  }
  if (nthin <= 0) {
    Rcpp::stop("nthin must be positive.");
  }
  if ((hyperbolic && ndim < 3) || (!hyperbolic && ndim < 1)) {
    Rcpp::stop("The paper engine requires H^d with d>=2 or E^d with d>=1.");
  }
  if (!(jump_beta > 0.0) || !(jump_theta > 0.0) || !(jump_z > 0.0) ||
      !(jump_w > 0.0)) {
    Rcpp::stop("All initial proposal scales must be positive.");
  }
  if (!(pr_sd_beta > 0.0) || !(pr_sd_theta > 0.0) || !(pr_sd_z > 0.0) ||
      !(pr_sd_w > 0.0) || !(pr_a_beta > 0.0) || !(pr_b_beta > 0.0) ||
      !(pr_a_theta > 0.0) || !(pr_b_theta > 0.0)) {
    Rcpp::stop("All prior scales, shapes, and rates must be positive.");
  }
  if (!(max_hyperbolic_radius > 0.0) || !std::isfinite(max_hyperbolic_radius)) {
    Rcpp::stop("max_hyperbolic_radius must be finite and positive.");
  }
  if (hyperbolic &&
      max_hyperbolic_radius >
          hyperboloid_utils::default_max_hyperbolic_radius &&
      !stable_intrinsic_proposal) {
    Rcpp::warning(
        "A radius cap above 13 is numerically unsafe with the legacy "
        "Lorentz-projection proposal; use stable_intrinsic_proposal=true.");
  }
  if (adaptation) {
    if (adapt_interval <= 0 || adapt_interval >= nburn) {
      Rcpp::stop("adapt_interval must be in [1, nburn).");
    }
    if (!(target_accept_1d > 0.0 && target_accept_1d < 1.0)) {
      Rcpp::stop("target_accept_1d must be in (0,1).");
    }
    if (!(adapt_exponent > 0.5 && adapt_exponent <= 1.0)) {
      Rcpp::stop("adapt_exponent must be in (0.5,1].");
    }
  }
  for (int i = 0; i < nrow; ++i) {
    for (int j = 0; j < ncol; ++j) {
      const double y = data(i, j);
      if (std::isnan(y)) {
        continue;
      }
      if (!std::isfinite(y) || !observed_binary(y)) {
        Rcpp::stop("data must contain only 0, 1, or NA/NaN.");
      }
    }
  }

  bool use_radial_block = false;
  double radial_jump = 0.15;
  bool has_init_state = false;
  Rcpp::List init_state;
  if (radial_block_control.isNotNull()) {
    Rcpp::List ctl(radial_block_control.get());
    use_radial_block = ctl.containsElementNamed("enabled")
                           ? Rcpp::as<bool>(ctl["enabled"])
                           : true;
    if (ctl.containsElementNamed("jump")) {
      radial_jump = Rcpp::as<double>(ctl["jump"]);
    }
    if (ctl.containsElementNamed("init_state")) {
      init_state = Rcpp::as<Rcpp::List>(ctl["init_state"]);
      has_init_state = true;
    }
  }
  if (use_radial_block && !hyperbolic) {
    Rcpp::stop("The radial block is defined only for the hyperbolic model.");
  }
  if (polar_split && (!hyperbolic || !use_radial_block)) {
    Rcpp::stop("polar_split requires a hyperbolic model with the radial block "
               "enabled.");
  }
  if (!(radial_jump > 0.0) || !std::isfinite(radial_jump)) {
    Rcpp::stop("The radial-block jump must be finite and positive.");
  }

  const double hn_scale_z = pr_sd_z;
  const double hn_scale_w = pr_sd_w;

  arma::vec beta(ncol, arma::fill::randu);
  arma::vec theta(nrow, arma::fill::randu);
  arma::mat z(nrow, ndim);
  arma::mat w(ncol, ndim);

  beta = 4.0 * beta - 2.0;
  theta = 4.0 * theta - 2.0;

  if (hyperbolic) {
    arma::rowvec origin(ndim, arma::fill::zeros);
    arma::rowvec initial_tangent(ndim, arma::fill::zeros);
    hyperboloid_utils::ExpMapResult initial_position;

    origin(0) = 1.0;
    for (int i = 0; i < nrow; ++i) {
      initial_tangent.zeros();
      for (int q = 1; q < ndim; ++q) {
        initial_tangent(q) = 0.1 * R::norm_rand();
      }
      initial_position = hyperboloid_utils::exp_map(origin, initial_tangent,
                                                    max_hyperbolic_radius);
      z.row(i) = initial_position.isValid ? initial_position.point : origin;
    }
    for (int j = 0; j < ncol; ++j) {
      initial_tangent.zeros();
      for (int q = 1; q < ndim; ++q) {
        initial_tangent(q) = 0.1 * R::norm_rand();
      }
      initial_position = hyperboloid_utils::exp_map(origin, initial_tangent,
                                                    max_hyperbolic_radius);
      w.row(j) = initial_position.isValid ? initial_position.point : origin;
    }
  } else {
    z.randn();
    z *= 0.1;
    w.randn();
    w *= 0.1;
  }

  if (has_init_state) {
    const char *required[] = {"beta", "theta", "z", "w"};
    for (const char *nm : required) {
      if (!init_state.containsElementNamed(nm)) {
        Rcpp::stop("init_state requires beta, theta, z, and w.");
      }
    }
    arma::vec b0 = Rcpp::as<arma::vec>(init_state["beta"]);
    arma::vec t0 = Rcpp::as<arma::vec>(init_state["theta"]);
    arma::mat z0 = Rcpp::as<arma::mat>(init_state["z"]);
    arma::mat w0 = Rcpp::as<arma::mat>(init_state["w"]);
    if ((int)b0.n_elem != ncol || (int)t0.n_elem != nrow ||
        (int)z0.n_rows != nrow || (int)z0.n_cols != ndim ||
        (int)w0.n_rows != ncol || (int)w0.n_cols != ndim) {
      Rcpp::stop("init_state dimensions do not match the model.");
    }
    if (!b0.is_finite() || !t0.is_finite() || !z0.is_finite() ||
        !w0.is_finite()) {
      Rcpp::stop("init_state contains non-finite values.");
    }
    if (hyperbolic) {
      auto valid_sheet = [&](const arma::mat &x) {
        for (arma::uword i = 0; i < x.n_rows; ++i) {
          double spatial_sum_of_squares = 0.0;
          for (arma::uword q = 1; q < x.n_cols; ++q) {
            spatial_sum_of_squares += x(i, q) * x(i, q);
          }
          const double expected_time_coordinate =
              std::sqrt(1.0 + spatial_sum_of_squares);
          const double relative_sheet_error =
              std::abs(x(i, 0) - expected_time_coordinate) /
              std::max(1.0, expected_time_coordinate);
          const double r = hyperboloid_utils::radius_hyperboloid(x.row(i));
          if (!(x(i, 0) > 0.0) || !std::isfinite(relative_sheet_error) ||
              relative_sheet_error > 1e-12 || !std::isfinite(r) ||
              r > max_hyperbolic_radius) {
            return false;
          }
        }
        return true;
      };
      if (!valid_sheet(z0) || !valid_sheet(w0)) {
        Rcpp::stop(
            "init_state positions are not on the upper hyperboloid sheet.");
      }
    }
    beta = b0;
    theta = t0;
    z = z0;
    w = w0;
    if (init_state.containsElementNamed("sigma_beta")) {
      pr_sd_beta = Rcpp::as<double>(init_state["sigma_beta"]);
    }
    if (init_state.containsElementNamed("sigma_theta")) {
      pr_sd_theta = Rcpp::as<double>(init_state["sigma_theta"]);
    }
    if (init_state.containsElementNamed("sigma_z")) {
      pr_sd_z = Rcpp::as<double>(init_state["sigma_z"]);
    }
    if (init_state.containsElementNamed("sigma_w")) {
      pr_sd_w = Rcpp::as<double>(init_state["sigma_w"]);
    }
    if (!(pr_sd_beta > 0.0) || !(pr_sd_theta > 0.0) || !(pr_sd_z > 0.0) ||
        !(pr_sd_w > 0.0) || !std::isfinite(pr_sd_beta) ||
        !std::isfinite(pr_sd_theta) || !std::isfinite(pr_sd_z) ||
        !std::isfinite(pr_sd_w)) {
      Rcpp::stop("init_state scales must be finite and positive.");
    }
  }

  auto distance = [&](const auto &a, const auto &b) {
    return hyperbolic ? hyperboloid_utils::distance_hyperboloid(a, b)
                      : euclidean_distance(a, b);
  };
  arma::mat dist(nrow, ncol);
  auto refresh_dist = [&]() {
    for (int i = 0; i < nrow; ++i) {
      for (int j = 0; j < ncol; ++j) {
        dist(i, j) = distance(z.row(i), w.row(j));
      }
    }
  };
  refresh_dist();

  arma::mat current_loglik(nrow, ncol, arma::fill::zeros);
  for (int i = 0; i < nrow; ++i) {
    for (int j = 0; j < ncol; ++j) {
      if (observed_binary(data(i, j))) {
        current_loglik(i, j) =
            bernoulli_loglik(data(i, j), theta(i) + beta(j) - dist(i, j));
      }
    }
  }

  arma::vec jump_beta_node(ncol);
  arma::vec jump_theta_node(nrow);
  arma::vec jump_z_node(nrow);
  arma::vec jump_w_node(ncol);
  arma::vec jump_radial_z_node(nrow);
  arma::vec jump_radial_w_node(ncol);

  jump_beta_node.fill(jump_beta);
  jump_theta_node.fill(jump_theta);
  jump_z_node.fill(jump_z);
  jump_w_node.fill(jump_w);
  jump_radial_z_node.fill(radial_jump);
  jump_radial_w_node.fill(radial_jump);

  arma::vec accept_beta(ncol, fill::zeros);
  arma::vec accept_theta(nrow, fill::zeros);
  arma::vec accept_z(nrow, fill::zeros);
  arma::vec accept_w(ncol, fill::zeros);
  arma::vec accept_radial_z(nrow, fill::zeros);
  arma::vec accept_radial_w(ncol, fill::zeros);

  arma::vec batch_accept_beta(ncol, fill::zeros);
  arma::vec batch_accept_theta(nrow, fill::zeros);
  arma::vec batch_accept_z(nrow, fill::zeros);
  arma::vec batch_accept_w(ncol, fill::zeros);
  arma::vec batch_accept_radial_z(nrow, fill::zeros);
  arma::vec batch_accept_radial_w(ncol, fill::zeros);

  arma::vec batch_count_beta(ncol, fill::zeros);
  arma::vec batch_count_theta(nrow, fill::zeros);
  arma::vec batch_count_z(nrow, fill::zeros);
  arma::vec batch_count_w(ncol, fill::zeros);
  arma::vec batch_count_radial_z(nrow, fill::zeros);
  arma::vec batch_count_radial_w(ncol, fill::zeros);

  double jump_sz = 2.4 / std::sqrt(2.0 * nrow * intrinsic_dim);
  double jump_sw = 2.4 / std::sqrt(2.0 * ncol * intrinsic_dim);
  double accept_sigma_z = 0.0;
  double accept_sigma_w = 0.0;
  double count_sigma_z = 0.0;
  double count_sigma_w = 0.0;
  double batch_accept_sigma_z = 0.0;
  double batch_accept_sigma_w = 0.0;
  double batch_count_sigma_z = 0.0;
  double batch_count_sigma_w = 0.0;
  double invalid_z = 0.0;
  double invalid_w = 0.0;
  double capped_z = 0.0;
  double capped_w = 0.0;
  double capped_z_postburn = 0.0;
  double capped_w_postburn = 0.0;
  double invalid_z_postburn = 0.0;
  double invalid_w_postburn = 0.0;
  double radial_boundary_z = 0.0;
  double radial_boundary_w = 0.0;
  double radial_proposal_z_postburn = 0.0;
  double radial_proposal_w_postburn = 0.0;
  double radial_accept_z_postburn = 0.0;
  double radial_accept_w_postburn = 0.0;
  double radial_boundary_z_postburn = 0.0;
  double radial_boundary_w_postburn = 0.0;

  const int first_saved = ((nburn + nthin - 1) / nthin) * nthin;
  if (first_saved >= niter) {
    Rcpp::stop("The schedule stores no draws.");
  }
  const int nsave = 1 + (niter - 1 - first_saved) / nthin;
  arma::mat out_beta(nsave, ncol);
  arma::mat out_theta(nsave, nrow);
  arma::cube out_z(nsave, nrow, ndim);
  arma::cube out_w(nsave, ncol, ndim);
  arma::vec out_sigma_theta(nsave);
  arma::vec out_sigma_beta(nsave);
  arma::vec out_sigma_z(nsave);
  arma::vec out_sigma_w(nsave);
  arma::vec out_log_posterior(nsave);
  arma::vec out_log_likelihood(nsave);
  arma::vec out_max_radius_z(nsave);
  arma::vec out_max_radius_w(nsave);
  arma::ivec out_iter(nsave);
  int save_index = 0;

  auto adapt_vec = [](arma::vec &jump, const arma::vec &acc,
                      const arma::vec &den, double weight, double target) {
    for (arma::uword i = 0; i < jump.n_elem; ++i) {
      if (den(i) <= 0.0) {
        continue;
      }
      const double log_jump =
          std::log(jump(i)) + weight * (acc(i) / den(i) - target);
      jump(i) = std::exp(std::max(-8.0, std::min(1.5, log_jump)));
    }
  };

  auto update_scale = [&](double &sigma, double step, double sum_of_squares,
                          double parameter_count, double halfnormal_scale,
                          double &acceptance_count, double &proposal_count) {
    const double log_sigma_current = std::log(sigma);
    const double log_sigma_candidate =
        log_sigma_current + step * R::norm_rand();
    const double sigma_candidate = std::exp(log_sigma_candidate);

    auto log_conditional = [&](double candidate) {
      return -parameter_count * std::log(candidate) -
             sum_of_squares / (2.0 * candidate * candidate) -
             candidate * candidate /
                 (2.0 * halfnormal_scale * halfnormal_scale);
    };

    proposal_count += 1.0;
    const double log_acceptance_ratio = log_conditional(sigma_candidate) -
                                        log_conditional(sigma) +
                                        log_sigma_candidate - log_sigma_current;
    if (std::log(R::unif_rand()) < log_acceptance_ratio) {
      sigma = sigma_candidate;
      acceptance_count += 1.0;
      return true;
    }
    return false;
  };

  // Reusable workspaces avoid repeated dynamic allocation in nodewise loops.
  arma::rowvec candidate_z(ndim);
  arma::rowvec candidate_w(ndim);
  arma::rowvec candidate_distance_z(ncol);
  arma::vec candidate_distance_w(nrow);
  arma::rowvec radial_tangent_z(intrinsic_dim);
  arma::rowvec radial_tangent_w(intrinsic_dim);
  arma::rowvec radial_position_z(ndim);
  arma::rowvec radial_position_w(ndim);
  arma::rowvec radial_distance_z(ncol);
  arma::vec radial_distance_w(nrow);
  arma::rowvec radial_direction_current(intrinsic_dim);
  arma::vec radial_partner_radius_w(ncol);
  arma::vec radial_partner_sinh_w(ncol);
  arma::mat radial_partner_direction_w(ncol, intrinsic_dim);
  arma::vec radial_partner_radius_z(nrow);
  arma::vec radial_partner_sinh_z(nrow);
  arma::mat radial_partner_direction_z(nrow, intrinsic_dim);
  arma::rowvec candidate_loglik_row(ncol, arma::fill::zeros);
  arma::vec candidate_loglik_column(nrow, arma::fill::zeros);
  arma::mat proposal_basis(ndim, intrinsic_dim);
  arma::rowvec proposal_coefficients(intrinsic_dim);
  arma::rowvec proposal_ambient_tangent(ndim);
  arma::rowvec proposal_snapped_point(ndim);
  arma::rowvec proposal_projected_tangent(ndim);

  for (int iter = 0; iter < niter; ++iter) {
    if (iter % 10 == 0) {
      Rcpp::checkUserInterrupt();
    }

    for (int j = 0; j < ncol; ++j) {
      const double candidate = R::rnorm(beta(j), jump_beta_node(j));
      double log_acceptance_ratio = R::dnorm4(candidate, 0.0, pr_sd_beta, 1) -
                                    R::dnorm4(beta(j), 0.0, pr_sd_beta, 1);

      for (int i = 0; i < nrow; ++i) {
        if (!observed_binary(data(i, j))) {
          candidate_loglik_column(i) = 0.0;
          continue;
        }
        candidate_loglik_column(i) =
            bernoulli_loglik(data(i, j), theta(i) + candidate - dist(i, j));
        log_acceptance_ratio +=
            candidate_loglik_column(i) - current_loglik(i, j);
      }

      if (std::log(R::unif_rand()) < log_acceptance_ratio) {
        beta(j) = candidate;
        current_loglik.col(j) = candidate_loglik_column;
        accept_beta(j) += 1.0 / niter;
        if (iter < nburn) {
          batch_accept_beta(j) += 1.0;
        }
      }
      if (iter < nburn) {
        batch_count_beta(j) += 1.0;
      }
    }

    for (int i = 0; i < nrow; ++i) {
      const double candidate = R::rnorm(theta(i), jump_theta_node(i));
      double log_acceptance_ratio = R::dnorm4(candidate, 0.0, pr_sd_theta, 1) -
                                    R::dnorm4(theta(i), 0.0, pr_sd_theta, 1);

      for (int j = 0; j < ncol; ++j) {
        if (!observed_binary(data(i, j))) {
          candidate_loglik_row(j) = 0.0;
          continue;
        }
        candidate_loglik_row(j) =
            bernoulli_loglik(data(i, j), candidate + beta(j) - dist(i, j));
        log_acceptance_ratio += candidate_loglik_row(j) - current_loglik(i, j);
      }

      if (std::log(R::unif_rand()) < log_acceptance_ratio) {
        theta(i) = candidate;
        current_loglik.row(i) = candidate_loglik_row;
        accept_theta(i) += 1.0 / niter;
        if (iter < nburn) {
          batch_accept_theta(i) += 1.0;
        }
      }
      if (iter < nburn) {
        batch_count_theta(i) += 1.0;
      }
    }

    for (int i = 0; i < nrow; ++i) {
      bool valid = true;
      bool proposal_hit_radius_cap = false;

      if (hyperbolic) {
        if (polar_split &&
            hyperboloid_utils::radius_hyperboloid(z.row(i)) > 0.0) {
          valid = hyperboloid_utils::propose_angular_fixed_radius_into(
              z.row(i), jump_z_node(i), proposal_coefficients, candidate_z);
        } else if (stable_intrinsic_proposal) {
          valid = hyperboloid_utils::propose_intrinsic_wrapped_normal_into(
              z.row(i), jump_z_node(i), proposal_coefficients, candidate_z,
              proposal_hit_radius_cap, max_hyperbolic_radius);
        } else {
          valid = hyperboloid_utils::propose_wrapped_normal_into(
              z.row(i), jump_z_node(i), proposal_basis, proposal_coefficients,
              proposal_ambient_tangent, proposal_snapped_point,
              proposal_projected_tangent, candidate_z, proposal_hit_radius_cap,
              max_hyperbolic_radius);
        }
      } else {
        candidate_z = z.row(i);
        for (int q = 0; q < ndim; ++q) {
          candidate_z(q) += jump_z_node(i) * R::norm_rand();
        }
      }

      if (valid) {
        double log_acceptance_ratio = 0.0;
        if (hyperbolic) {
          log_acceptance_ratio +=
              hyperboloid_utils::log_wrapped_normal_prior_density(
                  candidate_z, pr_sd_z * pr_sd_z) -
              hyperboloid_utils::log_wrapped_normal_prior_density(
                  z.row(i), pr_sd_z * pr_sd_z);
        } else {
          for (int q = 0; q < ndim; ++q) {
            log_acceptance_ratio += R::dnorm4(candidate_z(q), 0.0, pr_sd_z, 1) -
                                    R::dnorm4(z(i, q), 0.0, pr_sd_z, 1);
          }
        }

        for (int j = 0; j < ncol; ++j) {
          candidate_distance_z(j) = distance(candidate_z, w.row(j));
          if (!observed_binary(data(i, j))) {
            candidate_loglik_row(j) = 0.0;
            continue;
          }
          candidate_loglik_row(j) = bernoulli_loglik(
              data(i, j), theta(i) + beta(j) - candidate_distance_z(j));
          log_acceptance_ratio +=
              candidate_loglik_row(j) - current_loglik(i, j);
        }

        if (std::log(R::unif_rand()) < log_acceptance_ratio) {
          z.row(i) = candidate_z;
          dist.row(i) = candidate_distance_z;
          current_loglik.row(i) = candidate_loglik_row;
          accept_z(i) += 1.0 / niter;
          if (iter < nburn) {
            batch_accept_z(i) += 1.0;
          }
        }
      } else {
        invalid_z += 1.0;
        if (iter >= nburn) {
          invalid_z_postburn += 1.0;
        }
        if (proposal_hit_radius_cap) {
          capped_z += 1.0;
          if (iter >= nburn) {
            capped_z_postburn += 1.0;
          }
        }
      }
      if (iter < nburn) {
        batch_count_z(i) += 1.0;
      }
    }

    for (int j = 0; j < ncol; ++j) {
      bool valid = true;
      bool proposal_hit_radius_cap = false;

      if (hyperbolic) {
        if (polar_split &&
            hyperboloid_utils::radius_hyperboloid(w.row(j)) > 0.0) {
          valid = hyperboloid_utils::propose_angular_fixed_radius_into(
              w.row(j), jump_w_node(j), proposal_coefficients, candidate_w);
        } else if (stable_intrinsic_proposal) {
          valid = hyperboloid_utils::propose_intrinsic_wrapped_normal_into(
              w.row(j), jump_w_node(j), proposal_coefficients, candidate_w,
              proposal_hit_radius_cap, max_hyperbolic_radius);
        } else {
          valid = hyperboloid_utils::propose_wrapped_normal_into(
              w.row(j), jump_w_node(j), proposal_basis, proposal_coefficients,
              proposal_ambient_tangent, proposal_snapped_point,
              proposal_projected_tangent, candidate_w, proposal_hit_radius_cap,
              max_hyperbolic_radius);
        }
      } else {
        candidate_w = w.row(j);
        for (int q = 0; q < ndim; ++q) {
          candidate_w(q) += jump_w_node(j) * R::norm_rand();
        }
      }

      if (valid) {
        double log_acceptance_ratio = 0.0;
        if (hyperbolic) {
          log_acceptance_ratio +=
              hyperboloid_utils::log_wrapped_normal_prior_density(
                  candidate_w, pr_sd_w * pr_sd_w) -
              hyperboloid_utils::log_wrapped_normal_prior_density(
                  w.row(j), pr_sd_w * pr_sd_w);
        } else {
          for (int q = 0; q < ndim; ++q) {
            log_acceptance_ratio += R::dnorm4(candidate_w(q), 0.0, pr_sd_w, 1) -
                                    R::dnorm4(w(j, q), 0.0, pr_sd_w, 1);
          }
        }

        for (int i = 0; i < nrow; ++i) {
          candidate_distance_w(i) = distance(z.row(i), candidate_w);
          if (!observed_binary(data(i, j))) {
            candidate_loglik_column(i) = 0.0;
            continue;
          }
          candidate_loglik_column(i) = bernoulli_loglik(
              data(i, j), theta(i) + beta(j) - candidate_distance_w(i));
          log_acceptance_ratio +=
              candidate_loglik_column(i) - current_loglik(i, j);
        }

        if (std::log(R::unif_rand()) < log_acceptance_ratio) {
          w.row(j) = candidate_w;
          dist.col(j) = candidate_distance_w;
          current_loglik.col(j) = candidate_loglik_column;
          accept_w(j) += 1.0 / niter;
          if (iter < nburn) {
            batch_accept_w(j) += 1.0;
          }
        }
      } else {
        invalid_w += 1.0;
        if (iter >= nburn) {
          invalid_w_postburn += 1.0;
        }
        if (proposal_hit_radius_cap) {
          capped_w += 1.0;
          if (iter >= nburn) {
            capped_w_postburn += 1.0;
          }
        }
      }
      if (iter < nburn) {
        batch_count_w(j) += 1.0;
      }
    }

    if (use_radial_block) {
      constexpr double minimum_radius = 0.0;

      for (int j = 0; j < ncol; ++j) {
        double spatial_norm_squared = 0.0;
        for (int q = 1; q < ndim; ++q) {
          spatial_norm_squared += w(j, q) * w(j, q);
        }
        const double spatial_norm = std::sqrt(spatial_norm_squared);
        radial_partner_sinh_w(j) = spatial_norm;
        radial_partner_radius_w(j) = std::asinh(spatial_norm);
        if (spatial_norm > minimum_radius) {
          for (int q = 0; q < intrinsic_dim; ++q) {
            radial_partner_direction_w(j, q) = w(j, q + 1) / spatial_norm;
          }
        } else {
          radial_partner_direction_w.row(j).zeros();
        }
      }

      for (int i = 0; i < nrow; ++i) {
        radial_tangent_z = hyperboloid_utils::log_origin(z.row(i));
        const double radius_current = arma::norm(radial_tangent_z, 2);
        if (!(radius_current > minimum_radius) ||
            radius_current > max_hyperbolic_radius ||
            !std::isfinite(radius_current)) {
          continue;
        }
        if (iter < nburn) {
          batch_count_radial_z(i) += 1.0;
        } else {
          radial_proposal_z_postburn += 1.0;
        }

        const double log_radius_current = std::log(radius_current);
        const double log_radius_candidate =
            log_radius_current + jump_radial_z_node(i) * R::norm_rand();
        const double radius_candidate = std::exp(log_radius_candidate);
        if (!(radius_candidate > minimum_radius) ||
            radius_candidate > max_hyperbolic_radius ||
            !std::isfinite(radius_candidate)) {
          radial_boundary_z += 1.0;
          if (iter >= nburn) {
            radial_boundary_z_postburn += 1.0;
          }
          continue;
        }

        radial_position_z = hyperboloid_utils::exp_origin(
            (radius_candidate / radius_current) * radial_tangent_z,
            max_hyperbolic_radius);
        if (!radial_position_z.is_finite()) {
          continue;
        }
        radial_direction_current = radial_tangent_z / radius_current;
        const double sinh_radius_candidate = std::sinh(radius_candidate);

        const double theta_candidate =
            theta(i) - radius_current + radius_candidate;
        double log_acceptance_ratio =
            R::dnorm4(theta_candidate, 0.0, pr_sd_theta, 1) -
            R::dnorm4(theta(i), 0.0, pr_sd_theta, 1) -
            (radius_candidate * radius_candidate -
             radius_current * radius_current) /
                (2.0 * pr_sd_z * pr_sd_z) +
            intrinsic_dim * (log_radius_candidate - log_radius_current);

        for (int j = 0; j < ncol; ++j) {
          radial_distance_z(j) =
              hyperboloid_utils::geodesic_distance_cached_polar(
                  radius_candidate, sinh_radius_candidate,
                  radial_direction_current, radial_partner_radius_w(j),
                  radial_partner_sinh_w(j), radial_partner_direction_w.row(j));
          if (!observed_binary(data(i, j))) {
            candidate_loglik_row(j) = 0.0;
            continue;
          }
          candidate_loglik_row(j) = bernoulli_loglik(
              data(i, j), theta_candidate + beta(j) - radial_distance_z(j));
          log_acceptance_ratio +=
              candidate_loglik_row(j) - current_loglik(i, j);
        }

        if (std::log(R::unif_rand()) < log_acceptance_ratio) {
          z.row(i) = radial_position_z;
          theta(i) = theta_candidate;
          dist.row(i) = radial_distance_z;
          current_loglik.row(i) = candidate_loglik_row;
          accept_radial_z(i) += 1.0 / niter;
          if (iter < nburn) {
            batch_accept_radial_z(i) += 1.0;
          } else {
            radial_accept_z_postburn += 1.0;
          }
        }
      }

      for (int i = 0; i < nrow; ++i) {
        double spatial_norm_squared = 0.0;
        for (int q = 1; q < ndim; ++q) {
          spatial_norm_squared += z(i, q) * z(i, q);
        }
        const double spatial_norm = std::sqrt(spatial_norm_squared);
        radial_partner_sinh_z(i) = spatial_norm;
        radial_partner_radius_z(i) = std::asinh(spatial_norm);
        if (spatial_norm > minimum_radius) {
          for (int q = 0; q < intrinsic_dim; ++q) {
            radial_partner_direction_z(i, q) = z(i, q + 1) / spatial_norm;
          }
        } else {
          radial_partner_direction_z.row(i).zeros();
        }
      }

      for (int j = 0; j < ncol; ++j) {
        radial_tangent_w = hyperboloid_utils::log_origin(w.row(j));
        const double radius_current = arma::norm(radial_tangent_w, 2);
        if (!(radius_current > minimum_radius) ||
            radius_current > max_hyperbolic_radius ||
            !std::isfinite(radius_current)) {
          continue;
        }
        if (iter < nburn) {
          batch_count_radial_w(j) += 1.0;
        } else {
          radial_proposal_w_postburn += 1.0;
        }

        const double log_radius_current = std::log(radius_current);
        const double log_radius_candidate =
            log_radius_current + jump_radial_w_node(j) * R::norm_rand();
        const double radius_candidate = std::exp(log_radius_candidate);
        if (!(radius_candidate > minimum_radius) ||
            radius_candidate > max_hyperbolic_radius ||
            !std::isfinite(radius_candidate)) {
          radial_boundary_w += 1.0;
          if (iter >= nburn) {
            radial_boundary_w_postburn += 1.0;
          }
          continue;
        }

        radial_position_w = hyperboloid_utils::exp_origin(
            (radius_candidate / radius_current) * radial_tangent_w,
            max_hyperbolic_radius);
        if (!radial_position_w.is_finite()) {
          continue;
        }
        radial_direction_current = radial_tangent_w / radius_current;
        const double sinh_radius_candidate = std::sinh(radius_candidate);

        const double beta_candidate =
            beta(j) - radius_current + radius_candidate;
        double log_acceptance_ratio =
            R::dnorm4(beta_candidate, 0.0, pr_sd_beta, 1) -
            R::dnorm4(beta(j), 0.0, pr_sd_beta, 1) -
            (radius_candidate * radius_candidate -
             radius_current * radius_current) /
                (2.0 * pr_sd_w * pr_sd_w) +
            intrinsic_dim * (log_radius_candidate - log_radius_current);

        for (int i = 0; i < nrow; ++i) {
          radial_distance_w(i) =
              hyperboloid_utils::geodesic_distance_cached_polar(
                  radius_candidate, sinh_radius_candidate,
                  radial_direction_current, radial_partner_radius_z(i),
                  radial_partner_sinh_z(i), radial_partner_direction_z.row(i));
          if (!observed_binary(data(i, j))) {
            candidate_loglik_column(i) = 0.0;
            continue;
          }
          candidate_loglik_column(i) = bernoulli_loglik(
              data(i, j), theta(i) + beta_candidate - radial_distance_w(i));
          log_acceptance_ratio +=
              candidate_loglik_column(i) - current_loglik(i, j);
        }

        if (std::log(R::unif_rand()) < log_acceptance_ratio) {
          w.row(j) = radial_position_w;
          beta(j) = beta_candidate;
          dist.col(j) = radial_distance_w;
          current_loglik.col(j) = candidate_loglik_column;
          accept_radial_w(j) += 1.0 / niter;
          if (iter < nburn) {
            batch_accept_radial_w(j) += 1.0;
          } else {
            radial_accept_w_postburn += 1.0;
          }
        }
      }
    }

    // Inverse-gamma Gibbs updates for additive-effect variances.
    double sum_of_squares = arma::dot(theta, theta);
    pr_sd_theta = std::sqrt(2.0 * (pr_b_theta + 0.5 * sum_of_squares) /
                            R::rchisq(2.0 * pr_a_theta + nrow));
    if (estimate_beta_sd) {
      sum_of_squares = arma::dot(beta, beta);
      pr_sd_beta = std::sqrt(2.0 * (pr_b_beta + 0.5 * sum_of_squares) /
                             R::rchisq(2.0 * pr_a_beta + ncol));
    }

    double sum_of_squares_z = 0.0;
    double sum_of_squares_w = 0.0;
    for (int i = 0; i < nrow; ++i) {
      const double radius = latent_radius(z.row(i), hyperbolic);
      sum_of_squares_z += radius * radius;
    }
    for (int j = 0; j < ncol; ++j) {
      const double radius = latent_radius(w.row(j), hyperbolic);
      sum_of_squares_w += radius * radius;
    }

    const bool scale_z_accepted =
        update_scale(pr_sd_z, jump_sz, sum_of_squares_z, nrow * intrinsic_dim,
                     hn_scale_z, accept_sigma_z, count_sigma_z);
    const bool scale_w_accepted =
        update_scale(pr_sd_w, jump_sw, sum_of_squares_w, ncol * intrinsic_dim,
                     hn_scale_w, accept_sigma_w, count_sigma_w);

    if (iter < nburn) {
      batch_count_sigma_z += 1.0;
      batch_count_sigma_w += 1.0;
      if (scale_z_accepted) {
        batch_accept_sigma_z += 1.0;
      }
      if (scale_w_accepted) {
        batch_accept_sigma_w += 1.0;
      }
    }

    if (adaptation && iter < nburn && (iter + 1) % adapt_interval == 0) {
      const double batch = (iter + 1.0) / adapt_interval;
      const double weight = 1.0 / std::pow(batch, adapt_exponent);

      adapt_vec(jump_beta_node, batch_accept_beta, batch_count_beta, weight,
                target_accept_1d);
      adapt_vec(jump_theta_node, batch_accept_theta, batch_count_theta, weight,
                target_accept_1d);

      const double latent_target =
          intrinsic_dim >= 5
              ? 0.234
              : 0.44 - (0.44 - 0.234) * (intrinsic_dim - 1.0) / 4.0;
      adapt_vec(jump_z_node, batch_accept_z, batch_count_z, weight,
                latent_target);
      adapt_vec(jump_w_node, batch_accept_w, batch_count_w, weight,
                latent_target);

      if (use_radial_block) {
        adapt_vec(jump_radial_z_node, batch_accept_radial_z,
                  batch_count_radial_z, weight, target_accept_1d);
        adapt_vec(jump_radial_w_node, batch_accept_radial_w,
                  batch_count_radial_w, weight, target_accept_1d);
      }

      if (batch_count_sigma_z > 0.0) {
        const double log_jump =
            std::log(jump_sz) +
            weight *
                (batch_accept_sigma_z / batch_count_sigma_z - target_accept_1d);
        jump_sz = std::exp(std::max(-8.0, std::min(3.0, log_jump)));
      }
      if (batch_count_sigma_w > 0.0) {
        const double log_jump =
            std::log(jump_sw) +
            weight *
                (batch_accept_sigma_w / batch_count_sigma_w - target_accept_1d);
        jump_sw = std::exp(std::max(-8.0, std::min(3.0, log_jump)));
      }

      batch_accept_beta.zeros();
      batch_accept_theta.zeros();
      batch_accept_z.zeros();
      batch_accept_w.zeros();
      batch_accept_radial_z.zeros();
      batch_accept_radial_w.zeros();
      batch_count_beta.zeros();
      batch_count_theta.zeros();
      batch_count_z.zeros();
      batch_count_w.zeros();
      batch_count_radial_z.zeros();
      batch_count_radial_w.zeros();
      batch_accept_sigma_z = 0.0;
      batch_accept_sigma_w = 0.0;
      batch_count_sigma_z = 0.0;
      batch_count_sigma_w = 0.0;
    }

    if (iter >= nburn && iter % nthin == 0) {
      out_beta.row(save_index) = beta.t();
      out_theta.row(save_index) = theta.t();

      for (int i = 0; i < nrow; ++i) {
        for (int q = 0; q < ndim; ++q) {
          out_z(save_index, i, q) = z(i, q);
        }
      }
      for (int j = 0; j < ncol; ++j) {
        for (int q = 0; q < ndim; ++q) {
          out_w(save_index, j, q) = w(j, q);
        }
      }

      out_sigma_theta(save_index) = pr_sd_theta;
      out_sigma_beta(save_index) = pr_sd_beta;
      out_sigma_z(save_index) = pr_sd_z;
      out_sigma_w(save_index) = pr_sd_w;
      out_iter(save_index) = iter;
      out_max_radius_z(save_index) = 0.0;
      out_max_radius_w(save_index) = 0.0;

      for (int i = 0; i < nrow; ++i) {
        out_max_radius_z(save_index) = std::max(
            out_max_radius_z(save_index), latent_radius(z.row(i), hyperbolic));
      }
      for (int j = 0; j < ncol; ++j) {
        out_max_radius_w(save_index) = std::max(
            out_max_radius_w(save_index), latent_radius(w.row(j), hyperbolic));
      }

      double log_posterior = 0.0;
      for (int i = 0; i < nrow; ++i) {
        log_posterior += R::dnorm4(theta(i), 0.0, pr_sd_theta, 1);
      }
      for (int j = 0; j < ncol; ++j) {
        log_posterior += R::dnorm4(beta(j), 0.0, pr_sd_beta, 1);
      }

      if (hyperbolic) {
        for (int i = 0; i < nrow; ++i) {
          log_posterior += hyperboloid_utils::log_wrapped_normal_prior_density(
              z.row(i), pr_sd_z * pr_sd_z);
        }
        for (int j = 0; j < ncol; ++j) {
          log_posterior += hyperboloid_utils::log_wrapped_normal_prior_density(
              w.row(j), pr_sd_w * pr_sd_w);
        }
      } else {
        for (int i = 0; i < nrow; ++i) {
          for (int q = 0; q < ndim; ++q) {
            log_posterior += R::dnorm4(z(i, q), 0.0, pr_sd_z, 1);
          }
        }
        for (int j = 0; j < ncol; ++j) {
          for (int q = 0; q < ndim; ++q) {
            log_posterior += R::dnorm4(w(j, q), 0.0, pr_sd_w, 1);
          }
        }
      }

      log_posterior +=
          log_invgamma_variance_as_sd(pr_sd_theta, pr_a_theta, pr_b_theta);
      if (estimate_beta_sd) {
        log_posterior +=
            log_invgamma_variance_as_sd(pr_sd_beta, pr_a_beta, pr_b_beta);
      }
      log_posterior += log_halfnormal_sd(pr_sd_z, hn_scale_z);
      log_posterior += log_halfnormal_sd(pr_sd_w, hn_scale_w);

      log_posterior += arma::accu(current_loglik);

      out_log_posterior(save_index) = log_posterior;
      out_log_likelihood(save_index) = arma::accu(current_loglik);
      ++save_index;
    }
  }

  const double accept_rate_sigma_z =
      count_sigma_z > 0.0 ? accept_sigma_z / count_sigma_z : NA_REAL;
  const double accept_rate_sigma_w =
      count_sigma_w > 0.0 ? accept_sigma_w / count_sigma_w : NA_REAL;

  Rcpp::List out = Rcpp::List::create(
      Rcpp::Named("beta") = out_beta, Rcpp::Named("theta") = out_theta,
      Rcpp::Named("z") = out_z, Rcpp::Named("w") = out_w,
      Rcpp::Named("sigma_theta") = out_sigma_theta,
      Rcpp::Named("sigma_beta") = out_sigma_beta,
      Rcpp::Named("sigma_z") = out_sigma_z,
      Rcpp::Named("sigma_w") = out_sigma_w,
      Rcpp::Named("gamma") = arma::vec(nsave, arma::fill::ones),
      Rcpp::Named("log_posterior") = out_log_posterior,
      Rcpp::Named("log_likelihood") = out_log_likelihood,
      Rcpp::Named("retained_iteration") = out_iter,
      Rcpp::Named("max_radius_z") = out_max_radius_z,
      Rcpp::Named("max_radius_w") = out_max_radius_w,
      Rcpp::Named("accept_beta") = accept_beta,
      Rcpp::Named("accept_theta") = accept_theta,
      Rcpp::Named("accept_z") = accept_z, Rcpp::Named("accept_w") = accept_w,
      Rcpp::Named("jump_node_beta") = jump_beta_node,
      Rcpp::Named("jump_node_theta") = jump_theta_node,
      Rcpp::Named("jump_node_z") = jump_z_node,
      Rcpp::Named("jump_node_w") = jump_w_node,
      Rcpp::Named("jump_sigma_z") = jump_sz,
      Rcpp::Named("jump_sigma_w") = jump_sw,
      Rcpp::Named("accept_sigma_z") = accept_rate_sigma_z,
      Rcpp::Named("accept_sigma_w") = accept_rate_sigma_w,
      Rcpp::Named("invalid_position_rate_z") =
          invalid_z / (niter * static_cast<double>(nrow)),
      Rcpp::Named("invalid_position_rate_w") =
          invalid_w / (niter * static_cast<double>(ncol)),
      Rcpp::Named("capped_proposal_rate_z") =
          capped_z / (niter * static_cast<double>(nrow)),
      Rcpp::Named("capped_proposal_rate_w") =
          capped_w / (niter * static_cast<double>(ncol)),
      Rcpp::Named("postburn_capped_proposal_rate_z") =
          capped_z_postburn / ((niter - nburn) * static_cast<double>(nrow)),
      Rcpp::Named("postburn_capped_proposal_rate_w") =
          capped_w_postburn / ((niter - nburn) * static_cast<double>(ncol)),
      Rcpp::Named("numerical_failure_rate_z") =
          (invalid_z - capped_z) / (niter * static_cast<double>(nrow)),
      Rcpp::Named("numerical_failure_rate_w") =
          (invalid_w - capped_w) / (niter * static_cast<double>(ncol)),
      Rcpp::Named("postburn_numerical_failure_rate_z") =
          (invalid_z_postburn - capped_z_postburn) /
          ((niter - nburn) * static_cast<double>(nrow)),
      Rcpp::Named("postburn_numerical_failure_rate_w") =
          (invalid_w_postburn - capped_w_postburn) /
          ((niter - nburn) * static_cast<double>(ncol)),
      Rcpp::Named("radial_boundary_rate_z") =
          use_radial_block
              ? radial_boundary_z / (niter * static_cast<double>(nrow))
              : 0.0,
      Rcpp::Named("radial_boundary_rate_w") =
          use_radial_block
              ? radial_boundary_w / (niter * static_cast<double>(ncol))
              : 0.0,
      Rcpp::Named("postburn_radial_acceptance_rate_z") =
          use_radial_block && radial_proposal_z_postburn > 0.0
              ? radial_accept_z_postburn / radial_proposal_z_postburn
              : 0.0,
      Rcpp::Named("postburn_radial_acceptance_rate_w") =
          use_radial_block && radial_proposal_w_postburn > 0.0
              ? radial_accept_w_postburn / radial_proposal_w_postburn
              : 0.0,
      Rcpp::Named("postburn_radial_boundary_rate_z") =
          use_radial_block && radial_proposal_z_postburn > 0.0
              ? radial_boundary_z_postburn / radial_proposal_z_postburn
              : 0.0,
      Rcpp::Named("postburn_radial_boundary_rate_w") =
          use_radial_block && radial_proposal_w_postburn > 0.0
              ? radial_boundary_w_postburn / radial_proposal_w_postburn
              : 0.0,
      Rcpp::Named("hyperbolic") = hyperbolic, Rcpp::Named("ndim") = ndim,
      Rcpp::Named("latent_dim") = intrinsic_dim,
      Rcpp::Named("fix_gamma") = true,
      Rcpp::Named("fix_beta_sd") = !estimate_beta_sd,
      Rcpp::Named("fix_theta_sd") = false, Rcpp::Named("sample_zw_sd") = 2,
      Rcpp::Named("use_radial_block") = use_radial_block,
      Rcpp::Named("polar_split") = polar_split,
      Rcpp::Named("stable_intrinsic_proposal") = stable_intrinsic_proposal,
      Rcpp::Named("max_hyperbolic_radius") = max_hyperbolic_radius,
      Rcpp::Named("continued_from_state") = has_init_state,
      Rcpp::Named("adaptation") = adaptation,
      Rcpp::Named("adapt_interval") = adapt_interval,
      Rcpp::Named("target_accept_1d") = target_accept_1d,
      Rcpp::Named("adapt_exponent") = adapt_exponent);

  if (use_radial_block) {
    out["accept_radial_z"] = accept_radial_z;
    out["accept_radial_w"] = accept_radial_w;
    out["jump_radial_z_final"] = jump_radial_z_node;
    out["jump_radial_w_final"] = jump_radial_w_node;
  }
  return out;
}
