#ifndef LSIRM_PAPER_HYPERBOLOID_UTILS_H
#define LSIRM_PAPER_HYPERBOLOID_UTILS_H

#include <RcppArmadillo.h>
#include <cmath>
#include <limits>

namespace hyperboloid_utils {

constexpr double default_max_hyperbolic_radius = 13.0;

template <typename VectorX, typename VectorY>
inline double lorentz_dot(const VectorX &x, const VectorY &y) {
  if (x.n_elem != y.n_elem || x.n_elem < 2) {
    Rcpp::stop("lorentz_dot requires equal vectors of length at least two.");
  }
  double ans = -x(0) * y(0);
  for (arma::uword j = 1; j < x.n_elem; ++j) {
    ans += x(j) * y(j);
  }
  return ans;
}

template <typename VectorX, typename VectorY>
inline double geodesic_distance_polar(const VectorX &x, const VectorY &y) {
  if (x.n_elem != y.n_elem || x.n_elem < 2) {
    Rcpp::stop("Polar distance requires equal vectors of length at least two.");
  }
  const arma::uword d = x.n_elem - 1;
  double spatial_norm_x_squared = 0.0;
  double spatial_norm_y_squared = 0.0;
  for (arma::uword q = 1; q <= d; ++q) {
    spatial_norm_x_squared += x(q) * x(q);
    spatial_norm_y_squared += y(q) * y(q);
  }
  const double spatial_norm_x = std::sqrt(spatial_norm_x_squared);
  const double spatial_norm_y = std::sqrt(spatial_norm_y_squared);
  const double radius_x = std::asinh(spatial_norm_x);
  const double radius_y = std::asinh(spatial_norm_y);
  if (!(spatial_norm_x > 0.0) || !(spatial_norm_y > 0.0)) {
    return std::abs(radius_x - radius_y);
  }
  double direction_difference_squared = 0.0;
  for (arma::uword q = 1; q <= d; ++q) {
    const double difference = x(q) / spatial_norm_x - y(q) / spatial_norm_y;
    direction_difference_squared += difference * difference;
  }
  double half = 0.5 * std::sqrt(direction_difference_squared);
  half = std::max(0.0, std::min(1.0, half));
  const double radial_term = std::sinh(0.5 * (radius_x - radius_y));
  return 2.0 *
         std::asinh(std::sqrt(radial_term * radial_term +
                              spatial_norm_x * spatial_norm_y * half * half));
}

// Exact half-distance identity when polar quantities are already available.
// This is useful for a radial proposal: the proposed direction is unchanged,
// and every partner radius and direction can be reused throughout the sweep.
template <typename DirectionX, typename DirectionY>
inline double geodesic_distance_cached_polar(
    double radius_x, double sinh_radius_x, const DirectionX &direction_x,
    double radius_y, double sinh_radius_y, const DirectionY &direction_y) {
  if (direction_x.n_elem != direction_y.n_elem || direction_x.n_elem < 1) {
    Rcpp::stop("Cached polar distance requires equal nonempty directions.");
  }
  if (!(radius_x >= 0.0) || !(radius_y >= 0.0) || !(sinh_radius_x >= 0.0) ||
      !(sinh_radius_y >= 0.0) || !std::isfinite(radius_x) ||
      !std::isfinite(radius_y) || !std::isfinite(sinh_radius_x) ||
      !std::isfinite(sinh_radius_y)) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  if (!(sinh_radius_x > 0.0) || !(sinh_radius_y > 0.0)) {
    return std::abs(radius_x - radius_y);
  }
  double direction_difference_squared = 0.0;
  for (arma::uword q = 0; q < direction_x.n_elem; ++q) {
    const double difference = direction_x(q) - direction_y(q);
    direction_difference_squared += difference * difference;
  }
  double half_squared = 0.25 * direction_difference_squared;
  half_squared = std::max(0.0, std::min(1.0, half_squared));
  const double radial_term = std::sinh(0.5 * (radius_x - radius_y));
  return 2.0 *
         std::asinh(std::sqrt(radial_term * radial_term +
                              sinh_radius_x * sinh_radius_y * half_squared));
}

template <typename Vector>
inline void fill_lorentz_orthonormal_basis(const Vector &x, arma::mat &basis) {
  if (x.n_elem < 2) {
    Rcpp::stop("A hyperboloid point needs at least two coordinates.");
  }
  const arma::uword d = x.n_elem - 1;
  if (basis.n_rows != x.n_elem || basis.n_cols != d) {
    basis.set_size(x.n_elem, d);
  }
  basis.zeros();

  double spatial_norm_squared = 0.0;
  for (arma::uword q = 1; q < x.n_elem; ++q) {
    spatial_norm_squared += x(q) * x(q);
  }
  const double spatial_norm = std::sqrt(spatial_norm_squared);
  if (!(spatial_norm > 1e-12)) {
    for (arma::uword j = 0; j < d; ++j) {
      basis(j + 1, j) = 1.0;
    }
    return;
  }
  for (arma::uword j = 0; j < d; ++j) {
    const double direction_j = x(j + 1) / spatial_norm;
    basis(0, j) = x(j + 1);
    for (arma::uword k = 0; k < d; ++k) {
      const double direction_k = x(k + 1) / spatial_norm;
      basis(k + 1, j) =
          (j == k ? 1.0 : 0.0) + (x(0) - 1.0) * direction_j * direction_k;
    }
  }
}

template <typename Vector>
inline arma::mat lorentz_orthonormal_basis(const Vector &x) {
  arma::mat basis;
  fill_lorentz_orthonormal_basis(x, basis);
  return basis;
}

inline arma::rowvec project_to_tangent_space(const arma::rowvec &x,
                                             const arma::rowvec &u) {
  if (x.n_elem != u.n_elem) {
    Rcpp::stop("Tangent projection dimension mismatch.");
  }
  return u + lorentz_dot(u, x) * x;
}

struct ExpMapResult {
  arma::rowvec point;
  bool isValid;
  bool hitRadiusCap = false;
};

inline ExpMapResult
exp_map(const arma::rowvec &x, const arma::rowvec &vin,
        double max_hyperbolic_radius = default_max_hyperbolic_radius) {
  ExpMapResult out{x, false, false};
  if (x.n_elem != vin.n_elem || x.n_elem < 2) {
    Rcpp::stop("exp_map dimension mismatch.");
  }
  arma::rowvec xs = x;
  double sp2 = 0.0;
  for (arma::uword j = 1; j < xs.n_elem; ++j) {
    sp2 += xs(j) * xs(j);
  }
  xs(0) = std::sqrt(1.0 + sp2);
  arma::rowvec v = project_to_tangent_space(xs, vin);
  double norm2 = lorentz_dot(v, v);
  if (norm2 < 0.0) {
    norm2 = 0.0;
  }
  const double r = std::sqrt(norm2);
  arma::rowvec y;
  if (r < 1e-10) {
    y = xs;
  } else {
    const double c = std::cosh(r), s = std::sinh(r), a = s / r;
    if (!std::isfinite(c) || !std::isfinite(s) || !std::isfinite(a)) {
      return out;
    }
    y = c * xs + a * v;
  }
  sp2 = 0.0;
  for (arma::uword j = 1; j < y.n_elem; ++j) {
    sp2 += y(j) * y(j);
  }
  const double y0 = std::sqrt(1.0 + sp2);
  if (!std::isfinite(y0)) {
    return out;
  }
  if (std::asinh(std::sqrt(sp2)) > max_hyperbolic_radius) {
    out.hitRadiusCap = true;
    return out;
  }
  y(0) = y0;
  out.point = y;
  out.isValid = true;
  return out;
}

template <typename Vector>
inline bool
exp_map_into(const Vector &x, const arma::rowvec &input_tangent,
             arma::rowvec &snapped_point, arma::rowvec &projected_tangent,
             arma::rowvec &output_point, bool &hit_radius_cap,
             double max_hyperbolic_radius = default_max_hyperbolic_radius) {
  if (x.n_elem != input_tangent.n_elem || x.n_elem < 2) {
    Rcpp::stop("exp_map_into dimension mismatch.");
  }
  const arma::uword ambient_dimension = x.n_elem;
  if (snapped_point.n_elem != ambient_dimension) {
    snapped_point.set_size(ambient_dimension);
  }
  if (projected_tangent.n_elem != ambient_dimension) {
    projected_tangent.set_size(ambient_dimension);
  }
  if (output_point.n_elem != ambient_dimension) {
    output_point.set_size(ambient_dimension);
  }
  hit_radius_cap = false;

  double spatial_norm_squared = 0.0;
  for (arma::uword q = 1; q < ambient_dimension; ++q) {
    snapped_point(q) = x(q);
    spatial_norm_squared += x(q) * x(q);
  }
  snapped_point(0) = std::sqrt(1.0 + spatial_norm_squared);

  const double tangent_component = lorentz_dot(input_tangent, snapped_point);
  for (arma::uword q = 0; q < ambient_dimension; ++q) {
    projected_tangent(q) =
        input_tangent(q) + tangent_component * snapped_point(q);
  }

  double tangent_norm_squared =
      lorentz_dot(projected_tangent, projected_tangent);
  if (tangent_norm_squared < 0.0) {
    tangent_norm_squared = 0.0;
  }
  const double tangent_norm = std::sqrt(tangent_norm_squared);

  if (tangent_norm < 1e-10) {
    output_point = snapped_point;
  } else {
    const double hyperbolic_cosine = std::cosh(tangent_norm);
    const double hyperbolic_sine = std::sinh(tangent_norm);
    const double tangent_multiplier = hyperbolic_sine / tangent_norm;
    if (!std::isfinite(hyperbolic_cosine) || !std::isfinite(hyperbolic_sine) ||
        !std::isfinite(tangent_multiplier)) {
      return false;
    }
    for (arma::uword q = 0; q < ambient_dimension; ++q) {
      output_point(q) = hyperbolic_cosine * snapped_point(q) +
                        tangent_multiplier * projected_tangent(q);
    }
  }

  spatial_norm_squared = 0.0;
  for (arma::uword q = 1; q < ambient_dimension; ++q) {
    spatial_norm_squared += output_point(q) * output_point(q);
  }
  const double time_coordinate = std::sqrt(1.0 + spatial_norm_squared);
  if (!std::isfinite(time_coordinate)) {
    return false;
  }
  if (std::asinh(std::sqrt(spatial_norm_squared)) > max_hyperbolic_radius) {
    hit_radius_cap = true;
    return false;
  }
  output_point(0) = time_coordinate;
  return true;
}

struct LogMapResult {
  arma::rowvec vector;
  bool isValid;
  double distance;
};

inline LogMapResult log_map(const arma::rowvec &from, const arma::rowvec &to) {
  LogMapResult out{arma::rowvec(), false, 0.0};
  if (from.n_elem != to.n_elem || from.n_elem < 2) {
    Rcpp::stop("log_map dimension mismatch.");
  }
  const double r = geodesic_distance_polar(from, to);
  out.distance = r;
  if (r < 1e-10) {
    out.vector = arma::zeros<arma::rowvec>(from.n_elem);
    out.isValid = true;
    return out;
  }
  const double sh = std::sinh(r), ch = std::cosh(r);
  if (!std::isfinite(sh) || !std::isfinite(ch) || !(sh > 0.0)) {
    return out;
  }
  const double factor = r < 1e-4 ? 1.0 - r * r / 6.0 : r / sh;
  out.vector = factor * (to - ch * from);
  out.isValid = out.vector.is_finite();
  return out;
}

inline double log_volume_correction(double r, double d) {
  if (r < 1e-10) {
    return -(d - 1.0) * std::log1p(r * r / 6.0);
  }
  const double sh = std::sinh(r);
  if (!(sh > 0.0) || !std::isfinite(sh)) {
    return -std::numeric_limits<double>::infinity();
  }
  return (d - 1.0) * (std::log(r) - std::log(sh));
}

inline double log_wrapped_normal_proposal_density(const arma::rowvec &from,
                                                  const arma::rowvec &to,
                                                  double sigma_sq) {
  if (!(sigma_sq > 0.0) || !std::isfinite(sigma_sq)) {
    return -std::numeric_limits<double>::infinity();
  }
  const LogMapResult lm = log_map(from, to);
  if (!lm.isValid) {
    return -std::numeric_limits<double>::infinity();
  }
  const double d = from.n_elem - 1.0;
  return -lm.distance * lm.distance / (2.0 * sigma_sq) -
         0.5 * d * std::log(2.0 * M_PI * sigma_sq) +
         log_volume_correction(lm.distance, d);
}

template <typename Vector>
inline double log_wrapped_normal_prior_density(const Vector &x,
                                               double sigma_sq) {
  if (!(sigma_sq > 0.0) || !std::isfinite(sigma_sq) || x.n_elem < 2) {
    return -std::numeric_limits<double>::infinity();
  }
  double radius_squared = 0.0;
  for (arma::uword q = 1; q < x.n_elem; ++q) {
    radius_squared += x(q) * x(q);
  }
  const double radius = std::asinh(std::sqrt(radius_squared));
  const double dimension = x.n_elem - 1.0;
  return -radius * radius / (2.0 * sigma_sq) -
         0.5 * dimension * std::log(2.0 * M_PI * sigma_sq) +
         log_volume_correction(radius, dimension);
}

struct WrappedNormalProposalResult {
  arma::rowvec point;
  bool isValid;
  bool hitRadiusCap = false;
};

template <typename Vector>
inline bool propose_wrapped_normal_into(
    const Vector &x, double sigma, arma::mat &basis, arma::rowvec &coefficients,
    arma::rowvec &ambient_tangent, arma::rowvec &snapped_point,
    arma::rowvec &projected_tangent, arma::rowvec &candidate,
    bool &hit_radius_cap,
    double max_hyperbolic_radius = default_max_hyperbolic_radius) {
  if (!(sigma > 0.0) || !std::isfinite(sigma)) {
    hit_radius_cap = false;
    return false;
  }
  fill_lorentz_orthonormal_basis(x, basis);
  if (coefficients.n_elem != basis.n_cols) {
    coefficients.set_size(basis.n_cols);
  }
  if (ambient_tangent.n_elem != basis.n_rows) {
    ambient_tangent.set_size(basis.n_rows);
  }
  coefficients.randn();
  coefficients *= sigma;
  ambient_tangent.zeros();
  for (arma::uword q = 0; q < basis.n_rows; ++q) {
    for (arma::uword k = 0; k < basis.n_cols; ++k) {
      ambient_tangent(q) += basis(q, k) * coefficients(k);
    }
  }
  return exp_map_into(x, ambient_tangent, snapped_point, projected_tangent,
                      candidate, hit_radius_cap, max_hyperbolic_radius);
}

// Numerically stable equivalent of the wrapped-normal proposal above.  The
// Gaussian coefficients are coordinates in the canonical boosted orthonormal
// frame, so their Euclidean norm is the tangent norm exactly.  Carrying that
// invariant forward avoids forming a Lorentz basis, reprojecting an already
// tangent vector, and recovering its norm by cancellation of ambient terms.
template <typename Vector>
inline bool intrinsic_wrapped_normal_from_coefficients_into(
    const Vector &x, const arma::rowvec &coefficients, arma::rowvec &candidate,
    bool &hit_radius_cap,
    double max_hyperbolic_radius = default_max_hyperbolic_radius) {
  hit_radius_cap = false;
  if (x.n_elem < 2 || coefficients.n_elem != x.n_elem - 1 ||
      !coefficients.is_finite() || !(max_hyperbolic_radius > 0.0) ||
      !std::isfinite(max_hyperbolic_radius)) {
    return false;
  }
  if (candidate.n_elem != x.n_elem) {
    candidate.set_size(x.n_elem);
  }
  const double tangent_norm = arma::norm(coefficients, 2);
  if (!std::isfinite(tangent_norm)) {
    return false;
  }

  double spatial_norm_squared = 0.0;
  for (arma::uword q = 1; q < x.n_elem; ++q) {
    spatial_norm_squared += x(q) * x(q);
  }
  const double spatial_norm = std::sqrt(spatial_norm_squared);
  const double snapped_time = std::sqrt(1.0 + spatial_norm_squared);

  if (tangent_norm < 1e-10) {
    candidate(0) = snapped_time;
    for (arma::uword q = 1; q < x.n_elem; ++q) {
      candidate(q) = x(q);
    }
    return true;
  }

  double radial_coefficient = 0.0;
  if (spatial_norm > 1e-12) {
    for (arma::uword q = 1; q < x.n_elem; ++q) {
      radial_coefficient += (x(q) / spatial_norm) * coefficients(q - 1);
    }
  }
  const double hyperbolic_cosine = std::cosh(tangent_norm);
  const double tangent_multiplier = std::sinh(tangent_norm) / tangent_norm;
  if (!std::isfinite(hyperbolic_cosine) || !std::isfinite(tangent_multiplier)) {
    return false;
  }

  double candidate_spatial_norm_squared = 0.0;
  for (arma::uword q = 1; q < x.n_elem; ++q) {
    const double direction = spatial_norm > 1e-12 ? x(q) / spatial_norm : 0.0;
    const double tangent_spatial =
        coefficients(q - 1) +
        (snapped_time - 1.0) * radial_coefficient * direction;
    candidate(q) =
        hyperbolic_cosine * x(q) + tangent_multiplier * tangent_spatial;
    candidate_spatial_norm_squared += candidate(q) * candidate(q);
  }
  const double candidate_radius =
      std::asinh(std::sqrt(candidate_spatial_norm_squared));
  if (!std::isfinite(candidate_radius)) {
    return false;
  }
  if (candidate_radius > max_hyperbolic_radius) {
    hit_radius_cap = true;
    return false;
  }
  candidate(0) = std::sqrt(1.0 + candidate_spatial_norm_squared);
  return std::isfinite(candidate(0));
}

template <typename Vector>
inline bool propose_intrinsic_wrapped_normal_into(
    const Vector &x, double sigma, arma::rowvec &coefficients,
    arma::rowvec &candidate, bool &hit_radius_cap,
    double max_hyperbolic_radius = default_max_hyperbolic_radius) {
  if (!(sigma > 0.0) || !std::isfinite(sigma)) {
    hit_radius_cap = false;
    return false;
  }
  const arma::uword intrinsic_dimension = x.n_elem - 1;
  if (coefficients.n_elem != intrinsic_dimension) {
    coefficients.set_size(intrinsic_dimension);
  }
  coefficients.randn();
  coefficients *= sigma;
  return intrinsic_wrapped_normal_from_coefficients_into(
      x, coefficients, candidate, hit_radius_cap, max_hyperbolic_radius);
}

inline WrappedNormalProposalResult propose_wrapped_normal(const arma::rowvec &x,
                                                          double sigma) {
  WrappedNormalProposalResult out{x, false, false};
  if (!(sigma > 0.0) || !std::isfinite(sigma)) {
    return out;
  }
  const arma::mat E = lorentz_orthonormal_basis(x);
  arma::rowvec b(E.n_cols, arma::fill::randn);
  b *= sigma;
  const ExpMapResult p = exp_map(x, (E * b.t()).t());
  if (!p.isValid) {
    out.hitRadiusCap = p.hitRadiusCap;
    return out;
  }
  out.point = p.point;
  out.isValid = true;
  return out;
}

template <typename VectorX, typename VectorY>
inline double distance_hyperboloid(const VectorX &x, const VectorY &y) {
  if (x.n_elem != y.n_elem || x.n_elem < 2) {
    Rcpp::stop("Distance dimension mismatch.");
  }
  const double a = -lorentz_dot(x, y);
  const double guard =
      1e4 * std::numeric_limits<double>::epsilon() * x(0) * y(0);
  if (std::isfinite(a) && a - 1.0 > guard) {
    return std::acosh(a);
  }
  return geodesic_distance_polar(x, y);
}

template <typename Vector> inline double radius_hyperboloid(const Vector &x) {
  double spatial_norm_squared = 0.0;
  for (arma::uword q = 1; q < x.n_elem; ++q) {
    spatial_norm_squared += x(q) * x(q);
  }
  return std::asinh(std::sqrt(spatial_norm_squared));
}

template <typename Vector> inline arma::rowvec log_origin(const Vector &x) {
  if (x.n_elem < 2) {
    Rcpp::stop("log_origin needs at least two coordinates.");
  }
  double spatial_norm_squared = 0.0;
  for (arma::uword q = 1; q < x.n_elem; ++q) {
    spatial_norm_squared += x(q) * x(q);
  }
  const double sv = std::sqrt(spatial_norm_squared);
  const double r = std::asinh(sv);
  arma::rowvec q(x.n_elem - 1, arma::fill::zeros);
  if (sv == 0.0) {
    return q;
  }
  for (arma::uword index = 0; index < q.n_elem; ++index) {
    q(index) = (r / sv) * x(index + 1);
  }
  return q;
}

inline arma::rowvec
exp_origin(const arma::rowvec &q,
           double max_hyperbolic_radius = default_max_hyperbolic_radius) {
  if (q.n_elem < 1) {
    Rcpp::stop("exp_origin needs a nonempty tangent vector.");
  }
  arma::rowvec x(q.n_elem + 1, arma::fill::zeros);
  const double r = arma::norm(q, 2);
  if (r == 0.0) {
    x(0) = 1.0;
    return x;
  }
  if (r > max_hyperbolic_radius || !std::isfinite(r)) {
    x.fill(std::numeric_limits<double>::quiet_NaN());
    return x;
  }
  x(0) = std::cosh(r);
  x.subvec(1, q.n_elem) = (std::sinh(r) / r) * q;
  return x;
}

// Symmetric geodesic random walk on the unit direction sphere at fixed
// hyperbolic radius.  The step angle has a centered Normal law and the tangent
// direction is uniform, so q(u'|u) depends only on the spherical distance and
// equals q(u|u').  This kernel is paired with a separate radial update.
template <typename Vector>
inline bool propose_angular_fixed_radius_into(const Vector &x,
                                              double angular_scale,
                                              arma::rowvec &gaussian,
                                              arma::rowvec &candidate,
                                              double minimum_radius = 0.0) {
  if (x.n_elem < 3 || !(angular_scale > 0.0) || !std::isfinite(angular_scale)) {
    return false;
  }
  double spatial_norm_squared = 0.0;
  for (arma::uword index = 1; index < x.n_elem; ++index) {
    spatial_norm_squared += x(index) * x(index);
  }
  const double spatial_norm = std::sqrt(spatial_norm_squared);
  const double radius = std::asinh(spatial_norm);
  if (!(radius > minimum_radius) || !std::isfinite(radius)) {
    return false;
  }
  arma::rowvec direction(x.n_elem - 1);
  for (arma::uword index = 0; index < direction.n_elem; ++index) {
    direction(index) = x(index + 1) / spatial_norm;
  }
  if (gaussian.n_elem != direction.n_elem) {
    gaussian.set_size(direction.n_elem);
  }
  gaussian.randn();
  gaussian -= arma::dot(gaussian, direction) * direction;
  const double tangent_norm = arma::norm(gaussian, 2);
  if (!(tangent_norm > 1e-12) || !std::isfinite(tangent_norm)) {
    return false;
  }
  gaussian /= tangent_norm;
  const double angle = angular_scale * R::norm_rand();
  arma::rowvec proposed_direction =
      std::cos(angle) * direction + std::sin(angle) * gaussian;
  const double proposed_direction_norm = arma::norm(proposed_direction, 2);
  if (!(proposed_direction_norm > 0.0) ||
      !std::isfinite(proposed_direction_norm)) {
    return false;
  }
  proposed_direction /= proposed_direction_norm;

  // This is a fixed-radius move.  Reconstruct directly from the current
  // spatial norm instead of calling exp_origin with an artificial cap only a
  // few machine epsilons above the radius.  The latter can reject a valid
  // angular proposal when roundoff makes ||radius * proposed_direction||
  // infinitesimally larger than radius, particularly in H^2.
  if (candidate.n_elem != x.n_elem) {
    candidate.set_size(x.n_elem);
  }
  double candidate_spatial_norm_squared = 0.0;
  for (arma::uword index = 0; index < proposed_direction.n_elem; ++index) {
    candidate(index + 1) = spatial_norm * proposed_direction(index);
    candidate_spatial_norm_squared +=
        candidate(index + 1) * candidate(index + 1);
  }
  candidate(0) = std::sqrt(1.0 + candidate_spatial_norm_squared);
  return candidate.is_finite();
}

} // namespace hyperboloid_utils

// Compatibility for older auxiliary sources that used the result type globally.
using ExpMapResult = hyperboloid_utils::ExpMapResult;

#endif
