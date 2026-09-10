// align_hyper.cpp — fast hyperbolic Procrustes alignment (apex_rotation) +
// per-node Fréchet mean, 1:1 port of R/procrustes.R, OpenMP-parallel over draws.
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace arma;

static inline double ldot(const rowvec& x, const rowvec& y) {
  double s = -x[0] * y[0];
  for (uword i = 1; i < x.n_elem; ++i) s += x[i] * y[i];
  return s;
}
static inline double lnorm_sq(const rowvec& v) { return -ldot(v, v); }

static inline double ldist(const rowvec& x, const rowvec& y) {
  double val = -ldot(x, y);
  if (val < 1.0) val = 1.0;
  return std::acosh(val);
}

// Fréchet mean of rows of `data` (each a Lorentz point), iterative.
static rowvec frechet_mean(const mat& data) {
  uword d = data.n_cols, n = data.n_rows;
  rowvec y(d, fill::zeros); y[0] = 1.0;
  double eps = 1.0; int iter = 0;
  // eps is the squared geodesic update; 1e-12 means a 1e-6 distance
  // tolerance.  The earlier 1e-6 threshold only guaranteed 1e-3 accuracy and
  // can create visible drift in generalized Procrustes iterations.
  while (eps > 1e-12 && iter < 1000) {
    ++iter;
    rowvec tmp(d, fill::zeros);
    for (uword i = 0; i < n; ++i) {
      double val = -ldot(data.row(i), y);
      if (val < 1.0 && val > 1.0 - 1e-9) val = 1.0;
      if (val < 1.0) continue;
      double sinh_sq = val * val - 1.0;
      if (sinh_sq >= 1e-12)
        tmp += (2.0 * std::acosh(val) / std::sqrt(sinh_sq)) * data.row(i);
    }
    double nsq = lnorm_sq(tmp);
    if (nsq > 1e-12) {
      rowvec y_new = tmp / std::sqrt(nsq);
      double dv = ldist(y_new, y);
      eps = dv * dv;
      y = y_new;
    } else break;
  }
  return y;
}

// Weighted counterpart used for bipartite configurations.  Weights are
// normalized by the caller; a common multiplicative constant is immaterial.
static rowvec frechet_mean_weighted(const mat& data, const vec& weight) {
  uword d = data.n_cols, n = data.n_rows;
  if (weight.n_elem != n) Rcpp::stop("weight length does not match rows");
  rowvec y(d, fill::zeros); y[0] = 1.0;
  double eps = 1.0; int iter = 0;
  while (eps > 1e-12 && iter < 1000) {
    ++iter;
    rowvec tmp(d, fill::zeros);
    for (uword i = 0; i < n; ++i) {
      double val = -ldot(data.row(i), y);
      if (val < 1.0 && val > 1.0 - 1e-9) val = 1.0;
      if (val < 1.0) continue;
      double sinh_sq = val * val - 1.0;
      if (sinh_sq >= 1e-12)
        tmp += weight[i] * (2.0 * std::acosh(val) / std::sqrt(sinh_sq)) * data.row(i);
    }
    double nsq = lnorm_sq(tmp);
    if (nsq > 1e-12) {
      rowvec y_new = tmp / std::sqrt(nsq);
      double dv = ldist(y_new, y);
      eps = dv * dv; y = y_new;
    } else break;
  }
  return y;
}

static inline rowvec log_map(const rowvec& base, const rowvec& target) {
  double lambda = -ldot(base, target);
  if (lambda < 1.0 && lambda > 1.0 - 1e-9) lambda = 1.0;
  double sinh_sq = lambda * lambda - 1.0;
  if (sinh_sq < 1e-12) return rowvec(base.n_elem, fill::zeros);
  double d = std::acosh(lambda), sinh_d = std::sqrt(sinh_sq);
  return (d / sinh_d) * (target - lambda * base);
}

static inline rowvec ptransport(const rowvec& x, const rowvec& y, const rowvec& v) {
  double lambda = -ldot(x, y);
  double num = ldot(y - lambda * x, v);
  double den = 1.0 + lambda;
  return v + (num / den) * (x + y);
}

static inline rowvec exp_point(const rowvec& v, const rowvec& base) {
  double nsq = std::abs(ldot(v, v));
  double nv = std::sqrt(nsq);
  if (nv < 1e-12) return base;
  return base * std::cosh(nv) + (v / nv) * std::sinh(nv);
}

// riemannian_translation_hpa for all rows of `data`.
static mat translate(const mat& data, const rowvec& mfrom, const rowvec& mto) {
  mat out(data.n_rows, data.n_cols);
  for (uword i = 0; i < data.n_rows; ++i) {
    rowvec v1 = log_map(mfrom, data.row(i));
    rowvec v2 = ptransport(mfrom, mto, v1);
    out.row(i) = exp_point(v2, mto);
  }
  return out;
}

// Tangent coords at apex o=(1,0,..): drop x0 of log_map(o,row).
static mat tangent_at_apex(const mat& cfg_o) {
  uword n = cfg_o.n_rows, dt = cfg_o.n_cols - 1, o0 = 0;
  mat X(n, dt);
  rowvec o(cfg_o.n_cols, fill::zeros); o[0] = 1.0;
  for (uword i = 0; i < n; ++i) {
    rowvec v = log_map(o, cfg_o.row(i));
    X.row(i) = v.subvec(1, cfg_o.n_cols - 1);
  }
  (void)o0;
  return X;
}

// One apex_rotation alignment of `target` to precomputed reference side
// (mu_ref, Xr). Returns aligned config (n x d). allow_reflection always TRUE.
// allow_boost=false keeps the origin the prior names.  The likelihood is invariant under the whole
// of SO+(d,1), but the prior is a function of the radius and so is invariant under the rotations
// only, which makes the boost a softly identified direction rather than a gauge.  Since the radius is
// reported and is not isometry invariant, carrying each draw's Frechet mean to a common point moves a
// reported quantity.  Matching Frechet means is also not the minimiser of the Procrustes objective
// here, the composition of two boosts being a boost times a rotation, so the step is a heuristic in
// this geometry where its Euclidean counterpart is exact.
static mat procrustes_apex(const mat& target, const rowvec& mu_ref,
                           const mat& Xr, const rowvec& o, bool allow_boost = true) {
  rowvec mu_tar = allow_boost ? frechet_mean(target) : o;
  mat tar_o = allow_boost ? translate(target, mu_tar, o) : target;
  mat Xt = tangent_at_apex(tar_o);
  mat M = Xt.t() * Xr;                 // dt x dt
  mat U, V; vec s;
  svd(U, s, V, M);
  mat R = U * V.t();                   // O(d) Procrustes (Xt -> Xr)
  mat Xt_rot = Xt * R;
  // exp of tangent (0, v) at apex
  mat tar_o_rot(target.n_rows, target.n_cols);
  for (uword i = 0; i < target.n_rows; ++i) {
    rowvec v(target.n_cols, fill::zeros);
    v.subvec(1, target.n_cols - 1) = Xt_rot.row(i);
    tar_o_rot.row(i) = exp_point(v, o);
  }
  return allow_boost ? translate(tar_o_rot, o, mu_ref) : tar_o_rot;
}

static mat procrustes_apex_weighted(const mat& target, const rowvec& mu_ref,
                                    const mat& Xr, const rowvec& o,
                                    const vec& weight, bool allow_boost) {
  rowvec mu_tar = allow_boost ? frechet_mean_weighted(target, weight) : o;
  mat tar_o = allow_boost ? translate(target, mu_tar, o) : target;
  mat Xt = tangent_at_apex(tar_o);
  mat weighted_ref = Xr.each_col() % weight;
  mat M = Xt.t() * weighted_ref; // t(Xt) diag(weight) Xr
  mat U, V; vec s; svd(U, s, V, M);
  mat R = U * V.t();
  mat Xt_rot = Xt * R;
  mat out(target.n_rows, target.n_cols);
  for (uword i = 0; i < target.n_rows; ++i) {
    rowvec v(target.n_cols, fill::zeros);
    v.subvec(1, target.n_cols - 1) = Xt_rot.row(i);
    out.row(i) = exp_point(v, o);
  }
  return allow_boost ? translate(out, o, mu_ref) : out;
}

// Euclidean Procrustes (translation + reflection, no dilation), 1:1 with
// MCMCpack::procrustes(X=target, Xstar=ref, translation=TRUE, dilation=FALSE).
// Returns aligned target = X R + tt.
static mat procrustes_euc(const mat& X, const mat& Xc, const rowvec& Xc_mean) {
  // Xc = ref centered, Xc_mean = colMeans(ref) (precomputed, constant)
  rowvec Xm = mean(X, 0);
  mat Xcen = X.each_row() - Xm;
  mat C = Xc.t() * Xcen;              // t(Xstar)*J*X  (m x m)
  mat U, V; vec s;
  svd(U, s, V, C);
  mat R = V * U.t();
  mat XR = X * R;
  rowvec tt = Xc_mean - mean(XR, 0);   // colMeans(Xstar - X R)
  return XR.each_row() + tt;
}

// [[Rcpp::export]]
Rcpp::List align_euc_cpp(arma::cube draws, arma::mat ref,
                         int n_cores = 4, bool verbose = true) {
  uword n = draws.n_rows, d = draws.n_cols, S = draws.n_slices;
  rowvec ref_mean = mean(ref, 0);
  mat ref_c = ref.each_row() - ref_mean;     // J*Xstar, computed once

  cube aligned(n, d, S);
  int chunk = std::max<int>(1, (int)S / 60);
  for (int base = 0; base < (int)S; base += chunk) {
    int hi = std::min<int>(base + chunk, (int)S);
#ifdef _OPENMP
#pragma omp parallel for num_threads(n_cores) schedule(dynamic)
#endif
    for (int s = base; s < hi; ++s)
      aligned.slice(s) = procrustes_euc(draws.slice(s), ref_c, ref_mean);
    if (verbose) {
      int pct = (int)(100.0 * hi / S), bars = pct / 2;
      Rcpp::Rcout << "\r  align [" << std::string(bars, '=')
                  << std::string(50 - bars, ' ') << "] " << pct << "%  "
                  << hi << "/" << S << std::flush;
    }
    Rcpp::checkUserInterrupt();
  }
  if (verbose) Rcpp::Rcout << std::endl;

  mat Z = mean(aligned, 2);                  // per-node mean over slices
  vec rad(n), ang(n);
  for (uword i = 0; i < n; ++i) {
    rad[i] = std::sqrt(Z(i, 0) * Z(i, 0) + Z(i, 1) * Z(i, 1));
    ang[i] = std::atan2(Z(i, 1), Z(i, 0));
  }
  return Rcpp::List::create(Rcpp::Named("radius") = rad,
                            Rcpp::Named("angle")  = ang,
                            Rcpp::Named("Z")      = Z);
}

// [[Rcpp::export]]
Rcpp::List align_frechet_cpp(arma::cube draws, arma::mat ref,
                             int n_cores = 4, bool verbose = true,
                             bool allow_boost = true) {
  // draws: n x d x S (slice s = one config).  ref: n x d (the MAP draw).
  uword n = draws.n_rows, d = draws.n_cols, S = draws.n_slices;
  rowvec o(d, fill::zeros); o[0] = 1.0;

  // reference side computed ONCE (constant across draws)
  rowvec mu_ref = allow_boost ? frechet_mean(ref) : o;
  mat ref_o = allow_boost ? translate(ref, mu_ref, o) : ref;
  mat Xr = tangent_at_apex(ref_o);

  cube aligned(n, d, S);
  int done = 0;
  int chunk = std::max<int>(1, (int)S / 60);   // ~60 progress ticks

  for (int base = 0; base < (int)S; base += chunk) {
    int hi = std::min<int>(base + chunk, (int)S);
#ifdef _OPENMP
#pragma omp parallel for num_threads(n_cores) schedule(dynamic)
#endif
    for (int s = base; s < hi; ++s)
      aligned.slice(s) = procrustes_apex(draws.slice(s), mu_ref, Xr, o, allow_boost);

    done = hi;
    if (verbose) {
      int pct = (int)(100.0 * done / S);
      int bars = pct / 2;
      Rcpp::Rcout << "\r  align [" << std::string(bars, '=')
                  << std::string(50 - bars, ' ') << "] " << pct << "%  "
                  << done << "/" << S << std::flush;
    }
    Rcpp::checkUserInterrupt();
  }
  if (verbose) Rcpp::Rcout << std::endl;

  // per-node Fréchet mean over aligned slices
  mat Z(n, d);
#ifdef _OPENMP
#pragma omp parallel for num_threads(n_cores) schedule(dynamic)
#endif
  for (int i = 0; i < (int)n; ++i) {
    mat nodi(S, d);
    for (uword s = 0; s < S; ++s) nodi.row(s) = aligned.slice(s).row(i);
    Z.row(i) = frechet_mean(nodi);
  }

  vec rad(n), ang(n);
  for (uword i = 0; i < n; ++i) {
    double x0 = std::max(1.0, Z(i, 0));
    rad[i] = std::acosh(x0);
    ang[i] = std::atan2(Z(i, 2), Z(i, 1));
  }
  return Rcpp::List::create(Rcpp::Named("radius") = rad,
                            Rcpp::Named("angle")  = ang,
                            Rcpp::Named("Z")      = Z);
}

// Weighted hyperbolic generalized-Procrustes building block.  A single common
// isometry is applied to all rows in a draw, while landmark weights determine
// the alignment only.  Nodewise posterior Fréchet means remain equally
// weighted over posterior draws.
// [[Rcpp::export]]
Rcpp::List align_frechet_weighted_cpp(arma::cube draws, arma::mat ref,
                                      arma::vec weight, int n_cores = 4,
                                      bool verbose = true,
                                      bool allow_boost = false) {
  uword n = draws.n_rows, d = draws.n_cols, S = draws.n_slices;
  if (ref.n_rows != n || ref.n_cols != d) Rcpp::stop("reference dimension mismatch");
  if (weight.n_elem != n || !weight.is_finite() || any(weight <= 0.0))
    Rcpp::stop("weights must be finite, positive, and match the number of rows");
  weight /= accu(weight);
  rowvec o(d, fill::zeros); o[0] = 1.0;
  rowvec mu_ref = allow_boost ? frechet_mean_weighted(ref, weight) : o;
  mat ref_o = allow_boost ? translate(ref, mu_ref, o) : ref;
  mat Xr = tangent_at_apex(ref_o);

  cube aligned(n, d, S);
  int chunk = std::max<int>(1, (int)S / 60);
  for (int base = 0; base < (int)S; base += chunk) {
    int hi = std::min<int>(base + chunk, (int)S);
#ifdef _OPENMP
#pragma omp parallel for num_threads(n_cores) schedule(dynamic)
#endif
    for (int s = base; s < hi; ++s)
      aligned.slice(s) = procrustes_apex_weighted(
        draws.slice(s), mu_ref, Xr, o, weight, allow_boost);
    if (verbose) {
      int pct = (int)(100.0 * hi / S), bars = pct / 2;
      Rcpp::Rcout << "\r  weighted align [" << std::string(bars, '=')
                  << std::string(50 - bars, ' ') << "] " << pct << "%  "
                  << hi << "/" << S << std::flush;
    }
    Rcpp::checkUserInterrupt();
  }
  if (verbose) Rcpp::Rcout << std::endl;

  mat Z(n, d);
#ifdef _OPENMP
#pragma omp parallel for num_threads(n_cores) schedule(dynamic)
#endif
  for (int i = 0; i < (int)n; ++i) {
    mat nodi(S, d);
    for (uword s = 0; s < S; ++s) nodi.row(s) = aligned.slice(s).row(i);
    Z.row(i) = frechet_mean(nodi);
  }
  vec rad(n), ang(n);
  for (uword i = 0; i < n; ++i) {
    rad[i] = std::acosh(std::max(1.0, Z(i, 0)));
    ang[i] = std::atan2(Z(i, 2), Z(i, 1));
  }
  return Rcpp::List::create(Rcpp::Named("radius") = rad,
                            Rcpp::Named("angle") = ang,
                            Rcpp::Named("Z") = Z);
}

// [[Rcpp::export]]
arma::vec weighted_alignment_rms_cpp(arma::cube draws, arma::mat ref,
                                     arma::vec weight, int n_cores = 4,
                                     bool allow_boost = false) {
  uword n = draws.n_rows, d = draws.n_cols, S = draws.n_slices;
  if (ref.n_rows != n || ref.n_cols != d || weight.n_elem != n)
    Rcpp::stop("dimension mismatch");
  if (!weight.is_finite() || any(weight <= 0.0))
    Rcpp::stop("weights must be finite and positive");
  weight /= accu(weight);
  rowvec o(d, fill::zeros); o[0] = 1.0;
  rowvec mu_ref = allow_boost ? frechet_mean_weighted(ref, weight) : o;
  mat ref_o = allow_boost ? translate(ref, mu_ref, o) : ref;
  mat Xr = tangent_at_apex(ref_o);
  vec rms(S);
#ifdef _OPENMP
#pragma omp parallel for num_threads(n_cores) schedule(dynamic)
#endif
  for (int s = 0; s < (int)S; ++s) {
    mat a = procrustes_apex_weighted(
      draws.slice(s), mu_ref, Xr, o, weight, allow_boost);
    double loss = 0.0;
    for (uword i = 0; i < n; ++i) {
      double di = ldist(a.row(i), ref.row(i));
      loss += weight[i] * di * di;
    }
    rms[s] = std::sqrt(loss);
  }
  return rms;
}

// [[Rcpp::export]]
Rcpp::List diagnose_hyper_alignment_cpp(arma::cube draws, arma::mat ref,
                                        Rcpp::IntegerVector pair_i,
                                        Rcpp::IntegerVector pair_j,
                                        int n_cores = 4) {
  // Diagnostics after one common joint alignment per draw. Pair indices arrive
  // from R as one-based indices.  Rotation-only diagnostics preserve the apex
  // and therefore every radius; full diagnostics additionally allow a boost.
  uword n = draws.n_rows, d = draws.n_cols, S = draws.n_slices;
  int P = pair_i.size();
  if (pair_j.size() != P) Rcpp::stop("pair_i and pair_j must have equal length");
  int nt = std::max(1, n_cores);

  rowvec o(d, fill::zeros); o[0] = 1.0;
  mat Xr_rot = tangent_at_apex(ref);
  rowvec mu_ref_full = frechet_mean(ref);
  mat ref_o_full = translate(ref, mu_ref_full, o);
  mat Xr_full = tangent_at_apex(ref_o_full);

  vec pair_ref(P);
  for (int p = 0; p < P; ++p) {
    int i = pair_i[p] - 1, j = pair_j[p] - 1;
    if (i < 0 || j < 0 || i >= (int)n || j >= (int)n || i == j)
      Rcpp::stop("invalid sampled pair index");
    pair_ref[p] = ldist(ref.row(i), ref.row(j));
  }

  mat csum(n, nt, fill::zeros), ssum(n, nt, fill::zeros);
  mat rsum(n, nt, fill::zeros), r2sum(n, nt, fill::zeros);
  mat node_ref_sq(n, nt, fill::zeros);
  vec rms_rot(S, fill::zeros), rms_full(S, fill::zeros);
  vec pair_rmse(S, fill::zeros), pair_cor(S, fill::zeros);

#ifdef _OPENMP
#pragma omp parallel for num_threads(n_cores) schedule(dynamic)
#endif
  for (int s = 0; s < (int)S; ++s) {
    int tid = 0;
#ifdef _OPENMP
    tid = omp_get_thread_num();
#endif
    mat target = draws.slice(s);
    mat A = procrustes_apex(target, o, Xr_rot, o, false);
    mat B = procrustes_apex(target, mu_ref_full, Xr_full, o, true);
    double er = 0.0, ef = 0.0;
    for (uword i = 0; i < n; ++i) {
      double rr = std::acosh(std::max(1.0, A(i, 0)));
      double aa = std::atan2(A(i, 2), A(i, 1));
      csum(i, tid) += std::cos(aa); ssum(i, tid) += std::sin(aa);
      rsum(i, tid) += rr; r2sum(i, tid) += rr * rr;
      double dr = ldist(A.row(i), ref.row(i));
      double df = ldist(B.row(i), ref.row(i));
      er += dr * dr; ef += df * df;
      node_ref_sq(i, tid) += dr * dr;
    }
    rms_rot[s] = std::sqrt(er / n);
    rms_full[s] = std::sqrt(ef / n);

    double se = 0.0, sx = 0.0, sy = 0.0, sxx = 0.0, syy = 0.0, sxy = 0.0;
    for (int p = 0; p < P; ++p) {
      int i = pair_i[p] - 1, j = pair_j[p] - 1;
      double x = pair_ref[p], y = ldist(target.row(i), target.row(j));
      double e = y - x; se += e * e;
      sx += x; sy += y; sxx += x * x; syy += y * y; sxy += x * y;
    }
    pair_rmse[s] = std::sqrt(se / P);
    double den = std::sqrt(std::max(0.0, (P * sxx - sx * sx) * (P * syy - sy * sy)));
    pair_cor[s] = den > 0.0 ? (P * sxy - sx * sy) / den : NA_REAL;
  }

  vec cs = sum(csum, 1), ss = sum(ssum, 1);
  vec rs = sum(rsum, 1), r2s = sum(r2sum, 1);
  vec resultant = sqrt(square(cs) + square(ss)) / (double)S;
  vec mean_angle = atan2(ss, cs);
  vec mean_radius = rs / (double)S;
  vec var_radius = clamp(r2s / (double)S - square(mean_radius), 0.0, datum::inf);
  vec node_rms = sqrt(sum(node_ref_sq, 1) / (double)S);

  return Rcpp::List::create(
    Rcpp::Named("angle_resultant") = resultant,
    Rcpp::Named("angle_mean") = mean_angle,
    Rcpp::Named("radius_mean") = mean_radius,
    Rcpp::Named("radius_sd") = sqrt(var_radius),
    Rcpp::Named("node_rms_to_map") = node_rms,
    Rcpp::Named("draw_rms_rotation") = rms_rot,
    Rcpp::Named("draw_rms_full") = rms_full,
    Rcpp::Named("pair_rmse_to_map") = pair_rmse,
    Rcpp::Named("pair_cor_to_map") = pair_cor,
    Rcpp::Named("pair_ref_distance") = pair_ref
  );
}
