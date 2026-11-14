#ifndef JOINTVTX_COMMON_H
#define JOINTVTX_COMMON_H

/*
  Shared types, constants and helpers of the joint primary + secondary vertex
  fitter (AlephJointVtx).

  Eigen and the standard library only, so that the fitter core also builds
  outside the stage1 environment. Units: cm, GeV, rad.

  The helpers reproduce the semantics of the Python / NumPy reference
  implementation where results depend on them (floor modulo, NaN ordering,
  argmax, sorting, summation and LAPACK operation order). No function throws:
  a failure gives a false or non-zero status and NaN outputs.
*/

#include <array>
#include <cmath>
#include <limits>
#include <numeric>
#include <vector>
#include <algorithm>

#include <Eigen/Dense>

#include "aleph_units.h"

// fixed-size vectorizable Eigen types (Vec2, Vec4, Mat2, ...) inside std::vector
// rely on the over-aligned operator new of C++17
#if !defined(__cpp_aligned_new)
#error "jointvtx_common.h needs C++17 aligned new"
#endif

namespace FCCAnalyses {
namespace AlephJointVtx {

// ---------------------------------------------------------------- types
template <int N> using Vec = Eigen::Matrix<double, N, 1>;
template <int R, int C> using Mat = Eigen::Matrix<double, R, C>;
using Vec2 = Vec<2>;
using Vec3 = Vec<3>;
using Vec5 = Vec<5>;
using Mat2 = Mat<2, 2>;
using Mat3 = Mat<3, 3>;
using Mat5 = Mat<5, 5>;
using Mat35 = Mat<3, 5>;
using VecX = Eigen::VectorXd;
using MatX = Eigen::MatrixXd;

// element (i, j) = p[i * C + j], as a C-ordered NumPy array
template <int R, int C>
inline Mat<R, C> fromRowMajor(const double* p) {
  return Eigen::Map<const Eigen::Matrix<double, R, C, (C == 1 ? Eigen::ColMajor : Eigen::RowMajor)>>(p);
}

// ---------------------------------------------------------------- constants
// momentum constant [GeV cm]: pT = JV_PT_CONST / |omega|
constexpr double JV_PT_CONST = AlephUnits::kPtPerTeslaCm * AlephUnits::kBz;
// pT = JV_PT_HALF / |C| for the fitter curvature C = omega / 2
constexpr double JV_PT_HALF = 0.5 * JV_PT_CONST;

constexpr double JV_M_PI = AlephMasses::kPiCh;
constexpr double JV_M_P = AlephMasses::kProton;
constexpr double JV_M_KS = AlephMasses::kKs;
constexpr double JV_M_LAM = AlephMasses::kLambda;
constexpr double JV_M_E = AlephMasses::kElectron;

constexpr double JV_PI = 3.141592653589793;   // np.pi
constexpr double JV_TWO_PI = 2.0 * JV_PI;     // 2 * np.pi, exact doubling
constexpr double JV_NAN = std::numeric_limits<double>::quiet_NaN();
constexpr double JV_INF = std::numeric_limits<double>::infinity();

// V0 hypothesis code; Lambda: the proton is the positive daughter, anti-Lambda: the antiproton the negative one
enum V0Hyp : int { kV0K0S = 1, kV0Lambda = 2, kV0AntiLambda = 3, kV0Conversion = 4 };

// ---------------------------------------------------------------- non-finite convention
template <class D>
inline bool allFinite(const Eigen::DenseBase<D>& m) {
  for (Eigen::Index j = 0; j < m.cols(); ++j)
    for (Eigen::Index i = 0; i < m.rows(); ++i)
      if (!std::isfinite(m(i, j))) return false;
  return true;
}

template <class D>
inline void setNaN(Eigen::DenseBase<D>& m) { m.setConstant(JV_NAN); }

// ---------------------------------------------------------------- scalar semantics
// a % m of Python floats and np.remainder: result has the sign of m (it can equal m
// when a tiny negative a is rounded up); NaN for m == 0 or non-finite a
inline double pyMod(double a, double m) {
  double r = std::fmod(a, m);
  if (m == 0.0) return r;
  if (r != 0.0) {
    if ((m < 0.0) != (r < 0.0)) r += m;
  } else {
    r = std::copysign(0.0, m);
  }
  return r;
}

// (d + np.pi) % (2 * np.pi) - np.pi: angle difference wrapped to [-pi, pi]
inline double wrapPi(double d) { return pyMod(d + JV_PI, JV_TWO_PI) - JV_PI; }

// np.maximum / np.minimum: NaN if either argument is NaN, else the second on ties (signed zeros)
inline double npMaximum(double a, double b) { return (a > b || std::isnan(a)) ? a : b; }
inline double npMinimum(double a, double b) { return (a < b || std::isnan(a)) ? a : b; }

// np.clip(x, lo, hi) for floats: NaN propagates, max(x, lo) first, then min(., hi)
inline double npClip(double x, double lo, double hi) {
  const double t = std::isnan(x) ? x : (x > lo ? x : lo);
  return std::isnan(t) ? t : (t < hi ? t : hi);
}
// Python's built-in max(a, b) / min(a, b) (also inside numba) are std::max / std::min:
// the first argument is kept unless the second compares strictly greater / smaller.

// NumPy sort order: a before b; NaN after every number, NaNs equal
inline bool npLess(double a, double b) { return a < b || (b != b && a == a); }

// ---------------------------------------------------------------- array semantics
// np.argmax: index of the first NaN if any, else of the first maximum; -1 if empty
template <class V>
inline int argmaxNaN(const V& v) {
  const int n = static_cast<int>(v.size());
  if (n == 0) return -1;
  int k = 0;
  auto m = v[0];
  if (m != m) return 0;
  for (int i = 1; i < n; ++i) {
    if (!(v[i] <= m)) {
      m = v[i];
      k = i;
      if (m != m) break;
    }
  }
  return k;
}

// np.argsort(x, kind='stable') and list.sort: ascending, ties in index order, NaN last.
// For np.argsort(-x) pass the negated values: negation is exact.
template <class V>
inline std::vector<int> stableArgsort(const V& v) {
  std::vector<int> idx(static_cast<size_t>(v.size()));
  std::iota(idx.begin(), idx.end(), 0);
  std::stable_sort(idx.begin(), idx.end(), [&v](int a, int b) { return npLess(v[a], v[b]); });
  return idx;
}

// int(np.clip(np.searchsorted(edges, x, side='right') - 1, 0, len(edges) - 2)):
// bin [edges[k], edges[k+1]) clipped to the first / last bin; NaN goes to the last bin;
// -1 if fewer than two edges
template <class V>
inline int binRight(const V& edges, double x) {
  const int n = static_cast<int>(edges.size());
  if (n < 2) return -1;
  int lo = 0, hi = n;
  while (lo < hi) {
    const int mid = lo + (hi - lo) / 2;
    if (npLess(x, static_cast<double>(edges[mid]))) hi = mid;
    else lo = mid + 1;
  }
  int k = lo - 1;
  if (k < 0) k = 0;
  if (k > n - 2) k = n - 2;
  return k;
}

// np.sum of a 1-D float64 array: pairwise summation (blocks of 128, 8 accumulators)
inline double npSumRange(const double* a, int n) {
  if (n < 8) {
    double r = 0.0;
    for (int i = 0; i < n; ++i) r += a[i];
    return r;
  }
  if (n <= 128) {
    double r[8];
    for (int j = 0; j < 8; ++j) r[j] = a[j];
    int i = 8;
    for (; i < n - (n % 8); i += 8)
      for (int j = 0; j < 8; ++j) r[j] += a[i + j];
    double res = ((r[0] + r[1]) + (r[2] + r[3])) + ((r[4] + r[5]) + (r[6] + r[7]));
    for (; i < n; ++i) res += a[i];
    return res;
  }
  int n2 = n / 2;
  n2 -= n2 % 8;
  return npSumRange(a, n2) + npSumRange(a + n2, n - n2);
}
template <class V>
inline double npSum(const V& v) {
  std::vector<double> t(v.size());
  for (size_t i = 0; i < t.size(); ++i) t[i] = v[i];
  return npSumRange(t.data(), static_cast<int>(t.size()));
}

// left-to-right sum from 0: Python 3.11 sum() of floats, numba's arr.sum()
template <class V>
inline double seqSum(const V& v) {
  double r = 0.0;
  for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(v.size()); ++i) r += v[i];
  return r;
}

// a @ b of two 1-D float64 arrays with up to 13 elements (NumPy and numba, BLAS ddot):
// left to right from 0; longer vectors: detail::dotBlasN
template <class A, class B>
inline double dotSeq(const A& a, const B& b) {
  double r = 0.0;
  for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(a.size()); ++i) r += a[i] * b[i];
  return r;
}

// M @ x for a C-ordered matrix M with N <= 4 columns (NumPy and numba, BLAS gemv on the
// transpose): per row, products of even and odd columns summed separately, then added
template <int R, int N>
inline Vec<R> matVecC(const Mat<R, N>& M, const Vec<N>& x) {
  static_assert(N <= 4, "operation order verified for up to 4 columns");
  Vec<R> y;
  for (int i = 0; i < R; ++i) {
    double e = M(i, 0) * x[0];
    double o = N > 1 ? M(i, 1) * x[1] : 0.0;
    if (N > 2) e += M(i, 2) * x[2];
    if (N > 3) o += M(i, 3) * x[3];
    y[i] = N > 1 ? e + o : e;
  }
  return y;
}

// M @ x for an F-ordered matrix M (for instance a numba linalg result), and x @ M for
// a C-ordered one with M passed transposed: per row, left to right from 0
template <int R, int N>
inline Vec<R> matVecSeq(const Mat<R, N>& M, const Vec<N>& x) {
  Vec<R> y;
  for (int i = 0; i < R; ++i) {
    double s = 0.0;
    for (int k = 0; k < N; ++k) s += M(i, k) * x[k];
    y[i] = s;
  }
  return y;
}

// ---------------------------------------------------------------- linear algebra
// Operation order of OpenBLAS getrf / getrs / getri as called by np.linalg (and by
// numba's np.linalg): bit-identical for the sizes stated per function.
enum LinStatus : int {
  kLinOk = 0,
  kLinSingular = 1,   // exactly zero pivot: np.linalg raises LinAlgError (numpy and numba)
  kLinNonFinite = 2   // NaN / inf in the input: numba raises; numpy does not raise and returns non-finite
                      // values (finite ones for an inf on the diagonal)
};

// call sites where a failed fit or solve makes the reference raise (its run ends); the event is then marked failed.
// JOINTVTX_FAIL_AT(site) is false unless a test build defines it to force a failure at that site
enum FailSite : int {
  kFsTeardown = 1, kFsRefPvChi2, kFsSignif, kFsVdSupply, kFsCandFit, kFsCandMassFit, kFsMergeFit, kFsMergeRefit,
  kFsV0RefDen, kFsV0PairFit, kFsV0HfChi2, kFsWeightInv, kFsWeightChi2, kFsPvRefit, kFsSvRefit, kFsFinalChi2,
  kFsSeedPair, kFsGrowChi2Pv, kFsGrowChi2Seed, kFsGrowFit, kFsSplitFit, kFsSplitHalfFit, kFsOutSignif
};
#ifndef JOINTVTX_FAIL_AT
#define JOINTVTX_FAIL_AT(site) false
#endif

// partial-pivot LU, left-looking (column by column, dot products, multipliers scaled by
// the reciprocal pivot); pivot = first entry of largest magnitude. a holds L (unit, below
// the diagonal) and U; row k was swapped with row piv[k] at step k.
template <int N>
struct LU {
  Mat<N, N> a;
  std::array<int, N> piv{};
  int nswap = 0;
  LinStatus status = kLinNonFinite;
};

template <int N>
inline LU<N> luFactor(const Mat<N, N>& A) {
  LU<N> f;
  f.a = A;
  if (!allFinite(A)) {
    f.status = kLinNonFinite;
    return f;
  }
  f.status = kLinOk;
  Mat<N, N>& a = f.a;
  for (int j = 0; j < N; ++j) {
    for (int i = 1; i < j; ++i) {
      double s = 0.0;
      for (int k = 0; k < i; ++k) s += a(i, k) * a(k, j);
      a(i, j) -= s;
    }
    for (int i = j; i < N; ++i) {
      double s = 0.0;
      for (int k = 0; k < j; ++k) s += a(i, k) * a(k, j);
      a(i, j) -= s;
    }
    int p = j;
    double m = std::fabs(a(j, j));
    for (int i = j + 1; i < N; ++i)
      if (std::fabs(a(i, j)) > m) {
        m = std::fabs(a(i, j));
        p = i;
      }
    f.piv[j] = p;
    if (p != j) {
      a.row(j).swap(a.row(p));
      ++f.nswap;
    }
    if (a(j, j) != 0.0) {
      const double r = 1.0 / a(j, j);
      for (int i = j + 1; i < N; ++i) a(i, j) *= r;
    } else {
      f.status = kLinSingular;
    }
  }
  return f;
}

// getrs with one right-hand side: column-oriented substitutions, division by the pivots
template <int N>
inline Vec<N> luSolve(const LU<N>& f, const Vec<N>& b) {
  Vec<N> x = b;
  for (int k = 0; k < N; ++k)
    if (f.piv[k] != k) std::swap(x[k], x[f.piv[k]]);
  for (int j = 0; j < N; ++j)
    for (int i = j + 1; i < N; ++i) x[i] -= x[j] * f.a(i, j);
  for (int j = N - 1; j >= 0; --j) {
    x[j] /= f.a(j, j);
    for (int i = 0; i < j; ++i) x[i] -= x[j] * f.a(i, j);
  }
  return x;
}

// np.linalg.solve(A, b) with 1-D b, in Python and in numba code (exact for N <= 5)
template <int N>
inline LinStatus npSolve(const Mat<N, N>& A, const Vec<N>& b, Vec<N>& x) {
  const LU<N> f = luFactor(A);
  if (f.status != kLinOk || !allFinite(b)) {
    setNaN(x);
    return f.status != kLinOk ? f.status : kLinNonFinite;
  }
  x = luSolve(f, b);
  return kLinOk;
}

// np.linalg.inv(A) in Python code: gesv against the identity, the substitutions multiply
// by the reciprocal pivots (exact for N <= 3)
template <int N>
inline LinStatus npInv(const Mat<N, N>& A, Mat<N, N>& Ai) {
  const LU<N> f = luFactor(A);
  if (f.status != kLinOk) {
    setNaN(Ai);
    return f.status;
  }
  for (int c = 0; c < N; ++c) {
    Vec<N> x = Vec<N>::Zero();
    x[c] = 1.0;
    for (int k = 0; k < N; ++k)
      if (f.piv[k] != k) std::swap(x[k], x[f.piv[k]]);
    for (int j = 0; j < N; ++j)
      for (int i = j + 1; i < N; ++i) x[i] -= x[j] * f.a(i, j);
    for (int j = N - 1; j >= 0; --j) {
      x[j] *= 1.0 / f.a(j, j);
      for (int i = 0; i < j; ++i) x[i] -= x[j] * f.a(i, j);
    }
    Ai.col(c) = x;
  }
  return kLinOk;
}

// np.linalg.inv(A) inside numba code: getrf + getri (inverse of U, then the L solve with
// dot products, then the column swaps); differs from npInv in the last bits (exact for N <= 3)
template <int N>
inline LinStatus nbInv(const Mat<N, N>& A, Mat<N, N>& Ai) {
  const LU<N> f = luFactor(A);
  if (f.status != kLinOk) {
    setNaN(Ai);
    return f.status;
  }
  Mat<N, N> a = f.a;
  for (int j = 0; j < N; ++j) {
    a(j, j) = 1.0 / a(j, j);
    const double ajj = -a(j, j);
    for (int jj = 0; jj < j; ++jj) {
      const double t = a(jj, j);
      if (t != 0.0) {
        for (int i = 0; i < jj; ++i) a(i, j) += t * a(i, jj);
        a(jj, j) = t * a(jj, jj);
      }
    }
    for (int i = 0; i < j; ++i) a(i, j) *= ajj;
  }
  std::array<double, N> work{};
  for (int j = N - 1; j >= 0; --j) {
    for (int i = j + 1; i < N; ++i) {
      work[i] = a(i, j);
      a(i, j) = 0.0;
    }
    if (j < N - 1) {
      for (int i = 0; i < N; ++i) {
        double s = 0.0;
        for (int jj = j + 1; jj < N; ++jj) s += a(i, jj) * work[jj];
        a(i, j) += -1.0 * s;
      }
    }
  }
  for (int j = N - 2; j >= 0; --j) {
    const int jp = f.piv[j];
    if (jp != j) a.col(j).swap(a.col(jp));
  }
  Ai = a;
  return kLinOk;
}

// np.linalg.det(A), numpy and numba: sign * exp(sum of log |U_kk|); 0 if singular,
// NaN for non-finite input (exact for N <= 5)
template <int N>
inline double npDet(const Mat<N, N>& A) {
  const LU<N> f = luFactor(A);
  if (f.status == kLinNonFinite) return JV_NAN;
  if (f.status == kLinSingular) return 0.0;
  double sgn = (f.nswap % 2) ? -1.0 : 1.0;
  double acc = 0.0;
  for (int k = 0; k < N; ++k) {
    double v = f.a(k, k);
    if (v < 0.0) {
      sgn = -sgn;
      v = -v;
    }
    acc += std::log(v);
  }
  return sgn * std::exp(acc);
}

// ---------------------------------------------------------------- reference operation orders
// Products, powers, linear algebra and the 3x3 eigen-decomposition in the operation order of the reference
// (NumPy 1.26 and numba with OpenBLAS 0.3.20 built for NEHALEM: SSE2 kernels, no fused multiply-add).
namespace detail {

// numba's u ** 4 and u ** 3 for an integer literal exponent: exponentiation by squaring
inline double pow4(double u) {
  const double u2 = u * u;
  return u2 * u2;
}
inline double pow3(double u) { return u * (u * u); }

// x ** 2 of a Python / NumPy float: libm pow(x, 2.0). The exponent is read through a volatile so that the compiler
// does not turn the call into x * x, which differs in the last bit for about 1 in 1000 arguments.
inline double powSquare(double x) {
  volatile double e = 2.0;
  return std::pow(x, static_cast<double>(e));
}

// np.dot of two contiguous float64 vectors of any length (OpenBLAS ddot, SSE2 kernel): blocks of 16 elements in
// 8 lane accumulators, lanes added as ((0 + 2) + (4 + 6)) + ((1 + 3) + (5 + 7)), then the rest left to right
inline double dotBlasN(const double* a, const double* b, int n) {
  const int nb = n & -16;
  double s = 0.;
  if (nb > 0) {
    double L[8] = {0., 0., 0., 0., 0., 0., 0., 0.};
    for (int j = 0; j < nb; ++j) L[j % 8] += a[j] * b[j];
    s = ((L[0] + L[2]) + (L[4] + L[6])) + ((L[1] + L[3]) + (L[5] + L[7]));
  }
  for (int j = nb; j < n; ++j) s += a[j] * b[j];
  return s;
}

// element (i, j) of an R x M matrix product with inner size 5 (numba and NumPy `@`, any memory order of the
// operands; R, M >= 2), from its five products p[k] = A(i, k) B(k, j). Last column of an odd M: rows inside full
// blocks of 8 left to right, the other rows ((p0 + p2) + p4) + (p1 + p3); other columns: rows inside full blocks of 4
// left to right, the other rows (((p0 + p2) + p1) + p3) + p4.
inline double sum5Blas(const double* p, int i, int j, int R, int M) {
  const double seq = (((p[0] + p[1]) + p[2]) + p[3]) + p[4];
  if (M % 2 == 1 && j == M - 1) return i < (R / 8) * 8 ? seq : ((p[0] + p[2]) + p[4]) + (p[1] + p[3]);
  if (i < (R / 4) * 4) return seq;
  return (((p[0] + p[2]) + p[1]) + p[3]) + p[4];
}

// A @ cov @ A.T (A 3x5, cov 5x5; any memory order), evaluated as (A @ cov) @ A.T
inline Mat3 sandwich(const Mat35& A, const Mat5& cov) {
  double a[3][5], c[5][5], t[3][5], p[5];
  for (int i = 0; i < 3; ++i)
    for (int k = 0; k < 5; ++k) a[i][k] = A(i, k);
  for (int i = 0; i < 5; ++i)
    for (int k = 0; k < 5; ++k) c[i][k] = cov(i, k);
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 5; ++j) {
      for (int k = 0; k < 5; ++k) p[k] = a[i][k] * c[k][j];
      t[i][j] = sum5Blas(p, i, j, 3, 5);
    }
  Mat3 M;
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) {
      for (int k = 0; k < 5; ++k) p[k] = t[i][k] * a[j][k];
      M(i, j) = sum5Blas(p, i, j, 3, 3);
    }
  return M;
}

// u @ C @ u for a 3x3 C whose memory order in the reference is given: (u @ C) sums left to right per column for a
// C-ordered C and in even / odd terms for a Fortran-ordered one (numba linalg results, sums of Fortran-ordered
// arrays); then the dot product left to right
inline double quadForm(const Vec3& u, const Mat3& C, bool fortran) {
  const Vec3 y = fortran ? matVecC<3, 3>(C.transpose(), u) : matVecSeq<3, 3>(C.transpose(), u);
  return dotSeq(y, u);
}

// sum of the products p[k], k < m, of element (i, j) of the 3 x m @ m x 3 product (left operand Fortran-ordered,
// right C-ordered; OpenBLAS dgemm, SSE2 kernels; exact for m <= 200): nb = 4 (m / 4) terms in lanes, then the rest
// in order
inline double gemm3mBlas(const double* p, int m, int i, int j) {
  const int nb = 4 * (m / 4);
  if (j < 2) {
    double s = 0.;
    if (nb > 0) {
      double l0 = 0., l1 = 0.;
      for (int k = 0; k < nb; ++k) {
        if (k % 2 == 0)
          l0 += p[k];
        else if (k < nb - 1)
          l1 += p[k];
      }
      s = l0 + l1;
      s += p[nb - 1];
    }
    for (int k = nb; k < m; ++k) s += p[k];
    return s;
  }
  if (i < 2) {
    double l0 = 0., l1 = 0.;
    for (int k = 0; k < nb; ++k) {
      if (k % 2 == 0)
        l0 += p[k];
      else
        l1 += p[k];
    }
    for (int k = nb; k < m; ++k) l0 += p[k];
    return l0 + l1;
  }
  double a[4] = {0., 0., 0., 0.};
  for (int k = 0; k < nb; ++k) a[k % 4] += p[k];
  double s0 = a[0] + a[2];
  const double s1 = a[1] + a[3];
  for (int k = nb; k < m; ++k) s0 += p[k];
  return s0 + s1;
}

// row i of a Fortran-ordered 3 x m matrix times an m-vector (OpenBLAS dgemv_n, three-row tail): products p[k]
// added in pairs within blocks of 4, then one by one
inline double gemv3mBlas(const double* p, int m) {
  double t = 0.;
  int k = 0;
  for (; k + 4 <= m; k += 4) {
    t += p[k] + p[k + 1];
    t += p[k + 2] + p[k + 3];
  }
  for (; k < m; ++k) t += p[k];
  return t;
}

// (n x 3, C-ordered) @ 3-vector rows: matVecC order, or left to right for a single row (dot)
inline double rowDot3(const double* h, const double* x, bool single) {
  if (single) return (h[0] * x[0] + h[1] * x[1]) + h[2] * x[2];
  return (h[0] * x[0] + h[2] * x[2]) + h[1] * x[1];
}

// v @ J for a C-ordered 3x5 J (BLAS gemv on the transpose): per column, left to right
inline Vec5 vecMat35(const Vec3& v, const Mat35& J) {
  Vec5 y;
  for (int j = 0; j < 5; ++j) y[j] = (v[0] * J(0, j) + v[1] * J(1, j)) + v[2] * J(2, j);
  return y;
}

// v @ C @ v for a C-ordered 5x5 C: (v @ C) per column left to right, then the dot product left to right
inline double quadForm5C(const Vec5& v, const Mat5& C) {
  Vec5 y;
  for (int j = 0; j < 5; ++j) {
    double s = 0.;
    for (int k = 0; k < 5; ++k) s += v[k] * C(k, j);
    y[j] = s;
  }
  return dotSeq(y, v);
}

// np.cross(a, b) of two 3-vectors (NumPy's order of the products)
inline Vec3 npCross(const Vec3& a, const Vec3& b) {
  Vec3 c;
  c[0] = a[1] * b[2] - a[2] * b[1];
  c[1] = a[2] * b[0] - a[0] * b[2];
  c[2] = a[0] * b[1] - a[1] * b[0];
  return c;
}

// Euclidean norm of a 3-vector: sqrt(v @ v)
inline double norm3(const Vec3& v) { return std::sqrt(dotSeq(v, v)); }

// luFactor / npSolve / nbInv on plain arrays, same operation order
template <int N>
inline LinStatus luFactorA(double a[N][N], int piv[N]) {
  for (int i = 0; i < N; ++i)
    for (int j = 0; j < N; ++j)
      if (!std::isfinite(a[i][j])) return kLinNonFinite;
  LinStatus st = kLinOk;
  for (int j = 0; j < N; ++j) {
    for (int i = 1; i < j; ++i) {
      double s = 0.0;
      for (int k = 0; k < i; ++k) s += a[i][k] * a[k][j];
      a[i][j] -= s;
    }
    for (int i = j; i < N; ++i) {
      double s = 0.0;
      for (int k = 0; k < j; ++k) s += a[i][k] * a[k][j];
      a[i][j] -= s;
    }
    int p = j;
    double m = std::fabs(a[j][j]);
    for (int i = j + 1; i < N; ++i)
      if (std::fabs(a[i][j]) > m) {
        m = std::fabs(a[i][j]);
        p = i;
      }
    piv[j] = p;
    if (p != j)
      for (int k = 0; k < N; ++k) std::swap(a[j][k], a[p][k]);
    if (a[j][j] != 0.0) {
      const double r = 1.0 / a[j][j];
      for (int i = j + 1; i < N; ++i) a[i][j] *= r;
    } else {
      st = kLinSingular;
    }
  }
  return st;
}

// np.linalg.solve(A, b) (numba); x = NaN on failure
template <int N>
inline LinStatus npSolveA(const double A[N][N], const double* b, double* x) {
  double a[N][N];
  int piv[N];
  for (int i = 0; i < N; ++i)
    for (int j = 0; j < N; ++j) a[i][j] = A[i][j];
  LinStatus st = luFactorA<N>(a, piv);
  if (st == kLinOk)
    for (int i = 0; i < N; ++i)
      if (!std::isfinite(b[i])) st = kLinNonFinite;
  if (st != kLinOk) {
    for (int i = 0; i < N; ++i) x[i] = JV_NAN;
    return st;
  }
  for (int i = 0; i < N; ++i) x[i] = b[i];
  for (int k = 0; k < N; ++k)
    if (piv[k] != k) std::swap(x[k], x[piv[k]]);
  for (int j = 0; j < N; ++j)
    for (int i = j + 1; i < N; ++i) x[i] -= x[j] * a[i][j];
  for (int j = N - 1; j >= 0; --j) {
    x[j] /= a[j][j];
    for (int i = 0; i < j; ++i) x[i] -= x[j] * a[i][j];
  }
  return kLinOk;
}

// np.linalg.inv(A) inside numba code; Ai = NaN on failure
template <int N>
inline LinStatus nbInvA(const double A[N][N], double Ai[N][N]) {
  double a[N][N];
  int piv[N];
  for (int i = 0; i < N; ++i)
    for (int j = 0; j < N; ++j) a[i][j] = A[i][j];
  const LinStatus st = luFactorA<N>(a, piv);
  if (st != kLinOk) {
    for (int i = 0; i < N; ++i)
      for (int j = 0; j < N; ++j) Ai[i][j] = JV_NAN;
    return st;
  }
  for (int j = 0; j < N; ++j) {
    a[j][j] = 1.0 / a[j][j];
    const double ajj = -a[j][j];
    for (int jj = 0; jj < j; ++jj) {
      const double t = a[jj][j];
      if (t != 0.0) {
        for (int i = 0; i < jj; ++i) a[i][j] += t * a[i][jj];
        a[jj][j] = t * a[jj][jj];
      }
    }
    for (int i = 0; i < j; ++i) a[i][j] *= ajj;
  }
  double work[N];
  for (int i = 0; i < N; ++i) work[i] = 0.0;
  for (int j = N - 1; j >= 0; --j) {
    for (int i = j + 1; i < N; ++i) {
      work[i] = a[i][j];
      a[i][j] = 0.0;
    }
    if (j < N - 1) {
      for (int i = 0; i < N; ++i) {
        double s = 0.0;
        for (int jj = j + 1; jj < N; ++jj) s += a[i][jj] * work[jj];
        a[i][j] += -1.0 * s;
      }
    }
  }
  for (int j = N - 2; j >= 0; --j) {
    const int jp = piv[j];
    if (jp != j)
      for (int i = 0; i < N; ++i) std::swap(a[i][j], a[i][jp]);
  }
  for (int i = 0; i < N; ++i)
    for (int j = 0; j < N; ++j) Ai[i][j] = a[i][j];
  return kLinOk;
}

// ---- symmetric 3x3 eigen-decomposition in the operation order of LAPACK dsyevd (jobz 'V',
// uplo 'L'; reference LAPACK 3.9 with the OpenBLAS level-1/2 kernels): dsytd2, dsteqr, dorm2r.
// Arrays are 1-based as in LAPACK. Scaling of matrices with norm below 1e-120 or above 1e140
// is not reproduced.
constexpr double kLaSafmin = std::numeric_limits<double>::min();
constexpr double kLaEps = 0.5 * std::numeric_limits<double>::epsilon();

inline double fsign(double a, double b) { return b >= 0 ? std::fabs(a) : -std::fabs(a); }

inline double dlapy2(double x, double y) {
  const double xa = std::fabs(x), ya = std::fabs(y);
  const double w = std::max(xa, ya), z = std::min(xa, ya);
  if (z == 0.0) return w;
  const double q = z / w;
  return w * std::sqrt(1.0 + q * q);
}

inline void dlartg(double f, double g, double& cs, double& sn, double& r) {
  const double safmn2 = std::ldexp(1.0, -484), safmx2 = 1.0 / safmn2;
  if (g == 0.0) {
    cs = 1.0;
    sn = 0.0;
    r = f;
    return;
  }
  if (f == 0.0) {
    cs = 0.0;
    sn = 1.0;
    r = g;
    return;
  }
  double f1 = f, g1 = g;
  double scale = std::max(std::fabs(f1), std::fabs(g1));
  if (scale >= safmx2) {
    int cnt = 0;
    do {
      ++cnt;
      f1 *= safmn2;
      g1 *= safmn2;
      scale = std::max(std::fabs(f1), std::fabs(g1));
    } while (scale >= safmx2 && cnt < 20);
    r = std::sqrt(f1 * f1 + g1 * g1);
    cs = f1 / r;
    sn = g1 / r;
    for (int i = 0; i < cnt; ++i) r *= safmx2;
  } else if (scale <= safmn2) {
    int cnt = 0;
    do {
      ++cnt;
      f1 *= safmx2;
      g1 *= safmx2;
      scale = std::max(std::fabs(f1), std::fabs(g1));
    } while (scale <= safmn2);
    r = std::sqrt(f1 * f1 + g1 * g1);
    cs = f1 / r;
    sn = g1 / r;
    for (int i = 0; i < cnt; ++i) r *= safmn2;
  } else {
    r = std::sqrt(f1 * f1 + g1 * g1);
    cs = f1 / r;
    sn = g1 / r;
  }
  if (std::fabs(f) > std::fabs(g) && cs < 0) {
    cs = -cs;
    sn = -sn;
    r = -r;
  }
}

inline void dlaev2(double a, double b, double c, double& rt1, double& rt2, double& cs1, double& sn1) {
  const double sm = a + c, df = a - c, adf = std::fabs(df), tb = b + b, ab = std::fabs(tb);
  double acmx, acmn;
  if (std::fabs(a) > std::fabs(c)) {
    acmx = a;
    acmn = c;
  } else {
    acmx = c;
    acmn = a;
  }
  double rt;
  if (adf > ab) {
    const double q = ab / adf;
    rt = adf * std::sqrt(1.0 + q * q);
  } else if (adf < ab) {
    const double q = adf / ab;
    rt = ab * std::sqrt(1.0 + q * q);
  } else {
    rt = ab * std::sqrt(2.0);
  }
  int sgn1;
  if (sm < 0) {
    rt1 = 0.5 * (sm - rt);
    sgn1 = -1;
    rt2 = (acmx / rt1) * acmn - (b / rt1) * b;
  } else if (sm > 0) {
    rt1 = 0.5 * (sm + rt);
    sgn1 = 1;
    rt2 = (acmx / rt1) * acmn - (b / rt1) * b;
  } else {
    rt1 = 0.5 * rt;
    rt2 = -0.5 * rt;
    sgn1 = 1;
  }
  int sgn2;
  double cs;
  if (df >= 0) {
    cs = df + rt;
    sgn2 = 1;
  } else {
    cs = df - rt;
    sgn2 = -1;
  }
  if (std::fabs(cs) > ab) {
    const double ct = -tb / cs;
    sn1 = 1.0 / std::sqrt(1.0 + ct * ct);
    cs1 = ct * sn1;
  } else if (ab == 0) {
    cs1 = 1.0;
    sn1 = 0.0;
  } else {
    const double tn = -cs / tb;
    cs1 = 1.0 / std::sqrt(1.0 + tn * tn);
    sn1 = tn * cs1;
  }
  if (sgn1 == sgn2) {
    const double tn = cs1;
    cs1 = -sn1;
    sn1 = tn;
  }
}

// dlasr side 'R', pivot 'V' on columns col0 .. col0 + mm - 1 of Z; rotation j = (c[j], s[j]), j = 1 .. mm - 1
inline void dlasrRV(bool backward, int mm, const double* c, const double* s, double Z[4][4], int col0) {
  for (int jj = 1; jj < mm; ++jj) {
    const int j = backward ? mm - jj : jj;
    const double ct = c[j], st = s[j];
    if (ct != 1.0 || st != 0.0) {
      for (int i = 1; i <= 3; ++i) {
        const double temp = Z[i][col0 + j];
        Z[i][col0 + j] = ct * temp - st * Z[i][col0 + j - 1];
        Z[i][col0 + j - 1] = st * temp + ct * Z[i][col0 + j - 1];
      }
    }
  }
}

// dsteqr (compz 'I') for n = 3: d[1..3], e[1..2], Z starts as the identity; 0 = converged
inline int dsteqr3(double d[4], double e[4], double Z[4][4]) {
  const int n = 3, nm1 = 2, nmaxit = n * 30;
  const double eps2 = kLaEps * kLaEps;
  double wc[4] = {0., 0., 0., 0.}, ws[4] = {0., 0., 0., 0.};
  int jtot = 0, l1 = 1;
  while (l1 <= n) {
    if (l1 > 1) e[l1 - 1] = 0.0;
    int m = n;
    for (int q = l1; q <= nm1; ++q) {
      const double tst = std::fabs(e[q]);
      if (tst == 0.0) {
        m = q;
        break;
      }
      if (tst <= (std::sqrt(std::fabs(d[q])) * std::sqrt(std::fabs(d[q + 1]))) * kLaEps) {
        e[q] = 0.0;
        m = q;
        break;
      }
    }
    int l = l1;
    const int lsv = l;
    int lend = m;
    const int lendsv = lend;
    l1 = m + 1;
    if (lend == l) continue;
    double anorm = 0.0;
    for (int i = l; i <= lend; ++i) anorm = std::max(anorm, std::fabs(d[i]));
    for (int i = l; i < lend; ++i) anorm = std::max(anorm, std::fabs(e[i]));
    if (anorm == 0.0) continue;
    if (std::fabs(d[lend]) < std::fabs(d[l])) {
      lend = lsv;
      l = lendsv;
    }
    if (lend > l) {
      while (true) {
        m = lend;
        for (int q = l; q < lend; ++q) {
          const double tst = std::fabs(e[q]) * std::fabs(e[q]);
          if (tst <= (eps2 * std::fabs(d[q])) * std::fabs(d[q + 1]) + kLaSafmin) {
            m = q;
            break;
          }
        }
        if (m < lend) e[m] = 0.0;
        double p = d[l];
        if (m == l) {
          d[l] = p;
          ++l;
          if (l <= lend) continue;
          break;
        }
        if (m == l + 1) {
          double rt1, rt2, c, s;
          dlaev2(d[l], e[l], d[l + 1], rt1, rt2, c, s);
          const double cc[2] = {0., c}, ss[2] = {0., s};
          dlasrRV(true, 2, cc, ss, Z, l);
          d[l] = rt1;
          d[l + 1] = rt2;
          e[l] = 0.0;
          l += 2;
          if (l <= lend) continue;
          break;
        }
        if (jtot == nmaxit) break;
        ++jtot;
        double g = (d[l + 1] - p) / (2.0 * e[l]);
        double r = dlapy2(g, 1.0);
        g = d[m] - p + (e[l] / (g + fsign(r, g)));
        double s = 1.0, c = 1.0;
        p = 0.0;
        for (int i = m - 1; i >= l; --i) {
          const double f = s * e[i], b = c * e[i];
          dlartg(g, f, c, s, r);
          if (i != m - 1) e[i + 1] = r;
          g = d[i + 1] - p;
          r = (d[i] - g) * s + 2.0 * c * b;
          p = s * r;
          d[i + 1] = g + p;
          g = c * r - b;
          wc[i - l + 1] = c;
          ws[i - l + 1] = -s;
        }
        dlasrRV(true, m - l + 1, wc, ws, Z, l);
        d[l] = d[l] - p;
        e[l] = g;
      }
    } else {
      while (true) {
        m = lend;
        for (int q = l; q > lend; --q) {
          const double tst = std::fabs(e[q - 1]) * std::fabs(e[q - 1]);
          if (tst <= (eps2 * std::fabs(d[q])) * std::fabs(d[q - 1]) + kLaSafmin) {
            m = q;
            break;
          }
        }
        if (m > lend) e[m - 1] = 0.0;
        double p = d[l];
        if (m == l) {
          d[l] = p;
          --l;
          if (l >= lend) continue;
          break;
        }
        if (m == l - 1) {
          double rt1, rt2, c, s;
          dlaev2(d[l - 1], e[l - 1], d[l], rt1, rt2, c, s);
          const double cc[2] = {0., c}, ss[2] = {0., s};
          dlasrRV(false, 2, cc, ss, Z, l - 1);
          d[l - 1] = rt1;
          d[l] = rt2;
          e[l - 1] = 0.0;
          l -= 2;
          if (l >= lend) continue;
          break;
        }
        if (jtot == nmaxit) break;
        ++jtot;
        double g = (d[l - 1] - p) / (2.0 * e[l - 1]);
        double r = dlapy2(g, 1.0);
        g = d[m] - p + (e[l - 1] / (g + fsign(r, g)));
        double s = 1.0, c = 1.0;
        p = 0.0;
        for (int i = m; i < l; ++i) {
          const double f = s * e[i], b = c * e[i];
          dlartg(g, f, c, s, r);
          if (i != m) e[i - 1] = r;
          g = d[i] - p;
          r = (d[i + 1] - g) * s + 2.0 * c * b;
          p = s * r;
          d[i] = g + p;
          g = c * r - b;
          wc[i - m + 1] = c;
          ws[i - m + 1] = s;
        }
        dlasrRV(false, l - m + 1, wc, ws, Z, m);
        d[l] = d[l] - p;
        e[l - 1] = g;
      }
    }
    if (jtot >= nmaxit) return 1;
  }
  for (int ii = 2; ii <= n; ++ii) {
    const int i = ii - 1;
    int k = i;
    double p = d[i];
    for (int j = ii; j <= n; ++j)
      if (d[j] < p) {
        k = j;
        p = d[j];
      }
    if (k != i) {
      d[k] = d[i];
      d[i] = p;
      for (int r = 1; r <= n; ++r) std::swap(Z[r][i], Z[r][k]);
    }
  }
  return 0;
}

// np.linalg.eigh of a symmetric 3x3 (lower triangle read): ascending eigenvalues w, eigenvectors
// in the columns of U; 0 = converged
inline int eigh3(const Mat3& S, Vec3& w, Mat3& U) {
  const int n = 3;
  double a[4][4], d[4] = {0., 0., 0., 0.}, e[4] = {0., 0., 0., 0.}, tau[4] = {0., 0., 0., 0.};
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) a[i + 1][j + 1] = S(i, j);
  // dsytd2, lower: one reflector for n = 3 (the second has length 1)
  for (int i = 1; i < n; ++i) {
    double taui = 0.0;
    if (n - i > 1) {
      double alpha = a[i + 1][i];
      const double xnorm = std::fabs(a[i + 2][i]);
      if (xnorm != 0.0) {
        const double beta = -fsign(dlapy2(alpha, xnorm), alpha);
        taui = (beta - alpha) / beta;
        a[i + 2][i] = a[i + 2][i] * (1.0 / (alpha - beta));
        alpha = beta;
      }
      a[i + 1][i] = alpha;
    }
    e[i] = a[i + 1][i];
    if (taui != 0.0) {
      a[i + 1][i] = 1.0;
      const int m = n - i;
      double v[3], y[3] = {0., 0., 0.}, Am[3][3];
      for (int k = 0; k < m; ++k) v[k] = a[i + 1 + k][i];
      for (int r = 0; r < m; ++r)
        for (int c = 0; c < m; ++c) Am[r][c] = a[i + 1 + r][i + 1 + c];
      for (int j = 0; j < m; ++j) {
        const double t1 = taui * v[j];
        double t2 = 0.0;
        y[j] = y[j] + t1 * Am[j][j];
        for (int q = j + 1; q < m; ++q) {
          y[q] = y[q] + t1 * Am[q][j];
          t2 = t2 + Am[q][j] * v[q];
        }
        y[j] = y[j] + taui * t2;
      }
      double dt = 0.0;
      for (int k = 0; k < m; ++k) dt += y[k] * v[k];
      const double alpha2 = -(0.5 * taui * dt);
      for (int k = 0; k < m; ++k) y[k] = y[k] + alpha2 * v[k];
      for (int j = 0; j < m; ++j) {
        for (int k = j; k < m; ++k) Am[k][j] = Am[k][j] + (-v[j]) * y[k];
        for (int k = j; k < m; ++k) Am[k][j] = Am[k][j] + (-y[j]) * v[k];
      }
      for (int r = 0; r < m; ++r)
        for (int c = 0; c < m; ++c) a[i + 1 + r][i + 1 + c] = Am[r][c];
      a[i + 1][i] = e[i];
    }
    d[i] = a[i][i];
    tau[i] = taui;
  }
  d[n] = a[n][n];
  double Z[4][4];
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 4; ++j) Z[i][j] = (i == j) ? 1.0 : 0.0;
  const int info = dsteqr3(d, e, Z);
  // dormtr ('L', 'L', 'N') through dorm2r: reflector k acts on rows k + 1 .. n, applied for k = n - 1 .. 1
  for (int k = n - 1; k >= 1; --k) {
    if (tau[k] == 0.0) continue;
    const int m = n - k;
    double v[3] = {0., 0., 0.};
    v[0] = 1.0;
    for (int q = 1; q < m; ++q) v[q] = a[k + 1 + q][k];
    for (int j = 1; j <= n; ++j) {
      double ev = Z[k + 1][j] * v[0], od = 0.0;
      for (int q = 1; q < m; ++q) {
        if (q % 2 == 0) ev = ev + Z[k + 1 + q][j] * v[q];
        else od = (q == 1) ? Z[k + 1 + q][j] * v[q] : od + Z[k + 1 + q][j] * v[q];
      }
      const double wj = m > 1 ? ev + od : ev;
      const double al = (-tau[k]) * wj;
      for (int q = 0; q < m; ++q) Z[k + 1 + q][j] = Z[k + 1 + q][j] + al * v[q];
    }
  }
  for (int i = 0; i < 3; ++i) {
    w[i] = d[i + 1];
    for (int j = 0; j < 3; ++j) U(i, j) = Z[i + 1][j + 1];
  }
  return info;
}

}  // namespace detail

}  // namespace AlephJointVtx
}  // namespace FCCAnalyses

#endif
