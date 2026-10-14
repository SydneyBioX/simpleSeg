#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <cstdint>
using namespace Rcpp;

// Normalised 1D Gaussian of EBImage::gblur(): 2 * ceiling(3 * sigma) + 1 taps.
static std::vector<double> gaussianKernel(double sigma, int &r) {
  r = (int) std::ceil(3 * sigma);
  std::vector<double> g(2 * r + 1);
  double s = 0;
  for (int k = 0; k <= 2 * r; k++) {
    double x = k - r;
    g[k] = std::exp(-(x * x) / (2 * sigma * sigma));
    s += g[k];
  }
  for (double &v : g) v /= s;
  return g;
}

static inline int wrap(int i, int n) {
  i %= n;
  return i < 0 ? i + n : i;
}

// Separable Gaussian blur of a matrix.
// circular = true: same size as x, wrapping around the edges (as EBImage's
// filter2() does). circular = false: only the pixels at least r from every
// edge, an (nr - 2r) x (nc - 2r) matrix. A pixel's value is computed from the
// same inputs in the same order in both modes, so a window read with an r
// pixel margin gives exactly the values of the whole image.
// [[Rcpp::export(name = ".gaussianBlurCpp")]]
NumericMatrix gaussianBlur(NumericMatrix x, double sigma, bool circular) {
  int r;
  std::vector<double> g = gaussianKernel(sigma, r);
  int K = 2 * r + 1, nr = x.nrow(), nc = x.ncol();
  int onr = circular ? nr : nr - 2 * r, onc = circular ? nc : nc - 2 * r;
  if (onr < 1 || onc < 1) stop("image is smaller than the blur kernel");

  // pass 1: along the first dimension, for every column
  std::vector<double> tmp((size_t) onr * nc), buf(onr + 2 * r);
  for (int j = 0; j < nc; j++) {
    const double *col = &x[(size_t) j * nr];
    if (circular) {
      for (int i = 0; i < nr + 2 * r; i++) buf[i] = col[wrap(i - r, nr)];
    } else {
      for (int i = 0; i < nr; i++) buf[i] = col[i];
    }
    double *out = &tmp[(size_t) j * onr];
    for (int i = 0; i < onr; i++) {
      double s = 0;
      for (int k = 0; k < K; k++) s += g[k] * buf[i + k];
      out[i] = s;
    }
  }

  // pass 2: along the second dimension
  NumericMatrix res(onr, onc);
  std::vector<double> acc(onr);
  for (int j = 0; j < onc; j++) {
    std::fill(acc.begin(), acc.end(), 0.0);
    for (int k = 0; k < K; k++) {
      int jj = circular ? wrap(j + k - r, nc) : j + k;
      const double *src = &tmp[(size_t) jj * onr];
      double w = g[k];
      for (int i = 0; i < onr; i++) acc[i] += w * src[i];
    }
    std::copy(acc.begin(), acc.end(), &res[(size_t) j * onr]);
  }
  return res;
}

// Pixel count and sums of the pixel positions for labels 1..n of a label
// matrix; label 0 and labels above n are ignored. Positions are those in the
// whole image, label[1, 1] being pixel (x0, y0), so the sums are exact
// integers that identify a cell in any tile containing it. Returns n x 3.
// [[Rcpp::export(name = ".labelMomentsCpp")]]
NumericMatrix labelMoments(IntegerMatrix label, int n, double x0, double y0) {
  NumericMatrix out(n, 3);
  int nr = label.nrow(), nc = label.ncol();
  for (int j = 0; j < nc; j++) {
    for (int i = 0; i < nr; i++) {
      int l = label[(size_t) j * nr + i];
      if (l > 0 && l <= n) {
        out(l - 1, 0) += 1;
        out(l - 1, 1) += i + x0;
        out(l - 1, 2) += j + y0;
      }
    }
  }
  return out;
}

// Sums of `values` for labels 1..n of `label` (same length), in pixel order.
// [[Rcpp::export(name = ".labelSumsCpp")]]
NumericVector labelSums(IntegerVector label, NumericVector values, int n) {
  NumericVector out(n);
  R_xlen_t len = label.size();
  for (R_xlen_t p = 0; p < len; p++) {
    int l = label[p];
    if (l > 0 && l <= n) out[l - 1] += values[p];
  }
  return out;
}

// A number in [0, 1) from a pixel's position in the whole image.
static inline double pixelHash(uint32_t x, uint32_t y) {
  uint64_t z = ((uint64_t) x << 32) ^ (uint64_t) y;
  z += 0x9e3779b97f4a7c15ULL;
  z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
  z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
  z ^= z >> 31;
  return (double) (z >> 11) * 0x1.0p-53;
}

// Adds eps * hash(x, y) to the positive pixels of h, where (x, y) is the
// pixel's position in the whole image (x0, y0: position of h[1, 1]). Equal
// heights then differ, so the watershed does not depend on how ties are
// ordered, and a tile gives the same result as the whole image.
// [[Rcpp::export(name = ".breakTiesCpp")]]
NumericMatrix breakTies(NumericMatrix h, int x0, int y0, double eps) {
  NumericMatrix out = clone(h);
  int nr = h.nrow(), nc = h.ncol();
  for (int j = 0; j < nc; j++) {
    for (int i = 0; i < nr; i++) {
      double &v = out[(size_t) j * nr + i];
      if (v > 0) v += eps * pixelHash(x0 + i, y0 + j);
    }
  }
  return out;
}
