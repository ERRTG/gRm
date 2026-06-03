#include <algorithm>
#include <cmath>
#include <cstdint>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

#include <R.h>
#include <R_ext/Rdynload.h>
#include <Rinternals.h>

namespace {

struct ExactError : public std::runtime_error {
  explicit ExactError(const std::string &message) : std::runtime_error(message) {}
};

struct SourceLcg {
  explicit SourceLcg(int seed) : state(static_cast<std::uint32_t>(clamp_seed(seed))), draws(0) {}

  double next() {
    state = static_cast<std::uint32_t>(static_cast<std::uint64_t>(state) * 134775813ULL + 1ULL);
    ++draws;
    return static_cast<double>(state) / 4294967296.0;
  }

  static int clamp_seed(int seed) {
    if (seed < 0) return 0;
    if (seed > 255) return 255;
    return seed;
  }

  std::uint32_t state;
  std::uint64_t draws;
};

struct PreparedSlice {
  int nrow = 0;
  int ncol = 0;
  int total = 0;
  std::vector<int> row_total;
  std::vector<int> col_total;
  std::vector<double> log_factorial;
  std::vector<double> expected;
};

int as_single_int(SEXP value, const char *name) {
  if (Rf_length(value) != 1) {
    throw ExactError(std::string("`") + name + "` must have length 1.");
  }
  if (TYPEOF(value) == INTSXP) {
    int out = INTEGER(value)[0];
    if (out == NA_INTEGER) throw ExactError(std::string("`") + name + "` must not be NA.");
    return out;
  }
  if (TYPEOF(value) == REALSXP) {
    double out = REAL(value)[0];
    if (!R_FINITE(out)) throw ExactError(std::string("`") + name + "` must be finite.");
    return static_cast<int>(out);
  }
  throw ExactError(std::string("`") + name + "` must be numeric.");
}

double as_single_double(SEXP value, const char *name) {
  if (Rf_length(value) != 1) {
    throw ExactError(std::string("`") + name + "` must have length 1.");
  }
  double out;
  if (TYPEOF(value) == INTSXP) {
    if (INTEGER(value)[0] == NA_INTEGER) throw ExactError(std::string("`") + name + "` must not be NA.");
    out = static_cast<double>(INTEGER(value)[0]);
  } else if (TYPEOF(value) == REALSXP) {
    out = REAL(value)[0];
  } else {
    throw ExactError(std::string("`") + name + "` must be numeric.");
  }
  if (!R_FINITE(out)) throw ExactError(std::string("`") + name + "` must be finite.");
  return out;
}

bool as_single_bool(SEXP value, const char *name) {
  if (Rf_length(value) != 1) {
    throw ExactError(std::string("`") + name + "` must have length 1.");
  }
  if (TYPEOF(value) != LGLSXP) {
    throw ExactError(std::string("`") + name + "` must be logical.");
  }
  int out = LOGICAL(value)[0];
  if (out == NA_LOGICAL) throw ExactError(std::string("`") + name + "` must not be NA.");
  return out == TRUE;
}

SEXP list_element(SEXP list, const char *name) {
  SEXP names = Rf_getAttrib(list, R_NamesSymbol);
  if (Rf_length(names) != Rf_length(list)) {
    throw ExactError("Prepared slices must be named lists.");
  }
  for (R_xlen_t i = 0; i < XLENGTH(list); ++i) {
    if (std::strcmp(CHAR(STRING_ELT(names, i)), name) == 0) {
      return VECTOR_ELT(list, i);
    }
  }
  throw ExactError(std::string("Prepared slice is missing `") + name + "`.");
}

SEXP optional_list_element(SEXP list, const char *name) {
  SEXP names = Rf_getAttrib(list, R_NamesSymbol);
  if (Rf_length(names) != Rf_length(list)) {
    throw ExactError("Prepared slices must be named lists.");
  }
  for (R_xlen_t i = 0; i < XLENGTH(list); ++i) {
    if (std::strcmp(CHAR(STRING_ELT(names, i)), name) == 0) {
      return VECTOR_ELT(list, i);
    }
  }
  return R_NilValue;
}

std::vector<int> integer_vector(SEXP value, int expected_length, const char *name) {
  if (Rf_length(value) != expected_length) {
    throw ExactError(std::string("Prepared `") + name + "` has the wrong length.");
  }
  std::vector<int> out(static_cast<std::size_t>(expected_length));
  if (TYPEOF(value) == INTSXP) {
    for (int i = 0; i < expected_length; ++i) {
      int cell = INTEGER(value)[i];
      if (cell == NA_INTEGER) throw ExactError(std::string("Prepared `") + name + "` must not contain NA.");
      out[static_cast<std::size_t>(i)] = cell;
    }
  } else if (TYPEOF(value) == REALSXP) {
    for (int i = 0; i < expected_length; ++i) {
      double cell = REAL(value)[i];
      if (!R_FINITE(cell)) throw ExactError(std::string("Prepared `") + name + "` must be finite.");
      out[static_cast<std::size_t>(i)] = static_cast<int>(cell);
    }
  } else {
    throw ExactError(std::string("Prepared `") + name + "` must be numeric.");
  }
  return out;
}

std::vector<double> double_vector(SEXP value, int expected_length, const char *name) {
  if (TYPEOF(value) != REALSXP || Rf_length(value) != expected_length) {
    throw ExactError(std::string("Prepared `") + name + "` must be a numeric vector with the expected length.");
  }
  std::vector<double> out(static_cast<std::size_t>(expected_length));
  for (int i = 0; i < expected_length; ++i) {
    double cell = REAL(value)[i];
    if (!R_FINITE(cell)) throw ExactError(std::string("Prepared `") + name + "` must be finite.");
    out[static_cast<std::size_t>(i)] = cell;
  }
  return out;
}

int rounded_cell(SEXP matrix, R_xlen_t index) {
  double value;
  if (TYPEOF(matrix) == INTSXP) {
    int cell = INTEGER(matrix)[index];
    if (cell == NA_INTEGER) throw ExactError("Slice tables must not contain NA values.");
    return cell;
  }
  if (TYPEOF(matrix) == REALSXP) {
    value = REAL(matrix)[index];
  } else {
    throw ExactError("Slice tables must be integer or numeric matrices.");
  }
  if (!R_FINITE(value)) throw ExactError("Slice tables must contain only finite values.");
  return static_cast<int>(std::nearbyint(value));
}

inline int index2(int row, int col, int nrow) {
  return row + nrow * col;
}

void fill_expected(PreparedSlice &prepared);

PreparedSlice prepare_slice(SEXP matrix) {
  SEXP dim = Rf_getAttrib(matrix, R_DimSymbol);
  if (Rf_length(dim) != 2) throw ExactError("Each slice must be a two-way matrix.");
  PreparedSlice prepared;
  prepared.nrow = INTEGER(dim)[0];
  prepared.ncol = INTEGER(dim)[1];
  if (prepared.nrow < 1 || prepared.ncol < 1) {
    throw ExactError("Slice matrices must have positive dimensions.");
  }
  prepared.row_total.assign(prepared.nrow, 0);
  prepared.col_total.assign(prepared.ncol, 0);

  for (int col = 0; col < prepared.ncol; ++col) {
    for (int row = 0; row < prepared.nrow; ++row) {
      int cell = rounded_cell(matrix, static_cast<R_xlen_t>(row) + static_cast<R_xlen_t>(prepared.nrow) * col);
      if (cell < 0) throw ExactError("Slice tables must not contain negative counts.");
      prepared.row_total[row] += cell;
      prepared.col_total[col] += cell;
      prepared.total += cell;
    }
  }
  prepared.log_factorial.resize(static_cast<std::size_t>(prepared.total) + 1U);
  if (prepared.total >= 0) prepared.log_factorial[0] = 0.0L;
  for (int i = 1; i <= prepared.total; ++i) {
    if (i < 1001) {
      prepared.log_factorial[static_cast<std::size_t>(i)] =
        prepared.log_factorial[static_cast<std::size_t>(i - 1)] + std::log(static_cast<double>(i));
    } else {
      double n = static_cast<double>(i);
      double ln_sqr_2pi = std::log((2.0 * 3.141592653589793238462643383279502884) *
                                   (2.0 * 3.141592653589793238462643383279502884));
      prepared.log_factorial[static_cast<std::size_t>(i)] =
        ln_sqr_2pi + (n + 0.5) * std::log(n) - n + 1.0 / (24.0 * n) - 2.756;
    }
  }
  fill_expected(prepared);
  return prepared;
}

void fill_log_factorial(PreparedSlice &prepared) {
  prepared.log_factorial.resize(static_cast<std::size_t>(prepared.total) + 1U);
  if (prepared.total >= 0) prepared.log_factorial[0] = 0.0L;
  for (int i = 1; i <= prepared.total; ++i) {
    if (i < 1001) {
      prepared.log_factorial[static_cast<std::size_t>(i)] =
        prepared.log_factorial[static_cast<std::size_t>(i - 1)] + std::log(static_cast<double>(i));
    } else {
      double n = static_cast<double>(i);
      double ln_sqr_2pi = std::log((2.0 * 3.141592653589793238462643383279502884) *
                                   (2.0 * 3.141592653589793238462643383279502884));
      prepared.log_factorial[static_cast<std::size_t>(i)] =
        ln_sqr_2pi + (n + 0.5) * std::log(n) - n + 1.0 / (24.0 * n) - 2.756;
    }
  }
}

void fill_expected(PreparedSlice &prepared) {
  prepared.expected.assign(static_cast<std::size_t>(prepared.nrow) * prepared.ncol, 0.0);
  if (prepared.total <= 0) return;
  for (int col = 0; col < prepared.ncol; ++col) {
    double col_share = static_cast<double>(prepared.col_total[col]) / prepared.total;
    for (int row = 0; row < prepared.nrow; ++row) {
      // Source trace: SKbigtab.Transfer_BT_to_XYZ_TABLE fills RTAB2 once as
      // row margin * (column margin / total). SKrandom.GENTAB1 then passes
      // that stored expected table to SkStat.RCCHI for every generated table.
      prepared.expected[index2(row, col, prepared.nrow)] =
        static_cast<double>(prepared.row_total[row]) * col_share;
    }
  }
}

PreparedSlice prepare_slice_counts(const std::vector<int> &tab, int nrow, int ncol) {
  PreparedSlice prepared;
  prepared.nrow = nrow;
  prepared.ncol = ncol;
  prepared.row_total.assign(nrow, 0);
  prepared.col_total.assign(ncol, 0);
  for (int col = 0; col < ncol; ++col) {
    for (int row = 0; row < nrow; ++row) {
      int cell = tab[index2(row, col, nrow)];
      prepared.row_total[row] += cell;
      prepared.col_total[col] += cell;
      prepared.total += cell;
    }
  }
  fill_log_factorial(prepared);
  fill_expected(prepared);
  return prepared;
}

std::vector<PreparedSlice> prepare_slices(SEXP slices) {
  if (!Rf_isNewList(slices)) throw ExactError("`slices` must be a list of matrices.");
  std::vector<PreparedSlice> prepared;
  R_xlen_t n = XLENGTH(slices);
  prepared.reserve(static_cast<std::size_t>(n));
  for (R_xlen_t i = 0; i < n; ++i) {
    PreparedSlice slice = prepare_slice(VECTOR_ELT(slices, i));
    if (slice.total > 0) prepared.push_back(std::move(slice));
  }
  return prepared;
}

PreparedSlice prepared_slice_from_list(SEXP value) {
  if (!Rf_isNewList(value)) {
    throw ExactError("Each prepared slice must be a list.");
  }
  PreparedSlice prepared;
  prepared.nrow = as_single_int(list_element(value, "cdim"), "cdim");
  prepared.ncol = as_single_int(list_element(value, "rdim"), "rdim");
  prepared.total = as_single_int(list_element(value, "grand_total"), "grand_total");
  if (prepared.nrow < 1 || prepared.ncol < 1 || prepared.total < 0) {
    throw ExactError("Prepared slice dimensions and total must be valid.");
  }
  prepared.row_total = integer_vector(list_element(value, "row_total"), prepared.nrow, "row_total");
  prepared.col_total = integer_vector(list_element(value, "col_total"), prepared.ncol, "col_total");
  {
    std::vector<double> log_factorial = double_vector(list_element(value, "log_factorial"), prepared.total + 1, "log_factorial");
    prepared.log_factorial.assign(log_factorial.begin(), log_factorial.end());
  }
  SEXP expected = optional_list_element(value, "expected");
  if (expected != R_NilValue) {
    std::vector<double> expected_values = double_vector(expected, prepared.nrow * prepared.ncol, "expected");
    prepared.expected.assign(expected_values.begin(), expected_values.end());
  } else {
    fill_expected(prepared);
  }
  return prepared;
}

std::vector<PreparedSlice> prepared_slices_from_list(SEXP slices) {
  if (!Rf_isNewList(slices)) throw ExactError("`prepared_slices` must be a list.");
  std::vector<PreparedSlice> prepared;
  R_xlen_t n = XLENGTH(slices);
  prepared.reserve(static_cast<std::size_t>(n));
  for (R_xlen_t i = 0; i < n; ++i) {
    PreparedSlice slice = prepared_slice_from_list(VECTOR_ELT(slices, i));
    if (slice.total > 0) prepared.push_back(std::move(slice));
  }
  return prepared;
}

int source_round(double value) {
  return static_cast<int>(std::nearbyint(value));
}

double cell_probability(const PreparedSlice &prepared,
                        int free_n,
                        int column1,
                        int column2,
                        int row1,
                        int row2,
                        int t11) {
  int t12 = column1 - t11;
  int t21 = row1 - t11;
  int t22 = free_n - row1 - t12;
  double logp = prepared.log_factorial[static_cast<std::size_t>(column1)] +
    prepared.log_factorial[static_cast<std::size_t>(column2)] +
    prepared.log_factorial[static_cast<std::size_t>(row1)] +
    prepared.log_factorial[static_cast<std::size_t>(row2)] -
    prepared.log_factorial[static_cast<std::size_t>(t11)] -
    prepared.log_factorial[static_cast<std::size_t>(t12)] -
    prepared.log_factorial[static_cast<std::size_t>(t21)] -
    prepared.log_factorial[static_cast<std::size_t>(t22)] -
    prepared.log_factorial[static_cast<std::size_t>(free_n)];
  return std::exp(logp);
}

std::vector<int> gentab1(const PreparedSlice &prepared, SourceLcg &rng) {
  std::vector<int> generated(static_cast<std::size_t>(prepared.nrow) * prepared.ncol, 0);
  std::vector<int> generated_rows(prepared.ncol, 0);
  int generated_total = 0;

  if (prepared.total == 0) return generated;

  for (int ic = 0; ic < prepared.nrow - 1; ++ic) {
    if (prepared.row_total[ic] == 0) continue;

    int generated_column = 0;
    int free_n = prepared.total - generated_total;

    for (int ir = 0; ir < prepared.ncol - 1; ++ir) {
      int t11 = 0;
      if ((prepared.col_total[ir] - generated_rows[ir]) != 0 && free_n != 0) {
        int column1 = prepared.row_total[ic] - generated_column;
        int row1 = prepared.col_total[ir] - generated_rows[ir];
        int column2 = free_n - column1;
        int row2 = free_n - row1;
        int tmax = std::min(column1, row1);
        int tmin = std::max(column1 + row1 - free_n, 0);
        int expected = source_round(static_cast<double>(column1) * (static_cast<double>(row1) / free_n));
        double draw = rng.next();

        int step_down = expected - tmin;
        int step_up = tmax - expected;
        int step_min = std::min(step_down, step_up);
        double cumulative = 0.0;
        bool found = false;

        t11 = expected;
        cumulative += cell_probability(prepared, free_n, column1, column2, row1, row2, t11);
        if (cumulative >= draw) {
          found = true;
        }

        for (int step = 1; !found && step <= step_min; ++step) {
          t11 = expected + step;
          cumulative += cell_probability(prepared, free_n, column1, column2, row1, row2, t11);
          if (cumulative >= draw) {
            found = true;
            break;
          }

          t11 = expected - step;
          cumulative += cell_probability(prepared, free_n, column1, column2, row1, row2, t11);
          if (cumulative >= draw) {
            found = true;
            break;
          }
        }

        if (!found) {
          if (step_down > step_up) {
            for (int step = step_min + 1; step <= step_down; ++step) {
              t11 = expected - step;
              cumulative += cell_probability(prepared, free_n, column1, column2, row1, row2, t11);
              if (cumulative >= draw) break;
            }
          } else if (step_up >= step_min + 1) {
            for (int step = step_min + 1; step <= step_up; ++step) {
              t11 = expected + step;
              cumulative += cell_probability(prepared, free_n, column1, column2, row1, row2, t11);
              if (cumulative >= draw) break;
            }
          }
        }
      }

      generated[index2(ic, ir, prepared.nrow)] = t11;
      generated_total += t11;
      generated_column += t11;
      free_n = free_n - prepared.col_total[ir] + generated_rows[ir];
      generated_rows[ir] += t11;
    }

    int t11 = prepared.row_total[ic] - generated_column;
    generated[index2(ic, prepared.ncol - 1, prepared.nrow)] = t11;
    generated_total += t11;
    generated_rows[prepared.ncol - 1] += t11;
  }

  for (int ir = 0; ir < prepared.ncol; ++ir) {
    generated[index2(prepared.nrow - 1, ir, prepared.nrow)] = prepared.col_total[ir] - generated_rows[ir];
  }

  return generated;
}

double chi_square(const PreparedSlice &prepared, const std::vector<int> &tab) {
  if (prepared.total <= 0) return 0.0;
  double chi = 0.0;
  for (int col = 0; col < prepared.ncol; ++col) {
    for (int row = 0; row < prepared.nrow; ++row) {
      double expected = prepared.expected[index2(row, col, prepared.nrow)];
      if (expected <= 0.0) continue;
      double residual = static_cast<double>(tab[index2(row, col, prepared.nrow)]) - expected;
      chi += residual * residual / expected;
    }
  }
  return chi;
}

struct GammaCounts {
  double ppq = 0.0;
  double pmq = 0.0;
};

GammaCounts gamma_counts(const PreparedSlice &prepared, const std::vector<int> &tab) {
  std::vector<double> cumulative(static_cast<std::size_t>(prepared.nrow) * prepared.ncol, 0.0);
  for (int col = 0; col < prepared.ncol; ++col) {
    for (int row = 0; row < prepared.nrow; ++row) {
      double value = static_cast<double>(tab[index2(row, col, prepared.nrow)]);
      double up = row > 0 ? cumulative[index2(row - 1, col, prepared.nrow)] : 0.0;
      double left = col > 0 ? cumulative[index2(row, col - 1, prepared.nrow)] : 0.0;
      double diag = (row > 0 && col > 0) ? cumulative[index2(row - 1, col - 1, prepared.nrow)] : 0.0;
      cumulative[index2(row, col, prepared.nrow)] = value + up + left - diag;
    }
  }

  auto cell_sum = [&](int row_to, int col_to) -> double {
    if (row_to < 0 || col_to < 0) return 0.0;
    return cumulative[index2(row_to, col_to, prepared.nrow)];
  };

  float p = 0.0F;
  float q = 0.0F;
  for (int col = 0; col < prepared.ncol; ++col) {
    for (int row = 0; row < prepared.nrow; ++row) {
      double n = static_cast<double>(tab[index2(row, col, prepared.nrow)]);
      double less_less = cell_sum(row - 1, col - 1);
      double greater_greater = prepared.total -
        cell_sum(row, prepared.ncol - 1) -
        cell_sum(prepared.nrow - 1, col) +
        cell_sum(row, col);
      double less_greater = cell_sum(row - 1, prepared.ncol - 1) - cell_sum(row - 1, col);
      double greater_less = cell_sum(prepared.nrow - 1, col - 1) - cell_sum(row, col - 1);
      p = static_cast<float>(p + static_cast<float>(n * (less_less + greater_greater)));
      q = static_cast<float>(q + static_cast<float>(n * (less_greater + greater_less)));
    }
  }

  GammaCounts out;
  out.ppq = static_cast<double>(static_cast<float>(p + q));
  out.pmq = static_cast<double>(static_cast<float>(p - q));
  return out;
}

struct ObservedStats {
  double chi_square = 0.0;
  int df = 0;
  double ppq = 0.0;
  double pmq = 0.0;
  double s = 0.0;
};

ObservedStats observed_stats(const PreparedSlice &prepared, const std::vector<int> &tab) {
  ObservedStats out;
  if (prepared.total <= 0) return out;

  int positive_rows = 0;
  int positive_cols = 0;
  for (int row = 0; row < prepared.nrow; ++row) {
    if (prepared.row_total[row] > 0) ++positive_rows;
  }
  for (int col = 0; col < prepared.ncol; ++col) {
    if (prepared.col_total[col] > 0) ++positive_cols;
  }
  out.df = std::max(positive_rows - 1, 0) * std::max(positive_cols - 1, 0);

  for (int col = 0; col < prepared.ncol; ++col) {
    for (int row = 0; row < prepared.nrow; ++row) {
      double expected = prepared.expected[index2(row, col, prepared.nrow)];
      if (expected <= 0.0) continue;
      double residual = static_cast<double>(tab[index2(row, col, prepared.nrow)]) - expected;
      out.chi_square += residual * residual / expected;
    }
  }

  std::vector<double> cumulative(static_cast<std::size_t>(prepared.nrow) * prepared.ncol, 0.0);
  for (int col = 0; col < prepared.ncol; ++col) {
    for (int row = 0; row < prepared.nrow; ++row) {
      double value = static_cast<double>(tab[index2(row, col, prepared.nrow)]);
      double up = row > 0 ? cumulative[index2(row - 1, col, prepared.nrow)] : 0.0;
      double left = col > 0 ? cumulative[index2(row, col - 1, prepared.nrow)] : 0.0;
      double diag = (row > 0 && col > 0) ? cumulative[index2(row - 1, col - 1, prepared.nrow)] : 0.0;
      cumulative[index2(row, col, prepared.nrow)] = value + up + left - diag;
    }
  }
  auto cell_sum = [&](int row_to, int col_to) -> double {
    if (row_to < 0 || col_to < 0) return 0.0;
    return cumulative[index2(row_to, col_to, prepared.nrow)];
  };

  std::vector<double> m(static_cast<std::size_t>(prepared.nrow) * prepared.ncol, 0.0);
  double p = 0.0;
  double q = 0.0;
  for (int col = 0; col < prepared.ncol; ++col) {
    for (int row = 0; row < prepared.nrow; ++row) {
      double n = static_cast<double>(tab[index2(row, col, prepared.nrow)]);
      double less_less = cell_sum(row - 1, col - 1);
      double greater_greater = prepared.total -
        cell_sum(row, prepared.ncol - 1) -
        cell_sum(prepared.nrow - 1, col) +
        cell_sum(row, col);
      double less_greater = cell_sum(row - 1, prepared.ncol - 1) - cell_sum(row - 1, col);
      double greater_less = cell_sum(prepared.nrow - 1, col - 1) - cell_sum(row, col - 1);
      double aij = less_less + greater_greater;
      double dij = less_greater + greater_less;
      m[index2(row, col, prepared.nrow)] = aij - dij;
      p += n * aij;
      q += n * dij;
    }
  }
  out.ppq = p + q;
  out.pmq = p - q;
  if (out.ppq > 0.0) {
    double s = -out.pmq * (out.pmq / prepared.total);
    for (int col = 0; col < prepared.ncol; ++col) {
      for (int row = 0; row < prepared.nrow; ++row) {
        double value = static_cast<double>(tab[index2(row, col, prepared.nrow)]);
        double mm = m[index2(row, col, prepared.nrow)];
        s += value * mm * mm;
      }
    }
    out.s = 4.0 * s;
  }
  return out;
}

bool informative_slice(const PreparedSlice &prepared) {
  int positive_rows = 0;
  int positive_cols = 0;
  for (int row_total : prepared.row_total) {
    if (row_total > 0) ++positive_rows;
  }
  for (int col_total : prepared.col_total) {
    if (col_total > 0) ++positive_cols;
  }
  return positive_rows >= 2 && positive_cols >= 2;
}

struct SimulationResult {
  int chi_exceed = 0;
  int gamma_exceed = 0;
  int gamma_directional_exceed = 0;
  int nsim = 0;
  std::uint64_t draw_count = 0;
  std::uint32_t final_seed = 0;
};

double seq_t(int exceed, int sim, double p0) {
  if (sim < 21) return 0.0;
  double root = std::sqrt(static_cast<double>(sim));
  return static_cast<double>(exceed) / root - root * p0;
}

SimulationResult simulate_chi_gamma(const std::vector<PreparedSlice> &slices,
                                    double observed_chi,
                                    double observed_gamma,
                                    int nsim,
                                    int seed,
                                    bool do_chi,
                                    bool do_gamma,
                                    bool sequential,
                                    int seq_limit) {
  if (nsim < 1) throw ExactError("`nsim` must be a positive integer.");
  SourceLcg rng(seed);
  SimulationResult result;
  result.nsim = nsim;
  bool chi_status = !do_chi;
  bool gamma_status = !do_gamma;
  const double seq_p0 = 0.05;
  const double seq_boundary = 1.058;
  if (seq_limit < 1) seq_limit = nsim;
  const double observed_chi_source = observed_chi;
  const double observed_gamma_source = observed_gamma;

  for (int sim = 0; sim < nsim; ++sim) {
    double chi_total = 0.0;
    double ppq_total = 0.0;
    double pmq_total = 0.0;
    for (const PreparedSlice &slice : slices) {
      std::vector<int> generated = gentab1(slice, rng);
      if (do_chi) {
        chi_total += chi_square(slice, generated);
      }
      if (do_gamma) {
        GammaCounts counts = gamma_counts(slice, generated);
        ppq_total += counts.ppq;
        pmq_total += counts.pmq;
      }
    }
    if (do_chi && chi_total >= observed_chi_source) {
      ++result.chi_exceed;
    }
    if (do_gamma) {
      double simulated_gamma = ppq_total > 0.0 ? pmq_total / ppq_total : 0.0;
      if ((observed_gamma_source > 0.0 && simulated_gamma >= observed_gamma_source) ||
          (observed_gamma_source < 0.0 && simulated_gamma <= observed_gamma_source) ||
          observed_gamma_source == 0.0) {
        ++result.gamma_directional_exceed;
      }
      if (std::fabs(simulated_gamma) >= std::fabs(observed_gamma_source)) {
        ++result.gamma_exceed;
      }
    }
    int completed = sim + 1;
    if (sequential && do_chi) {
      if (seq_t(result.chi_exceed, completed, seq_p0) >= seq_boundary ||
          result.chi_exceed >= seq_limit) {
        chi_status = true;
      }
    }
    if (sequential && do_gamma) {
      if (seq_t(result.gamma_exceed, completed, seq_p0) >= seq_boundary ||
          result.gamma_exceed >= seq_limit) {
        gamma_status = true;
      }
    }
    if (sequential && chi_status && gamma_status) {
      result.nsim = completed;
      break;
    }
  }

  result.draw_count = rng.draws;
  result.final_seed = rng.state;
  return result;
}

SEXP scalar_with_attrs(double value, int exceed, const SimulationResult &result) {
  SEXP out = PROTECT(Rf_allocVector(REALSXP, 1));
  REAL(out)[0] = value;
  Rf_setAttrib(out, Rf_install("exceed"), Rf_ScalarInteger(exceed));
  Rf_setAttrib(out, Rf_install("nsim"), Rf_ScalarInteger(result.nsim));
  Rf_setAttrib(out, Rf_install("draw_count"), Rf_ScalarReal(static_cast<double>(result.draw_count)));
  Rf_setAttrib(out, Rf_install("rng_draws"), Rf_ScalarReal(static_cast<double>(result.draw_count)));
  Rf_setAttrib(out, Rf_install("final_seed"), Rf_ScalarReal(static_cast<double>(result.final_seed)));
  UNPROTECT(1);
  return out;
}

double source_p_value(int exceed, int nsim) {
  return static_cast<double>(static_cast<float>(exceed) / static_cast<float>(nsim));
}

SEXP result_list(const SimulationResult &result) {
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 10));
  SEXP names = PROTECT(Rf_allocVector(STRSXP, 10));
  SET_STRING_ELT(names, 0, Rf_mkChar("p_chi"));
  SET_STRING_ELT(names, 1, Rf_mkChar("p_gamma"));
  SET_STRING_ELT(names, 2, Rf_mkChar("chi_exceed"));
  SET_STRING_ELT(names, 3, Rf_mkChar("gamma_exceed"));
  SET_STRING_ELT(names, 4, Rf_mkChar("nsim"));
  SET_STRING_ELT(names, 5, Rf_mkChar("draw_count"));
  SET_STRING_ELT(names, 6, Rf_mkChar("rng_draws"));
  SET_STRING_ELT(names, 7, Rf_mkChar("final_seed"));
  SET_STRING_ELT(names, 8, Rf_mkChar("p_gamma_directional"));
  SET_STRING_ELT(names, 9, Rf_mkChar("gamma_directional_exceed"));
  SET_VECTOR_ELT(out, 0, Rf_ScalarReal(source_p_value(result.chi_exceed, result.nsim)));
  SET_VECTOR_ELT(out, 1, Rf_ScalarReal(source_p_value(result.gamma_exceed, result.nsim)));
  SET_VECTOR_ELT(out, 2, Rf_ScalarInteger(result.chi_exceed));
  SET_VECTOR_ELT(out, 3, Rf_ScalarInteger(result.gamma_exceed));
  SET_VECTOR_ELT(out, 4, Rf_ScalarInteger(result.nsim));
  SET_VECTOR_ELT(out, 5, Rf_ScalarReal(static_cast<double>(result.draw_count)));
  SET_VECTOR_ELT(out, 6, Rf_ScalarReal(static_cast<double>(result.draw_count)));
  SET_VECTOR_ELT(out, 7, Rf_ScalarReal(static_cast<double>(result.final_seed)));
  SET_VECTOR_ELT(out, 8, Rf_ScalarReal(source_p_value(result.gamma_directional_exceed, result.nsim)));
  SET_VECTOR_ELT(out, 9, Rf_ScalarInteger(result.gamma_directional_exceed));
  Rf_setAttrib(out, R_NamesSymbol, names);
  UNPROTECT(2);
  return out;
}

}  // namespace

extern "C" SEXP gRm_screen_j_exact_kernel(SEXP prepared_slices,
                                              SEXP observed_chi,
                                              SEXP observed_gamma,
                                              SEXP nsim,
                                              SEXP seed,
                                              SEXP compute_chi,
                                              SEXP compute_gamma) {
  try {
    bool do_chi = as_single_bool(compute_chi, "compute_chi");
    bool do_gamma = as_single_bool(compute_gamma, "compute_gamma");
    if (!do_chi && !do_gamma) {
      throw ExactError("At least one of `compute_chi` or `compute_gamma` must be TRUE.");
    }
    std::vector<PreparedSlice> prepared = prepared_slices_from_list(prepared_slices);
    SimulationResult result = simulate_chi_gamma(
      prepared,
      as_single_double(observed_chi, "observed_chi"),
      as_single_double(observed_gamma, "observed_gamma"),
      as_single_int(nsim, "nsim"),
      as_single_int(seed, "seed"),
      do_chi,
      do_gamma,
      do_chi && do_gamma,
      as_single_int(nsim, "nsim")
    );
    return result_list(result);
  } catch (const std::exception &err) {
    Rf_error("%s", err.what());
  }
}

extern "C" SEXP gRm_screen_j_exact_chi_gamma_slices(SEXP slices,
                                                        SEXP observed_chi,
                                                        SEXP observed_gamma,
                                                        SEXP nsim,
                                                        SEXP seed,
                                                        SEXP sequential,
                                                        SEXP seq_limit) {
  try {
    std::vector<PreparedSlice> prepared = prepare_slices(slices);
    SimulationResult result = simulate_chi_gamma(
      prepared,
      as_single_double(observed_chi, "observed_chi"),
      as_single_double(observed_gamma, "observed_gamma"),
      as_single_int(nsim, "nsim"),
      as_single_int(seed, "seed"),
      true,
      true,
      as_single_bool(sequential, "sequential"),
      as_single_int(seq_limit, "seq_limit")
    );
    return result_list(result);
  } catch (const std::exception &err) {
    Rf_error("%s", err.what());
  }
}

extern "C" SEXP gRm_screen_j_exact_chi_slices(SEXP slices,
                                                  SEXP observed_chi,
                                                  SEXP nsim,
                                                  SEXP seed) {
  try {
    std::vector<PreparedSlice> prepared = prepare_slices(slices);
    SimulationResult result = simulate_chi_gamma(
      prepared,
      as_single_double(observed_chi, "observed_chi"),
      0.0,
      as_single_int(nsim, "nsim"),
      as_single_int(seed, "seed"),
      true,
      false,
      false,
      as_single_int(nsim, "nsim")
    );
    return scalar_with_attrs(static_cast<double>(result.chi_exceed) / result.nsim, result.chi_exceed, result);
  } catch (const std::exception &err) {
    Rf_error("%s", err.what());
  }
}

extern "C" SEXP gRm_screen_j_exact_gamma_slices(SEXP slices,
                                                    SEXP observed_gamma,
                                                    SEXP nsim,
                                                    SEXP seed) {
  try {
    std::vector<PreparedSlice> prepared = prepare_slices(slices);
    SimulationResult result = simulate_chi_gamma(
      prepared,
      0.0,
      as_single_double(observed_gamma, "observed_gamma"),
      as_single_int(nsim, "nsim"),
      as_single_int(seed, "seed"),
      false,
      true,
      false,
      as_single_int(nsim, "nsim")
    );
    return scalar_with_attrs(static_cast<double>(result.gamma_exceed) / result.nsim, result.gamma_exceed, result);
  } catch (const std::exception &err) {
    Rf_error("%s", err.what());
  }
}

extern "C" SEXP gRm_screen_j_conditional_bias_test(SEXP x,
                                                       SEXP y,
                                                       SEXP x_dim,
                                                       SEXP y_dim,
                                                       SEXP condition_values,
                                                       SEXP condition_dims,
                                                       SEXP valid,
                                                       SEXP exact,
                                                       SEXP nsim,
                                                       SEXP seed,
                                                       SEXP sequential,
                                                       SEXP seq_limit) {
  try {
    int n = Rf_length(x);
    if (Rf_length(y) != n || Rf_length(valid) != n) {
      throw ExactError("`x`, `y`, and `valid` must have the same length.");
    }
    if (TYPEOF(valid) != LGLSXP) {
      throw ExactError("`valid` must be logical.");
    }
    int nrow = as_single_int(x_dim, "x_dim");
    int ncol = as_single_int(y_dim, "y_dim");
    if (nrow < 1 || ncol < 1) {
      throw ExactError("`x_dim` and `y_dim` must be positive.");
    }

    SEXP condition_dim_attr = Rf_getAttrib(condition_values, R_DimSymbol);
    if (Rf_length(condition_dim_attr) != 2) {
      throw ExactError("`condition_values` must be a matrix.");
    }
    int condition_n = INTEGER(condition_dim_attr)[0];
    int condition_cols = INTEGER(condition_dim_attr)[1];
    if (condition_n != n) {
      throw ExactError("`condition_values` must have one row per observation.");
    }
    std::vector<int> cdims = integer_vector(condition_dims, condition_cols, "condition_dims");

    auto numeric_at = [](SEXP value, int index, const char *name) -> int {
      if (TYPEOF(value) == INTSXP) {
        int out = INTEGER(value)[index];
        return out == NA_INTEGER ? NA_INTEGER : out;
      }
      if (TYPEOF(value) == REALSXP) {
        double out = REAL(value)[index];
        if (!R_FINITE(out)) return NA_INTEGER;
        return static_cast<int>(out);
      }
      throw ExactError(std::string("`") + name + "` must be numeric.");
    };

    std::map<long long, std::vector<int>> tables;
    for (int obs = 0; obs < n; ++obs) {
      if (LOGICAL(valid)[obs] != TRUE) continue;
      int xv = numeric_at(x, obs, "x");
      int yv = numeric_at(y, obs, "y");
      if (xv < 1 || xv > nrow || yv < 1 || yv > ncol) continue;

      long long key = 1;
      long long multiplier = 1;
      bool keep = true;
      for (int condition_col = 0; condition_col < condition_cols; ++condition_col) {
        int cv = numeric_at(condition_values, obs + n * condition_col, "condition_values");
        int cdim = cdims[condition_col];
        if (cv < 1 || cv > cdim) {
          keep = false;
          break;
        }
        key += static_cast<long long>(cv - 1) * multiplier;
        multiplier *= static_cast<long long>(cdim);
      }
      if (!keep) continue;

      auto inserted = tables.emplace(key, std::vector<int>(static_cast<std::size_t>(nrow) * ncol, 0));
      std::vector<int> &tab = inserted.first->second;
      ++tab[index2(xv - 1, yv - 1, nrow)];
    }

    double chi_total = 0.0;
    int df_total = 0;
    double ppq_total = 0.0;
    double pmq_total = 0.0;
    double s_total = 0.0;
    std::vector<PreparedSlice> exact_slices;
    bool do_exact = as_single_bool(exact, "exact");

    for (const auto &entry : tables) {
      PreparedSlice prepared = prepare_slice_counts(entry.second, nrow, ncol);
      ObservedStats stats = observed_stats(prepared, entry.second);
      chi_total += stats.chi_square;
      df_total += stats.df;
      ppq_total += stats.ppq;
      pmq_total += stats.pmq;
      s_total += stats.s;
      if (!do_exact || informative_slice(prepared)) {
        if (prepared.total > 0) exact_slices.push_back(std::move(prepared));
      }
    }

    double gamma = ppq_total > 0.0 ? pmq_total / ppq_total : 0.0;
    double s = ppq_total > 0.0 ? s_total / ppq_total / ppq_total : 0.0;
    SimulationResult sim_result;
    if (do_exact) {
      sim_result = simulate_chi_gamma(
        exact_slices,
        chi_total,
        gamma,
        as_single_int(nsim, "nsim"),
        as_single_int(seed, "seed"),
        true,
        true,
        as_single_bool(sequential, "sequential"),
        as_single_int(seq_limit, "seq_limit")
      );
    }

    SEXP out = PROTECT(Rf_allocVector(VECSXP, 13));
    SEXP names = PROTECT(Rf_allocVector(STRSXP, 13));
    const char *out_names[] = {
      "chi_square", "df", "ppq", "pmq", "s", "gamma",
      "p_chi_exact", "p_gamma_exact", "exact_nsim",
      "chi_exceed", "gamma_exceed", "draw_count", "final_seed"
    };
    for (int i = 0; i < 13; ++i) SET_STRING_ELT(names, i, Rf_mkChar(out_names[i]));
    SET_VECTOR_ELT(out, 0, Rf_ScalarReal(chi_total));
    SET_VECTOR_ELT(out, 1, Rf_ScalarInteger(df_total));
    SET_VECTOR_ELT(out, 2, Rf_ScalarReal(ppq_total));
    SET_VECTOR_ELT(out, 3, Rf_ScalarReal(pmq_total));
    SET_VECTOR_ELT(out, 4, Rf_ScalarReal(s));
    SET_VECTOR_ELT(out, 5, Rf_ScalarReal(gamma));
    SET_VECTOR_ELT(out, 6, Rf_ScalarReal(do_exact ? source_p_value(sim_result.chi_exceed, sim_result.nsim) : NA_REAL));
    SET_VECTOR_ELT(out, 7, Rf_ScalarReal(do_exact ? source_p_value(sim_result.gamma_exceed, sim_result.nsim) : NA_REAL));
    SET_VECTOR_ELT(out, 8, Rf_ScalarInteger(do_exact ? sim_result.nsim : 0));
    SET_VECTOR_ELT(out, 9, Rf_ScalarInteger(do_exact ? sim_result.chi_exceed : NA_INTEGER));
    SET_VECTOR_ELT(out, 10, Rf_ScalarInteger(do_exact ? sim_result.gamma_exceed : NA_INTEGER));
    SET_VECTOR_ELT(out, 11, Rf_ScalarReal(do_exact ? static_cast<double>(sim_result.draw_count) : NA_REAL));
    SET_VECTOR_ELT(out, 12, Rf_ScalarReal(do_exact ? static_cast<double>(sim_result.final_seed) : NA_REAL));
    Rf_setAttrib(out, R_NamesSymbol, names);
    UNPROTECT(2);
    return out;
  } catch (const std::exception &err) {
    Rf_error("%s", err.what());
  }
}

extern "C" SEXP gRm_screen_j_source_random_draws(SEXP seed, SEXP n) {
  try {
    int n_draws = as_single_int(n, "n");
    if (n_draws < 0) throw ExactError("`n` must not be negative.");
    SourceLcg rng(as_single_int(seed, "seed"));
    SEXP out = PROTECT(Rf_allocVector(REALSXP, n_draws));
    for (int index = 0; index < n_draws; ++index) {
      REAL(out)[index] = rng.next();
    }
    UNPROTECT(1);
    return out;
  } catch (const std::exception &err) {
    Rf_error("%s", err.what());
  }
}

extern "C" SEXP gRm_item_parameters_top_ice_field(SEXP gamma_top, SEXP max_score);
extern "C" SEXP gRm_item_parameters_extended_top_ice_fields(SEXP item_counts,
                                                                SEXP score_counts,
                                                                SEXP raw_max,
                                                                SEXP max_step,
                                                                SEXP max_delta);
extern "C" SEXP gRm_item_parameters_extended_gamma(SEXP item_counts,
                                                       SEXP score_counts,
                                                       SEXP raw_max,
                                                       SEXP max_step,
                                                       SEXP max_delta);
extern "C" SEXP gRm_item_parameters_extended_ice_fields(SEXP item_counts,
                                                            SEXP score_counts,
                                                            SEXP raw_max,
                                                            SEXP max_step,
                                                            SEXP max_delta);
extern "C" SEXP gRm_item_parameters_ice_fields_from_gamma(SEXP item_gamma,
                                                              SEXP raw_max);
extern "C" SEXP gRm_global_homogeneity_expected_summary(SEXP item_gamma,
                                                            SEXP score_counts,
                                                            SEXP score_item_n,
                                                            SEXP raw_max,
                                                            SEXP least_score,
                                                            SEXP largest_score);
extern "C" SEXP gRm_global_homogeneity_expected_summary_from_fit(SEXP full_item_counts,
                                                                     SEXP full_score_counts,
                                                                     SEXP group_score_counts,
                                                                     SEXP score_item_n,
                                                                     SEXP raw_max,
                                                                     SEXP least_score,
                                                                     SEXP largest_score,
                                                                     SEXP max_step,
                                                                     SEXP max_delta);

static const R_CallMethodDef CallEntries[] = {
  {"_gRm_screen_j_exact_kernel", reinterpret_cast<DL_FUNC>(&gRm_screen_j_exact_kernel), 7},
  {"gRm_screen_j_exact_chi_gamma_slices", reinterpret_cast<DL_FUNC>(&gRm_screen_j_exact_chi_gamma_slices), 7},
  {"gRm_screen_j_exact_chi_slices", reinterpret_cast<DL_FUNC>(&gRm_screen_j_exact_chi_slices), 4},
  {"gRm_screen_j_exact_gamma_slices", reinterpret_cast<DL_FUNC>(&gRm_screen_j_exact_gamma_slices), 4},
  {"gRm_screen_j_conditional_bias_test", reinterpret_cast<DL_FUNC>(&gRm_screen_j_conditional_bias_test), 12},
  {"gRm_item_parameters_top_ice_field", reinterpret_cast<DL_FUNC>(&gRm_item_parameters_top_ice_field), 2},
  {"gRm_item_parameters_extended_top_ice_fields", reinterpret_cast<DL_FUNC>(&gRm_item_parameters_extended_top_ice_fields), 5},
  {"gRm_item_parameters_extended_gamma", reinterpret_cast<DL_FUNC>(&gRm_item_parameters_extended_gamma), 5},
  {"gRm_item_parameters_extended_ice_fields", reinterpret_cast<DL_FUNC>(&gRm_item_parameters_extended_ice_fields), 5},
  {"gRm_item_parameters_ice_fields_from_gamma", reinterpret_cast<DL_FUNC>(&gRm_item_parameters_ice_fields_from_gamma), 2},
  {"gRm_global_homogeneity_expected_summary", reinterpret_cast<DL_FUNC>(&gRm_global_homogeneity_expected_summary), 6},
  {"gRm_global_homogeneity_expected_summary_from_fit", reinterpret_cast<DL_FUNC>(&gRm_global_homogeneity_expected_summary_from_fit), 9},
  {"_gRm_screen_j_source_random_draws", reinterpret_cast<DL_FUNC>(&gRm_screen_j_source_random_draws), 2},
  {"gRm_screen_j_source_random_draws", reinterpret_cast<DL_FUNC>(&gRm_screen_j_source_random_draws), 2},
  {nullptr, nullptr, 0}
};

extern "C" void R_init_gRm(DllInfo *dll) {
  R_registerRoutines(dll, nullptr, CallEntries, nullptr, nullptr);
  R_useDynamicSymbols(dll, FALSE);
}
