#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstddef>
#include <cstring>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

#include <R.h>
#include <Rinternals.h>

namespace {

struct ExactError : public std::runtime_error {
  explicit ExactError(const std::string &message) : std::runtime_error(message) {}
};

struct SourceLcg {
  // Source trace: source/digram_source_20260817/skunits/SKrandom.pas::GENTAB1 consumes Pascal
  // Random values while traversing free cells. The Delphi-LCG recurrence is
  // the audited runtime compatibility stream pinned by
  // pascal_harness/SCREEN_EXACT_TRAJECTORY.pas::DelphiRandom.
  explicit SourceLcg(int seed) : state(static_cast<std::uint32_t>(clamp_seed(seed))), draws(0) {}

  double next() {
    state = static_cast<std::uint32_t>(static_cast<std::uint64_t>(state) * 134775813ULL + 1ULL);
    ++draws;
    const double value = static_cast<double>(state) / 4294967296.0;
    if (record_draws) draw_values.push_back(value);
    return value;
  }

  static int clamp_seed(int seed) {
    if (seed < 0) return 0;
    if (seed > 255) return 255;
    return seed;
  }

  std::uint32_t state;
  std::uint64_t draws;
  bool record_draws = false;
  std::vector<double> draw_values;
};

struct PreparedSlice {
  int nrow = 0;
  int ncol = 0;
  int ncell = 0;
  int total = 0;
  int positive_rows = 0;
  int positive_cols = 0;
  std::vector<int> row_total;
  std::vector<int> col_total;
  std::vector<int> positive_expected_index;
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

int numeric_at(SEXP value, R_xlen_t index, const char *name) {
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
}

inline int index2(int row, int col, int nrow) {
  return row + nrow * col;
}

void fill_expected(PreparedSlice &prepared);
void fill_slice_metadata(PreparedSlice &prepared);
void fill_log_factorial(PreparedSlice &prepared);

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
  fill_log_factorial(prepared);
  fill_expected(prepared);
  fill_slice_metadata(prepared);
  return prepared;
}

void fill_log_factorial(PreparedSlice &prepared) {
  // Source trace: source/digram_source_20260817/skunits/SKrandom.pas::GENTAB1 builds the
  // log-factorial table recursively through 1000 and then uses its preserved
  // large-n approximation. Both R-facing prepared-slice routes share this one
  // implementation so their probability trajectory cannot drift.
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
      // Source trace:
      // source/digram_source_20260817/skunits/SKbigtab.pas::Transfer_BT_to_XYZ_TABLE fills RTAB2
      // once as row margin * (column margin / total).
      // source/digram_source_20260817/skunits/SKrandom.pas::GENTAB1 then passes that stored
      // expected table to source/digram_source_20260817/skunits/SkStat.pas::RCCHI.
      prepared.expected[index2(row, col, prepared.nrow)] =
        static_cast<double>(prepared.row_total[row]) * col_share;
    }
  }
}

void fill_slice_metadata(PreparedSlice &prepared) {
  prepared.ncell = prepared.nrow * prepared.ncol;
  prepared.positive_rows = 0;
  prepared.positive_cols = 0;
  for (int row_total : prepared.row_total) {
    if (row_total > 0) ++prepared.positive_rows;
  }
  for (int col_total : prepared.col_total) {
    if (col_total > 0) ++prepared.positive_cols;
  }
  prepared.positive_expected_index.clear();
  prepared.positive_expected_index.reserve(static_cast<std::size_t>(prepared.ncell));
  // Source trace: source/digram_source_20260817/skunits/SkStat.pas::RCCHI accumulates cells as
  // `FOR I:=1 TO C DO FOR J:=1 TO R DO`, i.e. row/category of the first
  // variable first and then the second variable. Preserve that order because
  // exact tests compare generated chi-square totals with `>=`, so source-tied
  // tables can move by a few ulps under R/C column-major iteration.
  for (int row = 0; row < prepared.nrow; ++row) {
    for (int col = 0; col < prepared.ncol; ++col) {
      int index = index2(row, col, prepared.nrow);
      if (prepared.expected[static_cast<std::size_t>(index)] > 0.0) {
        prepared.positive_expected_index.push_back(index);
      }
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
  fill_slice_metadata(prepared);
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
  // Source trace: source/digram_source_20260817/skunits/SKrandom.pas::GENTAB1 evaluates each
  // feasible free-cell count from the same log-factorial hypergeometric term.
  // Mathematical step: exponentiate the four margin factorials minus the four
  // cell factorials and the free-total factorial in the preserved order.
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

struct SliceScratch {
  std::vector<int> generated;
  std::vector<int> generated_rows;
  std::vector<double> cumulative;

  void reset_for(const PreparedSlice &slice) {
    const std::size_t ncell = static_cast<std::size_t>(slice.ncell);
    if (generated.size() != ncell) generated.resize(ncell);
    std::fill(generated.begin(), generated.end(), 0);
    if (generated_rows.size() != static_cast<std::size_t>(slice.ncol)) {
      generated_rows.resize(static_cast<std::size_t>(slice.ncol));
    }
    std::fill(generated_rows.begin(), generated_rows.end(), 0);
    if (cumulative.size() != ncell) cumulative.resize(ncell);
    std::fill(cumulative.begin(), cumulative.end(), 0.0);
  }
};

void gentab1_into(const PreparedSlice &prepared, SourceLcg &rng, SliceScratch &scratch) {
  // Source trace: source/digram_source_20260817/skunits/SKrandom.pas::GENTAB1 visits every free
  // cell in source order, starts at the rounded expected count, alternates up
  // and down, then fills the final row/column deterministically.
  scratch.reset_for(prepared);
  std::vector<int> &generated = scratch.generated;
  std::vector<int> &generated_rows = scratch.generated_rows;
  int generated_total = 0;

  if (prepared.total == 0) return;

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
}

struct GeneratedStats {
  double chi_square = 0.0;
  double ppq = 0.0;
  double pmq = 0.0;
};

GeneratedStats generated_stats(const PreparedSlice &prepared,
                               const std::vector<int> &tab,
                               SliceScratch &scratch,
                               bool do_chi,
                               bool do_gamma) {
  GeneratedStats out;

  if (do_chi && prepared.total > 0) {
    for (int index : prepared.positive_expected_index) {
      double expected = prepared.expected[static_cast<std::size_t>(index)];
      double residual = static_cast<double>(tab[static_cast<std::size_t>(index)]) - expected;
      out.chi_square += residual * (residual / expected);
    }
  }

  if (!do_gamma) return out;

  std::vector<double> &cumulative = scratch.cumulative;
  const std::size_t ncell = static_cast<std::size_t>(prepared.ncell);
  std::fill(cumulative.begin(), cumulative.begin() + static_cast<std::ptrdiff_t>(ncell), 0.0);

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
      p += n * (less_less + greater_greater);
      q += n * (less_greater + greater_less);
    }
  }

  out.ppq = p + q;
  out.pmq = p - q;
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

  out.df = std::max(prepared.positive_rows - 1, 0) * std::max(prepared.positive_cols - 1, 0);

  for (int index : prepared.positive_expected_index) {
    double expected = prepared.expected[static_cast<std::size_t>(index)];
    double residual = static_cast<double>(tab[static_cast<std::size_t>(index)]) - expected;
    out.chi_square += residual * (residual / expected);
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
  return prepared.positive_rows >= 2 && prepared.positive_cols >= 2;
}

struct ExactSliceTrace {
  int simulation = 0;
  int slice = 0;
  int nrow = 0;
  int ncol = 0;
  std::vector<int> table;
  std::vector<double> random_draws;
  double chi_square = 0.0;
  double ppq = 0.0;
  double pmq = 0.0;
  std::uint64_t draw_count = 0;
  std::uint32_t final_seed = 0;
};

struct ExactSimulationTrace {
  int simulation = 0;
  double chi_square = 0.0;
  double ppq = 0.0;
  double pmq = 0.0;
  double gamma = 0.0;
  bool chi_ge_observed = false;
  bool gamma_ge_observed = false;
  int chi_exceed = 0;
  int gamma_exceed = 0;
  bool chi_status = false;
  bool gamma_status = false;
  bool stop = false;
  std::uint64_t draw_count = 0;
  std::uint32_t final_seed = 0;
};

struct SimulationResult {
  int chi_exceed = 0;
  int gamma_exceed = 0;
  int gamma_directional_exceed = 0;
  int nsim = 0;
  std::uint64_t draw_count = 0;
  std::uint32_t final_seed = 0;
  std::vector<ExactSliceTrace> slice_trace;
  std::vector<ExactSimulationTrace> simulation_trace;
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
                                    int seq_limit,
                                    bool record_trace = false) {
  if (nsim < 1) throw ExactError("`nsim` must be a positive integer.");
  SourceLcg rng(seed);
  rng.record_draws = record_trace;
  SimulationResult result;
  result.nsim = nsim;
  bool chi_status = !do_chi;
  bool gamma_status = !do_gamma;
  const double seq_p0 = 0.05;
  const double seq_boundary = 1.058;
  if (seq_limit < 1) seq_limit = nsim;
  // Source trace: source/digram_source_20260817/skunits/SKbias3.pas::XYZ_bias_ANALYSE stores
  // observed CHITOT in RESULTS[1,1], and source/digram_source_20260817/skunits/SKTypes.pas::
  // RESARRAY is explicitly SINGLE. GENTAB1's simulated CHI/PPQ/PMQ and
  // AbsGammaTot remain Pascal REAL (8-byte Double in the historical Delphi
  // target). Only the observed chi threshold is therefore quantized before
  // Evaluate_simulated_biasresults uses >=.
  const float observed_chi_source = static_cast<float>(observed_chi);
  const double observed_gamma_source = observed_gamma;
  std::vector<SliceScratch> scratch;
  scratch.resize(slices.size());
  for (std::size_t i = 0; i < slices.size(); ++i) {
    scratch[i].reset_for(slices[i]);
  }

  for (int sim = 0; sim < nsim; ++sim) {
    double chi_total = 0.0;
    double ppq_total = 0.0;
    double pmq_total = 0.0;
    for (std::size_t slice_index = 0; slice_index < slices.size(); ++slice_index) {
      const PreparedSlice &slice = slices[slice_index];
      SliceScratch &slice_scratch = scratch[slice_index];
      const std::size_t draw_start = rng.draw_values.size();
      gentab1_into(slice, rng, slice_scratch);
      const std::vector<int> &generated = slice_scratch.generated;
      GeneratedStats stats = generated_stats(slice, generated, slice_scratch, do_chi, do_gamma);
      if (do_chi) {
        chi_total += stats.chi_square;
      }
      if (do_gamma) {
        ppq_total += stats.ppq;
        pmq_total += stats.pmq;
      }
      if (record_trace) {
        ExactSliceTrace trace;
        trace.simulation = sim + 1;
        trace.slice = static_cast<int>(slice_index) + 1;
        trace.nrow = slice.nrow;
        trace.ncol = slice.ncol;
        trace.table = generated;
        trace.random_draws.assign(
          rng.draw_values.begin() + static_cast<std::ptrdiff_t>(draw_start),
          rng.draw_values.end()
        );
        trace.chi_square = stats.chi_square;
        trace.ppq = stats.ppq;
        trace.pmq = stats.pmq;
        trace.draw_count = rng.draws;
        trace.final_seed = rng.state;
        result.slice_trace.push_back(std::move(trace));
      }
    }
    const bool chi_ge_observed = do_chi && chi_total >= observed_chi_source;
    if (chi_ge_observed) {
      ++result.chi_exceed;
    }
    double simulated_gamma = 0.0;
    bool gamma_ge_observed = false;
    if (do_gamma) {
      simulated_gamma = ppq_total > 0.0 ? pmq_total / ppq_total : 0.0;
      if ((observed_gamma_source > 0.0 && simulated_gamma >= observed_gamma_source) ||
          (observed_gamma_source < 0.0 && simulated_gamma <= observed_gamma_source) ||
          observed_gamma_source == 0.0) {
        ++result.gamma_directional_exceed;
      }
      gamma_ge_observed = std::fabs(simulated_gamma) >= std::fabs(observed_gamma_source);
      if (gamma_ge_observed) {
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
    const bool stop_now = sequential && chi_status && gamma_status;
    if (record_trace) {
      ExactSimulationTrace trace;
      trace.simulation = completed;
      trace.chi_square = chi_total;
      trace.ppq = ppq_total;
      trace.pmq = pmq_total;
      trace.gamma = simulated_gamma;
      trace.chi_ge_observed = chi_ge_observed;
      trace.gamma_ge_observed = gamma_ge_observed;
      trace.chi_exceed = result.chi_exceed;
      trace.gamma_exceed = result.gamma_exceed;
      trace.chi_status = chi_status;
      trace.gamma_status = gamma_status;
      trace.stop = stop_now;
      trace.draw_count = rng.draws;
      trace.final_seed = rng.state;
      result.simulation_trace.push_back(trace);
    }
    if (stop_now) {
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
  SEXP exceed_attr = PROTECT(Rf_ScalarInteger(exceed));
  SEXP nsim_attr = PROTECT(Rf_ScalarInteger(result.nsim));
  SEXP draw_count_attr = PROTECT(Rf_ScalarReal(static_cast<double>(result.draw_count)));
  SEXP rng_draws_attr = PROTECT(Rf_ScalarReal(static_cast<double>(result.draw_count)));
  SEXP final_seed_attr = PROTECT(Rf_ScalarReal(static_cast<double>(result.final_seed)));
  Rf_setAttrib(out, Rf_install("exceed"), exceed_attr);
  Rf_setAttrib(out, Rf_install("nsim"), nsim_attr);
  Rf_setAttrib(out, Rf_install("draw_count"), draw_count_attr);
  Rf_setAttrib(out, Rf_install("rng_draws"), rng_draws_attr);
  Rf_setAttrib(out, Rf_install("final_seed"), final_seed_attr);
  UNPROTECT(6);
  return out;
}

double source_p_value(int exceed, int nsim) {
  // Source trace: source/digram_source_20260817/skunits/SKxyz1.PAS::XYZ_TEST and
  // source/digram_source_20260817/skunits/SKexa1.pas::EXA_SUMMARY1_2 store exact p-values as REAL
  // RESULTS entries from count / NSIM. Keep this in double precision so the
  // native path is bit-for-bit aligned with the R parity reference.
  return static_cast<double>(exceed) / static_cast<double>(nsim);
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

SEXP trace_tables(const SimulationResult &result) {
  const R_xlen_t n = static_cast<R_xlen_t>(result.slice_trace.size());
  SEXP out = PROTECT(Rf_allocVector(VECSXP, n));
  for (R_xlen_t index = 0; index < n; ++index) {
    const ExactSliceTrace &trace = result.slice_trace[static_cast<std::size_t>(index)];
    SEXP table = PROTECT(Rf_allocMatrix(INTSXP, trace.nrow, trace.ncol));
    for (R_xlen_t cell = 0; cell < XLENGTH(table); ++cell) {
      INTEGER(table)[cell] = trace.table[static_cast<std::size_t>(cell)];
    }
    SET_VECTOR_ELT(out, index, table);
    UNPROTECT(1);
  }
  UNPROTECT(1);
  return out;
}

SEXP trace_random_draws(const SimulationResult &result) {
  const R_xlen_t n = static_cast<R_xlen_t>(result.slice_trace.size());
  SEXP out = PROTECT(Rf_allocVector(VECSXP, n));
  for (R_xlen_t index = 0; index < n; ++index) {
    const std::vector<double> &draws =
      result.slice_trace[static_cast<std::size_t>(index)].random_draws;
    SEXP values = PROTECT(Rf_allocVector(REALSXP, static_cast<R_xlen_t>(draws.size())));
    for (R_xlen_t draw = 0; draw < XLENGTH(values); ++draw) {
      REAL(values)[draw] = draws[static_cast<std::size_t>(draw)];
    }
    SET_VECTOR_ELT(out, index, values);
    UNPROTECT(1);
  }
  UNPROTECT(1);
  return out;
}

SEXP trace_table_state(const SimulationResult &result) {
  const R_xlen_t n = static_cast<R_xlen_t>(result.slice_trace.size());
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 7));
  SEXP names = PROTECT(Rf_allocVector(STRSXP, 7));
  const char *column_names[] = {
    "sim", "slice", "chi_square", "ppq", "pmq", "draw_count", "final_seed"
  };
  for (int column = 0; column < 7; ++column) {
    SET_STRING_ELT(names, column, Rf_mkChar(column_names[column]));
  }
  SEXP simulation = PROTECT(Rf_allocVector(INTSXP, n));
  SEXP slice = PROTECT(Rf_allocVector(INTSXP, n));
  SEXP chi_square_values = PROTECT(Rf_allocVector(REALSXP, n));
  SEXP ppq = PROTECT(Rf_allocVector(REALSXP, n));
  SEXP pmq = PROTECT(Rf_allocVector(REALSXP, n));
  SEXP draw_count = PROTECT(Rf_allocVector(REALSXP, n));
  SEXP final_seed = PROTECT(Rf_allocVector(REALSXP, n));
  for (R_xlen_t index = 0; index < n; ++index) {
    const ExactSliceTrace &trace = result.slice_trace[static_cast<std::size_t>(index)];
    INTEGER(simulation)[index] = trace.simulation;
    INTEGER(slice)[index] = trace.slice;
    REAL(chi_square_values)[index] = trace.chi_square;
    REAL(ppq)[index] = trace.ppq;
    REAL(pmq)[index] = trace.pmq;
    REAL(draw_count)[index] = static_cast<double>(trace.draw_count);
    REAL(final_seed)[index] = static_cast<double>(trace.final_seed);
  }
  SET_VECTOR_ELT(out, 0, simulation);
  SET_VECTOR_ELT(out, 1, slice);
  SET_VECTOR_ELT(out, 2, chi_square_values);
  SET_VECTOR_ELT(out, 3, ppq);
  SET_VECTOR_ELT(out, 4, pmq);
  SET_VECTOR_ELT(out, 5, draw_count);
  SET_VECTOR_ELT(out, 6, final_seed);
  Rf_setAttrib(out, R_NamesSymbol, names);
  UNPROTECT(9);
  return out;
}

SEXP trace_simulation_state(const SimulationResult &result) {
  const R_xlen_t n = static_cast<R_xlen_t>(result.simulation_trace.size());
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 14));
  SEXP names = PROTECT(Rf_allocVector(STRSXP, 14));
  const char *column_names[] = {
    "sim", "chi_square", "ppq", "pmq", "gamma",
    "chi_ge_observed", "gamma_ge_observed", "chi_exceed", "gamma_exceed",
    "chi_status", "gamma_status", "stop", "draw_count", "final_seed"
  };
  for (int column = 0; column < 14; ++column) {
    SET_STRING_ELT(names, column, Rf_mkChar(column_names[column]));
  }
  SEXP simulation = PROTECT(Rf_allocVector(INTSXP, n));
  SEXP chi_square_values = PROTECT(Rf_allocVector(REALSXP, n));
  SEXP ppq = PROTECT(Rf_allocVector(REALSXP, n));
  SEXP pmq = PROTECT(Rf_allocVector(REALSXP, n));
  SEXP gamma = PROTECT(Rf_allocVector(REALSXP, n));
  SEXP chi_ge = PROTECT(Rf_allocVector(LGLSXP, n));
  SEXP gamma_ge = PROTECT(Rf_allocVector(LGLSXP, n));
  SEXP chi_exceed = PROTECT(Rf_allocVector(INTSXP, n));
  SEXP gamma_exceed = PROTECT(Rf_allocVector(INTSXP, n));
  SEXP chi_status = PROTECT(Rf_allocVector(LGLSXP, n));
  SEXP gamma_status = PROTECT(Rf_allocVector(LGLSXP, n));
  SEXP stop = PROTECT(Rf_allocVector(LGLSXP, n));
  SEXP draw_count = PROTECT(Rf_allocVector(REALSXP, n));
  SEXP final_seed = PROTECT(Rf_allocVector(REALSXP, n));
  for (R_xlen_t index = 0; index < n; ++index) {
    const ExactSimulationTrace &trace =
      result.simulation_trace[static_cast<std::size_t>(index)];
    INTEGER(simulation)[index] = trace.simulation;
    REAL(chi_square_values)[index] = trace.chi_square;
    REAL(ppq)[index] = trace.ppq;
    REAL(pmq)[index] = trace.pmq;
    REAL(gamma)[index] = trace.gamma;
    LOGICAL(chi_ge)[index] = trace.chi_ge_observed ? TRUE : FALSE;
    LOGICAL(gamma_ge)[index] = trace.gamma_ge_observed ? TRUE : FALSE;
    INTEGER(chi_exceed)[index] = trace.chi_exceed;
    INTEGER(gamma_exceed)[index] = trace.gamma_exceed;
    LOGICAL(chi_status)[index] = trace.chi_status ? TRUE : FALSE;
    LOGICAL(gamma_status)[index] = trace.gamma_status ? TRUE : FALSE;
    LOGICAL(stop)[index] = trace.stop ? TRUE : FALSE;
    REAL(draw_count)[index] = static_cast<double>(trace.draw_count);
    REAL(final_seed)[index] = static_cast<double>(trace.final_seed);
  }
  SET_VECTOR_ELT(out, 0, simulation);
  SET_VECTOR_ELT(out, 1, chi_square_values);
  SET_VECTOR_ELT(out, 2, ppq);
  SET_VECTOR_ELT(out, 3, pmq);
  SET_VECTOR_ELT(out, 4, gamma);
  SET_VECTOR_ELT(out, 5, chi_ge);
  SET_VECTOR_ELT(out, 6, gamma_ge);
  SET_VECTOR_ELT(out, 7, chi_exceed);
  SET_VECTOR_ELT(out, 8, gamma_exceed);
  SET_VECTOR_ELT(out, 9, chi_status);
  SET_VECTOR_ELT(out, 10, gamma_status);
  SET_VECTOR_ELT(out, 11, stop);
  SET_VECTOR_ELT(out, 12, draw_count);
  SET_VECTOR_ELT(out, 13, final_seed);
  Rf_setAttrib(out, R_NamesSymbol, names);
  UNPROTECT(16);
  return out;
}

SEXP trace_result_list(const SimulationResult &result) {
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 11));
  SEXP names = PROTECT(Rf_allocVector(STRSXP, 11));
  const char *out_names[] = {
    "p_chi", "p_gamma", "chi_exceed", "gamma_exceed", "nsim",
    "draw_count", "rng_draws", "final_seed", "p_gamma_directional",
    "gamma_directional_exceed", "trajectory"
  };
  for (int index = 0; index < 11; ++index) {
    SET_STRING_ELT(names, index, Rf_mkChar(out_names[index]));
  }
  SET_VECTOR_ELT(out, 0, Rf_ScalarReal(source_p_value(result.chi_exceed, result.nsim)));
  SET_VECTOR_ELT(out, 1, Rf_ScalarReal(source_p_value(result.gamma_exceed, result.nsim)));
  SET_VECTOR_ELT(out, 2, Rf_ScalarInteger(result.chi_exceed));
  SET_VECTOR_ELT(out, 3, Rf_ScalarInteger(result.gamma_exceed));
  SET_VECTOR_ELT(out, 4, Rf_ScalarInteger(result.nsim));
  SET_VECTOR_ELT(out, 5, Rf_ScalarReal(static_cast<double>(result.draw_count)));
  SET_VECTOR_ELT(out, 6, Rf_ScalarReal(static_cast<double>(result.draw_count)));
  SET_VECTOR_ELT(out, 7, Rf_ScalarReal(static_cast<double>(result.final_seed)));
  SET_VECTOR_ELT(out, 8, Rf_ScalarReal(source_p_value(
    result.gamma_directional_exceed, result.nsim
  )));
  SET_VECTOR_ELT(out, 9, Rf_ScalarInteger(result.gamma_directional_exceed));

  SEXP trajectory = PROTECT(Rf_allocVector(VECSXP, 4));
  SEXP trajectory_names = PROTECT(Rf_allocVector(STRSXP, 4));
  SET_STRING_ELT(trajectory_names, 0, Rf_mkChar("tables"));
  SET_STRING_ELT(trajectory_names, 1, Rf_mkChar("random_draws"));
  SET_STRING_ELT(trajectory_names, 2, Rf_mkChar("table_state"));
  SET_STRING_ELT(trajectory_names, 3, Rf_mkChar("simulations"));
  SEXP tables = PROTECT(trace_tables(result));
  SEXP random_draws = PROTECT(trace_random_draws(result));
  SEXP table_state = PROTECT(trace_table_state(result));
  SEXP simulations = PROTECT(trace_simulation_state(result));
  SET_VECTOR_ELT(trajectory, 0, tables);
  SET_VECTOR_ELT(trajectory, 1, random_draws);
  SET_VECTOR_ELT(trajectory, 2, table_state);
  SET_VECTOR_ELT(trajectory, 3, simulations);
  Rf_setAttrib(trajectory, R_NamesSymbol, trajectory_names);
  SET_VECTOR_ELT(out, 10, trajectory);
  Rf_setAttrib(out, R_NamesSymbol, names);
  UNPROTECT(8);
  return out;
}

struct ConditionalBiasBuild {
  std::vector<PreparedSlice> slices;
  double observed_chi = 0.0;
  int observed_df = 0;
  double observed_gamma = 0.0;
  double ppq = 0.0;
  double pmq = 0.0;
  double s = 0.0;
};

ConditionalBiasBuild build_conditional_bias_result(SEXP x,
                                                   SEXP y,
                                                   int nrow,
                                                   int ncol,
                                                   SEXP condition_values,
                                                   const std::vector<int> &condition_dims,
                                                   SEXP valid,
                                                   bool exact_requires_informative_slices) {
  R_xlen_t n = XLENGTH(x);
  if (XLENGTH(y) != n || XLENGTH(valid) != n) {
    throw ExactError("`x`, `y`, and `valid` must have the same length.");
  }
  if (TYPEOF(valid) != LGLSXP) {
    throw ExactError("`valid` must be logical.");
  }
  if (nrow < 1 || ncol < 1) {
    throw ExactError("`x_dim` and `y_dim` must be positive.");
  }

  SEXP condition_dim_attr = Rf_getAttrib(condition_values, R_DimSymbol);
  int condition_cols = 1;
  if (Rf_length(condition_dim_attr) == 2) {
    if (INTEGER(condition_dim_attr)[0] != n) {
      throw ExactError("`condition_values` must have one row per observation.");
    }
    condition_cols = INTEGER(condition_dim_attr)[1];
  } else {
    if (XLENGTH(condition_values) != n) {
      throw ExactError("`condition_values` must have one value per observation.");
    }
  }
  if (static_cast<int>(condition_dims.size()) != condition_cols) {
    throw ExactError("`condition_dims` must have one entry per conditioning column.");
  }

  std::map<long long, std::vector<int>> tables;
  for (R_xlen_t obs = 0; obs < n; ++obs) {
    if (LOGICAL(valid)[obs] != TRUE) continue;
    int xv = numeric_at(x, obs, "x");
    int yv = numeric_at(y, obs, "y");
    if (xv < 1 || xv > nrow || yv < 1 || yv > ncol) continue;

    long long key = 1;
    long long multiplier = 1;
    bool keep = true;
    for (int condition_col = 0; condition_col < condition_cols; ++condition_col) {
      R_xlen_t condition_index = Rf_length(condition_dim_attr) == 2 ?
        obs + n * static_cast<R_xlen_t>(condition_col) :
        obs;
      int cv = numeric_at(condition_values, condition_index, "condition_values");
      int cdim = condition_dims[static_cast<std::size_t>(condition_col)];
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

  ConditionalBiasBuild out;
  double s_total = 0.0;
  for (const auto &entry : tables) {
    PreparedSlice prepared = prepare_slice_counts(entry.second, nrow, ncol);
    ObservedStats stats = observed_stats(prepared, entry.second);
    out.observed_chi += stats.chi_square;
    out.observed_df += stats.df;
    out.ppq += stats.ppq;
    out.pmq += stats.pmq;
    s_total += stats.s;
    if (prepared.total > 0 &&
        (!exact_requires_informative_slices || informative_slice(prepared))) {
      out.slices.push_back(std::move(prepared));
    }
  }

  out.observed_gamma = out.ppq > 0.0 ? out.pmq / out.ppq : 0.0;
  out.s = out.ppq > 0.0 ? s_total / out.ppq / out.ppq : 0.0;
  return out;
}

}  // namespace

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

extern "C" SEXP gRm_screen_j_exact_chi_gamma_trace_slices(SEXP slices,
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
      as_single_int(seq_limit, "seq_limit"),
      true
    );
    return trace_result_list(result);
  } catch (const std::exception &err) {
    Rf_error("%s", err.what());
  }
}

extern "C" SEXP gRm_screen_j_exact_chi_slices(SEXP slices,
                                                  SEXP observed_chi,
                                                  SEXP nsim,
                                                  SEXP seed,
                                                  SEXP sequential,
                                                  SEXP seq_limit) {
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
      as_single_bool(sequential, "sequential"),
      as_single_int(seq_limit, "seq_limit")
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
    int nrow = as_single_int(x_dim, "x_dim");
    int ncol = as_single_int(y_dim, "y_dim");

    SEXP condition_dim_attr = Rf_getAttrib(condition_values, R_DimSymbol);
    if (Rf_length(condition_dim_attr) != 2) {
      throw ExactError("`condition_values` must be a matrix.");
    }
    int condition_cols = INTEGER(condition_dim_attr)[1];
    std::vector<int> cdims = integer_vector(condition_dims, condition_cols, "condition_dims");
    ConditionalBiasBuild built = build_conditional_bias_result(
      x, y, nrow, ncol, condition_values, cdims, valid, true
    );
    bool do_exact = as_single_bool(exact, "exact");
    SimulationResult sim_result;
    if (do_exact) {
      sim_result = simulate_chi_gamma(
        built.slices,
        built.observed_chi,
        built.observed_gamma,
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
    SET_VECTOR_ELT(out, 0, Rf_ScalarReal(built.observed_chi));
    SET_VECTOR_ELT(out, 1, Rf_ScalarInteger(built.observed_df));
    SET_VECTOR_ELT(out, 2, Rf_ScalarReal(built.ppq));
    SET_VECTOR_ELT(out, 3, Rf_ScalarReal(built.pmq));
    SET_VECTOR_ELT(out, 4, Rf_ScalarReal(built.s));
    SET_VECTOR_ELT(out, 5, Rf_ScalarReal(built.observed_gamma));
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

extern "C" SEXP gRm_screen_j_item_pair_conditional_exact(SEXP x,
                                                             SEXP y,
                                                             SEXP x_dim,
                                                             SEXP y_dim,
                                                             SEXP condition_values,
                                                             SEXP condition_dim,
                                                             SEXP valid,
                                                             SEXP nsim,
                                                             SEXP seed,
                                                             SEXP sequential,
                                                             SEXP seq_limit) {
  try {
    std::vector<int> condition_dims(1, as_single_int(condition_dim, "condition_dim"));
    ConditionalBiasBuild built = build_conditional_bias_result(
      x,
      y,
      as_single_int(x_dim, "x_dim"),
      as_single_int(y_dim, "y_dim"),
      condition_values,
      condition_dims,
      valid,
      false
    );

    SimulationResult exact = simulate_chi_gamma(
      built.slices,
      built.observed_chi,
      built.observed_gamma,
      as_single_int(nsim, "nsim"),
      as_single_int(seed, "seed"),
      true,
      true,
      as_single_bool(sequential, "sequential"),
      as_single_int(seq_limit, "seq_limit")
    );

    SEXP out = PROTECT(Rf_allocVector(VECSXP, 13));
    SEXP names = PROTECT(Rf_allocVector(STRSXP, 13));
    const char *out_names[] = {
      "chi", "df", "gamma", "ppq", "pmq", "s",
      "p_chi", "p_gamma", "nsim", "chi_exceed",
      "gamma_exceed", "rng_draws", "final_seed"
    };
    for (int i = 0; i < 13; ++i) SET_STRING_ELT(names, i, Rf_mkChar(out_names[i]));
    SET_VECTOR_ELT(out, 0, Rf_ScalarReal(built.observed_chi));
    SET_VECTOR_ELT(out, 1, Rf_ScalarInteger(built.observed_df));
    SET_VECTOR_ELT(out, 2, Rf_ScalarReal(built.observed_gamma));
    SET_VECTOR_ELT(out, 3, Rf_ScalarReal(built.ppq));
    SET_VECTOR_ELT(out, 4, Rf_ScalarReal(built.pmq));
    SET_VECTOR_ELT(out, 5, Rf_ScalarReal(built.s));
    SET_VECTOR_ELT(out, 6, Rf_ScalarReal(source_p_value(exact.chi_exceed, exact.nsim)));
    SET_VECTOR_ELT(out, 7, Rf_ScalarReal(source_p_value(exact.gamma_exceed, exact.nsim)));
    SET_VECTOR_ELT(out, 8, Rf_ScalarInteger(exact.nsim));
    SET_VECTOR_ELT(out, 9, Rf_ScalarInteger(exact.chi_exceed));
    SET_VECTOR_ELT(out, 10, Rf_ScalarInteger(exact.gamma_exceed));
    SET_VECTOR_ELT(out, 11, Rf_ScalarReal(static_cast<double>(exact.draw_count)));
    SET_VECTOR_ELT(out, 12, Rf_ScalarReal(static_cast<double>(exact.final_seed)));
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
