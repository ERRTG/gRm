#include <algorithm>
#include <climits>
#include <cmath>
#include <cstring>
#include <set>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include <R.h>
#include <Rinternals.h>

// Native backend for GLLRM expected margins.
//
// This file is a mechanical acceleration of
// R/gllrm_fit.R::calculate_gllrm_joint_expected_margins_r(), which is
// the source-readable reference for source/GLLRM_ESTIM.txt::
// CalculateBiasedGammaValues2, source/PAS_skunits/skbias12b.pas::
// Estimate_GLLRM, source/PAS_skunits/skbias22.pas::Gamma_calculation, and
// source/PAS_skunits/skbias22.pas::LD_Gamma_calculation. The native code only
// computes expected item, IJ/LD, and IX/DIF margins for the current state. R
// keeps ownership of the fitting loop, update equations, stopping rules,
// candidate ordering, likelihood calculation, and reporting gauge.

namespace {

struct GllrmExpectedError : public std::runtime_error {
  explicit GllrmExpectedError(const std::string &message) : std::runtime_error(message) {}
};

struct ComponentCache {
  std::vector<double> weights;
  double scale = 1.0;
  std::vector<double> gamma;
};

struct BackgroundCache {
  std::vector<ComponentCache> components;
  std::vector<std::vector<double>> rest;
  std::vector<double> full;
};

struct NativeComponent {
  std::vector<int> items;          // 1-based item indices from R.
  SEXP config_matrix = R_NilValue; // integer matrix, zero-based scores.
  std::vector<int> config_scores;
  SEXP ld_local = R_NilValue;      // integer matrix: ld_index, item1_pos, item2_pos.
};

struct NativeInput {
  int max_total_score = 0;
  std::vector<NativeComponent> components;
  std::vector<SEXP> dif_by_item;
  std::vector<int> item_raw_max;
  std::vector<int> group_scores;
  std::vector<double> group_counts;
  SEXP group_backgrounds = R_NilValue;
  std::vector<int> dif_backgrounds;
};

inline R_xlen_t index2(int row, int col, int nrow) {
  return static_cast<R_xlen_t>(row) +
    static_cast<R_xlen_t>(nrow) * static_cast<R_xlen_t>(col);
}

SEXP list_element(SEXP list, const char *name) {
  if (TYPEOF(list) != VECSXP) {
    throw GllrmExpectedError("Native GLLRM expected-margin input must be a list.");
  }
  SEXP names = Rf_getAttrib(list, R_NamesSymbol);
  if (Rf_length(names) != Rf_length(list)) {
    throw GllrmExpectedError("Native GLLRM expected-margin input must be a named list.");
  }
  for (R_xlen_t i = 0; i < XLENGTH(list); ++i) {
    if (std::strcmp(CHAR(STRING_ELT(names, i)), name) == 0) {
      return VECTOR_ELT(list, i);
    }
  }
  throw GllrmExpectedError(std::string("Native GLLRM expected-margin input is missing `") + name + "`.");
}

int as_single_int(SEXP value, const char *name) {
  if (Rf_length(value) != 1) {
    throw GllrmExpectedError(std::string("`") + name + "` must have length 1.");
  }
  if (TYPEOF(value) == INTSXP) {
    int out = INTEGER(value)[0];
    if (out == NA_INTEGER) throw GllrmExpectedError(std::string("`") + name + "` must not be NA.");
    return out;
  }
  if (TYPEOF(value) == REALSXP) {
    double out = REAL(value)[0];
    if (!R_FINITE(out)) throw GllrmExpectedError(std::string("`") + name + "` must be finite.");
    if (out != std::floor(out) || out < static_cast<double>(INT_MIN) ||
        out > static_cast<double>(INT_MAX)) {
      throw GllrmExpectedError(std::string("`") + name + "` must be integer-like.");
    }
    return static_cast<int>(out);
  }
  throw GllrmExpectedError(std::string("`") + name + "` must be numeric.");
}

std::vector<int> integer_vector(SEXP value, const char *name) {
  R_xlen_t n = XLENGTH(value);
  std::vector<int> out(static_cast<std::size_t>(n));
  if (TYPEOF(value) == INTSXP) {
    for (R_xlen_t i = 0; i < n; ++i) {
      int cell = INTEGER(value)[i];
      if (cell == NA_INTEGER) throw GllrmExpectedError(std::string("`") + name + "` must not contain NA.");
      out[static_cast<std::size_t>(i)] = cell;
    }
    return out;
  }
  if (TYPEOF(value) == REALSXP) {
    for (R_xlen_t i = 0; i < n; ++i) {
      double cell = REAL(value)[i];
      if (!R_FINITE(cell)) throw GllrmExpectedError(std::string("`") + name + "` must be finite.");
      if (cell != std::floor(cell) || cell < static_cast<double>(INT_MIN) ||
          cell > static_cast<double>(INT_MAX)) {
        throw GllrmExpectedError(std::string("`") + name + "` must contain integer-like values.");
      }
      out[static_cast<std::size_t>(i)] = static_cast<int>(cell);
    }
    return out;
  }
  throw GllrmExpectedError(std::string("`") + name + "` must be numeric.");
}

std::vector<double> double_vector(SEXP value, const char *name) {
  R_xlen_t n = XLENGTH(value);
  std::vector<double> out(static_cast<std::size_t>(n));
  if (TYPEOF(value) == REALSXP) {
    for (R_xlen_t i = 0; i < n; ++i) {
      double cell = REAL(value)[i];
      if (!R_FINITE(cell)) throw GllrmExpectedError(std::string("`") + name + "` must be finite.");
      out[static_cast<std::size_t>(i)] = cell;
    }
    return out;
  }
  if (TYPEOF(value) == INTSXP) {
    for (R_xlen_t i = 0; i < n; ++i) {
      int cell = INTEGER(value)[i];
      if (cell == NA_INTEGER) throw GllrmExpectedError(std::string("`") + name + "` must not contain NA.");
      out[static_cast<std::size_t>(i)] = static_cast<double>(cell);
    }
    return out;
  }
  throw GllrmExpectedError(std::string("`") + name + "` must be numeric.");
}

std::vector<SEXP> list_values(SEXP value, const char *name) {
  if (TYPEOF(value) != VECSXP) {
    throw GllrmExpectedError(std::string("`") + name + "` must be a list.");
  }
  std::vector<SEXP> out;
  out.reserve(static_cast<std::size_t>(XLENGTH(value)));
  for (R_xlen_t i = 0; i < XLENGTH(value); ++i) {
    out.push_back(VECTOR_ELT(value, i));
  }
  return out;
}

void require_matrix(SEXP value, const char *name) {
  SEXP dim = Rf_getAttrib(value, R_DimSymbol);
  if (TYPEOF(dim) != INTSXP || Rf_length(dim) != 2 ||
      INTEGER(dim)[0] < 0 || INTEGER(dim)[1] < 0) {
    throw GllrmExpectedError(std::string("`") + name + "` must be a matrix.");
  }
}

void require_numeric_matrix(SEXP value, const char *name) {
  require_matrix(value, name);
  if (TYPEOF(value) != REALSXP) {
    throw GllrmExpectedError(std::string("`") + name + "` must be a numeric matrix.");
  }
}

void require_numeric_matrix_list(SEXP value, const char *name) {
  if (TYPEOF(value) != VECSXP) {
    throw GllrmExpectedError(std::string("`") + name + "` must be a list.");
  }
  for (R_xlen_t i = 0; i < XLENGTH(value); ++i) {
    require_numeric_matrix(VECTOR_ELT(value, i), name);
  }
}

int matrix_nrows(SEXP matrix) {
  require_matrix(matrix, "matrix");
  return INTEGER(Rf_getAttrib(matrix, R_DimSymbol))[0];
}

int matrix_ncols(SEXP matrix) {
  require_matrix(matrix, "matrix");
  return INTEGER(Rf_getAttrib(matrix, R_DimSymbol))[1];
}

int integer_matrix_at(SEXP matrix, int row, int col, const char *name) {
  int nrow = matrix_nrows(matrix);
  int ncol = matrix_ncols(matrix);
  if (row < 0 || row >= nrow || col < 0 || col >= ncol) {
    throw GllrmExpectedError(std::string("`") + name + "` index is out of range.");
  }
  R_xlen_t idx = index2(row, col, nrow);
  if (TYPEOF(matrix) == INTSXP) {
    int out = INTEGER(matrix)[idx];
    if (out == NA_INTEGER) throw GllrmExpectedError(std::string("`") + name + "` must not contain NA.");
    return out;
  }
  if (TYPEOF(matrix) == REALSXP) {
    double out = REAL(matrix)[idx];
    if (!R_FINITE(out)) throw GllrmExpectedError(std::string("`") + name + "` must be finite.");
    if (out != std::floor(out) || out < static_cast<double>(INT_MIN) ||
        out > static_cast<double>(INT_MAX)) {
      throw GllrmExpectedError(std::string("`") + name + "` must contain integer-like values.");
    }
    return static_cast<int>(out);
  }
  throw GllrmExpectedError(std::string("`") + name + "` must be numeric.");
}

inline int integer_matrix_at_unchecked(SEXP matrix, int row, int col) {
  int nrow = INTEGER(Rf_getAttrib(matrix, R_DimSymbol))[0];
  R_xlen_t idx = index2(row, col, nrow);
  return TYPEOF(matrix) == INTSXP
    ? INTEGER(matrix)[idx]
    : static_cast<int>(REAL(matrix)[idx]);
}

inline double numeric_matrix_at_unchecked(SEXP matrix, int row, int col) {
  int nrow = INTEGER(Rf_getAttrib(matrix, R_DimSymbol))[0];
  return REAL(matrix)[index2(row, col, nrow)];
}

void require_integer_like_matrix(SEXP value, const char *name) {
  require_matrix(value, name);
  if (TYPEOF(value) != INTSXP && TYPEOF(value) != REALSXP) {
    throw GllrmExpectedError(std::string("`") + name + "` must be numeric.");
  }
  int nrow = matrix_nrows(value);
  int ncol = matrix_ncols(value);
  for (int col = 0; col < ncol; ++col) {
    for (int row = 0; row < nrow; ++row) {
      (void) integer_matrix_at(value, row, col, name);
    }
  }
}

void require_finite_nonnegative_matrix(SEXP value, const char *name) {
  require_numeric_matrix(value, name);
  for (R_xlen_t index = 0; index < XLENGTH(value); ++index) {
    double cell = REAL(value)[index];
    if (!R_FINITE(cell) || cell < 0.0) {
      throw GllrmExpectedError(
        std::string("`") + name + "` must contain non-negative finite values."
      );
    }
  }
}

NativeInput parse_native_input(SEXP native_input) {
  NativeInput input;
  SEXP components = list_element(native_input, "components");
  SEXP config_matrices = list_element(native_input, "config_matrices");
  SEXP config_scores = list_element(native_input, "config_scores");
  SEXP ld_local_matrices = list_element(native_input, "ld_local_matrices");
  std::vector<SEXP> component_values = list_values(components, "components");
  std::vector<SEXP> config_matrix_values = list_values(config_matrices, "config_matrices");
  std::vector<SEXP> config_score_values = list_values(config_scores, "config_scores");
  std::vector<SEXP> ld_local_values = list_values(ld_local_matrices, "ld_local_matrices");
  if (component_values.size() != config_matrix_values.size() ||
      component_values.size() != config_score_values.size() ||
      component_values.size() != ld_local_values.size()) {
    throw GllrmExpectedError("GLLRM expected-margin component metadata lists must have equal length.");
  }

  input.item_raw_max = integer_vector(list_element(native_input, "item_raw_max"), "item_raw_max");
  input.dif_by_item = list_values(list_element(native_input, "dif_by_item"), "dif_by_item");
  if (input.item_raw_max.size() != input.dif_by_item.size()) {
    throw GllrmExpectedError("`item_raw_max` must have one entry per DIF-by-item matrix.");
  }
  for (std::size_t item = 0; item < input.dif_by_item.size(); ++item) {
    require_matrix(input.dif_by_item[item], "dif_by_item");
    if (matrix_ncols(input.dif_by_item[item]) != 2) {
      throw GllrmExpectedError("Each DIF-by-item matrix must have two columns.");
    }
  }

  input.components.reserve(component_values.size());
  for (std::size_t i = 0; i < component_values.size(); ++i) {
    NativeComponent component;
    component.items = integer_vector(component_values[i], "components");
    component.config_matrix = config_matrix_values[i];
    component.config_scores = integer_vector(config_score_values[i], "config_scores");
    component.ld_local = ld_local_values[i];
    require_matrix(component.config_matrix, "config_matrices");
    require_matrix(component.ld_local, "ld_local_matrices");
    if (matrix_nrows(component.config_matrix) != static_cast<int>(component.config_scores.size())) {
      throw GllrmExpectedError("Each component config matrix must have one row per config score.");
    }
    if (matrix_ncols(component.config_matrix) != static_cast<int>(component.items.size())) {
      throw GllrmExpectedError("Each component config matrix must have one column per component item.");
    }
    if (matrix_ncols(component.ld_local) != 3) {
      throw GllrmExpectedError("Each local LD matrix must have three columns.");
    }
    input.components.push_back(component);
  }

  input.max_total_score = as_single_int(list_element(native_input, "max_total_score"), "max_total_score");
  input.group_scores = integer_vector(list_element(native_input, "group_scores"), "group_scores");
  input.group_counts = double_vector(list_element(native_input, "group_counts"), "group_counts");
  if (input.group_scores.size() != input.group_counts.size()) {
    throw GllrmExpectedError("`group_scores` and `group_counts` must have equal length.");
  }
  input.group_backgrounds = list_element(native_input, "group_backgrounds");
  require_matrix(input.group_backgrounds, "group_backgrounds");
  if (matrix_nrows(input.group_backgrounds) != static_cast<int>(input.group_scores.size())) {
    throw GllrmExpectedError("`group_backgrounds` must have one row per score/exogenous group.");
  }
  input.dif_backgrounds = integer_vector(list_element(native_input, "dif_backgrounds"), "dif_backgrounds");
  int background_count = matrix_ncols(input.group_backgrounds);
  for (int background : input.dif_backgrounds) {
    if (background < 1 || background > background_count) {
      throw GllrmExpectedError("`dif_backgrounds` contains an out-of-range background index.");
    }
  }
  return input;
}

// Validate the complete index graph once before any expected-margin loop.
// Source trace: source/PAS_skunits/skbias12b.pas::Initialize_GLLRMinfo and
// source/GLLRM_ESTIM.txt::CalculateBiasedGammaValues2 rely on globally bounded
// Pascal arrays. The R native boundary must establish those same bounds
// explicitly because malformed list metadata otherwise becomes unchecked
// pointer arithmetic below.
void validate_native_input_graph(const NativeInput &input,
                                 SEXP item_gamma,
                                 const std::vector<SEXP> &ld_parameters,
                                 const std::vector<SEXP> &dif_parameters) {
  const int n_items = static_cast<int>(input.item_raw_max.size());
  if (n_items < 1) {
    throw GllrmExpectedError("Native GLLRM input must contain at least one item.");
  }
  if (input.max_total_score < 0) {
    throw GllrmExpectedError("`max_total_score` must not be negative.");
  }

  require_finite_nonnegative_matrix(item_gamma, "item_gamma");
  if (matrix_nrows(item_gamma) != n_items) {
    throw GllrmExpectedError("`item_gamma` must have one row per item.");
  }
  const int item_gamma_columns = matrix_ncols(item_gamma);
  long long possible_total_score = 0;
  for (int item = 0; item < n_items; ++item) {
    int category_count = input.item_raw_max[static_cast<std::size_t>(item)];
    if (category_count < 1 || category_count > item_gamma_columns) {
      throw GllrmExpectedError(
        "`item_raw_max` entries must be positive and within `item_gamma` columns."
      );
    }
    possible_total_score += static_cast<long long>(category_count - 1);
    if (possible_total_score > INT_MAX) {
      throw GllrmExpectedError("The declared item score range is too large.");
    }
  }
  if (input.max_total_score != static_cast<int>(possible_total_score)) {
    throw GllrmExpectedError(
      "`max_total_score` must equal the sum of declared item maximum scores."
    );
  }

  if (input.components.empty()) {
    throw GllrmExpectedError("Native GLLRM input must contain item components.");
  }
  std::vector<int> item_membership(static_cast<std::size_t>(n_items), 0);
  std::vector<int> ld_reference_count(ld_parameters.size(), 0);
  for (std::size_t component_index = 0;
       component_index < input.components.size();
       ++component_index) {
    const NativeComponent &component = input.components[component_index];
    if (component.items.empty()) {
      throw GllrmExpectedError("GLLRM item components must not be empty.");
    }
    require_integer_like_matrix(component.config_matrix, "config_matrices");
    require_integer_like_matrix(component.ld_local, "ld_local_matrices");

    for (int item_one_based : component.items) {
      if (item_one_based < 1 || item_one_based > n_items) {
        throw GllrmExpectedError("`components` contains an out-of-range item index.");
      }
      int item = item_one_based - 1;
      ++item_membership[static_cast<std::size_t>(item)];
      if (item_membership[static_cast<std::size_t>(item)] > 1) {
        throw GllrmExpectedError("Each item must occur in exactly one GLLRM component.");
      }
    }

    int n_config = matrix_nrows(component.config_matrix);
    for (int config = 0; config < n_config; ++config) {
      int score_sum = 0;
      for (int local = 0;
           local < static_cast<int>(component.items.size());
           ++local) {
        int item = component.items[static_cast<std::size_t>(local)] - 1;
        int score = integer_matrix_at(
          component.config_matrix,
          config,
          local,
          "config_matrices"
        );
        int category_count = input.item_raw_max[static_cast<std::size_t>(item)];
        if (score < 0 || score >= category_count) {
          throw GllrmExpectedError(
            "`config_matrices` contains a score outside its item's category range."
          );
        }
        score_sum += score;
      }
      int declared_score = component.config_scores[static_cast<std::size_t>(config)];
      if (declared_score < 0 || declared_score > input.max_total_score) {
        throw GllrmExpectedError("`config_scores` contains an out-of-range score.");
      }
      if (declared_score != score_sum) {
        throw GllrmExpectedError(
          "Each `config_scores` value must equal its configuration row sum."
        );
      }
    }

    int n_ld_rows = matrix_nrows(component.ld_local);
    for (int row = 0; row < n_ld_rows; ++row) {
      int ld_index = integer_matrix_at(
        component.ld_local,
        row,
        0,
        "ld_local_matrices"
      );
      int item1_pos = integer_matrix_at(
        component.ld_local,
        row,
        1,
        "ld_local_matrices"
      );
      int item2_pos = integer_matrix_at(
        component.ld_local,
        row,
        2,
        "ld_local_matrices"
      );
      if (ld_index < 1 || ld_index > static_cast<int>(ld_parameters.size())) {
        throw GllrmExpectedError(
          "`ld_local_matrices` contains an out-of-range LD parameter index."
        );
      }
      if (item1_pos < 1 ||
          item1_pos > static_cast<int>(component.items.size()) ||
          item2_pos < 1 ||
          item2_pos > static_cast<int>(component.items.size()) ||
          item1_pos == item2_pos) {
        throw GllrmExpectedError(
          "`ld_local_matrices` contains invalid local item positions."
        );
      }
      ++ld_reference_count[static_cast<std::size_t>(ld_index - 1)];
      int item1 = component.items[static_cast<std::size_t>(item1_pos - 1)] - 1;
      int item2 = component.items[static_cast<std::size_t>(item2_pos - 1)] - 1;
      SEXP parameter = ld_parameters[static_cast<std::size_t>(ld_index - 1)];
      if (matrix_nrows(parameter) != input.item_raw_max[static_cast<std::size_t>(item1)] ||
          matrix_ncols(parameter) != input.item_raw_max[static_cast<std::size_t>(item2)]) {
        throw GllrmExpectedError(
          "Each LD parameter matrix must match its two item category counts."
        );
      }
    }
  }
  if (std::any_of(item_membership.begin(), item_membership.end(), [](int count) {
        return count != 1;
      })) {
    throw GllrmExpectedError("Each item must occur in exactly one GLLRM component.");
  }

  for (std::size_t index = 0; index < ld_parameters.size(); ++index) {
    require_finite_nonnegative_matrix(ld_parameters[index], "ld_parameters");
    if (ld_reference_count[index] != 1) {
      throw GllrmExpectedError(
        "Every LD parameter matrix must be referenced exactly once."
      );
    }
  }

  int background_count = matrix_ncols(input.group_backgrounds);
  require_integer_like_matrix(input.group_backgrounds, "group_backgrounds");
  for (int group = 0; group < matrix_nrows(input.group_backgrounds); ++group) {
    for (int background = 0; background < background_count; ++background) {
      if (integer_matrix_at(
            input.group_backgrounds,
            group,
            background,
            "group_backgrounds"
          ) < 1) {
        throw GllrmExpectedError(
          "`group_backgrounds` values must be positive one-based categories."
        );
      }
    }
  }
  for (std::size_t group = 0; group < input.group_scores.size(); ++group) {
    int score = input.group_scores[group];
    if (score < 0 || score > input.max_total_score) {
      throw GllrmExpectedError("`group_scores` contains an out-of-range score.");
    }
    if (input.group_counts[group] < 0.0) {
      throw GllrmExpectedError("`group_counts` must contain non-negative values.");
    }
  }

  std::vector<int> dif_reference_count(dif_parameters.size(), 0);
  std::set<int> referenced_backgrounds;
  for (int item = 0; item < n_items; ++item) {
    SEXP mappings = input.dif_by_item[static_cast<std::size_t>(item)];
    require_integer_like_matrix(mappings, "dif_by_item");
    int n_mapping = matrix_nrows(mappings);
    for (int row = 0; row < n_mapping; ++row) {
      int background = integer_matrix_at(mappings, row, 0, "dif_by_item");
      int dif_index = integer_matrix_at(mappings, row, 1, "dif_by_item");
      if (background < 1 || background > background_count) {
        throw GllrmExpectedError(
          "`dif_by_item` contains an out-of-range background index."
        );
      }
      if (dif_index < 1 || dif_index > static_cast<int>(dif_parameters.size())) {
        throw GllrmExpectedError(
          "`dif_by_item` contains an out-of-range DIF parameter index."
        );
      }
      referenced_backgrounds.insert(background);
      ++dif_reference_count[static_cast<std::size_t>(dif_index - 1)];
      SEXP parameter = dif_parameters[static_cast<std::size_t>(dif_index - 1)];
      if (matrix_nrows(parameter) != input.item_raw_max[static_cast<std::size_t>(item)]) {
        throw GllrmExpectedError(
          "Each DIF parameter matrix must match its item category count."
        );
      }
      int parameter_columns = matrix_ncols(parameter);
      for (int group = 0; group < matrix_nrows(input.group_backgrounds); ++group) {
        int category = integer_matrix_at(
          input.group_backgrounds,
          group,
          background - 1,
          "group_backgrounds"
        );
        if (category > parameter_columns) {
          throw GllrmExpectedError(
            "A `group_backgrounds` category exceeds its DIF parameter columns."
          );
        }
      }
    }
  }
  for (std::size_t index = 0; index < dif_parameters.size(); ++index) {
    require_finite_nonnegative_matrix(dif_parameters[index], "dif_parameters");
    if (dif_reference_count[index] != 1) {
      throw GllrmExpectedError(
        "Every DIF parameter matrix must be referenced exactly once."
      );
    }
  }

  std::set<int> declared_backgrounds;
  for (int background : input.dif_backgrounds) {
    if (!declared_backgrounds.insert(background).second) {
      throw GllrmExpectedError("`dif_backgrounds` must not contain duplicates.");
    }
  }
  if (declared_backgrounds != referenced_backgrounds) {
    throw GllrmExpectedError(
      "`dif_backgrounds` must exactly match backgrounds used by DIF mappings."
    );
  }
}

// Source trace: source/PAS_skunits/skbias22.pas::Gamma_calculation combines
// component score polynomials by truncated convolution over total score.
std::vector<double> convolve_truncated(const std::vector<double> &a,
                                       const std::vector<double> &b,
                                       int max_score) {
  std::vector<double> out(static_cast<std::size_t>(max_score) + 1U, 0.0);
  for (int i = 0; i <= max_score; ++i) {
    if (a[static_cast<std::size_t>(i)] == 0.0) continue;
    for (int j = 0; j + i <= max_score; ++j) {
      out[static_cast<std::size_t>(i + j)] +=
        a[static_cast<std::size_t>(i)] * b[static_cast<std::size_t>(j)];
    }
  }
  return out;
}

// Only included IX/DIF backgrounds affect component weights, so the cache key
// mirrors gllrm_background_cache_key(): GLLRM values joined by carriage return.
std::string background_key(SEXP group_backgrounds,
                           int group,
                           const std::vector<int> &dif_backgrounds) {
  if (dif_backgrounds.empty()) {
    return ".";
  }
  std::string out;
  for (std::size_t i = 0; i < dif_backgrounds.size(); ++i) {
    if (i > 0) out.push_back('\r');
    int background = dif_backgrounds[i] - 1;
    out += std::to_string(integer_matrix_at_unchecked(
      group_backgrounds,
      group,
      background
    ));
  }
  return out;
}

ComponentCache build_component_cache(const NativeInput &input,
                                     int component_index,
                                     int group_index,
                                     SEXP item_gamma,
                                     const std::vector<SEXP> &ld_parameters,
                                     const std::vector<SEXP> &dif_parameters) {
  const NativeComponent &component = input.components[static_cast<std::size_t>(component_index)];
  int n_config = matrix_nrows(component.config_matrix);
  ComponentCache cache;
  cache.weights.assign(static_cast<std::size_t>(n_config), 0.0);
  cache.gamma.assign(static_cast<std::size_t>(input.max_total_score) + 1U, 0.0);

  // Mirrors CalculateBiasedGammaValues2 / LD_Gamma_calculation: each component
  // configuration receives item, IX/DIF, and IJ/LD parameter factors.
  for (int config = 0; config < n_config; ++config) {
    double weight = 1.0;
    for (int local = 0; local < static_cast<int>(component.items.size()); ++local) {
      int item = component.items[static_cast<std::size_t>(local)] - 1;
      int score = integer_matrix_at_unchecked(
        component.config_matrix,
        config,
        local
      );
      // Item score parameters are the base component weight in the source loop.
      weight *= numeric_matrix_at_unchecked(item_gamma, item, score);

      SEXP dif_rows = input.dif_by_item[static_cast<std::size_t>(item)];
      int n_dif_rows = matrix_nrows(dif_rows);
      for (int dif_row = 0; dif_row < n_dif_rows; ++dif_row) {
        int background = integer_matrix_at_unchecked(dif_rows, dif_row, 0) - 1;
        int dif_index = integer_matrix_at_unchecked(dif_rows, dif_row, 1) - 1;
        int background_value = integer_matrix_at_unchecked(
          input.group_backgrounds,
          group_index,
          background
        );
        // IX/DIF parameters are selected by the relevant exogenous background
        // value for the current score/exogenous group.
        weight *= numeric_matrix_at_unchecked(
          dif_parameters[static_cast<std::size_t>(dif_index)],
          score,
          background_value - 1
        );
      }
    }

    int n_ld_rows = matrix_nrows(component.ld_local);
    for (int ld_row = 0; ld_row < n_ld_rows; ++ld_row) {
      int ld_index = integer_matrix_at_unchecked(component.ld_local, ld_row, 0) - 1;
      int item1_pos = integer_matrix_at_unchecked(component.ld_local, ld_row, 1) - 1;
      int item2_pos = integer_matrix_at_unchecked(component.ld_local, ld_row, 2) - 1;
      int score1 = integer_matrix_at_unchecked(
        component.config_matrix,
        config,
        item1_pos
      );
      int score2 = integer_matrix_at_unchecked(
        component.config_matrix,
        config,
        item2_pos
      );
      // IJ/LD parameters apply only to LD terms whose two items are inside the
      // current LD-connected component.
      weight *= numeric_matrix_at_unchecked(
        ld_parameters[static_cast<std::size_t>(ld_index)],
        score1,
        score2
      );
    }

    cache.weights[static_cast<std::size_t>(config)] = weight;
    int component_score = component.config_scores[static_cast<std::size_t>(config)];
    if (weight > 0.0 && component_score >= 0 && component_score <= input.max_total_score) {
      cache.gamma[static_cast<std::size_t>(component_score)] += weight;
    }
  }

  // Component gamma normalization follows the existing R source-stability
  // convention and preserves the downstream expected-margin ratios.
  cache.scale = *std::max_element(cache.gamma.begin(), cache.gamma.end());
  if (cache.scale <= 0.0) {
    cache.scale = 1.0;
  }
  for (double &value : cache.gamma) {
    value /= cache.scale;
  }
  return cache;
}

BackgroundCache build_background_cache(const NativeInput &input,
                                       int group_index,
                                       SEXP item_gamma,
                                       const std::vector<SEXP> &ld_parameters,
                                       const std::vector<SEXP> &dif_parameters) {
  BackgroundCache cache;
  int n_component = static_cast<int>(input.components.size());
  cache.components.reserve(static_cast<std::size_t>(n_component));
  for (int component = 0; component < n_component; ++component) {
    cache.components.push_back(
      build_component_cache(input, component, group_index, item_gamma, ld_parameters, dif_parameters)
    );
  }

  // Component weights and convolutions are cached by DIF background pattern,
  // matching calculate_gllrm_joint_expected_margins_r() and
  // gllrm_background_cache_key().
  std::vector<double> unit(static_cast<std::size_t>(input.max_total_score) + 1U, 0.0);
  unit[0] = 1.0;
  std::vector<std::vector<double>> prefix(static_cast<std::size_t>(n_component) + 1U, unit);
  std::vector<std::vector<double>> suffix(static_cast<std::size_t>(n_component) + 1U, unit);
  for (int i = 0; i < n_component; ++i) {
    prefix[static_cast<std::size_t>(i + 1)] = convolve_truncated(
      prefix[static_cast<std::size_t>(i)],
      cache.components[static_cast<std::size_t>(i)].gamma,
      input.max_total_score
    );
  }
  for (int i = n_component - 1; i >= 0; --i) {
    suffix[static_cast<std::size_t>(i)] = convolve_truncated(
      cache.components[static_cast<std::size_t>(i)].gamma,
      suffix[static_cast<std::size_t>(i + 1)],
      input.max_total_score
    );
  }
  cache.full = prefix[static_cast<std::size_t>(n_component)];
  cache.rest.reserve(static_cast<std::size_t>(n_component));
  for (int i = 0; i < n_component; ++i) {
    cache.rest.push_back(convolve_truncated(
      prefix[static_cast<std::size_t>(i)],
      suffix[static_cast<std::size_t>(i + 1)],
      input.max_total_score
    ));
  }
  return cache;
}

// Source trace: source/GLLRM_ESTIM.txt::CalculateBiasedGammaValues2 accumulates
// fitted item, IJ/LD, and IX/DIF cells conditional on total score.
void accumulate_group_expected(const NativeInput &input,
                               const BackgroundCache &cache,
                               int group_index,
                               SEXP expected_items,
                               std::vector<SEXP> &expected_ld,
                               std::vector<SEXP> &expected_dif) {
  int total_score = input.group_scores[static_cast<std::size_t>(group_index)];
  if (total_score < 0 || total_score > input.max_total_score) {
    return;
  }
  double denominator = cache.full[static_cast<std::size_t>(total_score)];
  if (denominator <= 0.0) {
    return;
  }
  double group_count = input.group_counts[static_cast<std::size_t>(group_index)];
  int expected_item_rows = matrix_nrows(expected_items);

  for (int component_index = 0; component_index < static_cast<int>(input.components.size()); ++component_index) {
    const NativeComponent &component = input.components[static_cast<std::size_t>(component_index)];
    const ComponentCache &component_cache = cache.components[static_cast<std::size_t>(component_index)];
    const std::vector<double> &rest_gamma = cache.rest[static_cast<std::size_t>(component_index)];
    int n_config = matrix_nrows(component.config_matrix);

    for (int config = 0; config < n_config; ++config) {
      int component_score = component.config_scores[static_cast<std::size_t>(config)];
      if (component_score > total_score) continue;
      double component_weight = component_cache.weights[static_cast<std::size_t>(config)];
      if (component_weight <= 0.0) continue;
      // expected = group_count * (component_weight / component_scale) *
      //   rest_gamma[total_score - component_score] / denominator
      double expected = group_count * (component_weight / component_cache.scale) *
        rest_gamma[static_cast<std::size_t>(total_score - component_score)] / denominator;
      if (expected <= 0.0) continue;

      for (int local = 0; local < static_cast<int>(component.items.size()); ++local) {
        int item = component.items[static_cast<std::size_t>(local)] - 1;
        int score = integer_matrix_at_unchecked(
          component.config_matrix,
          config,
          local
        );
        REAL(expected_items)[index2(item, score, expected_item_rows)] += expected;

        SEXP dif_rows = input.dif_by_item[static_cast<std::size_t>(item)];
        int n_dif_rows = matrix_nrows(dif_rows);
        for (int dif_row = 0; dif_row < n_dif_rows; ++dif_row) {
          int background = integer_matrix_at_unchecked(dif_rows, dif_row, 0) - 1;
          int dif_index = integer_matrix_at_unchecked(dif_rows, dif_row, 1) - 1;
          int background_value = integer_matrix_at_unchecked(
            input.group_backgrounds,
            group_index,
            background
          );
          SEXP matrix = expected_dif[static_cast<std::size_t>(dif_index)];
          REAL(matrix)[index2(score, background_value - 1, matrix_nrows(matrix))] += expected;
        }
      }

      int n_ld_rows = matrix_nrows(component.ld_local);
      for (int ld_row = 0; ld_row < n_ld_rows; ++ld_row) {
        int ld_index = integer_matrix_at_unchecked(component.ld_local, ld_row, 0) - 1;
        int item1_pos = integer_matrix_at_unchecked(component.ld_local, ld_row, 1) - 1;
        int item2_pos = integer_matrix_at_unchecked(component.ld_local, ld_row, 2) - 1;
        int score1 = integer_matrix_at_unchecked(
          component.config_matrix,
          config,
          item1_pos
        );
        int score2 = integer_matrix_at_unchecked(
          component.config_matrix,
          config,
          item2_pos
        );
        SEXP matrix = expected_ld[static_cast<std::size_t>(ld_index)];
        REAL(matrix)[index2(score1, score2, matrix_nrows(matrix))] += expected;
      }
    }
  }
}

std::vector<SEXP> duplicate_zero_list(SEXP source) {
  std::vector<SEXP> out;
  out.reserve(static_cast<std::size_t>(XLENGTH(source)));
  for (R_xlen_t i = 0; i < XLENGTH(source); ++i) {
    SEXP one = VECTOR_ELT(source, i);
    if (TYPEOF(one) != REALSXP) {
      throw GllrmExpectedError("GLLRM parameter matrices must be numeric.");
    }
    std::fill(REAL(one), REAL(one) + XLENGTH(one), 0.0);
    out.push_back(one);
  }
  return out;
}

}  // namespace

extern "C" SEXP gRm_gllrm_expected_margins(SEXP native_input,
                                            SEXP item_gamma,
                                            SEXP ld_parameters,
                                            SEXP dif_parameters) {
  try {
    NativeInput input = parse_native_input(native_input);
    require_numeric_matrix(item_gamma, "item_gamma");
    require_numeric_matrix_list(ld_parameters, "ld_parameters");
    require_numeric_matrix_list(dif_parameters, "dif_parameters");

    std::vector<SEXP> ld_parameter_values = list_values(ld_parameters, "ld_parameters");
    std::vector<SEXP> dif_parameter_values = list_values(dif_parameters, "dif_parameters");
    validate_native_input_graph(
      input,
      item_gamma,
      ld_parameter_values,
      dif_parameter_values
    );

    SEXP expected_items = PROTECT(Rf_duplicate(item_gamma));
    std::fill(REAL(expected_items), REAL(expected_items) + XLENGTH(expected_items), 0.0);

    SEXP expected_ld = PROTECT(Rf_duplicate(ld_parameters));
    std::vector<SEXP> expected_ld_values = duplicate_zero_list(expected_ld);

    SEXP expected_dif = PROTECT(Rf_duplicate(dif_parameters));
    std::vector<SEXP> expected_dif_values = duplicate_zero_list(expected_dif);

    std::unordered_map<std::string, BackgroundCache> background_cache;
    for (int group = 0; group < static_cast<int>(input.group_scores.size()); ++group) {
      std::string key = background_key(input.group_backgrounds, group, input.dif_backgrounds);
      auto found = background_cache.find(key);
      if (found == background_cache.end()) {
        found = background_cache.emplace(
          key,
          build_background_cache(input, group, item_gamma, ld_parameter_values, dif_parameter_values)
        ).first;
      }
      accumulate_group_expected(
        input,
        found->second,
        group,
        expected_items,
        expected_ld_values,
        expected_dif_values
      );
    }

    // R owns all fitting-loop state transitions. The native boundary returns
    // only expected margins for the current GLLRM state.
    SEXP out = PROTECT(Rf_allocVector(VECSXP, 3));
    SEXP names = PROTECT(Rf_allocVector(STRSXP, 3));
    SET_STRING_ELT(names, 0, Rf_mkChar("expected_items"));
    SET_STRING_ELT(names, 1, Rf_mkChar("expected_ld"));
    SET_STRING_ELT(names, 2, Rf_mkChar("expected_dif"));
    SET_VECTOR_ELT(out, 0, expected_items);
    SET_VECTOR_ELT(out, 1, expected_ld);
    SET_VECTOR_ELT(out, 2, expected_dif);
    Rf_setAttrib(out, R_NamesSymbol, names);
    UNPROTECT(5);
    return out;
  } catch (const std::exception &err) {
    Rf_error("%s", err.what());
  }
  return R_NilValue;
}
