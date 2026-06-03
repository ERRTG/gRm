#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <vector>

#include <R.h>
#include <Rinternals.h>

namespace {

#if defined(__GNUC__) && (defined(__i386__) || defined(__x86_64__))
class X87ExtendedPrecisionScope {
 public:
  X87ExtendedPrecisionScope() {
    __asm__ __volatile__("fnstcw %0" : "=m"(original_));
    unsigned short extended = static_cast<unsigned short>((original_ & ~0x0300U) | 0x0300U);
    __asm__ __volatile__("fldcw %0" : : "m"(extended));
  }
  ~X87ExtendedPrecisionScope() {
    __asm__ __volatile__("fldcw %0" : : "m"(original_));
  }

 private:
  unsigned short original_;
};
#else
class X87ExtendedPrecisionScope {
 public:
  X87ExtendedPrecisionScope() = default;
};
#endif

double as_single_double(SEXP value, const char *name) {
  if (Rf_length(value) != 1) {
    throw std::runtime_error(std::string("`") + name + "` must have length 1.");
  }
  double out;
  if (TYPEOF(value) == REALSXP) {
    out = REAL(value)[0];
  } else if (TYPEOF(value) == INTSXP) {
    if (INTEGER(value)[0] == NA_INTEGER) {
      throw std::runtime_error(std::string("`") + name + "` must not be NA.");
    }
    out = static_cast<double>(INTEGER(value)[0]);
  } else {
    throw std::runtime_error(std::string("`") + name + "` must be numeric.");
  }
  if (!R_FINITE(out)) {
    throw std::runtime_error(std::string("`") + name + "` must be finite.");
  }
  return out;
}

int as_single_int(SEXP value, const char *name) {
  if (Rf_length(value) != 1) {
    throw std::runtime_error(std::string("`") + name + "` must have length 1.");
  }
  int out;
  if (TYPEOF(value) == INTSXP) {
    out = INTEGER(value)[0];
    if (out == NA_INTEGER) {
      throw std::runtime_error(std::string("`") + name + "` must not be NA.");
    }
  } else if (TYPEOF(value) == REALSXP) {
    double numeric = REAL(value)[0];
    if (!R_FINITE(numeric)) {
      throw std::runtime_error(std::string("`") + name + "` must be finite.");
    }
    out = static_cast<int>(numeric);
  } else {
    throw std::runtime_error(std::string("`") + name + "` must be numeric.");
  }
  return out;
}

std::string field9_3_long_double(long double value) {
  char buffer[64];
  std::snprintf(buffer, sizeof(buffer), "%9.3Lf", value);
  return std::string(buffer);
}

inline int gamma_index(int item_index, int score, int max_category_count) {
  return item_index * max_category_count + score;
}

long double source_ln(long double value) {
#if defined(__GNUC__) && (defined(__i386__) || defined(__x86_64__))
  long double result;
  __asm__ __volatile__(
    "fldln2\n\t"
    "fldt %1\n\t"
    "fyl2x\n\t"
    "fstpt %0"
    : "=m"(result)
    : "m"(value)
    : "st"
  );
  return result;
#else
  return std::log(value);
#endif
}

std::string top_ice_field_from_gamma(long double gamma_top, int k) {
  volatile long double z = source_ln(gamma_top) / static_cast<long double>(k);
  long double top_ice = source_ln(gamma_top) - static_cast<long double>(k) * z;
  if (std::isnan(top_ice) || std::fabs(top_ice) >= 0.0005L) {
    return field9_3_long_double(top_ice);
  }
  return std::signbit(top_ice) ? "   -0.000" : "    0.000";
}

std::vector<std::string> skbias12_ice_fields_from_gamma(const std::vector<long double> &item_gamma,
                                                        int item_index,
                                                        int max_score,
                                                        int max_category_count) {
  std::vector<std::string> fields(static_cast<std::size_t>(max_category_count), "     ....");
  std::vector<long double> ice(static_cast<std::size_t>(max_category_count), 0.0L);
  if (max_score == 1) {
    long double gamma1 = item_gamma[gamma_index(item_index, 1, max_category_count)];
    if (gamma1 > 0.0L) {
      // skbias12 leaves ICE itself at zero in this branch and stores the
      // single-step value as ItemEffect.
      ice[0] = 0.0L;
      ice[1] = 0.0L;
    }
  } else if (max_score > 1) {
    long double location = 0.0L;
    for (int score = 1; score <= max_score; ++score) {
      long double current = item_gamma[gamma_index(item_index, score, max_category_count)];
      long double previous = item_gamma[gamma_index(item_index, score - 1, max_category_count)];
      long double pcm = 0.0L;
      if (current <= 0.000000000001L) {
        pcm = 9999999.0L;
      } else if (previous <= 0.000000000001L) {
        pcm = -9999999.0L;
      } else {
        pcm = source_ln(previous / current);
        location += pcm;
      }
    }
    location /= static_cast<long double>(max_score);

    long double theta = std::exp(location);
    long double thfactor = 1.0L / theta;
    std::vector<long double> probs(static_cast<std::size_t>(max_category_count), 0.0L);
    for (int score = 0; score <= max_score; ++score) {
      thfactor *= theta;
      probs[score] =
        item_gamma[gamma_index(item_index, score, max_category_count)] * thfactor;
    }

    if (probs[max_score] > 0.0L && probs[0] > 0.0L) {
      for (int score = 0; score <= max_score; ++score) {
        long double mice = probs[score] / probs[0];
        ice[score] = mice > 0.0L ? source_ln(mice) : -9.999L;
      }
    }
  }

  for (int score = 0; score <= max_score; ++score) {
    fields[score] = field9_3_long_double(ice[score]);
  }
  return fields;
}

struct ExtendedFitInput {
  int n_items = 0;
  int max_category_count = 0;
  int max_total_score = 0;
  std::vector<int> raw_max;
  std::vector<long double> item_counts;
  std::vector<long double> score_counts;
};

ExtendedFitInput parse_extended_fit_input(SEXP item_counts_sexp,
                                          SEXP score_counts_sexp,
                                          SEXP raw_max_sexp) {
  SEXP dim = Rf_getAttrib(item_counts_sexp, R_DimSymbol);
  if (Rf_length(dim) != 2) {
    throw std::runtime_error("`item_counts` must be a matrix.");
  }
  ExtendedFitInput input;
  input.n_items = INTEGER(dim)[0];
  input.max_category_count = INTEGER(dim)[1];
  if (input.n_items < 1 || input.max_category_count < 1) {
    throw std::runtime_error("`item_counts` must have positive dimensions.");
  }
  if (Rf_length(raw_max_sexp) != input.n_items) {
    throw std::runtime_error("`raw_max` must have one entry per item.");
  }

  input.raw_max.resize(input.n_items);
  int raw_max_total = 0;
  for (int item_index = 0; item_index < input.n_items; ++item_index) {
    int value;
    if (TYPEOF(raw_max_sexp) == INTSXP) {
      value = INTEGER(raw_max_sexp)[item_index];
      if (value == NA_INTEGER) throw std::runtime_error("`raw_max` must not contain NA.");
    } else if (TYPEOF(raw_max_sexp) == REALSXP) {
      double numeric = REAL(raw_max_sexp)[item_index];
      if (!R_FINITE(numeric)) throw std::runtime_error("`raw_max` must be finite.");
      value = static_cast<int>(numeric);
    } else {
      throw std::runtime_error("`raw_max` must be numeric.");
    }
    if (value < 1 || value > input.max_category_count) {
      throw std::runtime_error("`raw_max` entries must be valid category counts.");
    }
    input.raw_max[item_index] = value;
    raw_max_total += value - 1;
  }

  input.max_total_score = Rf_length(score_counts_sexp) - 1;
  if (input.max_total_score < 0 || input.max_total_score != raw_max_total) {
    throw std::runtime_error("`score_counts` has the wrong length.");
  }

  input.item_counts.assign(static_cast<std::size_t>(input.n_items) * input.max_category_count, 0.0L);
  for (int item_index = 0; item_index < input.n_items; ++item_index) {
    for (int score = 0; score < input.max_category_count; ++score) {
      R_xlen_t r_index = item_index + static_cast<R_xlen_t>(input.n_items) * score;
      double value;
      if (TYPEOF(item_counts_sexp) == INTSXP) {
        int cell = INTEGER(item_counts_sexp)[r_index];
        if (cell == NA_INTEGER) throw std::runtime_error("`item_counts` must not contain NA.");
        value = static_cast<double>(cell);
      } else if (TYPEOF(item_counts_sexp) == REALSXP) {
        value = REAL(item_counts_sexp)[r_index];
      } else {
        throw std::runtime_error("`item_counts` must be numeric.");
      }
      if (!R_FINITE(value) || value < 0) {
        throw std::runtime_error("`item_counts` must contain non-negative finite values.");
      }
      input.item_counts[gamma_index(item_index, score, input.max_category_count)] =
        static_cast<long double>(value);
    }
  }

  input.score_counts.assign(static_cast<std::size_t>(input.max_total_score) + 1U, 0.0L);
  for (int score = 0; score <= input.max_total_score; ++score) {
    double value;
    if (TYPEOF(score_counts_sexp) == INTSXP) {
      int cell = INTEGER(score_counts_sexp)[score];
      if (cell == NA_INTEGER) throw std::runtime_error("`score_counts` must not contain NA.");
      value = static_cast<double>(cell);
    } else if (TYPEOF(score_counts_sexp) == REALSXP) {
      value = REAL(score_counts_sexp)[score];
    } else {
      throw std::runtime_error("`score_counts` must be numeric.");
    }
    if (!R_FINITE(value) || value < 0) {
      throw std::runtime_error("`score_counts` must contain non-negative finite values.");
    }
    input.score_counts[score] = static_cast<long double>(value);
  }
  return input;
}

std::vector<long double> initial_extended_gamma(const ExtendedFitInput &input) {
  std::vector<long double> gamma(static_cast<std::size_t>(input.n_items) * input.max_category_count, 0.0L);
  for (int item_index = 0; item_index < input.n_items; ++item_index) {
    for (int score = 0; score < input.raw_max[item_index]; ++score) {
      if (input.item_counts[gamma_index(item_index, score, input.max_category_count)] > 0.0L) {
        gamma[gamma_index(item_index, score, input.max_category_count)] = 1.0L;
      }
    }
  }
  return gamma;
}

void item_score_bounds(const ExtendedFitInput &input,
                       const std::vector<long double> &item_gamma,
                       std::vector<int> &ifra,
                       std::vector<int> &itil,
                       std::vector<bool> &valid_item) {
  ifra.assign(input.n_items, 0);
  itil.assign(input.n_items, 0);
  valid_item.assign(input.n_items, false);
  for (int item_index = 0; item_index < input.n_items; ++item_index) {
    int first = input.raw_max[item_index];
    int last = 0;
    for (int score = 0; score < input.raw_max[item_index]; ++score) {
      long double value = item_gamma[gamma_index(item_index, score, input.max_category_count)];
      if (value > 0.0L) {
        if (score < first) first = score;
        if (score > last) last = score;
      }
    }
    ifra[item_index] = first;
    itil[item_index] = last;
    valid_item[item_index] = first < last;
  }
}

std::vector<long double> build_source_score_gamma_extended(const ExtendedFitInput &input,
                                                          const std::vector<long double> &item_gamma,
                                                          const std::vector<bool> &use_items) {
  std::vector<int> ifra;
  std::vector<int> itil;
  std::vector<bool> valid_item;
  item_score_bounds(input, item_gamma, ifra, itil, valid_item);

  int sfra = 0;
  int stil = 0;
  for (int item_index = 0; item_index < input.n_items; ++item_index) {
    if (use_items[item_index] && ifra[item_index] <= itil[item_index]) {
      sfra += ifra[item_index];
      stil += itil[item_index];
    }
  }

  std::vector<long double> gamma_values(static_cast<std::size_t>(input.max_total_score) + 1U, 0.0L);
  std::vector<long double> next_values(static_cast<std::size_t>(input.max_total_score) + 1U, 0.0L);
  if (sfra <= input.max_total_score) {
    gamma_values[sfra] = 1.0L;
    next_values[sfra] = 1.0L;
  }

  for (int item_index = 0; item_index < input.n_items; ++item_index) {
    if (!use_items[item_index] || ifra[item_index] > itil[item_index]) {
      continue;
    }
    for (int item_score = ifra[item_index]; item_score <= itil[item_index]; ++item_score) {
      int offset = item_score - ifra[item_index];
      for (int score = sfra; score <= stil - 1; ++score) {
        int target_score = score + offset;
        if (target_score <= input.max_total_score) {
          next_values[target_score] +=
            item_gamma[gamma_index(item_index, item_score, input.max_category_count)] *
            gamma_values[score];
        }
      }
    }
    for (int score = sfra; score <= stil; ++score) {
      if (score <= input.max_total_score) {
        next_values[score] -= gamma_values[score];
      }
    }
    gamma_values = next_values;
  }
  return gamma_values;
}

std::vector<long double> calculate_expected_extended(const ExtendedFitInput &input,
                                                     const std::vector<long double> &item_gamma) {
  std::vector<long double> expected(static_cast<std::size_t>(input.n_items) * input.max_category_count, 0.0L);
  std::vector<bool> use_items(input.n_items, true);
  std::vector<long double> full_gamma =
    build_source_score_gamma_extended(input, item_gamma, use_items);

  for (int item_index = 0; item_index < input.n_items; ++item_index) {
    use_items[item_index] = false;
    std::vector<long double> without_item =
      build_source_score_gamma_extended(input, item_gamma, use_items);
    use_items[item_index] = true;
    for (int item_score = 0; item_score < input.raw_max[item_index]; ++item_score) {
      long double item_parameter =
        item_gamma[gamma_index(item_index, item_score, input.max_category_count)];
      if (item_parameter <= 0.0L) continue;
      for (int score = item_score; score <= input.max_total_score; ++score) {
        long double score_count = input.score_counts[score];
        if (score_count == 0.0L) continue;
        long double denominator = full_gamma[score];
        if (denominator <= 0.0L) continue;
        expected[gamma_index(item_index, item_score, input.max_category_count)] +=
          item_parameter * without_item[score - item_score] * (score_count / denominator);
      }
    }
  }
  return expected;
}

long double update_extended(const ExtendedFitInput &input,
                            const std::vector<long double> &expected,
                            std::vector<long double> &item_gamma,
                            bool apply_update) {
  long double delta = 0.0L;
  for (int item_index = 0; item_index < input.n_items; ++item_index) {
    for (int item_score = 0; item_score < input.raw_max[item_index]; ++item_score) {
      int index = gamma_index(item_index, item_score, input.max_category_count);
      long double observed = input.item_counts[index];
      long double fitted = expected[index];
      long double ratio = 0.0L;
      if (observed > 0.0L) {
        ratio = fitted > 0.0L ? observed / fitted : 1.0L;
        long double cell_delta = std::fabs(fitted - observed);
        if (cell_delta > delta) delta = cell_delta;
      }
      if (apply_update) {
        item_gamma[index] *= ratio;
      }
    }
  }
  return delta;
}

void adjust_extended_gamma(const ExtendedFitInput &input, std::vector<long double> &item_gamma) {
  std::vector<int> ifra;
  std::vector<int> itil;
  std::vector<bool> valid_item;
  item_score_bounds(input, item_gamma, ifra, itil, valid_item);

  long double last_sgamma = 1.0L;
  int s_max = 0;
  int s_min = 0;
  for (int item_index = 0; item_index < input.n_items; ++item_index) {
    if (!valid_item[item_index]) continue;
    int fra = ifra[item_index];
    int til = itil[item_index];
    s_max += til;
    s_min += fra;
    long double alpha = item_gamma[gamma_index(item_index, fra, input.max_category_count)];
    if (alpha > 0.0L) {
      for (int score = fra; score <= til; ++score) {
        item_gamma[gamma_index(item_index, score, input.max_category_count)] /= alpha;
      }
    }
    long double top = item_gamma[gamma_index(item_index, til, input.max_category_count)];
    if (top > 0.0L) {
      last_sgamma *= top;
    }
  }

  long double alpha = 0.0L;
  if (s_max - s_min > 0) {
    alpha = -source_ln(last_sgamma) / static_cast<long double>(s_max - s_min);
  }

  for (int item_index = 0; item_index < input.n_items; ++item_index) {
    if (valid_item[item_index]) {
      int fra = ifra[item_index];
      int til = itil[item_index];
      for (int score = fra; score <= til; ++score) {
        item_gamma[gamma_index(item_index, score, input.max_category_count)] =
          std::exp(static_cast<long double>(score - fra) * alpha) *
          item_gamma[gamma_index(item_index, score, input.max_category_count)];
      }
    } else if (ifra[item_index] <= itil[item_index]) {
      item_gamma[gamma_index(item_index, itil[item_index], input.max_category_count)] = 1.0L;
    }
  }
}

std::vector<long double> fit_extended_gamma(const ExtendedFitInput &input,
                                            int max_step,
                                            long double max_delta) {
  std::vector<long double> item_gamma = initial_extended_gamma(input);
  long double delta = 0.0L;
  for (int n_step = 0; n_step < max_step; ++n_step) {
    std::vector<long double> expected = calculate_expected_extended(input, item_gamma);
    delta = update_extended(input, expected, item_gamma, false);
    update_extended(input, expected, item_gamma, true);
    adjust_extended_gamma(input, item_gamma);
    if (delta < max_delta) break;
  }
  return item_gamma;
}

SEXP extended_gamma_matrix(const ExtendedFitInput &input,
                           const std::vector<long double> &item_gamma,
                           SEXP item_counts) {
  SEXP out = PROTECT(Rf_allocMatrix(
    REALSXP,
    input.n_items,
    input.max_category_count
  ));
  for (int item_index = 0; item_index < input.n_items; ++item_index) {
    for (int score = 0; score < input.max_category_count; ++score) {
      R_xlen_t r_index = item_index + static_cast<R_xlen_t>(input.n_items) * score;
      REAL(out)[r_index] = static_cast<double>(
        item_gamma[gamma_index(item_index, score, input.max_category_count)]
      );
    }
  }

  SEXP dim_names = Rf_getAttrib(item_counts, R_DimNamesSymbol);
  if (Rf_length(dim_names) == 2) {
    Rf_setAttrib(out, R_DimNamesSymbol, dim_names);
  }
  UNPROTECT(1);
  return out;
}

std::vector<long double> build_gamma_from_matrix(const std::vector<long double> &item_gamma,
                                                 const std::vector<int> &raw_max,
                                                 int n_items,
                                                 int max_category_count,
                                                 int max_total_score,
                                                 int excluded_item) {
  std::vector<long double> gamma(static_cast<std::size_t>(max_total_score) + 1U, 0.0L);
  std::vector<long double> next(static_cast<std::size_t>(max_total_score) + 1U, 0.0L);
  gamma[0] = 1.0L;
  int current_max = 0;
  for (int item_index = 0; item_index < n_items; ++item_index) {
    if (item_index == excluded_item) continue;
    std::fill(next.begin(), next.end(), 0.0L);
    int next_max = current_max + raw_max[item_index] - 1;
    for (int score = 0; score <= current_max; ++score) {
      long double current = gamma[score];
      if (current == 0.0L) continue;
      for (int item_score = 0; item_score < raw_max[item_index]; ++item_score) {
        next[score + item_score] +=
          current * item_gamma[gamma_index(item_index, item_score, max_category_count)];
      }
    }
    gamma.swap(next);
    current_max = next_max;
  }
  return gamma;
}

struct HomogeneityTalSummary {
  long double n = 0.0L;
  long double mean = 0.0L;
  long double variance = 0.0L;
};

HomogeneityTalSummary summarize_homogeneity_tal(const std::vector<long double> &cells,
                                                int item_max) {
  HomogeneityTalSummary out;
  long double second = 0.0L;
  for (int item_score = 0; item_score <= item_max; ++item_score) {
    long double cell = cells[item_score];
    out.n += cell;
    out.mean += static_cast<long double>(item_score) * cell;
    second += static_cast<long double>(item_score * item_score) * cell;
  }
  if (out.n > 0.0L) {
    out.mean /= out.n;
    if (out.n > 1.0L) second /= out.n;
  }
  out.variance = out.n > 1.0L ? second - out.mean * out.mean : 0.0L;
  return out;
}

struct HomogeneityItemSummary {
  long double expected_n = 0.0L;
  long double expected_mean = 0.0L;
  long double expected_variance = 0.0L;
};

HomogeneityItemSummary global_homogeneity_item_summary_cpp(
    int item_index,
    int n_items,
    int max_category_count,
    int max_total_score,
    int source_item_max,
    int scfra,
    int sctil,
    const std::vector<int> &raw_max,
    const std::vector<long double> &item_gamma,
    const std::vector<long double> &score_counts,
    const std::vector<long double> &score_item_n,
    const std::vector<long double> &full_gamma,
    const std::vector<long double> &without_item) {
  int item_max = raw_max[item_index] - 1;
  std::vector<long double> expected_by_score(
    static_cast<std::size_t>(max_total_score + 1) * max_category_count,
    0.0L
  );
  std::vector<long double> expected_total(
    static_cast<std::size_t>(max_category_count),
    0.0L
  );

  for (int item_score = 0; item_score <= item_max; ++item_score) {
    long double item_weight =
      item_gamma[gamma_index(item_index, item_score, max_category_count)];
    if (item_weight <= 0.0L) continue;

    int first_score = item_score < scfra ? scfra : item_score;
    int last_score = item_score + (n_items - 1) * source_item_max;
    if (last_score > sctil) last_score = sctil;

    for (int score = first_score; score <= last_score; ++score) {
      int rest_score = score - item_score;
      if (rest_score < 0 || rest_score > max_total_score) continue;
      long double denominator = full_gamma[score];
      if (denominator <= 0.0L) continue;

      long double expected =
        item_weight * without_item[rest_score] * (score_counts[score] / denominator);
      expected_by_score[
        score + static_cast<std::size_t>(max_total_score + 1) * item_score
      ] += expected;
      if (score > 0 && score < max_total_score) expected_total[item_score] += expected;
    }
  }

  std::vector<long double> score_variance(
    static_cast<std::size_t>(max_total_score) + 1U,
    0.0L
  );
  for (int score = scfra; score <= sctil; ++score) {
    std::vector<long double> score_cells(static_cast<std::size_t>(item_max) + 1U, 0.0L);
    for (int item_score = 0; item_score <= item_max; ++item_score) {
      score_cells[item_score] = expected_by_score[
        score + static_cast<std::size_t>(max_total_score + 1) * item_score
      ];
    }
    score_variance[score] = summarize_homogeneity_tal(score_cells, item_max).variance;
  }

  HomogeneityTalSummary total = summarize_homogeneity_tal(expected_total, item_max);
  HomogeneityItemSummary out;
  out.expected_n = total.n;
  out.expected_mean = total.mean;

  double total_n = static_cast<double>(total.n);
  if (total_n > 0.0) {
    for (int score = scfra; score <= sctil; ++score) {
      double observed_n = static_cast<double>(
        score_item_n[score + static_cast<std::size_t>(max_total_score + 1) * item_index]
      );
      out.expected_variance +=
        static_cast<long double>(observed_n / total_n) * score_variance[score];
    }
  }

  return out;
}

}  // namespace

extern "C" SEXP gRm_item_parameters_top_ice_field(SEXP gamma_top, SEXP max_score) {
  try {
    X87ExtendedPrecisionScope precision_scope;
    double gamma_top_double = as_single_double(gamma_top, "gamma_top");
    int k = as_single_int(max_score, "max_score");
    if (k <= 0) {
      throw std::runtime_error("`max_score` must be positive.");
    }
    if (gamma_top_double <= 0) {
      throw std::runtime_error("`gamma_top` must be positive.");
    }

    // Source trace: source/GLLRM.txt output(4), ICE branch:
    // z := ln(x[i,k]) / k; write ln(x[i,k]) - k*z : 10:3.
    long double gamma_top_extended = static_cast<long double>(gamma_top_double);
    volatile long double z = source_ln(gamma_top_extended) / static_cast<long double>(k);
    long double top_ice = source_ln(gamma_top_extended) - static_cast<long double>(k) * z;

    std::string field;
    if (std::isnan(top_ice) || std::fabs(top_ice) >= 0.0005L) {
      field = field9_3_long_double(top_ice);
    } else if (std::signbit(top_ice)) {
      field = "   -0.000";
    } else {
      field = "    0.000";
    }

    return Rf_mkString(field.c_str());
  } catch (const std::exception &err) {
    Rf_error("%s", err.what());
  }
}

extern "C" SEXP gRm_item_parameters_extended_top_ice_fields(SEXP item_counts,
                                                                 SEXP score_counts,
                                                                 SEXP raw_max,
                                                                 SEXP max_step,
                                                                 SEXP max_delta) {
  try {
    X87ExtendedPrecisionScope precision_scope;
    ExtendedFitInput input = parse_extended_fit_input(item_counts, score_counts, raw_max);
    int max_step_value = as_single_int(max_step, "max_step");
    double max_delta_value = as_single_double(max_delta, "max_delta");
    if (max_step_value < 0) {
      throw std::runtime_error("`max_step` must not be negative.");
    }
    if (max_delta_value < 0) {
      throw std::runtime_error("`max_delta` must not be negative.");
    }

    std::vector<long double> item_gamma =
      fit_extended_gamma(input, max_step_value, static_cast<long double>(max_delta_value));

    SEXP out = PROTECT(Rf_allocVector(STRSXP, input.n_items));
    for (int item_index = 0; item_index < input.n_items; ++item_index) {
      int max_score = input.raw_max[item_index] - 1;
      std::string field = "    0.000";
      if (max_score > 0) {
        long double gamma_top =
          item_gamma[gamma_index(item_index, max_score, input.max_category_count)];
        if (gamma_top > 0.0L) {
          field = top_ice_field_from_gamma(gamma_top, max_score);
        }
      }
      SET_STRING_ELT(out, item_index, Rf_mkChar(field.c_str()));
    }

    SEXP row_names = Rf_getAttrib(item_counts, R_DimNamesSymbol);
    if (Rf_length(row_names) == 2) {
      SEXP item_names = VECTOR_ELT(row_names, 0);
      if (Rf_length(item_names) == input.n_items) {
        Rf_setAttrib(out, R_NamesSymbol, item_names);
      }
    }
    UNPROTECT(1);
    return out;
  } catch (const std::exception &err) {
    Rf_error("%s", err.what());
  }
}

extern "C" SEXP gRm_item_parameters_extended_gamma(SEXP item_counts,
                                                       SEXP score_counts,
                                                       SEXP raw_max,
                                                       SEXP max_step,
                                                       SEXP max_delta) {
  try {
    X87ExtendedPrecisionScope precision_scope;
    ExtendedFitInput input = parse_extended_fit_input(item_counts, score_counts, raw_max);
    int max_step_value = as_single_int(max_step, "max_step");
    double max_delta_value = as_single_double(max_delta, "max_delta");
    if (max_step_value < 0) {
      throw std::runtime_error("`max_step` must not be negative.");
    }
    if (max_delta_value < 0) {
      throw std::runtime_error("`max_delta` must not be negative.");
    }

    std::vector<long double> item_gamma =
      fit_extended_gamma(input, max_step_value, static_cast<long double>(max_delta_value));
    return extended_gamma_matrix(input, item_gamma, item_counts);
  } catch (const std::exception &err) {
    Rf_error("%s", err.what());
  }
}

extern "C" SEXP gRm_item_parameters_extended_ice_fields(SEXP item_counts,
                                                            SEXP score_counts,
                                                            SEXP raw_max,
                                                            SEXP max_step,
                                                            SEXP max_delta) {
  try {
    X87ExtendedPrecisionScope precision_scope;
    ExtendedFitInput input = parse_extended_fit_input(item_counts, score_counts, raw_max);
    int max_step_value = as_single_int(max_step, "max_step");
    double max_delta_value = as_single_double(max_delta, "max_delta");
    if (max_step_value < 0) {
      throw std::runtime_error("`max_step` must not be negative.");
    }
    if (max_delta_value < 0) {
      throw std::runtime_error("`max_delta` must not be negative.");
    }

    std::vector<long double> item_gamma =
      fit_extended_gamma(input, max_step_value, static_cast<long double>(max_delta_value));

    SEXP out = PROTECT(Rf_allocMatrix(
      STRSXP,
      input.n_items,
      input.max_category_count
    ));
    for (int item_index = 0; item_index < input.n_items; ++item_index) {
      int max_score = 0;
      for (int score = 0; score < input.raw_max[item_index]; ++score) {
        if (item_gamma[gamma_index(item_index, score, input.max_category_count)] > 0.0L) {
          max_score = score;
        }
      }
      std::vector<std::string> fields =
        skbias12_ice_fields_from_gamma(item_gamma, item_index, max_score, input.max_category_count);
      for (int score = 0; score < input.max_category_count; ++score) {
        R_xlen_t r_index = item_index + static_cast<R_xlen_t>(input.n_items) * score;
        SET_STRING_ELT(out, r_index, Rf_mkChar(fields[score].c_str()));
      }
    }

    SEXP dim_names = Rf_getAttrib(item_counts, R_DimNamesSymbol);
    if (Rf_length(dim_names) == 2) {
      Rf_setAttrib(out, R_DimNamesSymbol, dim_names);
    }
    UNPROTECT(1);
    return out;
  } catch (const std::exception &err) {
    Rf_error("%s", err.what());
  }
}

extern "C" SEXP gRm_item_parameters_ice_fields_from_gamma(SEXP item_gamma_sexp,
                                                              SEXP raw_max_sexp) {
  try {
    X87ExtendedPrecisionScope precision_scope;
    SEXP dim = Rf_getAttrib(item_gamma_sexp, R_DimSymbol);
    if (Rf_length(dim) != 2) {
      throw std::runtime_error("`item_gamma` must be a matrix.");
    }
    int n_items = INTEGER(dim)[0];
    int max_category_count = INTEGER(dim)[1];
    if (Rf_length(raw_max_sexp) != n_items) {
      throw std::runtime_error("`raw_max` must have one entry per item.");
    }

    std::vector<int> raw_max(n_items);
    for (int item_index = 0; item_index < n_items; ++item_index) {
      int value;
      if (TYPEOF(raw_max_sexp) == INTSXP) {
        value = INTEGER(raw_max_sexp)[item_index];
        if (value == NA_INTEGER) throw std::runtime_error("`raw_max` must not contain NA.");
      } else if (TYPEOF(raw_max_sexp) == REALSXP) {
        double numeric = REAL(raw_max_sexp)[item_index];
        if (!R_FINITE(numeric)) throw std::runtime_error("`raw_max` must be finite.");
        value = static_cast<int>(numeric);
      } else {
        throw std::runtime_error("`raw_max` must be numeric.");
      }
      if (value < 1 || value > max_category_count) {
        throw std::runtime_error("`raw_max` entries must be valid category counts.");
      }
      raw_max[item_index] = value;
    }

    std::vector<long double> item_gamma(static_cast<std::size_t>(n_items) * max_category_count, 0.0L);
    for (int item_index = 0; item_index < n_items; ++item_index) {
      for (int score = 0; score < max_category_count; ++score) {
        R_xlen_t r_index = item_index + static_cast<R_xlen_t>(n_items) * score;
        double value = REAL(item_gamma_sexp)[r_index];
        if (!R_FINITE(value)) {
          throw std::runtime_error("`item_gamma` must contain finite values.");
        }
        item_gamma[gamma_index(item_index, score, max_category_count)] =
          static_cast<long double>(value);
      }
    }

    SEXP out = PROTECT(Rf_allocMatrix(STRSXP, n_items, max_category_count));
    for (int item_index = 0; item_index < n_items; ++item_index) {
      int max_score = 0;
      for (int score = 0; score < raw_max[item_index]; ++score) {
        if (item_gamma[gamma_index(item_index, score, max_category_count)] > 0.0L) {
          max_score = score;
        }
      }
      std::vector<std::string> fields =
        skbias12_ice_fields_from_gamma(item_gamma, item_index, max_score, max_category_count);
      for (int score = 0; score < max_category_count; ++score) {
        R_xlen_t r_index = item_index + static_cast<R_xlen_t>(n_items) * score;
        SET_STRING_ELT(out, r_index, Rf_mkChar(fields[score].c_str()));
      }
    }

    SEXP dim_names = Rf_getAttrib(item_gamma_sexp, R_DimNamesSymbol);
    if (Rf_length(dim_names) == 2) {
      Rf_setAttrib(out, R_DimNamesSymbol, dim_names);
    }
    UNPROTECT(1);
    return out;
  } catch (const std::exception &err) {
    Rf_error("%s", err.what());
  }
}

extern "C" SEXP gRm_global_homogeneity_expected_summary(SEXP item_gamma_sexp,
                                                            SEXP score_counts_sexp,
                                                            SEXP score_item_n_sexp,
                                                            SEXP raw_max_sexp,
                                                            SEXP least_score_sexp,
                                                            SEXP largest_score_sexp) {
  try {
    X87ExtendedPrecisionScope precision_scope;
    SEXP dim = Rf_getAttrib(item_gamma_sexp, R_DimSymbol);
    if (Rf_length(dim) != 2) {
      throw std::runtime_error("`item_gamma` must be a matrix.");
    }
    int n_items = INTEGER(dim)[0];
    int max_category_count = INTEGER(dim)[1];
    int max_total_score = Rf_length(score_counts_sexp) - 1;
    if (max_total_score < 0 || Rf_length(raw_max_sexp) != n_items) {
      throw std::runtime_error("Invalid global homogeneity dimensions.");
    }

    SEXP score_item_dim = Rf_getAttrib(score_item_n_sexp, R_DimSymbol);
    if (Rf_length(score_item_dim) != 2 ||
        INTEGER(score_item_dim)[0] != max_total_score + 1 ||
        INTEGER(score_item_dim)[1] != n_items) {
      throw std::runtime_error("`score_item_n` must be a score-by-item matrix.");
    }

    std::vector<int> raw_max(n_items);
    int source_item_max = 0;
    for (int item_index = 0; item_index < n_items; ++item_index) {
      int value;
      if (TYPEOF(raw_max_sexp) == INTSXP) {
        value = INTEGER(raw_max_sexp)[item_index];
        if (value == NA_INTEGER) throw std::runtime_error("`raw_max` must not contain NA.");
      } else {
        double numeric = REAL(raw_max_sexp)[item_index];
        if (!R_FINITE(numeric)) throw std::runtime_error("`raw_max` must be finite.");
        value = static_cast<int>(numeric);
      }
      if (value < 1 || value > max_category_count) {
        throw std::runtime_error("`raw_max` entries must be valid category counts.");
      }
      raw_max[item_index] = value;
      if (value - 1 > source_item_max) source_item_max = value - 1;
    }

    std::vector<long double> item_gamma(
      static_cast<std::size_t>(n_items) * max_category_count,
      0.0L
    );
    for (int item_index = 0; item_index < n_items; ++item_index) {
      for (int item_score = 0; item_score < max_category_count; ++item_score) {
        R_xlen_t r_index = item_index + static_cast<R_xlen_t>(n_items) * item_score;
        double value = REAL(item_gamma_sexp)[r_index];
        if (!R_FINITE(value)) throw std::runtime_error("`item_gamma` must be finite.");
        item_gamma[gamma_index(item_index, item_score, max_category_count)] =
          static_cast<long double>(value);
      }
    }

    std::vector<long double> score_counts(
      static_cast<std::size_t>(max_total_score) + 1U,
      0.0L
    );
    for (int score = 0; score <= max_total_score; ++score) {
      double value = TYPEOF(score_counts_sexp) == INTSXP
        ? static_cast<double>(INTEGER(score_counts_sexp)[score])
        : REAL(score_counts_sexp)[score];
      if (!R_FINITE(value) || value < 0.0) {
        throw std::runtime_error("`score_counts` must be non-negative.");
      }
      score_counts[score] = static_cast<long double>(value);
    }

    std::vector<long double> score_item_n(
      static_cast<std::size_t>(max_total_score + 1) * n_items,
      0.0L
    );
    for (int score = 0; score <= max_total_score; ++score) {
      for (int item_index = 0; item_index < n_items; ++item_index) {
        R_xlen_t r_index = score + static_cast<R_xlen_t>(max_total_score + 1) * item_index;
        double value = TYPEOF(score_item_n_sexp) == INTSXP
          ? static_cast<double>(INTEGER(score_item_n_sexp)[r_index])
          : REAL(score_item_n_sexp)[r_index];
        if (!R_FINITE(value) || value < 0.0) {
          throw std::runtime_error("`score_item_n` must be non-negative.");
        }
        score_item_n[score + static_cast<std::size_t>(max_total_score + 1) * item_index] =
          static_cast<long double>(value);
      }
    }

    int least_score = as_single_int(least_score_sexp, "least_score");
    int largest_score = as_single_int(largest_score_sexp, "largest_score");
    int scfra = least_score > 1 ? least_score : 1;
    int sctil = largest_score < max_total_score - 1 ? largest_score : max_total_score - 1;

    ExtendedFitInput gamma_input;
    gamma_input.n_items = n_items;
    gamma_input.max_category_count = max_category_count;
    gamma_input.max_total_score = max_total_score;
    gamma_input.raw_max = raw_max;
    std::vector<bool> use_items(n_items, true);
    std::vector<long double> full_gamma =
      build_source_score_gamma_extended(gamma_input, item_gamma, use_items);

    SEXP expected_variance = PROTECT(Rf_allocVector(REALSXP, n_items));
    SEXP expected_n = PROTECT(Rf_allocVector(REALSXP, n_items));
    SEXP expected_mean = PROTECT(Rf_allocVector(REALSXP, n_items));

    for (int item_index = 0; item_index < n_items; ++item_index) {
      use_items[item_index] = false;
      std::vector<long double> without_item =
        build_source_score_gamma_extended(gamma_input, item_gamma, use_items);
      use_items[item_index] = true;

      HomogeneityItemSummary summary = global_homogeneity_item_summary_cpp(
        item_index,
        n_items,
        max_category_count,
        max_total_score,
        source_item_max,
        scfra,
        sctil,
        raw_max,
        item_gamma,
        score_counts,
        score_item_n,
        full_gamma,
        without_item
      );

      REAL(expected_variance)[item_index] = static_cast<double>(summary.expected_variance);
      REAL(expected_n)[item_index] = static_cast<double>(summary.expected_n);
      REAL(expected_mean)[item_index] = static_cast<double>(summary.expected_mean);
    }

    SEXP names = Rf_getAttrib(item_gamma_sexp, R_DimNamesSymbol);
    if (Rf_length(names) == 2) {
      SEXP row_names = VECTOR_ELT(names, 0);
      if (Rf_length(row_names) == n_items) {
        Rf_setAttrib(expected_variance, R_NamesSymbol, row_names);
        Rf_setAttrib(expected_n, R_NamesSymbol, row_names);
        Rf_setAttrib(expected_mean, R_NamesSymbol, row_names);
      }
    }

    SEXP out = PROTECT(Rf_allocVector(VECSXP, 3));
    SET_VECTOR_ELT(out, 0, expected_variance);
    SET_VECTOR_ELT(out, 1, expected_n);
    SET_VECTOR_ELT(out, 2, expected_mean);
    SEXP out_names = PROTECT(Rf_allocVector(STRSXP, 3));
    SET_STRING_ELT(out_names, 0, Rf_mkChar("expected_variance"));
    SET_STRING_ELT(out_names, 1, Rf_mkChar("expected_n"));
    SET_STRING_ELT(out_names, 2, Rf_mkChar("expected_mean"));
    Rf_setAttrib(out, R_NamesSymbol, out_names);
    UNPROTECT(5);
    return out;
  } catch (const std::exception &err) {
    Rf_error("%s", err.what());
  }
}

extern "C" SEXP gRm_global_homogeneity_expected_summary_from_fit(SEXP full_item_counts_sexp,
                                                                     SEXP full_score_counts_sexp,
                                                                     SEXP group_score_counts_sexp,
                                                                     SEXP score_item_n_sexp,
                                                                     SEXP raw_max_sexp,
                                                                     SEXP least_score_sexp,
                                                                     SEXP largest_score_sexp,
                                                                     SEXP max_step_sexp,
                                                                     SEXP max_delta_sexp) {
  try {
    X87ExtendedPrecisionScope precision_scope;
    ExtendedFitInput fit_input = parse_extended_fit_input(
      full_item_counts_sexp,
      full_score_counts_sexp,
      raw_max_sexp
    );
    int max_step_value = as_single_int(max_step_sexp, "max_step");
    double max_delta_value = as_single_double(max_delta_sexp, "max_delta");
    if (max_step_value < 0) {
      throw std::runtime_error("`max_step` must not be negative.");
    }
    if (max_delta_value < 0) {
      throw std::runtime_error("`max_delta` must not be negative.");
    }

    int n_items = fit_input.n_items;
    int max_category_count = fit_input.max_category_count;
    int max_total_score = fit_input.max_total_score;
    int source_item_max = 0;
    for (int item_index = 0; item_index < n_items; ++item_index) {
      if (fit_input.raw_max[item_index] - 1 > source_item_max) {
        source_item_max = fit_input.raw_max[item_index] - 1;
      }
    }

    if (Rf_length(group_score_counts_sexp) != max_total_score + 1) {
      throw std::runtime_error("`group_score_counts` has the wrong length.");
    }
    std::vector<long double> group_score_counts(
      static_cast<std::size_t>(max_total_score) + 1U,
      0.0L
    );
    for (int score = 0; score <= max_total_score; ++score) {
      double value = TYPEOF(group_score_counts_sexp) == INTSXP
        ? static_cast<double>(INTEGER(group_score_counts_sexp)[score])
        : REAL(group_score_counts_sexp)[score];
      if (!R_FINITE(value) || value < 0.0) {
        throw std::runtime_error("`group_score_counts` must be non-negative.");
      }
      group_score_counts[score] = static_cast<long double>(value);
    }

    SEXP score_item_dim = Rf_getAttrib(score_item_n_sexp, R_DimSymbol);
    if (Rf_length(score_item_dim) != 2 ||
        INTEGER(score_item_dim)[0] != max_total_score + 1 ||
        INTEGER(score_item_dim)[1] != n_items) {
      throw std::runtime_error("`score_item_n` must be a score-by-item matrix.");
    }
    std::vector<long double> score_item_n(
      static_cast<std::size_t>(max_total_score + 1) * n_items,
      0.0L
    );
    for (int score = 0; score <= max_total_score; ++score) {
      for (int item_index = 0; item_index < n_items; ++item_index) {
        R_xlen_t r_index = score + static_cast<R_xlen_t>(max_total_score + 1) * item_index;
        double value = TYPEOF(score_item_n_sexp) == INTSXP
          ? static_cast<double>(INTEGER(score_item_n_sexp)[r_index])
          : REAL(score_item_n_sexp)[r_index];
        if (!R_FINITE(value) || value < 0.0) {
          throw std::runtime_error("`score_item_n` must be non-negative.");
        }
        score_item_n[score + static_cast<std::size_t>(max_total_score + 1) * item_index] =
          static_cast<long double>(value);
      }
    }

    std::vector<long double> item_gamma = fit_extended_gamma(
      fit_input,
      max_step_value,
      static_cast<long double>(max_delta_value)
    );

    int least_score = as_single_int(least_score_sexp, "least_score");
    int largest_score = as_single_int(largest_score_sexp, "largest_score");
    int scfra = least_score > 1 ? least_score : 1;
    int sctil = largest_score < max_total_score - 1 ? largest_score : max_total_score - 1;

    std::vector<bool> use_items(n_items, true);
    std::vector<long double> full_gamma =
      build_source_score_gamma_extended(fit_input, item_gamma, use_items);

    SEXP expected_variance = PROTECT(Rf_allocVector(REALSXP, n_items));
    SEXP expected_n = PROTECT(Rf_allocVector(REALSXP, n_items));
    SEXP expected_mean = PROTECT(Rf_allocVector(REALSXP, n_items));

    for (int item_index = 0; item_index < n_items; ++item_index) {
      use_items[item_index] = false;
      std::vector<long double> without_item =
        build_source_score_gamma_extended(fit_input, item_gamma, use_items);
      use_items[item_index] = true;

      HomogeneityItemSummary summary = global_homogeneity_item_summary_cpp(
        item_index,
        n_items,
        max_category_count,
        max_total_score,
        source_item_max,
        scfra,
        sctil,
        fit_input.raw_max,
        item_gamma,
        group_score_counts,
        score_item_n,
        full_gamma,
        without_item
      );

      REAL(expected_variance)[item_index] = static_cast<double>(summary.expected_variance);
      REAL(expected_n)[item_index] = static_cast<double>(summary.expected_n);
      REAL(expected_mean)[item_index] = static_cast<double>(summary.expected_mean);
    }

    SEXP names = Rf_getAttrib(full_item_counts_sexp, R_DimNamesSymbol);
    if (Rf_length(names) == 2) {
      SEXP row_names = VECTOR_ELT(names, 0);
      if (Rf_length(row_names) == n_items) {
        Rf_setAttrib(expected_variance, R_NamesSymbol, row_names);
        Rf_setAttrib(expected_n, R_NamesSymbol, row_names);
        Rf_setAttrib(expected_mean, R_NamesSymbol, row_names);
      }
    }

    SEXP out = PROTECT(Rf_allocVector(VECSXP, 3));
    SET_VECTOR_ELT(out, 0, expected_variance);
    SET_VECTOR_ELT(out, 1, expected_n);
    SET_VECTOR_ELT(out, 2, expected_mean);
    SEXP out_names = PROTECT(Rf_allocVector(STRSXP, 3));
    SET_STRING_ELT(out_names, 0, Rf_mkChar("expected_variance"));
    SET_STRING_ELT(out_names, 1, Rf_mkChar("expected_n"));
    SET_STRING_ELT(out_names, 2, Rf_mkChar("expected_mean"));
    Rf_setAttrib(out, R_NamesSymbol, out_names);
    UNPROTECT(5);
    return out;
  } catch (const std::exception &err) {
    Rf_error("%s", err.what());
  }
}
