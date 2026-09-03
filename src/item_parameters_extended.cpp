#include <cfloat>
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <vector>

#include <R.h>
#include <Rinternals.h>

namespace {

#if defined(__GNUC__) && (defined(__i386__) || defined(__x86_64__))
constexpr bool kX87ExtendedAvailable = true;

class X87ExtendedPrecisionScope {
 public:
  X87ExtendedPrecisionScope() {
    __asm__ __volatile__("fnstcw %0" : "=m"(original_));
    unsigned short extended =
      static_cast<unsigned short>((original_ & ~0x0300U) | 0x0300U);
    __asm__ __volatile__("fldcw %0" : : "m"(extended));
  }

  ~X87ExtendedPrecisionScope() {
    __asm__ __volatile__("fldcw %0" : : "m"(original_));
  }

 private:
  unsigned short original_;
};
#else
constexpr bool kX87ExtendedAvailable = false;

class X87ExtendedPrecisionScope {
 public:
  X87ExtendedPrecisionScope() = default;
};
#endif

constexpr bool kPascalExtendedCompatible =
  kX87ExtendedAvailable && LDBL_MANT_DIG >= 64 && LDBL_MAX_EXP >= 16384;

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

int as_single_integer_like(SEXP value, const char *name) {
  double numeric = as_single_double(value, name);
  if (numeric != std::floor(numeric) || numeric < 0.0 || numeric > 2147483647.0) {
    throw std::runtime_error(
      std::string("`") + name + "` must be a non-negative integer-like value."
    );
  }
  return static_cast<int>(numeric);
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
  // Source trace: source/digram_source_20260817/skunits/skbias22.pas::write_iteminformation1 uses Pascal
  // Extended Ln. On supported GNU x86 builds, keep both the logarithm and its
  // cancellation operands in x87 Extended precision.
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
  // Compatibility boundary: this is only platform long-double arithmetic.
  // The capability routine reports false when that type cannot represent the
  // 64-bit-significand, 15-bit-exponent Pascal Extended format.
  return std::log(value);
#endif
}

struct TopIceDiagnostic {
  long double value;
  std::string field;
};

TopIceDiagnostic top_ice_diagnostic(long double gamma_top, int k) {
  // Source trace: source/digram_source_20260817/skunits/skbias22.pas::write_iteminformation1, output(4) ICE
  // branch. Preserve exact expression order: z := Ln(x[i,k]) / k, followed by
  // Ln(x[i,k]) - k*z. Direct fixed-field formatting derives the sign of a
  // displayed zero from the arithmetic result; no sign is imposed afterward.
  volatile long double z = source_ln(gamma_top) / static_cast<long double>(k);
  long double value = source_ln(gamma_top) - static_cast<long double>(k) * z;
  return TopIceDiagnostic{value, field9_3_long_double(value)};
}

std::vector<std::string> skbias12_ice_fields_from_gamma(
    const std::vector<long double> &item_gamma,
    int item_index,
    int max_score,
    int max_category_count) {
  std::vector<std::string> fields(
    static_cast<std::size_t>(max_category_count),
    "     ...."
  );
  std::vector<long double> ice(
    static_cast<std::size_t>(max_category_count),
    0.0L
  );

  if (max_score == 1) {
    long double gamma1 =
      item_gamma[gamma_index(item_index, 1, max_category_count)];
    if (gamma1 > 0.0L) {
      // Source trace: source/digram_source_20260817/skunits/skbias12.pas::CalculateICEandMICE.
      // The dichotomous branch retains zero ICE and stores the one-step
      // quantity as ItemEffect.
      ice[0] = 0.0L;
      ice[1] = 0.0L;
    }
  } else if (max_score > 1) {
    // Source trace: source/digram_source_20260817/skunits/skbias12.pas::CalculateICEandMICE.
    // Mathematical step: center adjacent-category log ratios at their mean,
    // then derive MICE and ICE in source category order.
    long double location = 0.0L;
    for (int score = 1; score <= max_score; ++score) {
      long double current =
        item_gamma[gamma_index(item_index, score, max_category_count)];
      long double previous =
        item_gamma[gamma_index(item_index, score - 1, max_category_count)];
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
    std::vector<long double> probs(
      static_cast<std::size_t>(max_category_count),
      0.0L
    );
    for (int score = 0; score <= max_score; ++score) {
      thfactor *= theta;
      probs[score] =
        item_gamma[gamma_index(item_index, score, max_category_count)] *
        thfactor;
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

SEXP make_named_capability_list() {
  const int n_fields = 8;
  SEXP out = PROTECT(Rf_allocVector(VECSXP, n_fields));
  SEXP names = PROTECT(Rf_allocVector(STRSXP, n_fields));
  const char *field_names[n_fields] = {
    "schema", "mantissa_bits", "max_exponent", "storage_bytes",
    "x87_available", "x87_control_scope", "x87_log_instruction",
    "pascal_extended_fixed_field_supported"
  };
  for (int index = 0; index < n_fields; ++index) {
    SET_STRING_ELT(names, index, Rf_mkChar(field_names[index]));
  }
  SET_VECTOR_ELT(out, 0, Rf_mkString("gRm-pascal-extended-capability-v1"));
  SET_VECTOR_ELT(out, 1, Rf_ScalarInteger(LDBL_MANT_DIG));
  SET_VECTOR_ELT(out, 2, Rf_ScalarInteger(LDBL_MAX_EXP));
  SET_VECTOR_ELT(out, 3, Rf_ScalarInteger(static_cast<int>(sizeof(long double))));
  SET_VECTOR_ELT(out, 4, Rf_ScalarLogical(kX87ExtendedAvailable));
  SET_VECTOR_ELT(out, 5, Rf_ScalarLogical(kX87ExtendedAvailable));
  SET_VECTOR_ELT(out, 6, Rf_ScalarLogical(kX87ExtendedAvailable));
  SET_VECTOR_ELT(out, 7, Rf_ScalarLogical(kPascalExtendedCompatible));
  Rf_setAttrib(out, R_NamesSymbol, names);
  UNPROTECT(2);
  return out;
}

}  // namespace

extern "C" SEXP gRm_item_parameters_extended_capabilities() {
  return make_named_capability_list();
}

extern "C" SEXP gRm_item_parameters_top_ice_field(SEXP gamma_top,
                                                      SEXP max_score) {
  try {
    X87ExtendedPrecisionScope precision_scope;
    double gamma_top_double = as_single_double(gamma_top, "gamma_top");
    int k = as_single_integer_like(max_score, "max_score");
    if (k <= 0) {
      throw std::runtime_error("`max_score` must be positive.");
    }
    if (gamma_top_double <= 0) {
      throw std::runtime_error("`gamma_top` must be positive.");
    }

    TopIceDiagnostic diagnostic = top_ice_diagnostic(
      static_cast<long double>(gamma_top_double),
      k
    );
    SEXP out = PROTECT(Rf_mkString(diagnostic.field.c_str()));
    SEXP value = PROTECT(Rf_ScalarReal(static_cast<double>(diagnostic.value)));
    const char *sign = diagnostic.value < 0.0L
      ? "negative"
      : (diagnostic.value > 0.0L ? "positive" : "zero");
    SEXP sign_value = PROTECT(Rf_mkString(sign));
    Rf_setAttrib(out, Rf_install("cancellation_value"), value);
    Rf_setAttrib(out, Rf_install("cancellation_sign"), sign_value);
    UNPROTECT(3);
    return out;
  } catch (const std::exception &err) {
    Rf_error("%s", err.what());
  }
}

extern "C" SEXP gRm_item_parameters_ice_fields_from_gamma(
    SEXP item_gamma_sexp,
    SEXP raw_max_sexp) {
  try {
    X87ExtendedPrecisionScope precision_scope;
    if (TYPEOF(item_gamma_sexp) != REALSXP && TYPEOF(item_gamma_sexp) != INTSXP) {
      throw std::runtime_error("`item_gamma` must be a numeric matrix.");
    }
    SEXP dim = Rf_getAttrib(item_gamma_sexp, R_DimSymbol);
    if (TYPEOF(dim) != INTSXP || Rf_length(dim) != 2) {
      throw std::runtime_error("`item_gamma` must be a matrix.");
    }
    int n_items = INTEGER(dim)[0];
    int max_category_count = INTEGER(dim)[1];
    if (n_items < 1 || max_category_count < 1) {
      throw std::runtime_error("`item_gamma` must have positive dimensions.");
    }
    if (Rf_length(raw_max_sexp) != n_items) {
      throw std::runtime_error("`raw_max` must have one entry per item.");
    }

    std::vector<int> raw_max(n_items);
    for (int item_index = 0; item_index < n_items; ++item_index) {
      SEXP value = PROTECT(Rf_ScalarReal(
        TYPEOF(raw_max_sexp) == INTSXP
          ? static_cast<double>(INTEGER(raw_max_sexp)[item_index])
          : (TYPEOF(raw_max_sexp) == REALSXP
            ? REAL(raw_max_sexp)[item_index]
            : NA_REAL)
      ));
      int category_count = as_single_integer_like(value, "raw_max");
      UNPROTECT(1);
      if (category_count < 1 || category_count > max_category_count) {
        throw std::runtime_error(
          "`raw_max` entries must be valid category counts."
        );
      }
      raw_max[item_index] = category_count;
    }

    std::vector<long double> item_gamma(
      static_cast<std::size_t>(n_items) * max_category_count,
      0.0L
    );
    for (int item_index = 0; item_index < n_items; ++item_index) {
      for (int score = 0; score < max_category_count; ++score) {
        R_xlen_t r_index =
          item_index + static_cast<R_xlen_t>(n_items) * score;
        double value;
        if (TYPEOF(item_gamma_sexp) == INTSXP) {
          int cell = INTEGER(item_gamma_sexp)[r_index];
          if (cell == NA_INTEGER) {
            throw std::runtime_error("`item_gamma` must not contain NA.");
          }
          value = static_cast<double>(cell);
        } else {
          value = REAL(item_gamma_sexp)[r_index];
        }
        if (!R_FINITE(value) || value < 0.0) {
          throw std::runtime_error(
            "`item_gamma` must contain non-negative finite values."
          );
        }
        item_gamma[gamma_index(item_index, score, max_category_count)] =
          static_cast<long double>(value);
      }
    }

    SEXP out = PROTECT(Rf_allocMatrix(
      STRSXP,
      n_items,
      max_category_count
    ));
    for (int item_index = 0; item_index < n_items; ++item_index) {
      int max_score = 0;
      for (int score = 0; score < raw_max[item_index]; ++score) {
        if (item_gamma[gamma_index(
              item_index,
              score,
              max_category_count
            )] > 0.0L) {
          max_score = score;
        }
      }
      std::vector<std::string> fields = skbias12_ice_fields_from_gamma(
        item_gamma,
        item_index,
        max_score,
        max_category_count
      );
      for (int score = 0; score < max_category_count; ++score) {
        R_xlen_t r_index =
          item_index + static_cast<R_xlen_t>(n_items) * score;
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
