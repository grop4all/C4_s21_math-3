#include "./s21_math.h"

#include <stdlib.h>

int s21_abs(int x) {
  if (x < 0) x *= -1;
  return x;
}

long double s21_fabs(double x) {
  if (x < 0) x *= -1;
  return x;
}

long double s21_ceil(double x) {
  long long unsigned u_x = s21_fabs(x);
  double sing = 1.0;
  if (x < 0) sing = -1.0;
  return sing < 0 ? u_x * sing : u_x * sing + (double)1;
}

long double s21_floor(double x) {
  if (x == s21_INFINITY) return s21_INFINITY;
  if (x == -s21_INFINITY) return -s21_INFINITY;
  long double y = s21_ceil(x);
  return y - 1;
}

long double s21_exp(double x) {
  long double ans = 1.0;
  int n = 1;
  long double eps = 1;
  int sing = 1;
  if (x < 0) {
    sing *= -1;
    x *= -1;
  }
  for (; s21_fabs(eps) > s21_E; n++) {
    eps *= x / n;
    ans += eps;
    if (ans > __DBL_MAX__) {
      ans = s21_INFINITY;
      break;
    }
  }
  return sing < 0 ? 1 / ans : ans;
}

long double s21_log(double x) {
  int ex_pow = 0;
  long double ans = 0;
  long double com = 0;
  long double exp = s21_exp(1);

  if (x == 0) return -s21_INFINITY;
  if (x == s21_INFINITY) return s21_INFINITY;
  if (x < 0) return s21_NAN;

  for (; x >= exp; x /= exp, ex_pow++)
    ;

  for (int i = 0; i < 100; i++) {
    com = ans;
    ans = com + 2 * (x - s21_exp(com)) / (x + s21_exp(com));
  }
  return ans + ex_pow;
}

long double s21_sqrt(double x) {
  long double left = 0;
  long double rigth = MAX(1, x);
  long double mid = (left + rigth) / 2;
  if (x < 0) return s21_NAN;
  if (x == s21_INFINITY) return s21_INFINITY;
  if (x == -s21_INFINITY) return s21_NAN;
  if (x == s21_NAN) return s21_NAN;
  for (; mid - left > s21_E;) {
    if (mid * mid > x) {
      rigth = mid;
    } else {
      left = mid;
    }
    mid = (left + rigth) / 2;
  }
  return mid;
}

long double s21_fmod(double x, double y) {
  if (y == 0) return s21_NAN;
  if (x == s21_INFINITY) return s21_NAN;
  if (x == -s21_INFINITY) return s21_NAN;
  if (x == s21_NAN) return s21_NAN;
  long long int mod;
  mod = x / y;
  long double res = (long double)x - mod * (long double)y;
  return res;
}

long double s21_pow(double base, double exp) {
  long double res = 1.;
  long double copy = base;
  long long int copy_exp_int = (long long int)exp;
  if (exp == 0) return 1;
  if ((exp == s21_INFINITY || exp == -s21_INFINITY) && base != 0 &&
      base != s21_INFINITY && base != -s21_INFINITY) {
    if (base == 1)
      return 1.;
    else if (base == -1) {
      if (exp == -s21_INFINITY)
        return 1.;
      else
        return -1.;
    } else
      return s21_NAN;
  }
  if (base == -1.) {
    if (copy_exp_int % 2 == 0)
      return 1.;
    else
      return -1.;
  }
  if (base == 0)
    if (exp == -1) return s21_INFINITY;

  if (exp == 1) return base;

  // for (long long int i = 0; i < copy_exp_int; i++) res *= base;
  // exp = exp - (long double)copy_exp_int;
  if (copy < 0) {
    copy = -copy;
    res *= s21_exp(exp * s21_log(copy));
    if (s21_fmod(exp, 2) != 0) res = -res;
  } else {
    res *= s21_exp(exp * s21_log(copy));
  }
  return res;
}

int s21_isnan(long double x) { return x == s21_NAN ? 1 : 0; }

long double s21_factor(int n) {
  if (n == 0) return 1;
  long double ans = 1.0;
  for (int i = 1; i <= n; ++i) ans *= (long double)i;
  return ans;
}

long double s21_sin(double x) {
  long double sum_sin = 0.;
  for (; x < -2 * s21_M_PI || 2 * s21_M_PI < x;) {
    if (x > 2 * s21_M_PI)
      x -= 2 * s21_M_PI;
    else
      x += 2 * s21_M_PI;
  }
  for (register int i = 0; i < 100; ++i) {
    long double temp = (i % 2 == 0 ? 1 : -1) * s21_pow(x, 2 * i + 1) /
                       (long double)s21_factor(i * 2 + 1);
    sum_sin += temp;
  }
  return sum_sin;
}

long double s21_atan(double x) {
  long double sum_atan = .0;
  if (s21_fabs(x) < 1.) {
    for (register int i = 0; i < 5000; i++) {
      sum_atan += (i % 2 == 0 ? 1 : -1) * s21_pow(x, 1 + (2 * i)) /
                  (long double)(1 + (2 * i));
    }
  } else {
    for (register int i = 0; i < 7000; i++) {
      sum_atan += (i % 2 == 0 ? 1 : -1) * s21_pow(x, -1 - (2 * i)) /
                  (long double)(1 + (2 * i));
    }
    sum_atan = s21_M_PI * s21_sqrt(x * x) / (2 * x) - sum_atan;
  }
  return sum_atan;
}

long double s21_cos(double x) { return s21_sin(x + s21_M_PI / 2); }

long double s21_tan(double x) {
  if (x == s21_M_PI / 2) return s21_NAN;
  return s21_sin(x) / s21_cos(x);
}

long double s21_acos(double x) {
  if (x == 1) return -s21_M_PI;
  if (x == -1) return s21_M_PI;
  return (long double)s21_M_PI / 2. - s21_asin(x);
}

long double s21_asin(double x) {
  if (x == s21_INFINITY || x == -s21_INFINITY) return s21_NAN;
  long double sum_asin = 0;
  if (x == -1) return -s21_M_PI / 2;
  if (x == 1) return s21_M_PI / 2;
  if (s21_fabs(x) < 1.) {
    for (register int i = 0; i < 500; ++i) {
      sum_asin += s21_factor(2 * i) * s21_pow(x, 2 * i + 1) /
                  (s21_pow(4, i) * s21_pow(s21_factor(i), 2) * (2 * i + 1));
    }
  }
  return sum_asin;
}
