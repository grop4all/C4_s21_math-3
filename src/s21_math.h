#include <limits.h>
#include <stdio.h>
#include <stdlib.h>

#define s21_E 1e-9
#define s21_INFINITY 1.0 / 0.0
#define s21_NAN 0.0 / 0.0
#define MAX(a, b) a > b ? a : b
#define s21_M_PI 3.1415926535897932384626433832795028

int s21_abs(int x);
long double s21_fabs(double x);
long double s21_ceil(double x);
long double s21_floor(double x);
long double s21_exp(double x);
long double s21_log(double x);
long double s21_pow(double base, double exp);
long double s21_fmod(double x, double y);
long double s21_sqrt(double x);

long double s21_acos(double x);
long double s21_asin(double x);
long double s21_atan(double x);

long double s21_cos(double x);
long double s21_sin(double x);
long double s21_tan(double x);
long double s21_factor(int n);

int s21_isnan(long double x);
