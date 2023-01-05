#include <math.h>
#define TRUE 1
#define FALSE 0

extern void __CRAB_assert(int);
extern void __CRAB_assume(int);
extern void __CRAB_get_range(double);
extern void __SEAHORN_error(int);

double logexp(double x) {
  /* Preconditions */
  __CRAB_assume(x >= -8.0);
  __CRAB_assume(x <=  8.0);

	double e = exp(x);
  double res = log((1.0 + e));
  __CRAB_get_range(res);
  return res;
}

double sphere(double x, double r, double lat, double lon) {
  /* Preconditions */
  __CRAB_assume(x >= -10.0);
  __CRAB_assume(x <=  10.0);
  __CRAB_assume(r >=  0.0);
  __CRAB_assume(r <=  10.0);
  __CRAB_assume(lat >= -1.570796);
  __CRAB_assume(lat <=  1.570796);
  __CRAB_assume(lon >= -3.14159265);
  __CRAB_assume(lon <=  3.14159265);

	double sinLat = sin(lat);
	double cosLon = cos(lon);
  double res = x + ((r * sinLat) * cosLon);
  __CRAB_get_range(res);
  return res;
}

double azimuth(double lat1, double lat2, double lon1, double lon2) {
  /* Preconditions */
  __CRAB_assume(lat1 >=  0.0);
  __CRAB_assume(lat1 <=  0.4);
  __CRAB_assume(lat2 >=  0.5);
  __CRAB_assume(lat2 <=  1.0);
  __CRAB_assume(lon1 >=  0.0);
  __CRAB_assume(lon1 <=  3.14159265);
  __CRAB_assume(lon2 >= -3.14159265);
  __CRAB_assume(lon2 <= -0.5);

	double dLon = lon2 - lon1;
	double s_lat1 = sin(lat1);
	double c_lat1 = cos(lat1);
	double s_lat2 = sin(lat2);
	double c_lat2 = cos(lat2);
	double s_dLon = sin(dLon);
	double c_dLon = cos(dLon);
  double res = atan(((c_lat2 * s_dLon) / ((c_lat1 * s_lat2) - ((s_lat1 * c_lat2) * c_dLon))));
  __CRAB_get_range(res);
  return res;
}

double floudas1(double x1, double x2, double x3, double x4, double x5, double x6) {
  double res =  (((((-25.0 * ((x1 - 2.0) * (x1 - 2.0))) - ((x2 - 2.0) * (x2 - 2.0))) - ((x3 - 1.0) * (x3 - 1.0))) - ((x4 - 4.0) * (x4 - 4.0))) - ((x5 - 1.0) * (x5 - 1.0))) - ((x6 - 4.0) * (x6 - 4.0));
  __CRAB_get_range(res);
  return res;
}

double floudas2(double x1, double x2) {
  double res = -x1 - x2;
  __CRAB_get_range(res);
  return res;
}

double floudas3(double x1, double x2) {
  double res = ((-12.0 * x1) - (7.0 * x2)) + (x2 * x2);
  __CRAB_get_range(res);
  return res;
}

double hartman3(double x1, double x2, double x3) {
  /* Preconditions */
  __CRAB_assume(x1 >=  0.0);
  __CRAB_assume(x1 <=  1.0);
  __CRAB_assume(x2 >=  0.0);
  __CRAB_assume(x2 <=  1.0);
  __CRAB_assume(x3 >=  0.0);
  __CRAB_assume(x3 <=  1.0);

	double e1 = ((3.0 * ((x1 - 0.3689) * (x1 - 0.3689))) + (10.0 * ((x2 - 0.117) * (x2 - 0.117)))) + (30.0 * ((x3 - 0.2673) * (x3 - 0.2673)));
	double e2 = ((0.1 * ((x1 - 0.4699) * (x1 - 0.4699))) + (10.0 * ((x2 - 0.4387) * (x2 - 0.4387)))) + (35.0 * ((x3 - 0.747) * (x3 - 0.747)));
	double e3 = ((3.0 * ((x1 - 0.1091) * (x1 - 0.1091))) + (10.0 * ((x2 - 0.8732) * (x2 - 0.8732)))) + (30.0 * ((x3 - 0.5547) * (x3 - 0.5547)));
	double e4 = ((0.1 * ((x1 - 0.03815) * (x1 - 0.03815))) + (10.0 * ((x2 - 0.5743) * (x2 - 0.5743)))) + (35.0 * ((x3 - 0.8828) * (x3 - 0.8828)));
	double exp1 = exp(-e1);
	double exp2 = exp(-e2);
	double exp3 = exp(-e3);
	double exp4 = exp(-e4);
  double res = -((((1.0 * exp1) + (1.2 * exp2)) + (3.0 * exp3)) + (3.2 * exp4));
  __CRAB_get_range(res);
  return res;
}

double hartman6(double x1, double x2, double x3, double x4, double x5, double x6) {
  /* Preconditions */
  __CRAB_assume(x1 >=  0.0);  __CRAB_assume(x4 >=  0.0);
  __CRAB_assume(x1 <=  1.0);  __CRAB_assume(x4 <=  1.0);
  __CRAB_assume(x2 >=  0.0);  __CRAB_assume(x5 >=  0.0);
  __CRAB_assume(x2 <=  1.0);  __CRAB_assume(x5 <=  1.0);
  __CRAB_assume(x3 >=  0.0);  __CRAB_assume(x6 >=  0.0);
  __CRAB_assume(x3 <=  1.0);  __CRAB_assume(x6 <=  1.0);

	double e1 = (((((10.0 * ((x1 - 0.1312) * (x1 - 0.1312))) + (3.0 * ((x2 - 0.1696) * (x2 - 0.1696)))) + (17.0 * ((x3 - 0.5569) * (x3 - 0.5569)))) + (3.5 * ((x4 - 0.0124) * (x4 - 0.0124)))) + (1.7 * ((x5 - 0.8283) * (x5 - 0.8283)))) + (8.0 * ((x6 - 0.5886) * (x6 - 0.5886)));
	double e2 = (((((0.05 * ((x1 - 0.2329) * (x1 - 0.2329))) + (10.0 * ((x2 - 0.4135) * (x2 - 0.4135)))) + (17.0 * ((x3 - 0.8307) * (x3 - 0.8307)))) + (0.1 * ((x4 - 0.3736) * (x4 - 0.3736)))) + (8.0 * ((x5 - 0.1004) * (x5 - 0.1004)))) + (14.0 * ((x6 - 0.9991) * (x6 - 0.9991)));
	double e3 = (((((3.0 * ((x1 - 0.2348) * (x1 - 0.2348))) + (3.5 * ((x2 - 0.1451) * (x2 - 0.1451)))) + (1.7 * ((x3 - 0.3522) * (x3 - 0.3522)))) + (10.0 * ((x4 - 0.2883) * (x4 - 0.2883)))) + (17.0 * ((x5 - 0.3047) * (x5 - 0.3047)))) + (8.0 * ((x6 - 0.665) * (x6 - 0.665)));
	double e4 = (((((17.0 * ((x1 - 0.4047) * (x1 - 0.4047))) + (8.0 * ((x2 - 0.8828) * (x2 - 0.8828)))) + (0.05 * ((x3 - 0.8732) * (x3 - 0.8732)))) + (10.0 * ((x4 - 0.5743) * (x4 - 0.5743)))) + (0.1 * ((x5 - 0.1091) * (x5 - 0.1091)))) + (14.0 * ((x6 - 0.0381) * (x6 - 0.0381)));
	double exp1 = exp(-e1);
	double exp2 = exp(-e2);
	double exp3 = exp(-e3);
	double exp4 = exp(-e4);
  double res = -((((1.0 * exp1) + (1.2 * exp2)) + (3.0 * exp3)) + (3.2 * exp4));
  __CRAB_get_range(res);
  return res;
}

double kepler0(double x1, double x2, double x3, double x4, double x5, double x6) {
  /* Preconditions */
  __CRAB_assume(x1 >=  4.00); __CRAB_assume(x4 >=  4.00);
  __CRAB_assume(x1 <=  6.36); __CRAB_assume(x4 <=  6.36);
  __CRAB_assume(x2 >=  4.00); __CRAB_assume(x5 >=  4.00);
  __CRAB_assume(x2 <=  6.36); __CRAB_assume(x5 <=  6.36);
  __CRAB_assume(x3 >=  4.00); __CRAB_assume(x6 >=  4.00);
  __CRAB_assume(x3 <=  6.36); __CRAB_assume(x6 <=  6.36);

	double res = ((((x2 * x5) + (x3 * x6)) - (x2 * x3)) - (x5 * x6)) + (x1 * (((((-x1 + x2) + x3) - x4) + x5) + x6));
  __CRAB_get_range(res);
  return res;
}

double kepler1(double x1, double x2, double x3, double x4) {
  /* Preconditions */
  __CRAB_assume(x1 >=  4.00);
  __CRAB_assume(x1 <=  6.36);
  __CRAB_assume(x2 >=  4.00);
  __CRAB_assume(x2 <=  6.36);
  __CRAB_assume(x3 >=  4.00);
  __CRAB_assume(x3 <=  6.36);
  __CRAB_assume(x4 >=  4.00);
  __CRAB_assume(x4 <=  6.36);

  double res = (((((((x1 * x4) * (((-x1 + x2) + x3) - x4)) + (x2 * (((x1 - x2) + x3) + x4))) + (x3 * (((x1 + x2) - x3) + x4))) - ((x2 * x3) * x4)) - (x1 * x3)) - (x1 * x2)) - x4;
  __CRAB_get_range(res);
  return res;
}

double kepler2(double x1, double x2, double x3, double x4, double x5, double x6) {
  /* Preconditions */
  __CRAB_assume(x1 >=  4.00); __CRAB_assume(x4 >=  4.00);
  __CRAB_assume(x1 <=  6.36); __CRAB_assume(x4 <=  6.36);
  __CRAB_assume(x2 >=  4.00); __CRAB_assume(x5 >=  4.00);
  __CRAB_assume(x2 <=  6.36); __CRAB_assume(x5 <=  6.36);
  __CRAB_assume(x3 >=  4.00); __CRAB_assume(x6 >=  4.00);
  __CRAB_assume(x3 <=  6.36); __CRAB_assume(x6 <=  6.36);

  double res = (((((((x1 * x4) * (((((-x1 + x2) + x3) - x4) + x5) + x6)) + ((x2 * x5) * (((((x1 - x2) + x3) + x4) - x5) + x6))) + ((x3 * x6) * (((((x1 + x2) - x3) + x4) + x5) - x6))) - ((x2 * x3) * x4)) - ((x1 * x3) * x5)) - ((x1 * x2) * x6)) - ((x4 * x5) * x6);
  __CRAB_get_range(res);
  return res;
}
