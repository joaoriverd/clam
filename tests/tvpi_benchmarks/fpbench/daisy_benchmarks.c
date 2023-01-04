#include <math.h>
#define TRUE 1
#define FALSE 0

double carthesianToPolar_radius(double x, double y) {
  /* Preconditions */
  __CRAB_assume(x >= 1.0);
  __CRAB_assume(x <= 100.0);
  __CRAB_assume(y >= 1.0);
  __CRAB_assume(y <= 100.0);

  double res = sqrt(((x * x) + (y * y)));
  __CRAB_get_range(res);
  return res;
}

double carthesianToPolar_theta(double x, double y) {
  /* Preconditions */
  __CRAB_assume(x >= 1.0);
  __CRAB_assume(x <= 100.0);
  __CRAB_assume(y >= 1.0);
  __CRAB_assume(y <= 100.0);

	double pi = 3.14159265359;
	double radiant = atan((y / x));
  double res = radiant * (180.0 / pi);
  __CRAB_get_range(res);
  return res;
}

double polarToCarthesian_X(double radius, double theta) {
  /* Preconditions */
  __CRAB_assume(radius >= 1.0);
  __CRAB_assume(radius <= 10.0);
  __CRAB_assume(theta  >= 0.0);
  __CRAB_assume(theta  <= 360.0);

	double pi = 3.14159265359;
	double radiant = theta * (pi / 180.0);
  double res = radius * cos(radiant);
  __CRAB_get_range(res);
  return res;
}

double polarToCarthesian_Y(double radius, double theta) {
  /* Preconditions */
  __CRAB_assume(radius >= 1.0);
  __CRAB_assume(radius <= 10.0);
  __CRAB_assume(theta  >= 0.0);
  __CRAB_assume(theta  <= 360.0);

  double pi = 3.14159265359;
	double radiant = theta * (pi / 180.0);
  double res = radius * sin(radiant);
  __CRAB_get_range(res);
  return res;
}

double instantaneousCurrent(double t, double resistance, double frequency, double inductance, double maxVoltage) {
  /* Preconditions */
  __CRAB_assume(t >= 0.0);
  __CRAB_assume(t <= 300.0);
  __CRAB_assume(resistance >= 1.0);
  __CRAB_assume(resistance <= 50.0);
  __CRAB_assume(frequency  >= 1.0);
  __CRAB_assume(frequency  <= 100.0);
  __CRAB_assume(inductance  >= 0.001);
  __CRAB_assume(inductance  <= 0.004);
  __CRAB_assume(maxVoltage  >= 1.0);
  __CRAB_assume(maxVoltage  <= 12.0);

	double pi = 3.14159265359;
	double impedance_re = resistance;
	double impedance_im = ((2.0 * pi) * frequency) * inductance;
	double denom = (impedance_re * impedance_re) + (impedance_im * impedance_im);
	double current_re = (maxVoltage * impedance_re) / denom;
	double current_im = -(maxVoltage * impedance_im) / denom;
	double maxCurrent = sqrt(((current_re * current_re) + (current_im * current_im)));
	double theta = atan((current_im / current_re));
  double res = maxCurrent * cos(((((2.0 * pi) * frequency) * t) + theta));
  __CRAB_get_range(res);
  return res;
}

double matrixDeterminant(double a, double b, double c, double d, double e, double f, double g, double h, double i) {
  /* Preconditions */
  __CRAB_assume(a >= -10.0); __CRAB_assume(b >= -10.0); __CRAB_assume(c >= -10.0);
  __CRAB_assume(a <=  10.0); __CRAB_assume(b <=  10.0); __CRAB_assume(c <=  10.0);
  __CRAB_assume(d >= -10.0); __CRAB_assume(e >= -10.0); __CRAB_assume(f >= -10.0);
  __CRAB_assume(d <=  10.0); __CRAB_assume(e <=  10.0); __CRAB_assume(f <=  10.0);
  __CRAB_assume(g >= -10.0); __CRAB_assume(h >= -10.0); __CRAB_assume(i >= -10.0);
  __CRAB_assume(g <=  10.0); __CRAB_assume(h <=  10.0); __CRAB_assume(i <=  10.0);

  double res = ((((a * e) * i) + ((b * f) * g)) + ((c * d) * h)) - ((((c * e) * g) + ((b * d) * i)) + ((a * f) * h));
  __CRAB_get_range(res);
  return res;
}

double matrixDeterminant2(double a, double b, double c, double d, double e, double f, double g, double h, double i) {
  /* Preconditions */
  __CRAB_assume(a >= -10.0); __CRAB_assume(b >= -10.0); __CRAB_assume(c >= -10.0);
  __CRAB_assume(a <=  10.0); __CRAB_assume(b <=  10.0); __CRAB_assume(c <=  10.0);
  __CRAB_assume(d >= -10.0); __CRAB_assume(e >= -10.0); __CRAB_assume(f >= -10.0);
  __CRAB_assume(d <=  10.0); __CRAB_assume(e <=  10.0); __CRAB_assume(f <=  10.0);
  __CRAB_assume(g >= -10.0); __CRAB_assume(h >= -10.0); __CRAB_assume(i >= -10.0);
  __CRAB_assume(g <=  10.0); __CRAB_assume(h <=  10.0); __CRAB_assume(i <=  10.0);

  double res = ((a * (e * i)) + ((g * (b * f)) + (c * (d * h)))) - ((e * (c * g)) + ((i * (b * d)) + (a * (f * h))));
  __CRAB_get_range(res);
  return res;
}

