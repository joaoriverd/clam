// RUN: %clam -O0 --crab-dom=oct --crab-check=assert --crab-sanity-checks "%s" 2>&1 | OutputCheck %s
// CHECK: ^1  Number of total safe checks$
// CHECK: ^0  Number of total error checks$
// CHECK: ^0  Number of total warning checks$

extern void __CRAB_assert(int);
extern void __CRAB_assume(int);
extern void __CRAB_get_range(double);
extern void __SEAHORN_error(int);

#include <math.h>
#define TRUE  1
#define FALSE 0

double doppler1(double u, double v, double T) {
  /* Preconditions */
  __CRAB_assume(u >= -100.0);
  __CRAB_assume(u <=  100.0);
  __CRAB_assume(v >=  20.0);
  __CRAB_assume(v <=  20000.0);
  __CRAB_assume(T >= -30.0);
  __CRAB_assume(T <=  50.0);

  double t1 = 331.4 + (0.6 * T);
  double res = (-t1 * v) / ((t1 + u) * (t1 + u));
  __CRAB_get_range(res);
  return res;
}

double doppler2(double u, double v, double T) {
  /* Preconditions */
  __CRAB_assume(u >= -125.0);
  __CRAB_assume(u <=  125.0);
  __CRAB_assume(v >=  15.0);
  __CRAB_assume(v <=  25000.0);
  __CRAB_assume(T >= -40.0);
  __CRAB_assume(T <=  60.0);

  double t1 = 331.4 + (0.6 * T);
  double res = (-t1 * v) / ((t1 + u) * (t1 + u));
  __CRAB_get_range(res);
  return res;
}

double doppler3(double u, double v, double T) {
  /* Preconditions */
  __CRAB_assume(u >= -30.0);
  __CRAB_assume(u <=  120.0);
  __CRAB_assume(v >=  320.0);
  __CRAB_assume(v <=  20300.0);
  __CRAB_assume(T >= -50.0);
  __CRAB_assume(T <=  30.0);

  double t1 = 331.4 + (0.6 * T);
  double res = (-t1 * v) / ((t1 + u) * (t1 + u));
  __CRAB_get_range(res);
  return res;
}

double rigidBody1(double x1, double x2, double x3) {
  /* Preconditions */
  __CRAB_assume(x1 >= -15.0);
  __CRAB_assume(x1 <=  15.0);
  __CRAB_assume(x2 >= -15.0);
  __CRAB_assume(x2 <=  15.0);
  __CRAB_assume(x3 >= -15.0);
  __CRAB_assume(x3 <=  15.0);

  double res = ((-(x1 * x2) - ((2.0 * x2) * x3)) - x1) - x3;
  __CRAB_get_range(res);
  return res;
}

double rigidBody2(double x1, double x2, double x3) {
  /* Preconditions */
  __CRAB_assume(x1 >= -15.0);
  __CRAB_assume(x1 <=  15.0);
  __CRAB_assume(x2 >= -15.0);
  __CRAB_assume(x2 <=  15.0);
  __CRAB_assume(x3 >= -15.0);
  __CRAB_assume(x3 <=  15.0);

  double res = ((((((2.0 * x1) * x2) * x3) + ((3.0 * x3) * x3)) - (((x2 * x1) * x2) * x3)) + ((3.0 * x3) * x3)) - x2;
  __CRAB_get_range(res);
  return res;
}

double jetEngine(double x1, double x2) {
  /* Preconditions */
  __CRAB_assume(x1 >= -5.0);
  __CRAB_assume(x1 <=  5.0);
  __CRAB_assume(x2 >= -20.0);
  __CRAB_assume(x2 <=  5.0);

  double t = (((3.0 * x1) * x1) + (2.0 * x2)) - x1;
  double t_42_ = (((3.0 * x1) * x1) - (2.0 * x2)) - x1;
  double d = (x1 * x1) + 1.0;
  double s = t / d;
  double s_42_ = t_42_ / d;
  double res = x1 + (((((((((2.0 * x1) * s) * (s - 3.0)) + ((x1 * x1) * ((4.0 * s) - 6.0))) * d) + (((3.0 * x1) * x1) * s)) + ((x1 * x1) * x1)) + x1) + (3.0 * s_42_));
  __CRAB_get_range(res);
  return res;
}

double turbine1(double v, double w, double r) {
  /* Preconditions */
  __CRAB_assume(v >= -4.5);
  __CRAB_assume(v <= -0.3);
  __CRAB_assume(w >=  0.4);
  __CRAB_assume(w <=  0.9);
  __CRAB_assume(r >=  3.8);
  __CRAB_assume(r <=  7.8);

  double res = ((3.0 + (2.0 / (r * r))) - (((0.125 * (3.0 - (2.0 * v))) * (((w * w) * r) * r)) / (1.0 - v))) - 4.5;
  __CRAB_get_range(res);
  return res;
}

double turbine2(double v, double w, double r) {
  /* Preconditions */
  __CRAB_assume(v >= -4.5);
  __CRAB_assume(v <= -0.3);
  __CRAB_assume(w >=  0.4);
  __CRAB_assume(w <=  0.9);
  __CRAB_assume(r >=  3.8);
  __CRAB_assume(r <=  7.8);


  double res = ((6.0 * v) - (((0.5 * v) * (((w * w) * r) * r)) / (1.0 - v))) - 2.5;
  __CRAB_get_range(res);
  return res;
}

double turbine3(double v, double w, double r) {
  /* Preconditions */
  __CRAB_assume(v >= -4.5);
  __CRAB_assume(v <= -0.3);
  __CRAB_assume(w >=  0.4);
  __CRAB_assume(w <=  0.9);
  __CRAB_assume(r >=  3.8);
  __CRAB_assume(r <=  7.8);


  double res = ((3.0 - (2.0 / (r * r))) - (((0.125 * (1.0 + (2.0 * v))) * (((w * w) * r) * r)) / (1.0 - v))) - 0.5;
  __CRAB_get_range(res);
  return res;
}

double verhulst(double x) {
  /* Preconditions */
  __CRAB_assume(x >= 0.1);
  __CRAB_assume(x <= 0.3);

  double r = 4.0;
  double K = 1.11;
  double res = (r * x) / (1.0 + (x / K));
  __CRAB_get_range(res);
  return res;
}

double predatorPrey(double x) {
  /* Preconditions */
  __CRAB_assume(x >= 0.1);
  __CRAB_assume(x <= 0.3);

  double r = 4.0;
  double K = 1.11;
  double res = ((r * x) * x) / (1.0 + ((x / K) * (x / K)));
  __CRAB_get_range(res);
  return res;
}

double carbonGas(double v) {
  /* Preconditions */
  __CRAB_assume(v >= 0.1);
  __CRAB_assume(v <= 0.5);

  double p = 35000000.0;
  double a = 0.401;
  double b = 4.27e-5;
  double t = 300.0;
  double n = 1000.0;
  double k = 1.3806503e-23;

  double res = ((p + ((a * (n / v)) * (n / v))) * (v - 4.27e-2)) - 4.1419509-18;
  __CRAB_get_range(res);
  return res;
}

double sine(double x) {
  /* Preconditions */
  __CRAB_assume(x >= -1.57079632679);
  __CRAB_assume(x <=  1.57079632679);

  double res = ((x - (((x * x) * x) / 6.0)) + (((((x * x) * x) * x) * x) / 120.0)) - (((((((x * x) * x) * x) * x) * x) * x) / 5040.0);
  __CRAB_get_range(res);
  return res;
}

double sqroot(double x) {
  /* Preconditions */
  __CRAB_assume(x >= 0);
  __CRAB_assume(x <= 1);

  double res = (((1.0 + (0.5 * x)) - ((0.125 * x) * x)) + (((0.0625 * x) * x) * x)) - ((((0.0390625 * x) * x) * x) * x);
  __CRAB_get_range(res);
  return res;
}

double sineOrder3(double x) {
  /* Preconditions */
  __CRAB_assume(x >= -2.0);
  __CRAB_assume(x <=  2.0);

  double res = (0.954929658551372 * x) - (0.12900613773279798 * ((x * x) * x));
  __CRAB_get_range(res);
  return res;
}

#if 0 //todo: need fix?
double smartRoot(double c) {
  /* Preconditions */
  __CRAB_assume(c >= -2.0);
  __CRAB_assume(c <=  2.0);

  double a = 3.0;
  double b = 3.5;
  double discr = (b * b) - ((a * c) * 4.0);
  __CRAB_assume(discr >= 0.1);

  double tmp_1;
  if (((b * b) - (a * c)) > 10.0) {
    double tmp_2;
    if (b > 0.0) {
      tmp_2 = (c * 2.0) / (-b - sqrt(discr));
    } else {
      tmp_2 = (-b + sqrt(discr)) / (a * 2.0);
    }
    tmp_1 = tmp_2;
  } else {
    tmp_1 = (-b + sqrt(discr)) / (a * 2.0);
  }

  __CRAB_get_range(tmp_1);
  return tmp_1;
}
#endif

double smartRoot_simplified(double c) {
  /* Preconditions */
  __CRAB_assume(c >= -2.0);
  __CRAB_assume(c <=  2.0);
  double discr = 12.25 - 12.0 * c;
  __CRAB_assume(discr > 0.1);
  double tmp_1;
  __CRAB_get_range(discr);
  __CRAB_get_range(c);
//  if (3.0 * c < 2.25) {
  if (c < 0.75) {
    tmp_1 = (c * 2.0) / (-3.5 - sqrt(discr));
  } else {
    tmp_1 = (-3.5 + sqrt(discr)) / 6.0;
  }
  __CRAB_get_range(tmp_1);
  return tmp_1;
}

double cav10(double x) {
  /* Preconditions */
  __CRAB_assume(x >= 0.0);
  __CRAB_assume(x <= 10.0);
  /* Expected result for res:
  IA:   [0, 102]
  TVPI: [0, 12]
  Rosa: [0, 3]
  */

  double tmp;
  double sqr_x = x*x;
  if (sqr_x - x >= 0.0) {
    tmp = x / 10.0;
  } else {
    tmp = sqr_x + 2.0;
  }
  __CRAB_get_range(tmp);
  return tmp;
}

double squareRoot3(double x) {
  __CRAB_assume(x >= 0.0);
  __CRAB_assume(x <= 10.0);
  /* Expected result for res:
     TVPI: [1, 3.31663]
  */

  double tmp;
  if (x < 1e-5) {
    tmp = 1.0 + (0.5 * x);
  } else {
    tmp = sqrt((1.0 + x));
  }
  __CRAB_get_range(tmp);
  return tmp;
}

double squareRoot3Invalid(double x) {
  __CRAB_assume(x >= 0.0);
  __CRAB_assume(x <= 10.0);
  /* Expected result for res:
     TVPI: [1, 3.31663]
  */
  double tmp;
  if (x < 0.0001) {
    tmp = 1.0 + (0.5 * x);
  } else {
    tmp = sqrt((1.0 + x));
  }
  __CRAB_get_range(tmp);
  return tmp;
}

double triangle(double a, double b, double c) {
  /* Preconditions */
  a = 9.0;
  __CRAB_assume(b >= 4.71);
  __CRAB_assume(b <= 4.89);
  __CRAB_assume(c >= 4.71);
  __CRAB_assume(c <= 4.89);
  /* Expected result for res:
     TVPI: [6.23682, 8.65856]
  */

  double s = ((a + b) + c) / 2.0;
  double res = sqrt(s * (s - a) * (s - b) * (s - c));
  __CRAB_get_range(res);
  return res;
}

double triangle1(double a, double b, double c) {
  /* Preconditions */
  __CRAB_assume(a >= 1.0);
  __CRAB_assume(a <= 9.0);
  __CRAB_assume(b >= 1.0);
  __CRAB_assume(b <= 9.0);
  __CRAB_assume(c >= 1.0);
  __CRAB_assume(c <= 9.0);
  /* Expected result for res:
     TVPI:
  */
  __CRAB_assume(a+b > c + 0.1);
  __CRAB_assume(a+c > b + 0.1);
  __CRAB_assume(b+c > a + 0.1);

  double s = (a + b + c) / 2.0;
  __CRAB_get_range(s);
  double temp = s * (s - a) * (s - b) * (s - c);
  __CRAB_get_range(temp);
  double res = sqrt(temp);
  __CRAB_get_range(res);
  return res;
}

#if 0
double ex21(double a, double b, double c) {
  double s = ((a + b) + c) / 2.0;
  return sqrt((((s * (s - a)) * (s - b)) * (s - c)));
}

double ex22(double a, double b, double c) {
  double s = ((a + b) + c) / 2.0;
  return sqrt((((s * (s - a)) * (s - b)) * (s - c)));
}

double ex23(double a, double b, double c) {
  double s = ((a + b) + c) / 2.0;
  return sqrt((((s * (s - a)) * (s - b)) * (s - c)));
}

double ex24(double a, double b, double c) {
  double s = ((a + b) + c) / 2.0;
  return sqrt((((s * (s - a)) * (s - b)) * (s - c)));
}

double ex25(double a, double b, double c) {
  double s = ((a + b) + c) / 2.0;
  return sqrt((((s * (s - a)) * (s - b)) * (s - c)));
}

double ex26(double a, double b, double c) {
  double s = ((a + b) + c) / 2.0;
  return sqrt((((s * (s - a)) * (s - b)) * (s - c)));
}

double ex27(double a, double b, double c) {
  double s = ((a + b) + c) / 2.0;
  return sqrt((((s * (s - a)) * (s - b)) * (s - c)));
}

double ex28(double a, double b, double c) {
  double s = ((a + b) + c) / 2.0;
  return sqrt((((s * (s - a)) * (s - b)) * (s - c)));
}

double ex29(double a, double b, double c) {
  double s = ((a + b) + c) / 2.0;
  return sqrt((((s * (s - a)) * (s - b)) * (s - c)));
}

double ex30(double a, double b, double c) {
  double s = ((a + b) + c) / 2.0;
  return sqrt((((s * (s - a)) * (s - b)) * (s - c)));
}

double ex31(double a, double b, double c) {
  double s = ((a + b) + c) / 2.0;
  return sqrt((((s * (s - a)) * (s - b)) * (s - c)));
}
#endif

double bspline3(double u) {
  __CRAB_assume(u >= 0);
  __CRAB_assume(u <= 1);
  double res = -((u * u) * u) / 6.0;
  __CRAB_get_range(res);
  return res;
}

#if 0
double triangleSorted(double a, double b, double c) {
  double tmp;
  if (a < b) {
    tmp = sqrt(((((c + (b + a)) * (a - (c - b))) * (a + (c - b))) * (c + (b - a)))) / 4.0;
  } else {
    tmp = sqrt(((((c + (a + b)) * (b - (c - a))) * (b + (c - a))) * (c + (a - b)))) / 4.0;
  }
  return tmp;
}

double n_body_simulation(double x0, double y0, double z0, double vx0, double vy0, double vz0) {
  double dt = 0.1;
  double solarMass = 39.47841760435743;
  double x = x0;
  double y = y0;
  double z = z0;
  double vx = vx0;
  double vy = vy0;
  double vz = vz0;
  double i = 0.0;
  int tmp = i < 100.0;
  while (tmp) {
    double distance = sqrt(x*x + y*y + z*z);
    double mag = dt / (distance * distance * distance);
    double vxNew = vx - (x * solarMass * mag);
    double x_1 = x + (dt * vxNew);

    double distance_2 = sqrt(x*x + y*y + z*z);
    double mag_3 = dt / (distance_2 * distance_2 * distance_2);
    double vyNew = vy - (y * solarMass * mag_3);
    double y_4 = y + (dt * vyNew);

    double distance_5 = sqrt(x*x + y*y + z*z);
    double mag_6 = dt / (distance_5 * distance_5 * distance_5);
    double vzNew = vz - (z * solarMass * mag_6);
    double z_7 = z + (dt * vzNew);

    double distance_8 = sqrt(x*x + y*y + z*z);
    double mag_9 = dt / (distance_8 * distance_8 * distance_8);
    double vx_10 = vx - (x * solarMass * mag_9);

    double distance_11 = sqrt(x*x + y*y + z*z);
    double mag_12 = dt / (distance_11 * distance_11 * distance_11);
    double vy_13 = vy - (y * solarMass * mag_12);

    double distance_14 = sqrt(x*x + y*y + z*z);
    double mag_15 = dt / (distance_14 * distance_14 * distance_14);
    double vz_16 = vz - (z * solarMass * mag_15);

    double i_17 = i + 1.0;

    x = x_1;
    y = y_4;
    z = z_7;
    vx = vx_10;
    vy = vy_13;
    vz = vz_16;
    i = i_17;
    tmp = i < 100.0;
  }
  __CRAB_get_range(x);
  return x;
}
#endif

double n_body_simulation_simplified(double x0, double y0, double z0, double vx0, double vy0, double vz0) {
  /* Notes:
   * Final result yields [-inf, inf]. This seems correct since a division by a possible zero will occur when
   * calculating "mag". */
  __CRAB_assume(x0 >= -6.0);
  __CRAB_assume(x0 <=  6.0);
  __CRAB_assume(y0 >= -6.0);
  __CRAB_assume(y0 <=  6.0);
  __CRAB_assume(z0 >= -0.2);
  __CRAB_assume(z0 <=  0.2);
  __CRAB_assume(vx0 >= -3.0);
  __CRAB_assume(vx0 <=  3.0);
  __CRAB_assume(vy0 >= -3.0);
  __CRAB_assume(vy0 <=  3.0);
  __CRAB_assume(vz0 >= -0.1);
  __CRAB_assume(vz0 <=  0.1);

  double dt = 0.1;
  double solarMass = 39.47841760435743;
  double x = x0;
  double y = y0;
  double z = z0;
  double vx = vx0;
  double vy = vy0;
  double vz = vz0;
  double i = 0.0;
  int tmp = i < 100.0;
  while (tmp) {
    double t1 = x*x + y*y + z*z;
    double distance = sqrt(t1);
    double mag = dt / (distance * distance * distance);

    double vx_10 = vx - (x * solarMass * mag);
    double vy_13 = vy - (y * solarMass * mag);
    double vz_16 = vz - (z * solarMass * mag);

    double x_1 = x + (dt * vx_10);
    double y_4 = y + (dt * vy_13);
    double z_7 = z + (dt * vz_16);

    double i_17 = i + 1.0;

    x = x_1;
    y = y_4;
    z = z_7;
    vx = vx_10;
    vy = vy_13;
    vz = vz_16;
    i = i_17;
    tmp = i < 100.0;
  }
  __CRAB_get_range(x);
  return x;
}

double Pendulum(double t0, double w0, double N) {
  /* Notes:
   * Final result yields [-inf, inf]. This seems correct since every iteration will increase the range of t and w.
   * The widening operator will make these variabels [-inf, inf]. */
  __CRAB_assume(t0 >= -2.0);
  __CRAB_assume(t0 <=  2.0);
  __CRAB_assume(w0 >= -5.0);
  __CRAB_assume(w0 <=  5.0);
  double t = t0;
  double w = w0;
  double n = 0.0;
  while (n < N) {
    double k1w = (-4.903325) * sin(t);
    double k2w = (-4.903325) * sin(t + 0.005*w);
    double k2t = w + 0.005*k1w;
    t = t + (0.01 * k2t);
    w = w + (0.01 * k2w);
    n = n + 1.0;
  }
  __CRAB_get_range(t);
  return t;
}

#if 0
double ex36(double x0) {
  double x = x0;
  double i = 0.0;
  int tmp = i < 10.0;
  while (tmp) {
    double x_1 = x - ((((x - (pow(x, 3.0) / 6.0)) + (pow(x, 5.0) / 120.0)) + (pow(x, 7.0) / 5040.0)) / (((1.0 - ((x * x) / 2.0)) + (pow(x, 4.0) / 24.0)) + (pow(x, 6.0) / 720.0)));
    double i_2 = i + 1.0;
    x = x_1;
    i = i_2;
    tmp = i < 10.0;
  }
  return x;
}
#endif